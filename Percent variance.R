
library(tidyverse); library(ggregplot); library(MCMCglmm)

setwd("~/Github/clover")

clo <- read_csv("./output/Clover_v1.0_NBCIreconciled_20201211.csv")

clo %>% 
  select(Database, Host_Original, Virus_Original) %>%
  unique %>%
  count(Database, Host_Original) %>%
  arrange(Host_Original) -> raw

clo %>% 
  select(Database, Host, Virus) %>%
  unique %>%
  count(Database, Host) -> clean

clean %>% pull(n) %>% qplot()
raw %>% pull(n) %>% qplot()

clean %>% 
  pivot_wider(names_from = Database, values_from = n) %>% 
  as.data.frame() %>% 
  pivot_longer(EID2:Shaw, names_to = "Database", values_to ="n") %>% 
  mutate(n=replace_na(n,0)) -> clean_pseudo

raw %>% 
  pivot_wider(names_from = Database, values_from = n) %>% 
  as.data.frame() %>% 
  pivot_longer(EID2:GMPD2, names_to = "Database", values_to ="n") %>% 
  mutate(n=replace_na(n,0)) -> raw_pseudo


# a simple model  ------------------------------------------------------------


resp<-"n"
fixed <- "Database"
mcmc <- as.formula(paste(resp, "~", paste(fixed)))


Prior2 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G2 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100)))

Prior3 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G2 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G3 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100)))


mf<-40


## no zeroes 
cleanMod <- MCMCglmm(fixed= n ~ 1, 
                     random= ~ Database + Host , 
                     prior=Prior2,
                     data=clean,
                     family = "poisson",
                     #pr=TRUE, 
                     nitt = 13000*mf,
                     thin = 10*mf,burnin=3000*mf)

rawMod <- MCMCglmm(fixed= n ~ 1, 
                   random= ~ Database + Host_Original , 
                   prior=Prior2,
                   data=raw,
                   family = "poisson",
                   #pr=TRUE, 
                   nitt = 13000*mf,
                   thin = 10*mf,burnin=3000*mf)


summary(cleanMod)
summary(rawMod)

plot(cleanMod)
plot(rawMod)

clean <- MCMCRep(cleanMod)  %>%  # database 5.06%  
  mutate(data="clean") 

raw <- MCMCRep(rawMod)  %>%   # database 5.4% 
  mutate(data="raw") %>% 
  mutate(Component = str_replace(Component, "Host_Original", "Host"))

propVar <- bind_rows(clean, raw)

saveRDS(cleanMod, "Outputs/cleanDataMCMC.rds")
saveRDS(rawMod, "Outputs/rawDataMCMC.rds")

## with pseudoabsences 
cleanModPseudo <- MCMCglmm(fixed= n ~ 1, 
                     random= ~ Database + Host , 
                     prior=Prior2,
                     data=clean_pseudo,
                     family = "poisson",
                     #pr=TRUE, 
                     nitt = 13000*mf,
                     thin = 10*mf,burnin=3000*mf)

rawModPseudo <- MCMCglmm(fixed= n ~ 1, 
                   random= ~ Database + Host_Original , 
                   prior=Prior2,
                   data=raw_pseudo,
                   family = "poisson",
                   #pr=TRUE, 
                   nitt = 13000*mf,
                   thin = 10*mf,burnin=3000*mf)


summary(cleanModPseudo)
summary(rawModPseudo)

plot(cleanModPseudo)
plot(rawModPseudo)

cleanPseudo <- MCMCRep(cleanModPseudo)  %>%  #8.98 database 
  mutate(data="clean") 

rawPseudo <- MCMCRep(rawModPseudo)  %>%   # 8.77 database 
  mutate(data="raw") %>% 
  mutate(Component = str_replace(Component, "Host_Original", "Host"))

propVarPseudo <- bind_rows(cleanPseudo, rawPseudo)

saveRDS(cleanModPseudo, "Outputs/cleanDataPseduoMCMC.rds")
saveRDS(rawModPseudo, "Outputs/rawDataPseudoMCMC.rds")

## multivariate models (bit overkill for the df structure)

MultiPriorPar<-list(R=list(V=diag(4), nu=4.002), 
                    G=list(G1=list(V=diag(4), nu=4, 
                                   alpha.mu=rep(0,4), alpha.V=diag(4)*100)))

clean_wide <- clean %>% 
  pivot_wider(names_from = Database, values_from = n) %>% 
  as.data.frame() 

raw_wide <- raw %>% 
  pivot_wider(names_from = Database, values_from = n) %>% 
  as.data.frame() 


rawModMulti <- MCMCglmm(cbind(EID2, GMPD2, HP3, Shaw) ~ trait -1, 
                         random= ~ us(trait):Host_Original , 
                         rcov = ~us(trait):units,  
                         prior=MultiPriorPar,
                         data=raw_wide,
                         family = c("poisson","poisson", "poisson", "poisson"),
                         #pr=TRUE, 
                         nitt = 13000*mf,
                         thin = 10*mf,burnin=3000*mf)

cleanModMulti <- MCMCglmm(cbind(EID2, GMPD2, HP3, Shaw) ~ trait -1, 
                        random= ~ us(trait):Host , 
                        rcov = ~us(trait):units,  
                        prior=MultiPriorPar,
                        data=clean_wide,
                        family = c("poisson","poisson", "poisson", "poisson"),
                        #pr=TRUE, 
                        nitt = 13000*mf,
                        thin = 10*mf,burnin=3000*mf)

summary(rawModMulti)
summary(cleanModMulti)

saveRDS(cleanModMulti, "Outputs/cleanDataMultiMCMC.rds")
saveRDS(rawModMulti, "Outputs/rawDataMultiMCMC.rds")

# prop var plots  ---------------------------------------------------------

propVar %>% 
  ggplot(aes(x=data, y=as.numeric(as.character(Mode)), fill=Component)) + 
  geom_bar(stat="identity", colour="transparent") + 
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 14) + 
  labs(x='Dataset', y='Proportion Variance') -> PropVarPlot 

propVarPseudo %>% 
  ggplot(aes(x=data, y=as.numeric(as.character(Mode)), fill=Component)) + 
  geom_bar(stat="identity", colour="transparent") + 
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 14) + 
  labs(x='Dataset', y='Proportion Variance') -> PropVarPlotPseudo  


# proportion not adding up to 1 is prop var explained by poisson distribution 
# presume there's just too much overdispersion / zeroinfl in pseudoabsence altered dataframes 


