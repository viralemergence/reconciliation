
library(tidyverse); library(ggregplot); library(MCMCglmm)

setwd("~/Github/clover")

#clo <- read_csv("./output/Clover_v1.0_NBCIreconciled_20201211.csv")
clo<- read.csv("./CLOVER_0.1_MammalViruses_AssociationsFlatFile.csv") #updated file 

clo %>% 
  select(Database, HostOriginal, VirusOriginal) %>%
  unique %>%
  count(Database, HostOriginal) %>%
  arrange(HostOriginal) -> raw

clo %>% 
  select(Database, Host, Virus) %>%
  unique %>%
  count(Database, Host) -> clean


# a simple model  ------------------------------------------------------------


resp<-"n"
fixed <- "Database"
mcmc <- as.formula(paste(resp, "~", paste(fixed)))


Prior2 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 0.02, alpha.mu = rep(0,1), alpha.V = diag(1)*1000),
                        G2 = list(V = diag(1), nu = 0.02, alpha.mu = rep(0,1), alpha.V = diag(1)*1000)))

Prior3 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G2 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G3 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100)))


mf<-200


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
                   random= ~ Database + HostOriginal , 
                   prior=Prior2,
                   data=raw,
                   family = "poisson",
                   #pr=TRUE, 
                   nitt = 13000*mf,
                   thin = 10*mf,burnin=3000*mf)

beepr::beep()


summary(cleanMod)
summary(rawMod)

plot(cleanMod)
plot(rawMod)

cleanVar <- MCMCRep(cleanMod)  %>%  # database 4.7 % 
  mutate(data="clean") 

rawVar <- MCMCRep(rawMod)  %>%   # database 8.8% 
  mutate(data="raw") %>% 
  mutate(Component = str_replace(Component, "Host_Original", "Host"))

propVar <- bind_rows(cleanVar, rawVar)

saveRDS(cleanMod, "Outputs/cleanDataMCMC.rds")
saveRDS(rawMod, "Outputs/rawDataMCMC.rds")



# prop var plots  ---------------------------------------------------------

propVar %>% 
  ggplot(aes(x=data, y=as.numeric(as.character(Mode)), fill=Component)) + 
  geom_bar(stat="identity", colour="transparent") + 
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 14) + 
  labs(x='Dataset', y='Proportion Variance') -> PropVarPlot 



