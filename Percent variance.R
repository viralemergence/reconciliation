
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


mf<-30

cleanMod <- MCMCglmm(fixed= n ~ 1, 
                     random= ~ Database + Host , 
                     prior=Prior2,
                     data=clean,
                     family = "poisson",
                     pr=TRUE, 
                     nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
                     thin = 10*mf,burnin=3000*mf)

rawMod <- MCMCglmm(fixed= n ~ 1, 
                   random= ~ Database + Host_Original , 
                   prior=Prior2,
                   data=raw,
                   family = "poisson",
                   pr=TRUE, 
                   nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
                   thin = 10*mf,burnin=3000*mf)


summary(cleanMod)
summary(rawMod)

plot(cleanMod)
plot(rawMod)

clean <- MCMCRep(cleanMod)  %>%  # database 5.66%  
  mutate(data="clean") 

raw <- MCMCRep(rawMod)  %>%   # database 6.98% 
  mutate(data="raw") %>% 
  mutate(Component = str_replace(Component, "Host_Original", "Host"))

propVar <- bind_rows(clean, raw)


# prop var plots  ---------------------------------------------------------

propVar %>% 
  ggplot(aes(x=data, y=as.numeric(as.character(Mode)), fill=Component)) + 
  geom_bar(stat="identity", colour="transparent") + 
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 14) + 
  labs(x='Dataset', y='Proportion Variance') -> PropVar  


# proportion not adding up to 1 is prop var explained by poisson distribution 


