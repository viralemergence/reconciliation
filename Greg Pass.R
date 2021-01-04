
# Greg Pass ####

library(tidyverse)

RClover <- Clover %>% dplyr::select(Host, Virus, Database) %>% unique

RClover %>% summarise(
  
  NHost = nunique(Host),
  NVirus = nunique(Virus),
  
) %>% 
  mutate(ShawImp1 = NHost/955 - 1,
         ShawImp2 = NVirus/733 - 1)

RClover %>% select(Host, Virus) %>% unique %>% nrow

Clover$DetectionMethod %>% table

Clover %>% 
  filter(DetectionMethod == "Antibodies") %>% 
  select(Host, Virus) %>% unique %>% nrow

RClover %>% select(Host, Virus) %>% unique %>% nrow

