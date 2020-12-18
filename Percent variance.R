
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
