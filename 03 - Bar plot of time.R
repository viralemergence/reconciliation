
library(ggthemr)
library(tidyverse)

setwd("~/Github/clover")

clo <- read_csv("./output/Clover_v1.0_NBCIreconciled_20201211.csv")

ggthemr("fresh")

clo %>% select(Host, Virus, Year, Database) %>%
  unique() %>%
  ggplot(aes(x = Year, fill = Database)) + 
  geom_histogram(binwidth = 1,
                 position = position_nudge(x = -0.5),
                 fill = 'grey25') + 
  theme_bw() + 
  scale_x_continuous(breaks = round(seq(1920, 2030, 10))) + 
  ylab("Number of associations") +
  theme(plot.margin = unit(c(0.9,0.9,0.9,0.9),"cm"),
        axis.title.x = element_text(vjust = -3),
        axis.title.y = element_text(vjust = 7),
        panel.border = element_rect(colour = 'black', fill = NA))
