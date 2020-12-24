
library(ggthemr)
library(tidyverse)

setwd("~/Github/clover")

clo <- read_csv("./output/Clover_v1.0_NBCIreconciled_20201211.csv")

ggthemr("fresh")

Clover %>% 
  select(Host, Virus, Year, Database) %>%
  mutate_at("Year", as.numeric) %>% 
  unique() %>%
  ggplot(aes(x = Year, fill = Database)) + 
  geom_histogram(binwidth = 1,
                 position = position_nudge(x = -0.5),
                 colour = 'grey25') +
  # geom_bar(position = "stack", colour = "grey25") +
  theme_bw() + 
  scale_x_continuous(breaks = round(seq(1920, 2030, 10))) + 
  ylab("Number of associations") +
  theme(plot.margin = unit(c(0.9,0.9,0.9,0.9),"cm"),
        axis.title.x = element_text(vjust = -3),
        axis.title.y = element_text(vjust = 7),
        panel.border = element_rect(colour = 'black', fill = NA)) +
  scale_fill_discrete_sequential(AlberPalettes[[3]], nmax = 8, order = 2:5+2)


