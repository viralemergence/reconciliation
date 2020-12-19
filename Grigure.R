
# Grigure ###

library(ggregplot); library(tidyverse); library(magrittr); library(igraph); library(colorspace); library(cowplot)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

ProportionalMatrix <- function(Matrix, Observations){
  
  N <- Matrix %>% nrow
  
  A <- rep(Observations, each = N)
  B <- rep(Observations, N)
  
  AMatrix <- matrix(A, ncol = N)
  BMatrix <- matrix(B, ncol = N)
  
  Matrix/(AMatrix + BMatrix - Matrix)
  
}

AsymmetricalMatrix <- function(Matrix, Observations){
  
  N <- Matrix %>% nrow
  
  A <- rep(Observations, each = N)
  
  AMatrix <- matrix(A, ncol = N)
  
  Matrix/(AMatrix)
  
}

list.files(full.names = T, pattern = "Clover") %>% read.csv -> Clover

# Hosts ####

Clover %>% dplyr::select(Host, Virus, Database) %>% 
  mutate_all(
    ~.x %>% str_replace_all(" ", "_")
  ) %>% unique -> ReducedClover

ReducedClover %>% dplyr::select(Host, Database) %>% unique %>% 
  table() -> M1

SocGraph <- graph_from_incidence_matrix(M1)

Proj <- bipartite.projection(SocGraph)$proj2

AM <- Proj %>% get.adjacency(attr = "weight") %>% as.matrix

Observations = colSums(M1)

Matrix <- AM

N <- Matrix %>% nrow

A <- rep(Observations, each = N)

AMatrix <- matrix(A, ncol = N)

AM <- (Matrix/(AMatrix))

AM %>% reshape2::melt() %>% 
  rename_all(~str_replace(.x, "Var", "DB")) %>% 
  mutate_at("value", 
            ~(.x*100) %>% round) %>% 
  mutate(Percent = paste0(value, "%")) %>% 
  filter(!DB1 == DB2) %>% 
  ggplot(aes(DB1, DB2)) + 
  geom_tile(aes(fill = value)) + 
  scale_y_discrete(limits = c(rev(colnames(AM)))) +
  coord_fixed() + 
  labs(x = NULL, y = NULL, fill = "% Shared") +
  geom_text(aes(label = Percent, colour = as.factor(value>60))) +
  scale_colour_manual(values = c("black", "white"), guide = F) +
  scale_fill_continuous_sequential(palette = AlberPalettes[[1]], limits = c(0, 100))

# Viruses ####