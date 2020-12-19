
# Grigure ###

library(ggregplot); library(tidyverse); library(magrittr); library(igraph); library(colorspace); library(cowplot)
library(patchwork)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

list.files(full.names = T, pattern = "Clover") %>% read.csv -> Clover

Clover %>% dplyr::select(Host, Virus, Database) %>% 
  mutate_all(
    ~.x %>% str_replace_all(" ", "_")
  ) %>% unique -> ReducedClover

# Percent in each group ####

# Hosts ####

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

Plot1 <- AM %>% reshape2::melt() %>% 
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

ReducedClover %>% dplyr::select(Virus, Database) %>% unique %>% 
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

Plot2 <- AM %>% reshape2::melt() %>% 
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

# Full associations ####

ReducedClover %>% dplyr::select(Host, Virus, Database) %>% unique %>% 
  unite("Assoc", Host, Virus, sep = "_") %>% 
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

Plot3 <- AM %>% reshape2::melt() %>% 
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


# Uncleaned ####

Clover %>% dplyr::select(Host = Host_Original, Virus = Virus_Original, Database) %>% 
  mutate_all(
    ~.x %>% str_replace_all(" ", "_")
  ) %>% unique -> ReducedClover

# Hosts ####

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

Plot1b <- AM %>% reshape2::melt() %>% 
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

ReducedClover %>% dplyr::select(Virus, Database) %>% unique %>% 
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

Plot2b <- AM %>% reshape2::melt() %>% 
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

# Full associations ####

ReducedClover %>% dplyr::select(Host, Virus, Database) %>% unique %>% 
  unite("Assoc", Host, Virus, sep = "_") %>% 
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

Plot3b <- AM %>% reshape2::melt() %>% 
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


(Plot1 + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Hosts") +
    Plot2 + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Viruses") +
    Plot3 + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Associations")) +
  plot_layout(guides = "collect") +
  ggsave("OverlapPercentFigure.jpeg", units = "mm", height = 120, width = 250)


(Plot1 + 
    labs(y = "Cleaned") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Hosts") +
    Plot2 + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Viruses") +
    Plot3 + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Associations"))/
  (Plot1b + 
     labs(y = "Uncleaned") +
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Hosts") +
     Plot2b + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Viruses") +
     Plot3b + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Associations")) +
  plot_layout(guides = "collect") +
  ggsave("OverlapPercentBothRows.jpeg", units = "mm", height = 200, width = 250)


