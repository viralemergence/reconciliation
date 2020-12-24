
list.files(full.names = T, pattern = "VIRION.csv") %>% read.csv -> Clover

Clover %>% dplyr::select(Host, Virus, Database) %>% 
  mutate_all(
    ~.x %>% str_replace_all(" ", "_")
  ) %>% unique -> ReducedClover

ReducedClover %<>% mutate_all(~.x %>% str_trim %>% str_replace_all(" ", "_"))

# Associations ####

ReducedClover %>% dplyr::select(Host, Virus, Database) %>% unique %>% 
  unite("Assoc", Host, Virus, sep = "_") %>% 
  table() -> M1

AssocFunction <- function(a){
  
  print(a)
  
  M1[M1[,a]==1,] %>% colSums %>% divide_by(sum(M1[,a])) %>% return
  
}

M1 %>% colnames %>% map(c(AssocFunction, as_tibble)) %>% bind_cols %>% as.matrix ->
  AM

dimnames(AM) <- list(colnames(M1), colnames(M1))

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
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_manual(values = c("black", "white"), guide = F) +
  scale_fill_continuous_sequential(palette = AlberPalettes[[1]], limits = c(0, 100))
