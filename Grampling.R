
# Grampling ####

Connectance <- 
  ReducedClover %>% dplyr::select(Host, Virus) %>% unique %>% 
  table() %>% unlist %>% Prev

1/Connectance

NVirus <- nunique(ReducedClover$Virus)
NHost <- nunique(ReducedClover$Host)

NAssocs <- NVirus*NHost
NAssocs <- ReducedClover %>% dplyr::select(Host, Virus) %>% unique %>% 
  nrow

MaximumMatrix <- matrix(1:NAssocs, ncol = NVirus, nrow = NHost)

MaximumMatrix %>% dim

Links <- ReducedClover$Database %>% table %>% as.list

NIterations <- 10

ProportionalMatrix <- function(Matrix, Observations){
  
  N <- Matrix %>% nrow
  
  A <- rep(Observations, each = N)
  B <- rep(Observations, N)
  
  AMatrix <- matrix(A, ncol = N)
  BMatrix <- matrix(B, ncol = N)
  
  Matrix/(AMatrix + BMatrix - Matrix)
  
}

1:NIterations %>% map(function(a){#
  
  print(a)
  
  Links %>% names %>% map(~data.frame(Association = sample(1:NAssocs, Links[[.x]]),
                                      Dataset = .x)) %>% 
    bind_rows %>% table %>% 
    
    graph_from_incidence_matrix %>% bipartite.projection %>% 
    extract2("proj2") %>% 
    get.adjacency(attr = "weight") %>% 
    as.matrix %>% 
    ProportionalMatrix(Observations = unlist(Links)) %>% 
    return
  
}) -> PropAssociationList

AverageOverlap <- PropAssociationList %>% reduce(add) %>% round(3) %>% divide_by(NIterations)

AverageOverlap
