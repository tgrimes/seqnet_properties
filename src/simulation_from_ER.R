library(tidyverse)
library(SeqNet)
library(igraph)
library(viridis)

# Read in each row from df_human_base and df_inetbio.
# set p and m based on its properties.
# run ~20 simulations and store topology in df.
# compare/plot results.

# Note: runtime is slow due to computing topological measures.

output_dir <- file.path("output", "simulations_ER")
if(!dir.exists(output_dir)) dir.create(output_dir)


run_sim <- function(params) {
  p <- params$p[1]
  s <- params$sparsity[1] 
  r <- params$max_degree[1] / p
  
  ###################################
  # Run simulation and save results.
  ###################################
  results <- cbind(params, 
                   m = NA,
                   alpha = NA)
  for(i in 1:n_sim) {
    cat("running simulation", i, "/", n_sim, "-", format(Sys.time()), "\n")
    
    pi <- s 
    graph_sim <- igraph::erdos.renyi.game(n = p, p.or.m = pi, type = "gnp",
                                          directed = FALSE)
    deg_sim <- igraph::degree(graph_sim)
    
    results <- rbind(results, 
                     c(tissue = i,
                       p = length(V(graph_sim)),
                       sparsity = 2 * length(E(graph_sim)) / (p * (p - 1)), # sparsity
                       max_deg = max(deg_sim), # maximum degree
                       diameter = igraph::diameter(graph_sim, directed = FALSE), # longest shortest path
                       avg_degree = mean(deg_sim), # average degree
                       avg_path_len = igraph::mean_distance(graph_sim), # avg path length
                       CC = igraph::transitivity(graph_sim, type = "global"), # CC
                       CC_alt = sum(igraph::transitivity(graph_sim, vids = V(graph_sim)[deg_sim > 1], type = "local")) / p, # CC alt.
                       betw_centr = igraph::centr_betw(graph_sim, directed = "FALSE")$centralization,
                       clo_centr = igraph::centr_clo(graph_sim, mode = "all")$centralization,
                       deg_centr = igraph::centr_degree(graph_sim)$centralization,
                       eig_centr = igraph::centr_eigen(graph_sim, directed = "FALSE")$centralization,
                       pi = pi))
    ###################################
  }
  
  return(results)
}

############################################################
# Run on tissue-specfic human networks
############################################################
for(k in 1:nrow(df_human_base)) {
  tissue <- df_human_base[k, 1]
  cat("Running simulation for", tissue, "tissue. ( k =", k, ")\n")
  sim <- run_sim(df_human_base[k, ])
  saveRDS(sim, file.path(output_dir, 
                         paste0(tissue, "_", format(Sys.time(), "%m_%d_%H%M"), 
                                ".rds")))
  gc()
}

############################################################
# Run on organism-specific networks
############################################################
for(k in 1:nrow(df_inetbio)) {
  tissue <- df_inetbio[k, 1]
  cat("Running simulation for", tissue, "tissue. ( k =", k, ")\n")
  sim <- run_sim(df_inetbio[k, ])
  saveRDS(sim, file.path(output_dir, 
                         paste0(tissue, "_", format(Sys.time(), "%m_%d_%H%M"), 
                                ".rds")))
  gc()
}
