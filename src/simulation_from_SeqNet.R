library(tidyverse)
library(SeqNet)
library(igraph)
library(viridis)

# Read in each row from df_human_base
# set p, s, alpha, and alpha_link_node based on its properties.
# run ~10 simulations and store topology in df.
# compare/plot results

# Notes: For higher avg_path_len, alpha_modules should be high (~50)

output_dir <- file.path("output", "simulations_SeqNet")
if(!dir.exists(output_dir)) dir.create(output_dir)

run_sim <- function(params) {
  p <- params$p[1]
  # Minimum sparsity is 1 / (p - 1). 
  min_sparsity <- 1 / (p - 1)
  s <- ifelse(params$sparsity[1] < min_sparsity, min_sparsity, params$sparsity[1])
  
  # NOTE: the following two functions will be implemented in newer versions 
  # of the SeqNet package.
  sample_link_nodes <- function(n, nodes, degree, alpha = 100, beta = 1, epsilon = 10^-5, ...) {
    prob <- pbeta(ecdf_cpp(degree), alpha_link_node, beta) + epsilon
    n <- 1
    link_nodes <- sample(nodes, n, prob = prob)
    return(link_nodes)
  }
  sample_module_nodes <- function(n, nodes, degree, nu = 0.01, ...) {
    if(nu <= 0) {
      stop("nu must be positive.")
    }
    m <- length(nodes)
    val <- 1 / m
    prob <- rep(val, m)
    prob[degree > 0] = val * nu * pbeta(ecdf_cpp(degree[degree > 0]), 
                                        alpha_modules, beta_modules)
    module_nodes <- sample(nodes, n, prob = prob)
    return(module_nodes)
  }
  ###################################
  
  
  ###################################
  # Run simulation and save results.
  ###################################
  results <- cbind(params,
                   alpha = NA,
                   alpha_link_node = NA,
                   inflection = NA,
                   nu = NA,
                   n_modules = NA,
                   neig_size = NA,
                   prob_remove = NA)
  for(i in 1:n_sim) {
    cat("running simulation", i, "/", n_sim, "-", format(Sys.time()), "\n")
    ###################################
    # SeqNet settings
    ###################################
    alpha <- 100
    beta <- 1 # ceiling(runif(1, 1, 100))
    # Note: Lower values for inflection decrease sparsity (i.e. more edges).
    inflection <- runif(1, 0.5, 1)
    alpha_modules <- 10 # ceiling(runif(1, 1, 50))
    beta_modules <- (alpha_modules - 1) * (1 - inflection) / inflection + 1 
    prob_rewire <- 0.5
    alpha_link_node <- ceiling(runif(1, 100, 1000))
    
    # Estimated number of modules that will be created:
    # Input: p, the network size; s, the desired sparsity.
    # Based on an average module size of m = 100.
    # NOTE: all of these calculations will be done automatically in newer
    # versions of the SeqNet package.
    m <- min(100, p)
    # Given nu, calculates the average number of modules in the network.
    n_modules <- function(nu) {
      if(p <= m) return(rep(1, length(nu)))
      log(0.0001) / log(1 - m / p) * nu + (1 - nu) * p / m
    }
    # Given nu, calculates the neigborhood size needed for desired sparsity s.
    neig_size <- function(nu) {
      pmax(1, ceiling(s * p * (p - 1) / (2 * m * n_modules(nu))))
    }
    sparsity <- function(nu) {
      (2 * m * n_modules(nu) * neig_size(nu)) / (p * (p - 1))
    }
    
    # Nu will be restricted to ensure neighborhood size is attainable.
    # Note: For nu = 0, there is no overlap of modules => increase neig size.
    #       For nu = 1, there is unrestricted overlap => flexible neig size,
    #       since removal step can trim excess edges.
    nu_set <- seq(0.01, 0.63, 0.01)^2 # Consider nu in (0.0001, 0.5).
    min_nu <- nu_set[which(neig_size(nu_set) < (m - 1) / 2)[1]]
    # Excess edges will be handled by edge removal step, but we don't want
    # nu to be too high, otherwise prob_removal will need to be set too high.
    max_nu <- tryCatch(
      max(nu_set[cummax(1 - s / sparsity(nu_set)) <= 0.5]),
      warning = function(w) return(1) # If 0.5 is too low, set max nu to 1.
    )
    nu <- 10^(runif(1, log10(min_nu), log10(max_nu)))
    
    n_modules_nu <- n_modules(nu)
    neig_size_nu <- neig_size(nu)
    
    s_est <- (2 * m * n_modules_nu * neig_size_nu) / (p * (p - 1)) # Estimated sparsity.
    prob_remove <- 1 - s / s_est
    min_module_size <- 2 * neig_size_nu + 1 # Minimum module size for desired neig.
    avg_module_size <- max(m, 1.1 * min_module_size) # Min size is 2 * neig + 1.
    
    network <- random_network(p = p,
                              nu = nu,
                              alpha = alpha,
                              beta = beta,
                              prob_rewire = prob_rewire,
                              prob_remove = prob_remove,
                              neig_size = neig_size_nu,
                              avg_module_size = avg_module_size,
                              sd_module_size = 30,
                              min_module_size = min_module_size,
                              max_module_size = 200,
                              sample_link_nodes_fn = sample_link_nodes,
                              sample_module_nodes_fn = sample_module_nodes)
    
    adj_sim <- get_adjacency_matrix(network)
    graph_sim <- igraph::graph_from_adjacency_matrix(adj_sim, mode = "undirected")
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
                       alpha = alpha,
                       alpha_link_node = alpha_link_node,
                       inflection = inflection,
                       nu = nu,
                       n_modules = n_modules_nu,
                       neig_size = neig_size_nu,
                       prob_remove = prob_remove))
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

