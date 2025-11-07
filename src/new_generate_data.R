### Load Packages 
library(ggplot2)
library(tidyr)
library(patchwork)
library(dplyr)
library(knitr)
library(kableExtra)
library(MASS)
library(nimble)
library(coda)
library(extraDistr)
library(igraph)
library(RColorBrewer)
library(nimbleNoBounds)
library(pracma)
library(extraDistr)

# Build J connected regions on a grid via graph-Voronoi (multi-source BFS)
graph_voronoi_regions <- function(g, J, avoid_singletons = TRUE, max_tries = 50) {
  n <- igraph::vcount(g)
  for (t in 1:max_tries) {
    seeds <- sample.int(n, J)                # pick J random seeds
    distJ <- igraph::distances(g, v = seeds, to = igraph::V(g))
    # Tie-break randomly to keep regions “random”
    region_id <- apply(distJ, 2, function(d) {
      m <- min(d)
      cands <- which(d == m)
      cands[sample.int(length(cands), 1)]
    })
    region_id <- as.integer(region_id)       # 1..J by seed order
    if (!avoid_singletons || all(tabulate(region_id, nbins = J) > 1)) {
      return(region_id)
    }
  }
  warning("Could not avoid singleton regions after max_tries; returning last assignment.")
  as.integer(region_id)
}

# Ensure regional contiguity by reassigning disconnected pieces ("stragglers")
enforce_contiguous_regions <- function(g, region_id, max_passes = 5) {
  pick_majority <- function(x) {
    if (length(x) == 0) return(NA_integer_)
    tab <- sort(table(x), decreasing = TRUE)
    winners <- names(tab)[tab == max(tab)]
    as.integer(sample(winners, 1))
  }
  for (pass in seq_len(max_passes)) {
    changed <- FALSE
    for (j in sort(unique(region_id))) {
      idx <- which(region_id == j)
      if (length(idx) <= 1) next
      comps <- igraph::components(igraph::induced_subgraph(g, vids = idx))
      if (comps$no <= 1) next
      keep_comp <- which.max(comps$csize)
      strag_local <- which(comps$membership != keep_comp)
      if (length(strag_local) == 0) next
      strag <- idx[strag_local]
      for (v in strag) {
        nb <- as.integer(igraph::neighbors(g, v))
        nb_regions <- region_id[nb]
        nb_regions <- nb_regions[nb_regions != j]
        new_r <- pick_majority(nb_regions)
        if (!is.na(new_r) && new_r != j) {
          region_id[v] <- new_r
          changed <- TRUE
        }
      }
    }
    # stop early if all regions are connected
    all_ok <- all(sapply(sort(unique(region_id)), function(j) {
      idx <- which(region_id == j)
      if (length(idx) <= 1) return(TRUE)
      igraph::components(igraph::induced_subgraph(g, vids = idx))$no == 1
    }))
    if (all_ok || !changed) break
  }
  region_id
}

new_generate_data <- function(root_n, rho, kappa, tau, J){
  
  n <- root_n^2
  
  # 1. Define a root_n x root_n grid of connections (4-neighbor)
  g <- make_lattice(dimvector = c(root_n, root_n), nei = 1, periodic = FALSE)
  W <- as_adjacency_matrix(g, sparse = FALSE)  
  storage.mode(W) <- "double"
  stopifnot(isSymmetric(W))                    # sanity check
  D <- diag(rowSums(W), nrow = n, ncol = n)    
  
  # 2. Simulate spatial field S with (proper) CAR structure
  Q <- kappa * (D - rho * W)
  Q_pd <- Q + diag(1e-6, n)                    
  
  Sigma <- solve(Q_pd)                          
  S <- mvrnorm(1, mu = rep(0, n), Sigma = Sigma)
  
  # 3. Simulate true counts P_i ~ Poisson(exp(S_i))
  base_pop <- rep(0,n)
  rural_urban <- rbinom(n ,1, exp(S)/(1 + exp(S)))
  rural_idx = which(rural_urban == 1)
  urban_idx = which(rural_urban == 0)
  base_pop[rural_idx] <- rlnorm(n=length(rural_idx), 9.9, 0.55)
  base_pop[urban_idx] <- rlnorm(n=length(rural_idx), 12.5, 0.85)
  
  lambda <- base_pop * exp(S)
  P <- rpois(n, lambda)
  
  # 4. Measurement noise 
  P_star <- round(P + rnorm(n, mean = 0, sd = tau))
  P_star <- pmax(P_star, 0)
  
  # 5. Random but connected regions on the lattice
  region_id <- graph_voronoi_regions(g, J, avoid_singletons = TRUE)
  region_id <- enforce_contiguous_regions(g, region_id)  # absorb stragglers

  # 6. Benchmark totals by region
  U <- tapply(P_star, region_id, sum)
  
  sim_data <- list(
    W = W,                       # dense 0/1 adjacency
    D = D,                       # dense diagonal degree matrix
    S = S,
    P = P,
    P_star = P_star,
    region_id = region_id,
    U = as.numeric(U),
    tau = tau,
    rho = rho
  )
  
  return ( sim_data )
}
