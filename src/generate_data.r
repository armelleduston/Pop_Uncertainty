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


generate_data <- function(n, rho, kappa, tau, J){
  
  # 1. Define a random web of connections, enforce connectivity
  repeat {
    g <- sample_gnp(n, p = 4/(n-1), directed = FALSE, loops = FALSE)
    if (is.connected(g)) break
  }
  W <- as.matrix(as_adjacency_matrix(g))
  D <- diag(rowSums(W))
  
  # 2. Simulate spatial field S with (proper) CAR structure
  
  # Construct precision matrix
  Q <- kappa * (D - rho * W)
  # Add a small jitter for numerical stability
  Q_pd <- Q + diag(1e-6, n)
  
  #  Sample S ~ N(0, Sigma = solve(Q_pd))
  Sigma <- solve(Q_pd)
  S <- mvrnorm(1, mu = rep(0, n), Sigma = Sigma)
  
  # 3. Simulate true counts P_i ~ Poisson(exp(S_i))
  lambda <- exp(S)
  P <- rpois(n, lambda)
  
  # 4. Add measurement noise: P*_i | P_i ~ Laplace(P_i, tau)
  P_star <- round(P + rlaplace(n, mu = 0, sigma=tau))
  P_star <- ifelse(P_star < 0, 0, P_star) 
  
  
  # 5. Define coarse regions and benchmark totals U_j
  region_id <- sample(1:J, n, replace = TRUE)
  U <- tapply(P_star, region_id, sum)
  
  # 6. Output simulated data
  sim_data <- list(
    W = W,
    D = D,
    S = S,
    P = P,
    P_star = P_star,
    region_id = region_id,
    U = as.numeric(U),
    tau = tau,
    rho = rho
  )
  
  return (sim_data)
}





