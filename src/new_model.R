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


new_model <- function(sim_data, bench = "none"){
  # Prepare inputs for NIMBLE
  n         <- length(sim_data$S)
  neighbors <- lapply(1:n, function(i) which(sim_data$W[i, ] == 1))
  num  <- sapply(neighbors, length)
  adj <- unlist(neighbors)
  weights <- rep(1, length(adj))
  CM <- as.carCM(adj, weights, num)
  C <- CM$C
  M <- CM$M
  mu <- rep(0, n)
  eta <- 0.1
  L <- length(adj)
  region_id <- sim_data$region_id
  J         <- length(unique(region_id))
  indicator <- matrix(0, J, n)
  for (j in 1:J) indicator[j, ] <- as.numeric(region_id == j)
  rho_max <- carMaxBound(C, adj, num, M)
  
  
  # Bundle constants
  constants <- list(
    n = n, # num of locations
    adj = adj, # adjacency matrix in vector form
    num = num, # number of neighbors
    tau = sim_data$tau, # privacy budget param
    mu = mu, # explicitly give mu = 0
    eta = eta, # discrepancy parameter
    L = L, # length of adj
    J = J, # length of region_id
    M = M, # pre-computed for dcar_proper
    C = C, # pre-computed for dcar_proper
    ind_mat = indicator, # to sum counts for benchmarking0
    rho_max = rho_max
  )
  
  data <- list(
    Pstar_obs = sim_data$P_star,   # length n
    U_obs     = sim_data$U,        # length J
    ones_p    = rep(1, n),         
    ones_u    = rep(1, J)          
  )
  
  
  # Initial values
  inits <- list(
    S           = rnorm(n, 0, 1),
    P           = pmax(sim_data$P_star, 1), 
    Pstar       = sim_data$P_star,
    rho         = 0.5,
    log_kappa   = 0
  )
  
  # Define the model 
  code <- nimbleCode({
    
    #### Hyperpriors ----------------------------------------------------------
    rho   ~ dunif(0, 1)
    log_kappa ~ dLogInvgamma(.001, .001)
    kappa <- exp(log_kappa)
    
    #### Spatial random effects (proper CAR) ----------------------------------
    S[1:n] ~ dcar_proper(
      mu    = mu[1:n],
      adj   = adj[1:L],
      num   = num[1:n],
      tau   = kappa,
      gamma = rho,
      M     = M[1:n],
      C     = C[1:L])
    
    #### Process model: true counts ------------------------------------------
    for (i in 1:n) {
      P[i] ~ dpois(exp(S[i]))
    }
    
    #### Measurement model: direct discrete Laplace noise --------------------
    for (i in 1:n) { 
      Pstar_obs[i] ~ droundnorm(P[i], tau)
    }
    
    #### Exact benchmarking to higher-level totals ---------------------------
    
    if (bench == "exact"){
      for (j in 1:J) {
        U_sum[j] <- inprod(P[1:n], ind_mat[j, 1:n])
        ones_u[j] ~ dconstraint(U_sum[j] == U_obs[j])
      }
    }
    if (bench == "inexact"){
      for (j in 1:J) {
        U_sum[j] <- inprod(P[1:n], ind_mat[j, 1:n])
        U_obs[j]  ~ dpois(eta * U_sum[j])
      }
    }
  })
  
  # Build and compile the model
  model <- nimbleModel(code, data = data, inits = inits, constants = constants)
  conf  <- configureMCMC(model)
  
  # Customize samplers
  # Use slice sampling for rho and kappa
  conf$removeSampler('rho')
  conf$addSampler(target = 'rho', type = 'slice')
  conf$removeSampler('log_kappa')
  conf$addSampler(target = 'log_kappa', type = 'slice')
  conf$addMonitors(c("S", "P"))
  
  
  # Build and compile MCMC
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(model)
  Cmcmc  <- compileNimble(Rmcmc, project = Cmodel)
  
  # Run MCMC
  niter   <- 30000
  nburnin <- 7500
  thin    <- 1
  samples <- runMCMC(Cmcmc, 
                      niter = niter, 
                      nburnin = nburnin, 
                      thin = thin, 
                      nchains = 3, 
                      samplesAsCodaMCMC = TRUE)
  
  return (samples)
}







