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

spatial_only <- function(sim_data){
  # Prepare inputs for NIMBLE
  n         <- length(sim_data$S)
  neighbors <- lapply(1:n, function(i) which(sim_data$W[i, ] == 1))
  num  <- sapply(neighbors, length)
  adj <- unlist(neighbors)
  mu <- rep(0, n)
  L <- length(adj)
  M <- CAR_calcM(num)
  C <- CAR_calcC(adj, num)
  
  # Bundle constants
  constants1 <- list(
    n = n,
    adj = adj,
    num = num,
    mu = mu,
    L = L,
    M = M,
    C = C
  )
  
  
  data1 <- list(
    Pstar = sim_data$P_star
  )
  
  # Initial values
  inits1 <- list(
    S = rnorm(n, 0, 1),
    rho = 0.5,
    log_kappa = 0
  )
  
  # Define the model 
  code1 <- nimbleCode({
    # Hyperpriors
    rho ~ dunif(0, 1)
    log_kappa ~ dLogInvgamma(0.001, 0.001)
    kappa <- exp(log_kappa)
    
    # Proper CAR prior on S
    S[1:n] ~ dcar_proper(mu = mu[1:n],
                         adj = adj[1:L], 
                         num = num[1:n],
                         tau = kappa, 
                         gamma = rho,
                         M = M[1:n],
                         C = C[1:L])
    
    # Process model: Poisson counts
    for(i in 1:n) {
      Pstar[i] ~ dpois(exp(S[i]))
    }
    
  })
  
  
  # Build and compile the model
  model1 <- nimbleModel(code1, data = data1, inits = inits1, constants = constants1)
  conf1  <- configureMCMC(model1)
  
  # Customize samplers
  # Use slice sampling for rho and kappa
  conf1$removeSampler('rho')
  conf1$addSampler(target = 'rho', type = 'slice')
  conf1$removeSampler('kappa')
  conf1$addSampler(target = 'kappa', type = 'slice', adaptInterval = 1000)
  conf1$addMonitors("S")
  
  # Build and compile MCMC
  Rmcmc1 <- buildMCMC(conf1)
  Cmodel1 <- compileNimble(model1)
  Cmcmc1  <- compileNimble(Rmcmc1, project = Cmodel1)
  
  # Run MCMC
  niter   <- 30000
  nburnin <- 7500
  thin <- 1
  samples1 <- runMCMC(Cmcmc1, 
                      niter = niter, 
                      nburnin = nburnin, 
                      thin = thin, 
                      nchains = 3, 
                      samplesAsCodaMCMC = TRUE)
  
  return(samples1)
}

