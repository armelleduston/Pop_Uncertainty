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


# Define the ddlaplace_nim distribution for NIMBLE
ddlaplace_nim <- nimbleFunction(
  run = function(x = double(0), 
                 mu = double(0), 
                 tau = double(0), 
                 log = integer(0, default = 0)) {
    returnType(double(0))
    
    # Parameter checks
    if(tau <= 0) return(NaN)
    
    # Calculate lambda parameter (lambda = exp(-1/tau))
    lambda <- exp(-1/tau)
    
    # Calculate normalization constant
    c <- (1 - lambda) / (1 + lambda)
    
    # Calculate PMF
    logProb <- log(c) + abs(x - mu) * log(lambda)
    
    # Return log or natural probability
    if(log) return(logProb)
    else return(exp(logProb))
  }
)

# Define the simulation function without using rgeom
rdlaplace_nim <- nimbleFunction(
  run = function(n = integer(0), mu = double(0), tau = double(0)) {
    returnType(double(0))
    
    # Parameter check
    if(tau <= 0) return(NaN)
    if(n != 1) stop("rdlaplace_nim only handles n = 1")
    
    # Calculate lambda
    lambda <- exp(-1/tau)
    
    # Probability of X <= mu
    p_leq_mu <- 1 / (1 + lambda)
    
    # Decide direction using uniform random number
    u <- runif(1)
    if(u < p_leq_mu) {
      # X <= mu (including equality)
      # Instead of using rgeom, implement geometric sampling directly
      v <- runif(1)
      # Convert uniform to geometric: k = floor(log(v) / log(lambda))
      k <- floor(log(v) / log(1 - lambda))
      return(mu - k)
    } else {
      # X > mu
      v <- runif(1)
      # Convert uniform to geometric: k = floor(log(v) / log(lambda))
      k <- floor(log(v) / log(1 - lambda))
      return(mu + k + 1)
    }
  }
)

# Register the distribution with NIMBLE
registerDistributions(list(
  ddlaplace_nim = list(
    BUGSdist = "ddlaplace_nim(mu, tau)",
    discrete = TRUE,
    types = c('value = double(0)', 'mu = double(0)', 'tau = double(0)'),
    pqAvail = FALSE,
    range = c(-Inf, Inf)
  )
))


# Prepare inputs for NIMBLE
n         <- length(sim_data$S)
neighbors <- lapply(1:n, function(i) which(sim_data$W[i, ] == 1))
num  <- sapply(neighbors, length)
adj <- unlist(neighbors)
mu <- rep(0, n)
L <- length(adj)
M <- CAR_calcM(num)
C <- CAR_calcC(adj, num)
region_id <- sim_data$region_id
J         <- length(unique(region_id))
indicator <- matrix(0, J, n)
for (j in 1:J) indicator[j, ] <- as.numeric(region_id == j)


# Bundle constants
constants2 <- list(
  n = n, # num of locations
  adj = adj, # adjacency matrix in vector form
  num = num, # number of neighbors
  tau = sim_data$tau, # privacy budget param
  mu = mu, # explicitely give mu = 0
  L = L, # length of adj
  J = J, # length of region_id
  M = M, # pre-computed for dcar_proper
  C = C, # pre-computed for dcar_proper
  ind_mat = indicator # to sum counts for benchmarking
)

data2 <- list(
  Pstar_obs = sim_data$P_star,   # length n
  U_obs     = sim_data$U,        # length J
  ones_p    = rep(1, n),         
  ones_u    = rep(1, J)          
)


# Initial values
inits2 <- list(
  S           = rnorm(n, 0, 1),
  P           = pmax(sim_data$P_star, 1), 
  Pstar       = sim_data$P_star,
  rho         = 0.5,
  log_kappa   = 0
)

# Define the model 
code2 <- nimbleCode({

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
    # Directly sample from discrete Laplace 
    Pstar[i] ~ ddlaplace_nim(P[i], tau)
    
    # Enforce non-negativity 
    Pstar_constrained[i] <- max(0, Pstar[i])
    
    # Enforce that the edited value equals the published count
    ones_p[i] ~ dconstraint(Pstar_constrained[i] == Pstar_obs[i])
  }

  #### Exact benchmarking to higher-level totals ---------------------------
  for (j in 1:J) {
    # Use constrained values for benchmarking
    U_sum[j] <- inprod(Pstar_constrained[1:n], ind_mat[j, 1:n])
    
    ones_u[j] ~ dconstraint(U_sum[j] == U_obs[j])
  }
})

# Build and compile the model
model2 <- nimbleModel(code2, data = data2, inits = inits2, constants = constants2)
conf2  <- configureMCMC(model2)

# Customize samplers
# Use slice sampling for rho and kappa
conf2$removeSampler('rho')
conf2$addSampler(target = 'rho', type = 'slice')
conf2$removeSampler('log_kappa')
conf2$addSampler(target = 'log_kappa', type = 'slice')
conf2$addMonitors(c("S", "P"))


# Build and compile MCMC
Rmcmc2 <- buildMCMC(conf2)
Cmodel2 <- compileNimble(model2)
Cmcmc2  <- compileNimble(Rmcmc2, project = Cmodel2)

# Run MCMC
niter   <- 30000
nburnin <- 7500
thin    <- 1
samples2 <- runMCMC(Cmcmc2, 
                   niter = niter, 
                   nburnin = nburnin, 
                   thin = thin, 
                   nchains = 3, 
                   samplesAsCodaMCMC = TRUE)

# Diagnostics
traceplot(as.mcmc.list(samples2)[, c("rho", "log_kappa", 
                                     "S[1]", "S[5]", "S[50]", 
                                     "P[1]", "P[5]", "P[50]")])

gelman.diag(as.mcmc.list(samples2)[, c("rho", "log_kappa", 
                                     "S[1]", "S[5]", "S[50]", 
                                     "P[1]", "P[5]", "P[50]")], 
            autoburnin = FALSE)
effectiveSize(as.mcmc.list(samples2)[, c("rho", "log_kappa", 
                                     "S[1]", "S[5]", "S[50]", 
                                     "P[1]", "P[5]", "P[50]")])

# Save samples for downstream analysis
save(samples2, file = 'model_fit_samples2.RData')
cat('Model fitting complete. Samples saved to model_fit_samples2.RData\n')



