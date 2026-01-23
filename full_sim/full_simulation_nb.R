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
library(parallel) # Required for parallelization

source("full_sim/full_model.R")
source("src/new_generate_data.R")

# define to generate data
root_n <- 15 
rho <- 0.7 
kappa <- 1.5 
tau <- 50 
J <- 6  
eta <- 0.3

R <- 1000

# MCMC settings
niter   <- 10000
nburnin <- 2000
thin    <- 1

# ------------------------------------------------------------------------------
# NIMBLE Custom Distributions (Must be defined in global scope)
# ------------------------------------------------------------------------------
## helper: log(exp(a) - exp(b)) for a >= b, numerically stable
logspace_sub <- nimbleFunction(
  run = function(a = double(0), b = double(0)) {
    returnType(double(0))
    if(b > a) return(logspace_sub(b, a))   # enforce a >= b
    if(b == a) return(-Inf)
    ## a + log(1 - exp(b-a))
    return(a + log1p(-exp(b - a)))
  }
)

ddnorm_nim <- nimbleFunction(
  run = function(x = double(0),
                 mean = double(0),
                 sd = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    
    if(sd <= 0) {
      if(log) return(-Inf) else return(0.0)
    }
    
    ## enforce integer support (as in extraDistr)
    if(x != floor(x)) {
      if(log) return(-Inf) else return(0.0)
    }
    
    ## compute log PMF stably: log( Phi(b) - Phi(a) )
    ## where a = (x-mean)/sd, b = (x+1-mean)/sd
    lp_b <- pnorm(x + 1.0, mean, sd, 1, 1)  # log Phi((x+1-mean)/sd)
    lp_a <- pnorm(x,       mean, sd, 1, 1)  # log Phi((x-mean)/sd)
    
    logProb <- logspace_sub(lp_b, lp_a)
    
    if(log) return(logProb) else return(exp(logProb))
  }
)

rdnorm_nim <- nimbleFunction(
  run = function(n = integer(0), mean = double(0), sd = double(0)) {
    returnType(double(0))
    if(sd <= 0) return(NaN)
    if(n != 1) stop("rdnorm_nim only handles n = 1")
    return(floor(rnorm(1, mean, sd)))
  }
)

registerDistributions(list(
  ddnorm_nim = list(
    BUGSdist = "ddnorm_nim(mean, sd)",
    discrete = TRUE,
    types = c('value = double(0)', 'mean = double(0)', 'sd = double(0)'),
    pqAvail = FALSE,
    range = c(-Inf, Inf)
  )
))

# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------

compile_model <- function(model){
  conf  <- configureMCMC(model)
  conf$removeSampler('rho')
  conf$addSampler(target = 'rho', type = 'slice')
  conf$removeSampler('log_kappa')
  conf$addSampler(target = 'log_kappa', type = 'slice')
  conf$addMonitors(c("S", "P"))
  
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(model)
  Cmcmc  <- compileNimble(Rmcmc, project = Cmodel)
  
  return(list(Cmodel = Cmodel, Cmcmc = Cmcmc))
}

get_results <- function(samples, data){
  true_Ps <- data$P
  modeled_Ps <- cbind(samples[[1]][,paste0("P[", 1:root_n^2, "]")])
  bias <- mean(apply(true_Ps - modeled_Ps, 2, mean))
  CIs <- apply(X=modeled_Ps, MARGIN=2 ,
               FUN= function(x) quantile(x, probs = c(0.025, 0.975)))
  coverage <- mean(ifelse(true_Ps >= CIs[1,] & true_Ps <= CIs[2,], 1, 0))
  
  gelman <- gelman.diag(as.mcmc.list(samples)[, c("rho", "log_kappa", 
                                         "S[1]", "P[1]")],
              autoburnin = FALSE)
  
  return (c(bias, coverage, c(gelman$psrf[,1])))
}

# ------------------------------------------------------------------------------
# Parallel Worker Function
# ------------------------------------------------------------------------------

run_chunk <- function(chunk_indices) {
  # Load libraries inside worker to ensure environment is correct
  library(nimble)
  library(coda)
  
  # 1. Compile model ONCE per worker
  # Generate dummy data just to define the model structure
  data_dummy <- new_generate_data(root_n, rho, kappa, tau, J)
  # Note: U_sd is now in 'data' in full_model_245.R, so we don't need to worry about it here
  no_bench_model <- DAS_model(data_dummy, bench = "none", eta = eta)
  
  # Compile
  objs <- compile_model(no_bench_model)
  no_bench_cmodel <- objs$Cmodel
  no_bench_cmcmc  <- objs$Cmcmc
  
  # Storage for this chunk
  chunk_results <- matrix(0, nrow = length(chunk_indices), ncol = 7)
  
  # 2. Loop over assigned iterations
  for (i in seq_along(chunk_indices)) {
    r <- chunk_indices[i]
    
    # Generate new data
    data_r <- new_generate_data(root_n, rho, kappa, tau, J)
    U_sd_r <- data_r$U * eta / 3
    
    # Update data in compiled model
    no_bench_cmodel$setData(Pstar_obs = data_r$P_star,
                            U_obs     = data_r$U,
                            U_sd      = U_sd_r)
    
    # Reset initial values 
    no_bench_cmodel$setInits(list(
      P         = pmax(data_r$P_star, 1),
      S         = log(pmax(data_r$P_star, 1)),
      rho       = 0.5,
      log_kappa = 0
    ))
    
    no_bench_samples <- NULL
    no_bench_time <- NA
    
    # Set 2 minute time limit
    setTimeLimit(elapsed = 120, transient = TRUE)
    
    try_no_bench <- try({
      no_bench_t <- system.time({
        no_bench_samples <- runMCMC(no_bench_cmcmc, niter = niter,
                                    nburnin = nburnin, 
                                    thin = thin, 
                                    nchains = 3, 
                                    samplesAsCodaMCMC = TRUE)
      })
      no_bench_time <- no_bench_t[3]
    }, silent = TRUE)
    
    # Reset time limit
    setTimeLimit(elapsed = Inf)
    
    if (inherits(try_no_bench, "try-error")) {
      chunk_results[i, ] <- NA
    } else {
      chunk_results[i, ] <- c(get_results(no_bench_samples, data_r), no_bench_time)
    }
  }
  
  return(chunk_results)
}

# ------------------------------------------------------------------------------
# Main Execution
# ------------------------------------------------------------------------------

# Detect cores (leave one free)
n_cores <- detectCores() - 1
print(paste("Running on", n_cores, "cores..."))

# Split indices 1:R into chunks for each core
chunks <- split(1:R, cut(1:R, n_cores, labels = FALSE))

# Run parallel simulation
# mclapply forks the process, so workers inherit functions/variables
results_list <- mclapply(chunks, run_chunk, mc.cores = n_cores)

# Combine results
no_bench_results <- do.call(rbind, results_list)

data.frame(rbind(apply(no_bench_results[,c(1:2, 7)], 2, mean))) |>
  kbl(caption = "Simulation Results for 1000 runs",
      booktabs = TRUE,
      digits = 4,
      col.names = c("Avg Bias", "Avg CI Coverage", "Avg Runtime (3 chains)")) |>
  kable_styling(full_width = FALSE,
                latex_options = c("hold_position", "striped"))

apply(no_bench_results[,3:6]>1.01, 2, sum)

