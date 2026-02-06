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
library(parallel)
library(R.utils)

# Set library path for custom package installation (if using setup_r_environment.R)
if (Sys.getenv("R_LIBS_USER") != "") {
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
}

source("full_sim/full_model.R")
source("src/new_generate_data.R")

# Define data generation parameters
root_n <- 15 
rho <- 0.7 
kappa <- 1.5 
tau <- 50 
J <- 6  
eta <- 0.3

R <- 500  # Total number of simulation runs

# Get SLURM array task ID and total tasks
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
n_tasks <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT", "1"))

# Calculate which simulations this task should run
runs_per_task <- ceiling(R / n_tasks)
start_idx <- (task_id - 1) * runs_per_task + 1
end_idx <- min(task_id * runs_per_task, R)
task_runs <- end_idx - start_idx + 1

cat(paste("Task", task_id, "of", n_tasks, ": Running simulations", start_idx, "to", end_idx, "\n"))

# MCMC settings
niter   <- 10000
nburnin <- 2000
thin    <- 5  # Thin more aggressively to reduce memory (stores 1600 samples instead of 8000)

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

compile_model <- function(model, monitor_S = TRUE){
  conf  <- configureMCMC(model)
  conf$removeSampler('rho')
  conf$addSampler(target = 'rho', type = 'slice')
  conf$removeSampler('log_kappa')
  conf$addSampler(target = 'log_kappa', type = 'slice')
  if (monitor_S) {
    conf$addMonitors(c("S", "P"))
  } else {
    conf$addMonitors(c("P"))
  }
  
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

get_results_bench <- function(samples, data){
  # For bench model, S is not monitored, so adapt the Gelman diagnostic
  true_Ps <- data$P
  modeled_Ps <- cbind(samples[[1]][,paste0("P[", 1:root_n^2, "]")])
  bias <- mean(apply(true_Ps - modeled_Ps, 2, mean))
  CIs <- apply(X=modeled_Ps, MARGIN=2 ,
               FUN= function(x) quantile(x, probs = c(0.025, 0.975)))
  coverage <- mean(ifelse(true_Ps >= CIs[1,] & true_Ps <= CIs[2,], 1, 0))
  
  gelman <- gelman.diag(as.mcmc.list(samples)[, c("rho", "log_kappa", "P[1]")],
              autoburnin = FALSE)
  
  return (c(bias, coverage, c(gelman$psrf[,1])))
}

# ------------------------------------------------------------------------------
# Compile Models ONCE (before parallel execution)
# ------------------------------------------------------------------------------

cat("Compiling NIMBLE models (this will take a few minutes)...\n")
data_dummy <- new_generate_data(root_n, rho, kappa, tau, J)

# Compile no_bench model
cat("  Compiling no_bench model...\n")
no_bench_model <- DAS_model(data_dummy, bench = "none", eta = eta)
objs_nb <- compile_model(no_bench_model, monitor_S = TRUE)

# Compile bench model
cat("  Compiling bench model...\n")
bench_model <- DAS_model(data_dummy, bench = "inexact", eta = eta)
objs_b <- compile_model(bench_model, monitor_S = FALSE)

cat("Model compilation complete!\n\n")

# ------------------------------------------------------------------------------
# Parallel Worker Function
# ------------------------------------------------------------------------------

run_chunk <- function(chunk_indices, no_bench_cmcmc, no_bench_cmodel, 
                      bench_cmcmc, bench_cmodel) {
  # Load libraries inside worker to ensure environment is correct
  library(nimble)
  library(coda)
  library(R.utils)
  
  # Storage for this chunk
  chunk_results_nobench <- matrix(0, nrow = length(chunk_indices), ncol = 7)
  chunk_results_bench <- matrix(0, nrow = length(chunk_indices), ncol = 6)
  
  # 2. Loop over assigned iterations
  for (i in seq_along(chunk_indices)) {
    r <- chunk_indices[i]
    
    # Generate new data ONCE for both models
    data_r <- new_generate_data(root_n, rho, kappa, tau, J)
    U_sd_r <- data_r$U * eta / 3
    
    # Create indicator matrix for bench model
    region_id <- data_r$region_id
    J_r <- length(unique(region_id))
    n <- root_n^2
    indicator <- matrix(0, J_r, n)
    for (j in 1:J_r) indicator[j, ] <- as.numeric(region_id == j)
    
    # ========================================================================
    # Run NO_BENCH model
    # ========================================================================
    no_bench_cmodel$setData(Pstar_obs = data_r$P_star,
                            U_obs     = data_r$U,
                            U_sd      = U_sd_r)
    
    no_bench_cmodel$setInits(list(
      P         = pmax(data_r$P_star, 1),
      S         = log(pmax(data_r$P_star, 1)),
      rho       = 0.5,
      log_kappa = 0
    ))
    
    no_bench_samples <- NULL
    no_bench_time <- NA
    
    try_no_bench <- try({
      withTimeout(
        {
          t_start <- proc.time()[3]
          no_bench_samples <- runMCMC(no_bench_cmcmc, niter = niter,
                                      nburnin = nburnin, 
                                      thin = thin, 
                                      nchains = 3, 
                                      samplesAsCodaMCMC = TRUE)
          no_bench_time <- proc.time()[3] - t_start
        },
        timeout = 180,
        onTimeout = "error"
      )
    }, silent = TRUE)
    
    if (inherits(try_no_bench, "try-error")) {
      chunk_results_nobench[i, ] <- NA
    } else {
      chunk_results_nobench[i, ] <- c(get_results(no_bench_samples, data_r), no_bench_time)
    }
    
    # ========================================================================
    # Run BENCH model
    # ========================================================================
    bench_cmodel$setData(Pstar_obs = data_r$P_star,
                         U_obs     = data_r$U,
                         U_sd      = U_sd_r,
                         ind_mat   = indicator)
    
    bench_cmodel$setInits(list(
      P         = pmax(data_r$P_star, 1),
      S         = log(pmax(data_r$P_star, 1)),
      rho       = 0.5,
      log_kappa = 0
    ))
    
    bench_samples <- NULL
    bench_time <- NA
    
    try_bench <- try({
      withTimeout(
        {
          t_start <- proc.time()[3]
          bench_samples <- runMCMC(bench_cmcmc, niter = niter,
                                   nburnin = nburnin, 
                                   thin = thin, 
                                   nchains = 3, 
                                   samplesAsCodaMCMC = TRUE)
          bench_time <- proc.time()[3] - t_start
        },
        timeout = 180,
        onTimeout = "error"
      )
    }, silent = TRUE)
    
    if (inherits(try_bench, "try-error")) {
      chunk_results_bench[i, ] <- NA
    } else {
      chunk_results_bench[i, ] <- c(get_results_bench(bench_samples, data_r), bench_time)
    }
    
    # Force garbage collection to free memory immediately
    rm(no_bench_samples, bench_samples, data_r)
    gc(verbose = FALSE)
  }
  
  return(list(nobench = chunk_results_nobench, bench = chunk_results_bench))
}

# ------------------------------------------------------------------------------
# Main Execution
# ------------------------------------------------------------------------------

# Use SLURM allocation instead of detectCores() to avoid over-parallelization
n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", detectCores() - 1))
cat(paste("Running on", n_cores, "cores within this task...\n"))

# Split the task's assigned indices into chunks for each core
task_indices <- start_idx:end_idx
chunks <- split(task_indices, cut(seq_along(task_indices), n_cores, labels = FALSE))

# Run parallel simulation for this task's subset
cat("Starting parallel simulation for task", task_id, "with", task_runs, "runs...\n")
results_list <- mclapply(chunks, run_chunk, 
                         no_bench_cmcmc = objs_nb$Cmcmc,
                         no_bench_cmodel = objs_nb$Cmodel,
                         bench_cmcmc = objs_b$Cmcmc,
                         bench_cmodel = objs_b$Cmodel,
                         mc.cores = n_cores)

# Combine results from this task
no_bench_results <- do.call(rbind, lapply(results_list, function(x) x$nobench))
bench_results <- do.call(rbind, lapply(results_list, function(x) x$bench))

# Add column names
colnames(no_bench_results) <- c("Bias", "CI_Coverage", "Gelman_rho", "Gelman_log_kappa", 
                                "Gelman_S1", "Gelman_P1", "Runtime")
colnames(bench_results) <- c("Bias", "CI_Coverage", "Gelman_rho", "Gelman_log_kappa", 
                             "Gelman_P1", "Runtime")

# Add run indices as first column
no_bench_results <- cbind(Run = start_idx:end_idx, no_bench_results)
bench_results <- cbind(Run = start_idx:end_idx, bench_results)

# Save results to task-specific CSV files
cat("Saving results to task-specific CSV files...\n")
write.csv(no_bench_results, 
          paste0("full_sim/results/no_bench_results_task_", task_id, ".csv"), 
          row.names = FALSE)
write.csv(bench_results, 
          paste0("full_sim/results/bench_results_task_", task_id, ".csv"), 
          row.names = FALSE)

# Print summary for this task
cat("\n=== TASK", task_id, "NO BENCH MODEL SUMMARY ===\n")
summary_nobench <- data.frame(
  Task = task_id,
  Avg_Bias = mean(no_bench_results[, "Bias"], na.rm = TRUE),
  Avg_Coverage = mean(no_bench_results[, "CI_Coverage"], na.rm = TRUE),
  Avg_Runtime = mean(no_bench_results[, "Runtime"], na.rm = TRUE),
  N_Failures = sum(is.na(no_bench_results[, "Runtime"]))
)
print(summary_nobench)

cat("\n=== TASK", task_id, "BENCH MODEL SUMMARY ===\n")
summary_bench <- data.frame(
  Task = task_id,
  Avg_Bias = mean(bench_results[, "Bias"], na.rm = TRUE),
  Avg_Coverage = mean(bench_results[, "CI_Coverage"], na.rm = TRUE),
  Avg_Runtime = mean(bench_results[, "Runtime"], na.rm = TRUE),
  N_Failures = sum(is.na(bench_results[, "Runtime"]))
)
print(summary_bench)

cat("\n=== GELMAN DIAGNOSTICS > 1.01 ===\n")
cat("No Bench Model:\n")
print(apply(no_bench_results[, 4:7] > 1.01, 2, sum, na.rm = TRUE))
cat("\nBench Model:\n")
print(apply(bench_results[, 4:6] > 1.01, 2, sum, na.rm = TRUE))

cat("\nTask", task_id, "results saved to:\n")
cat("  - full_sim/results/no_bench_results_task_", task_id, ".csv\n", sep = "")
cat("  - full_sim/results/bench_results_task_", task_id, ".csv\n", sep = "")
cat("\nTask", task_id, "completed successfully!\n")
