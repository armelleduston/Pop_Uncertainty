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

get_combined_results <- function(no_bench_samples, bench_samples, data, 
                                  no_bench_time, bench_time){
  true_Ps <- data$P
  true_Pstar <- data$P_star
  n <- length(true_Ps)
  
  # Extract posterior means (P_hat) for both models
  no_bench_Phat <- colMeans(no_bench_samples[[1]][, paste0("P[", 1:n, "]")])
  bench_Phat <- colMeans(bench_samples[[1]][, paste0("P[", 1:n, "]")])
  
  # Compute metrics for no_bench model
  no_bench_modeled_Ps <- cbind(no_bench_samples[[1]][, paste0("P[", 1:n, "]")])
  no_bench_bias <- mean(apply(true_Ps - no_bench_modeled_Ps, 2, mean))
  no_bench_CIs <- apply(X = no_bench_modeled_Ps, MARGIN = 2,
                        FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
  no_bench_coverage <- mean(ifelse(true_Ps >= no_bench_CIs[1,] & 
                                    true_Ps <= no_bench_CIs[2,], 1, 0))
  
  no_bench_gelman <- gelman.diag(as.mcmc.list(no_bench_samples)[, 
                                  c("rho", "log_kappa", "S[1]", "P[1]")],
                                 autoburnin = FALSE)
  
  # Compute metrics for bench model
  bench_modeled_Ps <- cbind(bench_samples[[1]][, paste0("P[", 1:n, "]")])
  bench_bias <- mean(apply(true_Ps - bench_modeled_Ps, 2, mean))
  bench_CIs <- apply(X = bench_modeled_Ps, MARGIN = 2,
                     FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
  bench_coverage <- mean(ifelse(true_Ps >= bench_CIs[1,] & 
                                 true_Ps <= bench_CIs[2,], 1, 0))
  
  bench_gelman <- gelman.diag(as.mcmc.list(bench_samples)[, 
                              c("rho", "log_kappa", "S[1]", "P[1]")],
                              autoburnin = FALSE)
  
  # Compute differences
  diff_bench_Pstar <- mean(abs(bench_Phat - true_Pstar))  # |P_hat_bench - P*|
  diff_nobench_Pstar <- mean(abs(no_bench_Phat - true_Pstar))  # |P_hat_nobench - P*|
  diff_bench_nobench <- mean(abs(bench_Phat - no_bench_Phat))  # |P_hat_bench - P_hat_nobench|
  
  # Return combined results as a vector
  return(c(
    no_bench_bias = no_bench_bias,
    no_bench_coverage = no_bench_coverage,
    no_bench_gelman_rho = no_bench_gelman$psrf["rho", 1],
    no_bench_gelman_log_kappa = no_bench_gelman$psrf["log_kappa", 1],
    no_bench_gelman_S1 = no_bench_gelman$psrf["S[1]", 1],
    no_bench_gelman_P1 = no_bench_gelman$psrf["P[1]", 1],
    no_bench_runtime = no_bench_time,
    bench_bias = bench_bias,
    bench_coverage = bench_coverage,
    bench_gelman_rho = bench_gelman$psrf["rho", 1],
    bench_gelman_log_kappa = bench_gelman$psrf["log_kappa", 1],
    bench_gelman_S1 = bench_gelman$psrf["S[1]", 1],
    bench_gelman_P1 = bench_gelman$psrf["P[1]", 1],
    bench_runtime = bench_time,
    diff_bench_Pstar = diff_bench_Pstar,
    diff_nobench_Pstar = diff_nobench_Pstar,
    diff_bench_nobench = diff_bench_nobench
  ))
}

# ------------------------------------------------------------------------------
# Compile Models ONCE (before parallel execution)
# ------------------------------------------------------------------------------

cat("Compiling NIMBLE models (this will take a few minutes)...\n")
data_dummy <- new_generate_data(root_n, rho, kappa, tau, J)

# Compile no_bench model
cat("  Compiling no_bench model...\n")
no_bench_model <- DAS_model(data_dummy, bench = "none", eta = eta)
objs_nb <- compile_model(no_bench_model)

# Compile bench model
cat("  Compiling bench model...\n")
bench_model <- DAS_model(data_dummy, bench = "inexact", eta = eta)
objs_b <- compile_model(bench_model)

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
  
  # Storage for this chunk - combined results
  chunk_results <- matrix(0, nrow = length(chunk_indices), ncol = 17)
  
  # Loop over assigned iterations
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
    
    # Compute combined results only if both models succeeded
    if (inherits(try_no_bench, "try-error") || inherits(try_bench, "try-error")) {
      chunk_results[i, ] <- NA
    } else {
      chunk_results[i, ] <- get_combined_results(no_bench_samples, bench_samples, 
                                                  data_r, no_bench_time, bench_time)
    }
    
    # Force garbage collection to free memory immediately
    rm(no_bench_samples, bench_samples, data_r)
    gc(verbose = FALSE)
  }
  
  return(chunk_results)
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
combined_results <- do.call(rbind, results_list)

# Add column names
colnames(combined_results) <- c(
  "NoBench_Bias", "NoBench_Coverage", 
  "NoBench_Gelman_rho", "NoBench_Gelman_log_kappa", 
  "NoBench_Gelman_S1", "NoBench_Gelman_P1", "NoBench_Runtime",
  "Bench_Bias", "Bench_Coverage", 
  "Bench_Gelman_rho", "Bench_Gelman_log_kappa", 
  "Bench_Gelman_S1", "Bench_Gelman_P1", "Bench_Runtime",
  "Diff_Bench_Pstar", "Diff_NoBench_Pstar", "Diff_Bench_NoBench"
)

# Add run indices as first column
combined_results <- cbind(Run = start_idx:end_idx, combined_results)

# Save results to task-specific CSV file
cat("Saving results to task-specific CSV file...\n")
write.csv(combined_results, 
          paste0("full_sim/results/simulation_results_task_", task_id, ".csv"), 
          row.names = FALSE)

# Print summary for this task
cat("\n=== TASK", task_id, "SUMMARY ===\n")
summary_df <- data.frame(
  Task = task_id,
  NoBench_Avg_Bias = mean(combined_results[, "NoBench_Bias"], na.rm = TRUE),
  NoBench_Avg_Coverage = mean(combined_results[, "NoBench_Coverage"], na.rm = TRUE),
  NoBench_Avg_Runtime = mean(combined_results[, "NoBench_Runtime"], na.rm = TRUE),
  Bench_Avg_Bias = mean(combined_results[, "Bench_Bias"], na.rm = TRUE),
  Bench_Avg_Coverage = mean(combined_results[, "Bench_Coverage"], na.rm = TRUE),
  Bench_Avg_Runtime = mean(combined_results[, "Bench_Runtime"], na.rm = TRUE),
  Avg_Diff_Bench_Pstar = mean(combined_results[, "Diff_Bench_Pstar"], na.rm = TRUE),
  Avg_Diff_NoBench_Pstar = mean(combined_results[, "Diff_NoBench_Pstar"], na.rm = TRUE),
  Avg_Diff_Bench_NoBench = mean(combined_results[, "Diff_Bench_NoBench"], na.rm = TRUE),
  N_Failures = sum(is.na(combined_results[, "NoBench_Runtime"]))
)
print(summary_df)

cat("\n=== GELMAN DIAGNOSTICS > 1.01 ===\n")
cat("No Bench Model:\n")
print(apply(combined_results[, c("NoBench_Gelman_rho", "NoBench_Gelman_log_kappa", 
                                  "NoBench_Gelman_S1", "NoBench_Gelman_P1")] > 1.01, 
            2, sum, na.rm = TRUE))
cat("\nBench Model:\n")
print(apply(combined_results[, c("Bench_Gelman_rho", "Bench_Gelman_log_kappa", 
                                  "Bench_Gelman_S1", "Bench_Gelman_P1")] > 1.01, 
            2, sum, na.rm = TRUE))

cat("\nTask", task_id, "results saved to:\n")
cat("  - full_sim/results/simulation_results_task_", task_id, ".csv\n", sep = "")
cat("\nTask", task_id, "completed successfully!\n")
