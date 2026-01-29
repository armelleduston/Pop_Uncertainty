#!/usr/bin/env Rscript
# combine_results.R
# Combines task-specific CSV files into final results files

library(dplyr)

cat("Combining simulation results...\n")

# Set working directory to project root
setwd("/Users/armelleduston/Documents/GitHub/Pop_Uncertainty")

# Find all task result files
no_bench_files <- list.files("full_sim/results", 
                             pattern = "^no_bench_results_task_\\d+\\.csv$", 
                             full.names = TRUE)
bench_files <- list.files("full_sim/results", 
                         pattern = "^bench_results_task_\\d+\\.csv$", 
                         full.names = TRUE)

cat(paste("Found", length(no_bench_files), "no_bench task files\n"))
cat(paste("Found", length(bench_files), "bench task files\n"))

# Combine no_bench results
if (length(no_bench_files) > 0) {
  no_bench_list <- lapply(no_bench_files, read.csv)
  no_bench_combined <- do.call(rbind, no_bench_list)
  no_bench_combined <- no_bench_combined[order(no_bench_combined$Run), ]
  
  write.csv(no_bench_combined, 
            "full_sim/results/no_bench_results_combined.csv", 
            row.names = FALSE)
  cat("No bench results combined and saved to: full_sim/results/no_bench_results_combined.csv\n")
  
  # Print summary
  cat("\n=== NO BENCH COMBINED SUMMARY ===\n")
  summary_nobench <- data.frame(
    Total_Runs = nrow(no_bench_combined),
    Avg_Bias = mean(no_bench_combined$Bias, na.rm = TRUE),
    Avg_Coverage = mean(no_bench_combined$CI_Coverage, na.rm = TRUE),
    Avg_Runtime = mean(no_bench_combined$Runtime, na.rm = TRUE),
    N_Failures = sum(is.na(no_bench_combined$Runtime))
  )
  print(summary_nobench)
  
  cat("\nGelman diagnostics > 1.01:\n")
  print(colSums(no_bench_combined[, c("Gelman_rho", "Gelman_log_kappa", 
                                       "Gelman_S1", "Gelman_P1")] > 1.01, na.rm = TRUE))
} else {
  cat("No no_bench task files found!\n")
}

# Combine bench results
if (length(bench_files) > 0) {
  bench_list <- lapply(bench_files, read.csv)
  bench_combined <- do.call(rbind, bench_list)
  bench_combined <- bench_combined[order(bench_combined$Run), ]
  
  write.csv(bench_combined, 
            "full_sim/results/bench_results_combined.csv", 
            row.names = FALSE)
  cat("\nBench results combined and saved to: full_sim/results/bench_results_combined.csv\n")
  
  # Print summary
  cat("\n=== BENCH COMBINED SUMMARY ===\n")
  summary_bench <- data.frame(
    Total_Runs = nrow(bench_combined),
    Avg_Bias = mean(bench_combined$Bias, na.rm = TRUE),
    Avg_Coverage = mean(bench_combined$CI_Coverage, na.rm = TRUE),
    Avg_Runtime = mean(bench_combined$Runtime, na.rm = TRUE),
    N_Failures = sum(is.na(bench_combined$Runtime))
  )
  print(summary_bench)
  
  cat("\nGelman diagnostics > 1.01:\n")
  print(colSums(bench_combined[, c("Gelman_rho", "Gelman_log_kappa", 
                                    "Gelman_P1")] > 1.01, na.rm = TRUE))
} else {
  cat("No bench task files found!\n")
}

cat("\nCombination complete!\n")
