# Script Name: Simulation_file
# Author: Armelle Duston
# Description: This file is used to run simulations of the methods in 
# `Simulation_study.rmd` using the cluster
#
# Dependencies:
#   - R packages: see below
#   - External files: generate_data_new.r, new_model.r
#


# Set working directory to this script's directory -----------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Auto-install missing packages ------------------------------------------------
pkgs <- c("ggplot2","tidyr","patchwork","dplyr","knitr","kableExtra",
           "MASS","nimble","coda","extraDistr","igraph","RColorBrewer",
           "nimbleNoBounds","pracma", "R.utils")
missing <- setdiff(pkgs, rownames(installed.packages()))
if(length(missing)) install.packages(missing, 
                                      repos = "https://cloud.r-project.org")

# Load necessary packages
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
library(R.utils)

# Clear workspace 
rm(list = ls())


droundnorm <- nimbleFunction(
  run = function(x = double(0),    
                 mu = double(0),   
                 sigma = double(0),   
                 log = double(0, default = 0)) {
    returnType(double(0))
    # Support: x must be (near) integer
    if (abs(x - round(x)) > 1e-8) {
      if (log == 1) return(-Inf)
      return(0.0)
    }
    if (sigma <= 0) {
      if (log == 1) return(-Inf)
      return(0.0)
    }
    a <- (x - 0.5 - mu) / sigma
    b <- (x + 0.5 - mu) / sigma
    p <- pnorm(b, 0, 1) - pnorm(a, 0, 1)   # CDF difference
    # Guard tiny probabilities to avoid log(0)
    if (p <= 0) {
      if (log == 1) return(-Inf)
      return(0.0)
    }
    if (log == 1) return(log(p))
    return(p)
  }
)

# RNG: sample normal then round to nearest integer
rroundnorm <- nimbleFunction(
  run = function(n = integer(0), mu = double(0), sigma = double(0)) {
    returnType(integer(0))
    if (n != 1) nimPrint("rroundnorm generates one value; ignoring n")
    if (sigma <= 0) nimStop("sigma must be > 0")
    z <- rnorm(1, mean = mu, sd = sigma)
    return(round(z))
  }
)

# Register for BUGS usage
registerDistributions(list(
  droundnorm = list(
    BUGSdist = "droundnorm(mu, sigma)",
    types    = c("value = double(0)", "mu = double(0)", "sigma = double(0)"),
    discrete = TRUE
  )
))


source("src/new_model.R")
source("src/new_generate_data.r")

etas=c(0.1, 0.2, 0.3)
kappas=c(0.01, 0.1, 1, 10)
taus = c(10, 25, 50, 100)
R = 3000

result_mat <- matrix(0, nrow = R, ncol = 4)

for (i in length(etas)){
  for (j in length(kappas)){
    for (k in length(taus)){
      for (r in 1:R){
        root_n <- 15 # root of number of points
        rho <- 0.6 # spatial correlation       
        kappa <- kappas[j] # precision scaling variable (for CAR)
        tau <- taus[k] # related to census privacy budget     
        J <- 10  # number of regions
        sim_data <- new_generate_data(root_n, 
                                      rho, 
                                      kappa, 
                                      tau, 
                                      J)

        
        no_bench <- withTimeout({
          new_model(sim_data, bench = "none")
        }, timeout = 300, onTimeout = "silent")
          
        bench <- withTimeout({
          new_model(sim_data, bench = "inexact", eta = eta[i])
        }, timeout = 300, onTimeout = "silent")
        
        true_Ps <- sim_data_new$P
        if (is.null(no_bench)==FALSE){
          # compute bias & coverage for no benchmarking
          no_bench_Ps <- cbind(no_bench$chain1[,paste0("P[", 1:root_n^2, "]")])
          no_bench_bias <- mean(apply(true_Ps - no_bench_Ps, 2, mean))
          no_bench_CIs <- apply(X=no_bench_Ps, MARGIN=2 , 
                                FUN= function(x) quantile(x, probs = c(0.025, 0.975))) 
          no_bench_coverage <- mean(ifelse(true_Ps >= no_bench_CIs[1,] & true_Ps <= no_bench_CIs[2,], 1, 0))
        } else{
          no_bench_bias <- NA
          bench_coverage <- NA
        }
        
        if(is.null(bench)==FALSE){
          # compute bias & coverage for inexact benchmarking
          bench_Ps <- cbind(no_bench$chain1[,paste0("P[", 1:root_n^2, "]")])
          bench_bias <- mean(apply(true_Ps - bench_Ps, 2, mean))
          bench_CIs <- apply(X=bench_Ps, MARGIN=2 , 
                             FUN= function(x) quantile(x, probs = c(0.025, 0.975))) 
          bench_coverage <- mean(ifelse(true_Ps >= bench_CIs[1,] & true_Ps <= bench_CIs[2,], 1, 0))
        } else{
          bench_bias <- NA
          bench_coverage <- NA
        }
        
        result_mat[r, 1:4] <- cbind(no_bench_bias,
                                    no_bench_coverage,
                                    bench_bias,
                                    bench_coverage)
      }
    }
  }
}


