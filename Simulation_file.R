# Script Name: Simulation_file
# Author: Armelle Duston
# Date: "`r Sys.Date()`"
# Description: This file is used to run simulations of the methods in 
# `Simulation_study.rmd` using the cluster
#
# Dependencies:
#   - R packages: see below
#   - External files: generate_data.r, complete_model.r
#
# Notes:
#   - This file is not meant to be run in isolation but using the slurm file
#     and work on the cluster

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

# Clear workspace 
rm(list = ls())

simulate <- function(rho, kappa, eta) {
  ## Generate Data ##############################################################
  
  source("src/new_generate_data.r")
  
  root_n <- 15 # root of number of points
  rho <- 0.6 # spatial correlation       
  kappa <- kappa # precision scaling variable (for CAR)
  tau <- tau # related to census privacy budget     
  J <- 10  # number of regions
  eta <- eta
  
  
  sim_data <- new_generate_data(root_n, 
                                0.6, 
                                kappa, 
                                tau, 
                                J)
  
  ###############################################################################
  
  
  ## Run model without benchmarking #############################################
  
  source("src/new_model.R")
  
  samples_no_bench <- new_model(sim_data,
                        bench = "none",
                        eta = eta)
  
  ###############################################################################
  
  ## Run model with inexact benchmarking ########################################
  
  source("src/new_model.R")
  
  samples_bench <- new_model(sim_data,
                       bench = "inexact",
                       eta = eta)
  
  ###############################################################################
  
  return ( list(samples_no_bench=samples_no_bench,
                samples_bench=samples_bench)) 
}

# Run via command line with args rho kappa eta and save results
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 3) {
    rho_val <- as.numeric(args[1])
    kappa_val <- as.numeric(args[2])
    eta_val <- as.numeric(args[3])
    dir.create("results", showWarnings = FALSE)
    out <- simulate(rho_val, kappa_val, eta_val)
    saveRDS(out, file = paste0("results/sim_rho", rho_val, "_kappa", kappa_val, "_eta", eta_val, ".rds"))
  }
}



