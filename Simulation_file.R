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


## Generate Data ##############################################################

source("src/generate_data.r")

n <- 200 # number of points
rho <- 0.7 # spatial correlation       
kappa <- 1.5 # precision scaling variable (for CAR)
tau <- 0.25 # related to census privacy budget     
J <- 10  # number of regions


sim_data <- generate_data(n=n, 
                          rho=rho, 
                          kappa=kappa, 
                          tau=tau,
                          J=J)

###############################################################################


## Define the ddlaplace_nim distribution for NIMBLE ###########################

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

###############################################################################

## Run model with exact benchmarking ##########################################

source("src/complete_model.r")

samples2 <- complete_model(sim_data, exact = TRUE)

###############################################################################

