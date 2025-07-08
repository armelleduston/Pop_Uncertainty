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



# 1. Define a spatial structure
n <- 100  
W <- matrix(0, n, n)
for(i in 1:n) {
  W[i, ifelse(i==1, n, i-1)] <- 1
  W[i, ifelse(i==n, 1, i+1)] <- 1
}
D <- diag(rowSums(W))

# 2. Simulate latent spatial field S with (proper) CAR structure
rho      <- 0.7        # spatial correlation parameter (0 < rho < 1)
kappa    <- 1.5          # precision scaling parameter

# Construct precision matrix
Q <- kappa * (D - rho * W)
# Add a small jitter for numerical stability
Q_pd <- Q + diag(1e-6, n)

#  Sample S ~ N(0, Sigma = solve(Q_pd))
Sigma <- solve(Q_pd)
S <- mvrnorm(1, mu = rep(0, n), Sigma = Sigma)

# 3. Simulate true counts P_i ~ Poisson(exp(S_i))
lambda <- exp(S)
P <- rpois(n, lambda)

# 4. Add measurement noise: P*_i | P_i ~ Laplace(P_i, tau)
tau <- 0.25      # known noise scale (related to "privacy budget")
P_star <- round(P + rlaplace(n, mu = 0, sigma=tau))
P_star <- ifelse(P_star < 0, 0, P_star) 


# 5. Define coarse regions and benchmark totals U_j
J <- 10  # number of regions
region_id <- sample(1:J, n, replace = TRUE)
U <- tapply(P_star, region_id, sum)

# 6. Output simulated data
sim_data <- list(
  W = W,
  D = D,
  S = S,
  P = P,
  P_star = P_star,
  region_id = region_id,
  U = as.numeric(U),
  tau = tau,
  rho = rho
)

