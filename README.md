# Population Uncertainty from Census Differential Privacy

## Project Overview

This project addresses the challenge of population count uncertainty resulting from differential privacy (DP) mechanisms implemented by the U.S. Census Bureau. These mechanisms intentionally inject noise into published population counts to protect individual privacy, creating uncertainty especially at fine geographic scales.

We develop hierarchical Bayesian models that leverage:
1. Spatial autocorrelation in population patterns
2. Known higher-level benchmark totals
3. The statistical properties of the noise-injection mechanism

## Models

We propose a hierarchical model to capture the noise-injection. The following is a model for exact benchmarking:

![model screenshot](model_screenshot.png)

Equations (3) and (4), which put a spatial prior on the unknown true population counts, will allow the model to leverage the spatial structure in population size/density to recover the latent true population counts. Equation (2) captures the DP noise injection into P_i from known distribution G. Finally, equation (1) encodes the application of the constraints in the DAS post-processing. We will embed the components of this formulation into the likelihood and prior for the disease mapping model.

The inexact benchmarking model uses the following benchmarking prior, with all of the other components as before: 

![inexact prior](inexact_prior.png)

For more details on the models, see main.tex

