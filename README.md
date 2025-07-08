# Population Uncertainty from Census Differential Privacy

## Project Overview

This project addresses the challenge of population count uncertainty resulting from differential privacy (DP) mechanisms implemented by the U.S. Census Bureau. These mechanisms intentionally inject noise into published population counts to protect individual privacy, creating uncertainty especially at fine geographic scales.

We develop hierarchical Bayesian models that leverage:
1. Spatial autocorrelation in population patterns
2. Known higher-level benchmark totals
3. The statistical properties of the noise-injection mechanism

## Models

### Exact Benchmarking Model

The model is formulated as follows (from main.tex):

1. **Data Layer**: Accounts for the noise in the observed data (census counts) in the form of differential privacy noise.
2. **Process Layer**: Describes the spatial process generating the population distribution, incorporating spatial autocorrelation.
3. **Parameter Layer**: Models the uncertainty in the parameters of the spatial process, allowing them to vary across different regions.

#### Model Equations

Let \( y_{i} \) be the noisy census count for area \( i \), \( \theta_{i} \) the true population count for area \( i \), and \( \epsilon_{i} \) the noise added to protect privacy.

**Data Layer**:
\[
y_{i} = \theta_{i} + \epsilon_{i}
\]
where \( \epsilon_{i} \sim \text{Laplace}(0, b) \) due to the Laplacian noise mechanism often used in DP.

**Process Layer**:
\[
\theta_{i} \sim \text{CAR}(\rho, \sigma^2)
\]
where CAR is the Conditional Autoregressive model for spatial autocorrelation, \( \rho \) is the spatial correlation parameter, and \( \sigma^2 \) is the process variance.

**Parameter Layer**:
\[
\rho \sim \text{Beta}(a, b)
\]
\[
\sigma^2 \sim \text{Inverse-Gamma}(c, d)
\]

### Hierarchical Model with Benchmarking

This model extends the Exact Benchmarking Model by incorporating external benchmark totals known for certain regions (e.g., state or county levels). It adjusts the estimates to be consistent with these benchmarks.

#### Model Equations

Let \( B_{j} \) be the known benchmark total for region \( j \).

**Data Layer**:
\[
y_{ij} = \theta_{ij} + \epsilon_{ij}
\]
where \( y_{ij} \) is the noisy count for area \( i \) in region \( j \), and \( \theta_{ij} \) is the true count.

**Process Layer**:
\[
\theta_{ij} \sim \text{CAR}(\rho, \sigma^2)
\]

**Parameter Layer**:
\[
\rho \sim \text{Beta}(a, b)
\]
\[
\sigma^2 \sim \text{Inverse-Gamma}(c, d)
\]
\[
\theta_{j}^* \sim \text{Normal}(\mu, \tau^2)
\]
\[
B_{j} \sim \text{Normal}(\theta_{j}^*, \sigma_B^2)
\]

**Benchmarking Adjustment**:
\[
\theta_{ij} = \theta_{ij} + (B_{j} - \bar{\theta}_{j})
\]
where \( \bar{\theta}_{j} \) is the average of \( \theta_{ij} \) in region \( j \).

## Implementation

The models are implemented in Python using PyMC3 for Bayesian inference. The code is structured as follows:

1. **Data Preparation**: Scripts to download, clean, and prepare the census data and any benchmark data.
2. **Model Definition**: PyMC3 scripts that define the hierarchical models described above.
3. **Inference**: Code to run the MCMC sampling and obtain posterior distributions of the parameters.
4. **Post-Processing**: Scripts to analyze the posterior samples, make predictions, and generate uncertainty intervals for the population estimates.

## Results

The models provide:
- **Adjusted Population Estimates**: More accurate estimates of the population counts with reduced bias.
- **Uncertainty Intervals**: Credible intervals for the population estimates, reflecting the uncertainty due to the differential privacy noise and the modeling process.
- **Parameter Estimates**: Posterior distributions of the model parameters, including the noise parameters and the spatial autocorrelation parameters.

## Conclusion

This project develops hierarchical Bayesian models that effectively adjust for the uncertainty introduced by differential privacy mechanisms in the U.S. Census data. The models leverage spatial autocorrelation and external benchmarks to provide more accurate and reliable population estimates.