# spatiocompo

An R package for Bayesian spatial compositional data modeling using Gaussian Markov Random Fields (GMRF) with a Dirichlet likelihood.

## Installation

```r
# From GitHub (the R package lives in the R-package/ subdirectory)
devtools::install_github("BehnazP/spatiocompo-r", subdir = "R-package")
```

## Overview

spatiocompo implements Bayesian inference for spatially referenced compositional data (e.g., land cover proportions). The model uses:

- **Observation model**: y_i | z_i, alpha ~ Dirichlet(alpha * z_i)
- **Latent field**: z = inverse-ALR(x)
- **Spatial prior**: x | kappa, rho ~ GMRF(0, rho x Q(kappa)^{-1})

Inference combines a Metropolis-adjusted Langevin (MALA) update for the latent field, regression coefficients, and Dirichlet precision with a Metropolis-within-Gibbs step for the spatial range (kappa) and cross-component covariance (rho).

## Quick Start

```r
library(spatiocompo)

# Simulate data
sim <- simulate_compo_data(n_side = 10, kappa_true = 0.5,
                            alpha_true = 20, seed = 123)

# Set up priors
priors <- create_priors(sim$G, d = 2)

# Initial values
theta0 <- list(alpha = 15, kappa = 0.3, rho = diag(2), x = sim$x_true)

# Run MCMC
fit <- mcmc_sampling(theta0, priors, sim$y, iter = 2000)

# Check results
plot_trace(fit)
```

## Built-in Data

The package includes European land cover data and a simulated dataset — no external files needed:

```r
# Simulated data
data(simulated_landcover)

# European pollen-based land cover (5 time periods)
data(europe_landcover)
result <- run_europe_model(time = "AD1900", model = "LPJ_KK10_ESM", iter = 5000)

# Plot the reconstructed composition maps, the observations, and the MCMC chains
plot_europe_map(result)        # reconstructed coniferous / broadleaved / unforested
plot_europe_obs("AD1900")      # observed compositions at the pollen sites
plot_trace(result$fit)         # MCMC trace plots
```

See the [R package README](R-package/README.md) for full documentation, covariate examples, and development instructions.

## Original Data and MATLAB Implementation

The original MATLAB implementation and the published land-cover map data of
Pirzamanbein et al. (2018) are maintained in the companion repository
[BehnazP/SpatioCompo](https://github.com/BehnazP/SpatioCompo). The
`europe_landcover` dataset bundled with this package is derived from those
reconstructions.

## Repository Structure

```
spatiocompo-r/
└── R-package/          # The R package (installable via devtools)
    ├── R/              # Source code
    ├── data/           # Built-in datasets (.rda)
    ├── man/            # Documentation
    ├── tests/          # Unit tests
    └── vignettes/      # Tutorials
```

## References

- Pirzamanbein, B., Lindstrom, J., Poska, A., & Gaillard, M.-J. (2018). Modelling Spatial Compositional Data: Reconstructions of past land cover and uncertainties. *Spatial Statistics*, 24, 14-31. [doi:10.1016/j.spasta.2018.03.005](https://doi.org/10.1016/j.spasta.2018.03.005)
- Pirzamanbein, B., Poska, A., & Lindstrom, J. (2020). Bayesian Reconstruction of Past Land Cover From Pollen Data: Model Robustness and Sensitivity to Auxiliary Variables. *Earth and Space Science*, 7(1), e2018EA000547. [doi:10.1029/2018EA000547](https://doi.org/10.1029/2018EA000547)

## License

GPL-3 (see [LICENSE.md](LICENSE.md)). The bundled `europe_landcover` data are derived from the reconstructions of Pirzamanbein et al. (2018); please cite that work when using the data.
