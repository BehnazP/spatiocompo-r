# spatiocompo

An R package for Bayesian inference on spatial compositional data using Gaussian Markov Random Fields (GMRF) and Dirichlet likelihood.

## Installation

```r
# From GitHub (R package is in the R-package subdirectory)
devtools::install_github("BehnazP/spatiocompo-r", subdir = "R-package")
```

## Overview

spatiocompo implements the methodology from:

> Pirzamanbein, B., Lindström, J., Poska, A., & Gaillard, M. J. (2018).
> Modelling spatial compositional data: Reconstructions of past land cover and uncertainties.
> Spatial statistics, 24, 14-31.
> [doi:10.1016/j.spasta.2018.03.005](https://doi.org/10.1016/j.spasta.2018.03.005)

## Quick Start

```r
library(spatiocompo)

# Simulate data
sim <- simulate_compo_data(n_side = 10, kappa_true = 0.5,
                            alpha_true = 20, seed = 123)

# Set up priors
priors <- create_priors(sim$G, d = 2)

# Initial values
theta0 <- list(
  alpha = 15,
  kappa = 0.3,
  rho = diag(2),
  x = sim$x_true
)

# Run MCMC
fit <- mcmc_sampling(theta0, priors, sim$y, iter = 2000)

# Check results
plot_trace(fit)
```

## Main Functions

| Function | Description |
|----------|-------------|
| `mcmc_sampling()` | Main MCMC sampler |
| `simulate_compo_data()` | Generate synthetic data |
| `create_priors()` | Set up prior distributions |
| `create_Q()` | Create GMRF precision matrix |
| `load_europe_data()` | Load built-in European land cover data |
| `run_europe_model()` | Run MCMC on European data |
| `plot_europe_map()` | Map the reconstructed compositions |
| `plot_europe_obs()` | Map the observed compositions at the sites |
| `plot_field()` | Map a latent field on the grid |
| `plot_trace()` | MCMC trace plots |
| `summary_mcmc()` | Posterior summaries |

## European Land Cover Data

The package includes built-in European land cover data from the paper. No external files needed:

```r
library(spatiocompo)

# Available models:
# - "intercept": Spatial field only (p=1)
# - "elevation": Intercept + standardized elevation (p=2)
# - "LPJ_KK10_ESM": Full model with LPJ-GUESS + elevation (p=4)

# Run with the full model
result <- run_europe_model(
  time = "AD1900",
  model = "LPJ_KK10_ESM",
  iter = 5000
)

# Check the MCMC chains and posterior summaries
plot_trace(result$fit)
summary_mcmc(result$fit)

# Map the reconstructed composition (coniferous / broadleaved / unforested)
plot_europe_map(result)

# Map the observed compositions at the pollen sites
plot_europe_obs("AD1900")

# Map a single latent ALR field on the grid
burnin <- ncol(result$fit$x) %/% 4
x_mean <- rowMeans(result$fit$x[, -(1:burnin)])
x_mat  <- matrix(x_mean, nrow = result$data$grid$N, ncol = result$data$d)
plot_field(x_mat, grid_size = result$data$grid$size,
           lon = result$data$grid$lon, lat = result$data$grid$lat,
           component = 1, main = "Latent field (alr coord 1)")

# Available time periods:
# "AD1900", "AD1725", "AD1425", "BC1000", "BC4000"
```

### Exploring the Data

```r
# Load the raw data
data(europe_landcover)

# See structure
str(europe_landcover, max.level = 2)

# Get AD1900 observations
obs <- europe_landcover$observations$AD1900
head(obs)

# Columns include:
# - lon, lat: coordinates
# - coniferous, broadleaved, unforested: pollen-based compositions
# - elevation: elevation at observation
# - adj_*: LPJ-GUESS predictions adjusted with KK10
# - lpj_*: raw LPJ-GUESS predictions
```

## Model with Custom Covariates

The package supports regression covariates. Use `setup_model()` for easy setup:

```r
library(spatiocompo)

# Your data. Here the bundled europe_landcover stands in for your own:
#   y         : n x D compositional observations (rows sum to 1)
#   coords    : n x 2 matrix of (lon, lat) on a 1-degree grid
#   elevation : numeric vector of length n
data(europe_landcover)
obs       <- europe_landcover$observations$AD1900
y         <- as.matrix(obs[, c("coniferous", "broadleaved", "unforested")])
coords    <- cbind(obs$lon, obs$lat)
elevation <- obs$elevation

# Easy setup with covariates
setup <- setup_model(
  y = y,
  obs_coords = coords,
  covariates = data.frame(elevation = elevation)
)

# Run MCMC
fit <- mcmc_sampling(
  theta0 = setup$theta0,
  priors = setup$priors,
  y = setup$y,
  A = setup$A,
  B = setup$B,
  iter = 5000
)

# Check results
summary_mcmc(fit)
```

### Manual Covariate Setup

For more control, build the covariate and observation matrices yourself.
This continues from the example above (reusing `y`, `coords`, `elevation`, `setup`):

```r
# Covariate matrix (auto-adds an intercept and standardizes)
B <- create_covariate_matrix(
  covariates = data.frame(elevation = elevation),
  d = setup$d,            # D - 1 components
  intercept = TRUE,
  standardize = TRUE
)

# Observation matrix mapping sites to grid cells
A <- create_observation_matrix(coords, setup$grid$coords, d = setup$d)

# Run MCMC with the manually built matrices
fit <- mcmc_sampling(setup$theta0, setup$priors, setup$y,
                     A = A, B = B, iter = 5000)
```

## Model

- **Observation**: y_i | z_i, α ~ Dirichlet(α · z_i)
- **Latent field**: z = inverse-ALR(x)
- **Spatial prior**: x | κ, ρ ~ GMRF(0, ρ ⊗ Q(κ)^{-1})

## Development

### Building from Source

```bash
# Navigate to the R-package directory
cd R-package

# Install pandoc for vignettes (macOS)
brew install pandoc

# Build the package
R CMD build .

# Check the package
R CMD check spatiocompo_*.tar.gz

# Install locally
R CMD INSTALL spatiocompo_*.tar.gz
```

### Running Tests

```r
# In R
devtools::test()

# Or from terminal
Rscript -e "testthat::test_local()"
```

### Generating Documentation

```r
# Requires roxygen2
roxygen2::roxygenise()
```

## Dependencies

**Required:**
- `Matrix` - Sparse matrix operations
- `stats`, `graphics`, `grDevices`, `utils` - Base R (modeling, plotting, progress output)

**Optional:**
- `testthat` - For running tests
- `knitr`, `rmarkdown` - For building vignettes

## License

GPL-3

## Citation

```bibtex
@article{pirzamanbein2018modelling,
  title   = {Modelling spatial compositional data: Reconstructions of past land cover and uncertainties},
  author  = {Pirzamanbein, Behnaz and Lindstr{\"o}m, Johan and Poska, Anneli and Gaillard, Marie-Jos{\'e}},
  journal = {Spatial Statistics},
  volume  = {24},
  pages   = {14--31},
  year    = {2018},
  doi     = {10.1016/j.spasta.2018.03.005}
}
```
