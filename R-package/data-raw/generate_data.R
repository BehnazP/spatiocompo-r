# Script to generate example datasets
# Run this from the package root directory

# Load all R functions
devtools::load_all()

# Generate simulated landcover data
set.seed(2024)
sim <- simulate_landcover(
  n_side = 20,
  pattern = "gradient",
  alpha_true = 80,
  kappa_true = 0.4,
  seed = 2024
)

# Create the dataset in the format we want
simulated_landcover <- list(
  y = sim$y,
  z_true = sim$z_true,
  x_true = sim$x_true,
  G = sim$G,
  coords = sim$grid$coords,
  params = sim$params
)

# Save to data directory
usethis::use_data(simulated_landcover, overwrite = TRUE)

cat("Dataset saved to data/simulated_landcover.rda\n")
