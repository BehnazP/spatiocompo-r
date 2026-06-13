#' Create Covariate Design Matrix
#'
#' Creates the block-diagonal covariate matrix B for use with mcmc_sampling().
#' This function handles the formatting required for the spatial compositional model.
#'
#' @param covariates A matrix or data frame of covariates (n x p), or a list of
#'   covariate vectors. An intercept column of 1s is automatically added unless
#'   `intercept = FALSE`.
#' @param d Number of compositional components minus 1 (D - 1).
#' @param intercept Logical; if TRUE (default), add an intercept column.
#' @param standardize Logical; if TRUE (default), standardize numeric covariates
#'   to have mean 0 and sd 1.
#'
#' @return A sparse block-diagonal matrix B of dimension (n*d) x (p*d), where
#'   p includes the intercept if added.
#'
#' @details
#' The spatial compositional model uses covariates in the linear predictor:
#' \deqn{\eta_i = B_i \beta + A_i x}
#' where \eqn{\eta_i} is the ALR-transformed composition at location i.
#'
#' This function creates the block-diagonal matrix B that applies the same
#' covariates to each compositional component. The block structure is:
#' \deqn{B = \begin{bmatrix} B_{single} & 0 \\ 0 & B_{single} \end{bmatrix}}
#'
#' @examples
#' # Create covariates from elevation
#' n <- 100
#' elevation <- rnorm(n)
#' B <- create_covariate_matrix(elevation, d = 2)
#' dim(B)  # (200, 4) for n*d x p*d where p=2 (intercept + elevation)
#'
#' # Multiple covariates
#' covs <- data.frame(
#'   elevation = rnorm(n),
#'   temperature = rnorm(n)
#' )
#' B <- create_covariate_matrix(covs, d = 2)
#' dim(B)  # (200, 6) for p=3
#'
#' # Without intercept
#' B <- create_covariate_matrix(elevation, d = 2, intercept = FALSE)
#' dim(B)  # (200, 2) for p=1
#'
#' @seealso [mcmc_sampling()], [create_observation_matrix()]
#'
#' @export
create_covariate_matrix <- function(covariates, d, intercept = TRUE,
                                    standardize = TRUE) {
  # Handle different input types
  if (is.vector(covariates) && !is.list(covariates)) {
    covariates <- matrix(covariates, ncol = 1)
  } else if (is.list(covariates) && !is.data.frame(covariates)) {
    covariates <- do.call(cbind, covariates)
  } else if (is.data.frame(covariates)) {
    covariates <- as.matrix(covariates)
  }

  n <- nrow(covariates)

  # Standardize numeric columns
  if (standardize) {
    for (j in seq_len(ncol(covariates))) {
      if (is.numeric(covariates[, j]) && sd(covariates[, j]) > 0) {
        covariates[, j] <- (covariates[, j] - mean(covariates[, j])) /
          sd(covariates[, j])
      }
    }
  }

  # Add intercept
  if (intercept) {
    covariates <- cbind(1, covariates)
  }

  p <- ncol(covariates)

  # Create block-diagonal matrix
  B <- Matrix::bdiag(replicate(d, covariates, simplify = FALSE))

  return(B)
}


#' Create Observation Matrix
#'
#' Creates the sparse observation matrix A that maps the spatial field to
#' observation locations.
#'
#' @param obs_coords Matrix (n x 2) of observation coordinates (lon, lat).
#' @param grid_coords Matrix (N x 2) of grid coordinates, or a list with
#'   components `lon` and `lat` vectors.
#' @param d Number of compositional components minus 1 (D - 1).
#'
#' @return A sparse block-diagonal matrix A of dimension (n*d) x (N*d).
#'
#' @details
#' Each observation is assigned to its nearest grid cell using Euclidean
#' distance. The matrix A has exactly one 1 per row, indicating which
#' grid cell each observation is associated with.
#'
#' @examples
#' # Create a 10x10 grid
#' grid_lon <- rep(1:10, each = 10)
#' grid_lat <- rep(1:10, 10)
#'
#' # Some observation locations
#' obs_lon <- runif(20, 1, 10)
#' obs_lat <- runif(20, 1, 10)
#'
#' A <- create_observation_matrix(
#'   obs_coords = cbind(obs_lon, obs_lat),
#'   grid_coords = cbind(grid_lon, grid_lat),
#'   d = 2
#' )
#' dim(A)  # (40, 200) for n*d x N*d
#'
#' @seealso [mcmc_sampling()], [create_covariate_matrix()]
#'
#' @export
create_observation_matrix <- function(obs_coords, grid_coords, d) {

  # Handle grid_coords as list
  if (is.list(grid_coords) && !is.matrix(grid_coords)) {
    grid_coords <- cbind(grid_coords$lon, grid_coords$lat)
  }

  n <- nrow(obs_coords)
  N <- nrow(grid_coords)

  # Find nearest grid cell for each observation
  A_idx <- numeric(n)
  for (i in seq_len(n)) {
    dists <- (grid_coords[, 1] - obs_coords[i, 1])^2 +
      (grid_coords[, 2] - obs_coords[i, 2])^2
    A_idx[i] <- which.min(dists)
  }

  # Create sparse matrix
  A_single <- Matrix::sparseMatrix(
    i = seq_len(n),
    j = A_idx,
    x = 1,
    dims = c(n, N)
  )

  # Block diagonal for d components
  A <- Matrix::bdiag(replicate(d, A_single, simplify = FALSE))

  return(A)
}


#' Set Up Spatial Compositional Model
#'
#' High-level function to set up all components needed for mcmc_sampling()
#' from raw data.
#'
#' @param y Observed compositional data (n x D matrix).
#' @param obs_coords Observation coordinates (n x 2 matrix with lon, lat).
#' @param grid_size Grid dimensions as c(n_lon, n_lat), or NULL to auto-detect.
#' @param covariates Optional covariates (n x p matrix, data frame, or list).
#' @param standardize Logical; standardize covariates (default TRUE).
#'
#' @return A list with components:
#'   \item{y}{Processed compositional data}
#'   \item{A}{Observation matrix}
#'   \item{B}{Covariate matrix (NULL if no covariates)}
#'   \item{G}{GMRF precision matrix}
#'   \item{priors}{Default priors}
#'   \item{theta0}{Initial values}
#'   \item{grid}{Grid information}
#'   \item{d}{Number of ALR components}
#'   \item{p}{Number of covariates (including intercept)}
#'
#' @examples
#' \dontrun{
#' # Set up model with elevation covariate
#' setup <- setup_model(
#'   y = observations,
#'   obs_coords = cbind(lon, lat),
#'   covariates = data.frame(elevation = elev)
#' )
#'
#' # Run MCMC
#' fit <- mcmc_sampling(
#'   theta0 = setup$theta0,
#'   priors = setup$priors,
#'   y = setup$y,
#'   A = setup$A,
#'   B = setup$B,
#'   iter = 5000
#' )
#' }
#'
#' @seealso [mcmc_sampling()], [create_covariate_matrix()],
#'   [create_observation_matrix()]
#'
#' @export
setup_model <- function(y, obs_coords, grid_size = NULL, covariates = NULL,
                        standardize = TRUE) {

  # Dimensions
  n <- nrow(y)
  D <- ncol(y)
  d <- D - 1

  # Process y
  y <- y / rowSums(y)
  y <- no_zero_one(y)

  # Determine grid
  lon <- obs_coords[, 1]
  lat <- obs_coords[, 2]

  if (is.null(grid_size)) {
    lon_range <- range(lon)
    lat_range <- range(lat)
    n_lon <- diff(lon_range) + 1
    n_lat <- diff(lat_range) + 1
    grid_size <- c(n_lon, n_lat)
  }

  N <- prod(grid_size)

  # Create grid coordinates
  # Q ordering: first dim (lon) varies fastest (Kronecker structure)
  grid_lon <- rep(seq(min(lon), by = 1, length.out = grid_size[1]),
                  grid_size[2])
  grid_lat <- rep(seq(min(lat), by = 1, length.out = grid_size[2]),
                  each = grid_size[1])
  grid_coords <- cbind(grid_lon, grid_lat)

  # Create observation matrix
  A <- create_observation_matrix(obs_coords, grid_coords, d)

  # Create covariate matrix
  if (!is.null(covariates)) {
    B <- create_covariate_matrix(covariates, d, intercept = TRUE,
                                 standardize = standardize)
    p <- ncol(B) / d
  } else {
    B <- NULL
    p <- 0
  }

  # Create GMRF precision matrix
  G <- create_Q(grid_size, alpha = 1)

  # Set up priors (MATLAB defaults)
  priors <- create_priors(G, d = d)
  priors$alpha.a <- 1.5
  priors$alpha.b <- 0.1
  priors$rho.df <- 10
  priors$kappa.range <- 1
  priors$kappa.sigma <- log(100) / sqrt(8)

  # Initial values
  set.seed(123)
  x0 <- matrix(rnorm(N * d, 0, 0.5), nrow = N, ncol = d)

  theta0 <- list(
    alpha = 10,
    kappa = 0.5,
    rho = diag(d),
    x = x0
  )

  if (p > 0) {
    theta0$beta <- matrix(0, nrow = p, ncol = d)
  }

  list(
    y = y,
    A = A,
    B = B,
    G = G,
    priors = priors,
    theta0 = theta0,
    grid = list(
      size = grid_size,
      N = N,
      lon = grid_lon,
      lat = grid_lat,
      coords = grid_coords
    ),
    coords = obs_coords,
    n = n,
    D = D,
    d = d,
    p = p
  )
}
