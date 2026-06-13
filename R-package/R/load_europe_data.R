# Suppress R CMD check NOTE for data() variables
utils::globalVariables("europe_landcover")

#' Load European Land Cover Data
#'
#' Loads the bundled European land cover dataset and prepares it for use with
#' mcmc_sampling().
#'
#' @param time Time period to load: "AD1900", "AD1725", "AD1425", "BC1000",
#'   or "BC4000".
#' @param model Model type: "intercept", "elevation", or "LPJ_KK10_ESM".
#'
#' @return A list with components:
#'   \item{y}{Observed compositional data (n x D)}
#'   \item{A}{Location matrix (sparse, n*d x N*d)}
#'   \item{B}{Covariate matrix (sparse, n*d x p*d), or NULL for intercept}
#'   \item{G}{GMRF precision matrix}
#'   \item{grid}{Grid information (size, coordinates)}
#'   \item{coords}{Observation coordinates}
#'
#' @details
#' The dataset contains pollen-based land cover reconstructions for Europe
#' at 5 time periods. Three model types are supported:
#'
#' \describe{
#'   \item{intercept}{Spatial field only, no covariates (p=1)}
#'   \item{elevation}{Intercept + standardized elevation (p=2)}
#'   \item{LPJ_KK10_ESM}{Full model with intercept + ALR(LPJ-adjusted) +
#'     elevation (p=4). This matches the model from the paper.}
#' }
#'
#' @examples
#' # Load AD1900 data with elevation model
#' data <- load_europe_data(time = "AD1900", model = "elevation")
#'
#' # Check dimensions
#' dim(data$y)  # observations
#' dim(data$B)  # covariates
#'
#' # Run MCMC
#' \dontrun{
#' result <- run_europe_model(time = "AD1900", model = "LPJ_KK10_ESM", iter = 5000)
#' }
#'
#' @seealso [run_europe_model()], [mcmc_sampling()]
#'
#' @export
load_europe_data <- function(time = "AD1900", model = "intercept") {

  # Load the data
  data("europe_landcover", package = "spatiocompo", envir = environment())

  if (!time %in% europe_landcover$time_periods) {
    stop("Time '", time, "' not available. Choose from: ",
         paste(europe_landcover$time_periods, collapse = ", "))
  }

  obs <- europe_landcover$observations[[time]]
  grid_data <- europe_landcover$grid[[time]]

  # Extract coordinates
  lon <- obs$lon
  lat <- obs$lat

  # Create compositional matrix
  y_raw <- cbind(obs$coniferous, obs$broadleaved, obs$unforested)
  y <- y_raw / rowSums(y_raw)
  y <- no_zero_one(y, eps = 1e-5)

  n <- nrow(y)
  D <- ncol(y)
  d <- D - 1

  # Create grid
  lon_range <- range(lon)
  lat_range <- range(lat)
  n_lon <- diff(lon_range) + 1
  n_lat <- diff(lat_range) + 1
  N <- n_lon * n_lat
  sz <- c(n_lon, n_lat)

  # Grid coordinates
  # Q ordering: first dim (lon) varies fastest (Kronecker structure)
  grid_lon <- rep(seq(lon_range[1], lon_range[2]), n_lat)
  grid_lat <- rep(seq(lat_range[1], lat_range[2]), each = n_lon)

  # Create observation matrix A (maps observations to nearest grid cell)
  A_idx <- numeric(n)
  for (i in seq_len(n)) {
    dists <- (grid_lon - lon[i])^2 + (grid_lat - lat[i])^2
    A_idx[i] <- which.min(dists)
  }
  A_single <- Matrix::sparseMatrix(i = 1:n, j = A_idx, x = 1, dims = c(n, N))
  A <- Matrix::bdiag(A_single, A_single)

  # Create covariate matrix B based on model
  B <- NULL
  p <- 0

  if (model == "intercept") {
    p <- 1
    B_single <- matrix(1, nrow = n, ncol = 1)
    B <- Matrix::bdiag(B_single, B_single)

  } else if (model == "elevation") {
    elev <- obs$elevation
    elev[elev < 0] <- 0
    tElev <- (elev - mean(elev)) / sd(elev)

    p <- 2
    B_single <- cbind(1, tElev)
    B <- Matrix::bdiag(B_single, B_single)

  } else if (model == "LPJ_KK10_ESM") {
    # LPJ-adjusted compositions
    adj_compo <- cbind(obs$adj_coniferous, obs$adj_broadleaved, obs$adj_unforested)
    adj_compo <- adj_compo / rowSums(adj_compo)
    adj_compo <- no_zero_one(adj_compo, eps = 1e-5)
    alradj <- alr(adj_compo)

    # Elevation
    elev <- obs$elevation
    elev[elev < 0] <- 0
    tElev <- (elev - mean(elev)) / sd(elev)

    p <- 4
    B_single <- cbind(1, alradj, tElev)
    B <- Matrix::bdiag(B_single, B_single)

  } else {
    stop("Model '", model, "' not implemented. ",
         "Available models: 'intercept', 'elevation', 'LPJ_KK10_ESM'")
  }

  # Create GMRF precision matrix
  G <- create_Q(sz, alpha = 1)

  list(
    y = y,
    A = A,
    B = B,
    G = G,
    grid = list(
      size = sz,
      N = N,
      lon = grid_lon,
      lat = grid_lat
    ),
    coords = list(
      lon = lon,
      lat = lat
    ),
    n = n,
    D = D,
    d = d,
    p = p,
    time = time,
    model = model
  )
}


#' Run MCMC on European Land Cover Data
#'
#' Convenience function to load the European land cover data and run
#' MCMC sampling.
#'
#' @param time Time period: "AD1900", "AD1725", "AD1425", "BC1000", or "BC4000".
#' @param model Model type: "intercept", "elevation", or "LPJ_KK10_ESM".
#' @param iter Number of MCMC iterations.
#' @param ... Additional arguments passed to mcmc_sampling().
#'
#' @return A list with MCMC results and data information.
#'
#' @examples
#' \dontrun{
#' # Run with full model
#' result <- run_europe_model(
#'   time = "AD1900",
#'   model = "LPJ_KK10_ESM",
#'   iter = 5000
#' )
#'
#' # Check results
#' plot_trace(result$fit)
#' summary_mcmc(result$fit)
#' }
#'
#' @seealso [load_europe_data()], [mcmc_sampling()]
#'
#' @export
run_europe_model <- function(time = "AD1900", model = "intercept",
                             iter = 2000, ...) {

  # Load data
  data <- load_europe_data(time = time, model = model)

  # Set up priors matching MATLAB defaults
  priors <- create_priors(data$G, d = data$d)
  priors$alpha.a <- 1.5
  priors$alpha.b <- 0.1
  priors$rho.df <- 10

  # Find good initial values via optimization (reduces burn-in)
  theta0 <- find_initial_values(
    y = data$y,
    priors = priors,
    A = data$A,
    B = data$B,
    w = 1
  )

  # Run MCMC
  fit <- mcmc_sampling(
    theta0 = theta0,
    priors = priors,
    y = data$y,
    A = data$A,
    B = data$B,
    iter = iter,
    ...
  )

  list(
    fit = fit,
    data = data,
    priors = priors,
    theta0 = theta0
  )
}
