#' Simulate Spatial Compositional Data
#'
#' Generates synthetic spatial compositional data from the generative model,
#' useful for testing and validating the MCMC sampler.
#'
#' @param n_side Grid size (creates n_side x n_side grid). Default is 20.
#' @param D Number of compositional components. Default is 3.
#' @param kappa_true True spatial range parameter. Default is 0.5.
#' @param rho_true True covariance matrix (d x d). Default is identity.
#' @param alpha_true True Dirichlet precision parameter. Default is 50.
#' @param neumann Use Neumann boundary conditions. Default is TRUE.
#' @param seed Random seed for reproducibility. Default is NULL.
#'
#' @return A list containing:
#'   \item{y}{Observed compositional data (N x D matrix)}
#'   \item{x_true}{True latent field (N x d matrix)}
#'   \item{z_true}{True compositional means (N x D matrix)}
#'   \item{G}{Base precision matrix}
#'   \item{Q}{Full precision matrix}
#'   \item{params}{List of true parameters (kappa, rho, alpha)}
#'   \item{grid}{List with n_side and coordinates}
#'
#' @details
#' The generative process is:
#' \enumerate{
#'   \item Generate latent field x from \eqn{GMRF(0, (\rho \otimes Q(\kappa))^{-1})}
#'   \item Transform to compositions: z = invalr(x)
#'   \item Sample observations: y_i ~ Dirichlet(alpha * z_i)
#' }
#'
#' @seealso [simulate_landcover()] for a more realistic example,
#'   [mcmc_sampling()] for inference.
#'
#' @examples
#' # Basic simulation
#' sim <- simulate_compo_data(n_side = 15, kappa_true = 0.3,
#'                            alpha_true = 100, seed = 42)
#'
#' # Check dimensions
#' dim(sim$y)  # 225 x 3
#' dim(sim$x_true)  # 225 x 2
#'
#' # Visualize one component
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   df <- data.frame(
#'     x = sim$grid$coords[, 1],
#'     y = sim$grid$coords[, 2],
#'     z1 = sim$z_true[, 1]
#'   )
#'   ggplot(df, aes(x, y, fill = z1)) +
#'     geom_tile() +
#'     scale_fill_viridis_c() +
#'     coord_fixed() +
#'     theme_minimal()
#' }
#'
#' @export
simulate_compo_data <- function(n_side = 20, D = 3, kappa_true = 0.5,
                                rho_true = NULL, alpha_true = 50,
                                neumann = TRUE, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  d <- D - 1
  N <- n_side^2

  # Default covariance
  if (is.null(rho_true)) {
    rho_true <- diag(d)
  }

  # Create grid coordinates
  coords <- expand.grid(x = 1:n_side, y = 1:n_side)

  # Create precision matrices
  G <- create_Q(c(n_side, n_side), alpha = 1, neumann = neumann)
  Q_result <- Q_rhoxQ(rho_true, kappa_true, G)
  Q <- Q_result$Q
  rhoxQ <- Q_result$rhoxQ

  # Sample latent field from GMRF
  # x ~ N(0, (rhoxQ)^{-1})
  # Use Cholesky of precision: if L = chol(Q), then x = L^{-T} z where z ~ N(0,I)

  # Add small diagonal for numerical stability
  rhoxQ_reg <- rhoxQ + 1e-6 * Matrix::Diagonal(N * d)

  # Cholesky factor
  L <- Matrix::chol(rhoxQ_reg)

  # Sample x
  z_noise <- rnorm(N * d)
  x_vec <- as.vector(Matrix::solve(L, z_noise))
  x_true <- matrix(x_vec, nrow = N, ncol = d)

  # Transform to compositional space
  z_true <- invalr(x_true)

  # Sample observations from Dirichlet
  y <- matrix(0, nrow = N, ncol = D)
  for (i in seq_len(N)) {
    concentration <- alpha_true * z_true[i, ]
    y[i, ] <- rdirichlet(1, concentration)
  }

  # Handle any numerical issues
  y <- no_zero_one(y)

  return(list(
    y = y,
    x_true = x_true,
    z_true = z_true,
    G = G,
    Q = Q,
    params = list(
      kappa = kappa_true,
      rho = rho_true,
      alpha = alpha_true
    ),
    grid = list(
      n_side = n_side,
      N = N,
      coords = as.matrix(coords)
    )
  ))
}


#' Simulate Land Cover Data
#'
#' Generates realistic spatial compositional data mimicking European land
#' cover patterns with forest, open land, and other categories.
#'
#' @param n_side Grid size. Default is 30.
#' @param pattern Type of spatial pattern:
#'   \itemize{
#'     \item "gradient": North-south gradient (more forest in north)
#'     \item "patchy": Clustered patches
#'     \item "random": Purely random (no trend)
#'   }
#' @param alpha_true Dirichlet precision (observation noise). Default is 80.
#' @param kappa_true Spatial range parameter. Default is 0.3.
#' @param seed Random seed.
#'
#' @return A list with the same structure as [simulate_compo_data()], plus:
#'   \item{categories}{Names of land cover categories}
#'   \item{trend}{Matrix of deterministic trend (if pattern != "random")}
#'
#' @details
#' This function creates more realistic simulated data with:
#' \itemize{
#'   \item Three categories: Forest, Open (grass/agriculture), Other (urban/water)
#'   \item Spatial correlation via GMRF
#'   \item Optional large-scale trend (e.g., latitude gradient)
#'   \item Cross-component correlation (forest and open are negatively correlated)
#' }
#'
#' @examples
#' # Simulate with north-south gradient
#' sim <- simulate_landcover(n_side = 25, pattern = "gradient", seed = 123)
#'
#' # Check categories
#' sim$categories
#'
#' # Compare true means vs observations
#' par(mfrow = c(2, 3))
#' for (i in 1:3) {
#'   image(matrix(sim$z_true[, i], 25, 25),
#'         main = paste("True", sim$categories[i]))
#' }
#' for (i in 1:3) {
#'   image(matrix(sim$y[, i], 25, 25),
#'         main = paste("Observed", sim$categories[i]))
#' }
#'
#' @export
simulate_landcover <- function(n_side = 30, pattern = c("gradient", "patchy", "random"),
                               alpha_true = 80, kappa_true = 0.3, seed = NULL) {

  pattern <- match.arg(pattern)

  if (!is.null(seed)) set.seed(seed)

  D <- 3  # Forest, Open, Other
  d <- D - 1
  N <- n_side^2

  # Covariance with negative correlation between forest and open
  rho_true <- matrix(c(1.0, -0.3,
                       -0.3, 1.0), nrow = 2, ncol = 2)

  # Create grid
  coords <- expand.grid(x = 1:n_side, y = 1:n_side)

  # Create base field
  G <- create_Q(c(n_side, n_side), alpha = 1, neumann = TRUE)
  Q_result <- Q_rhoxQ(rho_true, kappa_true, G)
  rhoxQ <- Q_result$rhoxQ

  # Add small diagonal for stability
  rhoxQ_reg <- rhoxQ + 1e-6 * Matrix::Diagonal(N * d)
  L <- Matrix::chol(rhoxQ_reg)

  # Sample random component
  z_noise <- rnorm(N * d)
  x_random <- as.vector(Matrix::solve(L, z_noise))
  x_random <- matrix(x_random, nrow = N, ncol = d)

  # Add trend based on pattern
  if (pattern == "gradient") {
    # North-south gradient: more forest in north
    y_coord <- coords[, 2] / n_side  # Scaled to [0, 1]
    trend <- cbind(
      2 * (y_coord - 0.5),  # Forest increases north
      -1 * (y_coord - 0.5)  # Open decreases north
    )
  } else if (pattern == "patchy") {
    # Create patchy landscape using sine waves
    x_coord <- coords[, 1] / n_side * 2 * pi
    y_coord <- coords[, 2] / n_side * 2 * pi
    trend <- cbind(
      sin(2 * x_coord) * cos(3 * y_coord),
      cos(2 * x_coord) * sin(2 * y_coord)
    )
  } else {
    trend <- matrix(0, nrow = N, ncol = d)
  }

  # Combine trend and random
  x_true <- trend + x_random

  # Transform to compositions
  z_true <- invalr(x_true)

  # Sample observations
  y <- matrix(0, nrow = N, ncol = D)
  for (i in seq_len(N)) {
    concentration <- alpha_true * z_true[i, ]
    y[i, ] <- rdirichlet(1, concentration)
  }
  y <- no_zero_one(y)

  categories <- c("Forest", "Open", "Other")

  return(list(
    y = y,
    x_true = x_true,
    z_true = z_true,
    G = G,
    Q = Q_result$Q,
    params = list(
      kappa = kappa_true,
      rho = rho_true,
      alpha = alpha_true
    ),
    grid = list(
      n_side = n_side,
      N = N,
      coords = as.matrix(coords)
    ),
    categories = categories,
    trend = trend
  ))
}
