#' Create Prior Specification
#'
#' Helper function to create a list of prior hyperparameters with sensible
#' defaults for the spatial compositional model.
#'
#' @param G Base precision matrix from [create_Q()] with alpha = 1.
#' @param d Number of ALR components (D - 1). Default is 2 (for D = 3).
#' @param alpha.a,alpha.b Shape and rate parameters for Gamma prior on alpha.
#'   Default is Gamma(1, 0.01) giving a weakly informative prior.
#' @param kappa.range,kappa.sigma Shape and rate parameters for Gamma prior
#'   on kappa. Default is Gamma(2, 1).
#' @param rho.df Degrees of freedom for Inverse-Wishart prior on rho.
#'   Default is d + 1 (minimally informative).
#' @param rho.Sigma Scale matrix for Inverse-Wishart prior. Default is
#'   identity matrix.
#' @param logstep.MALA Initial log step size for MALA. Default is 0.
#' @param logstep.Gibbs Initial log step size for Gibbs sampler. Default is 0.
#'
#' @return A list containing all prior hyperparameters needed for
#'   [mcmc_sampling()]:
#'   \itemize{
#'     \item alpha.a, alpha.b: Gamma prior for alpha
#'     \item kappa.range, kappa.sigma: Gamma prior for kappa
#'     \item rho.df, rho.Sigma: Inverse-Wishart prior for rho
#'     \item field.G: Base precision matrix
#'     \item field.lambda: Precomputed eigenvalues for log-determinant
#'     \item logstep.MALA, logstep.Gibbs: Initial step sizes
#'   }
#'
#' @details
#' Default priors are:
#' \itemize{
#'   \item alpha ~ Gamma(1, 0.01): Mean 100, weakly informative
#'   \item kappa ~ Gamma(2, 1): Mean 2, mode at 1
#'   \item rho ~ IW(d+1, I_d): Weakly informative, centered at identity
#' }
#'
#' The eigenvalues of G are precomputed for efficient log-determinant
#' calculation during sampling.
#'
#' @seealso [mcmc_sampling()], [find_initial_values()]
#'
#' @examples
#' # 20x20 grid with 3 components (d=2)
#' G <- create_Q(c(20, 20), alpha = 1)
#' priors <- create_priors(G, d = 2)
#' names(priors)
#'
#' # Custom priors for stronger spatial smoothing
#' priors <- create_priors(G, d = 2, kappa.range = 3, kappa.sigma = 2)
#'
#' @export
create_priors <- function(G, d = 2,
                          alpha.a = 1, alpha.b = 0.01,
                          kappa.range = 2, kappa.sigma = 1,
                          rho.df = NULL, rho.Sigma = NULL,
                          logstep.MALA = 0, logstep.Gibbs = 0) {

  # Default Inverse-Wishart parameters
  if (is.null(rho.df)) {
    rho.df <- d + 1
  }
  if (is.null(rho.Sigma)) {
    rho.Sigma <- diag(d)
  }

  # Precompute eigenvalues for efficient log-determinant
  lambda <- comp_lambda(G)

  priors <- list(
    alpha.a = alpha.a,
    alpha.b = alpha.b,
    kappa.range = kappa.range,
    kappa.sigma = kappa.sigma,
    rho.df = rho.df,
    rho.Sigma = rho.Sigma,
    field.G = G,
    field.lambda = lambda,
    logstep.MALA = logstep.MALA,
    logstep.Gibbs = logstep.Gibbs
  )

  return(priors)
}


#' Find Initial Values via Optimization
#'
#' Finds starting values for MCMC by optimizing the posterior, similar to
#' the MATLAB implementation. This significantly reduces burn-in time.
#'
#' @param y Observed compositional data (n x D matrix).
#' @param priors Prior specification from [create_priors()].
#' @param A Observation matrix (sparse, n*d x N*d). Default is identity.
#' @param B Covariate matrix (sparse, n*d x p*d). Default is NULL.
#' @param w Observation weights. Default is 1.
#' @param kappa_init Initial guess for kappa. Default is 0.5.
#' @param alpha_init Initial guess for alpha. Default is 10.
#' @param sigma_beta Prior variance for beta coefficients. Default is 1000.
#' @param verbose Logical; print optimization progress. Default is TRUE.
#'
#' @return A list suitable for the theta0 argument of [mcmc_sampling()]:
#'   \itemize{
#'     \item alpha: Optimized Dirichlet precision
#'     \item kappa: Optimized spatial range
#'     \item rho: Estimated d x d covariance matrix (from conditional posterior)
#'     \item x: Optimized latent field (N x d matrix)
#'     \item beta: Optimized regression coefficients (p x d matrix, if B provided)
#'   }
#'
#' @details
#' The optimization follows the MATLAB implementation:
#' \enumerate{
#'   \item Optimize (beta, x, alpha) jointly using L-BFGS-B with gradients
#'   \item Optimize kappa given x using bounded optimization
#'   \item Sample rho from its conditional inverse-Wishart posterior
#' }
#'
#' This provides much better starting values than naive initialization,
#' reducing the required burn-in period.
#'
#' @seealso [mcmc_sampling()], [create_priors()]
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' sim <- simulate_compo_data(n_side = 10, seed = 123)
#' priors <- create_priors(sim$G, d = 2)
#'
#' # Find optimized initial values
#' theta0 <- find_initial_values(sim$y, priors)
#'
#' # Run MCMC with optimized start
#' fit <- mcmc_sampling(theta0, priors, sim$y, iter = 2000)
#' }
#'
#' @export
find_initial_values <- function(y, priors, A = NULL, B = NULL, w = 1,
                                kappa_init = 0.5, alpha_init = 10,
                                sigma_beta = 1000, verbose = TRUE) {

  if (is.vector(y)) y <- matrix(y, nrow = 1)
  n <- nrow(y)
  D <- ncol(y)
  d <- D - 1

  N <- nrow(priors$field.G)

  # Ensure w is vector
  if (length(w) == 1) w <- rep(w, n)

  # Default A is identity
  if (is.null(A)) {
    A <- Matrix::Diagonal(N * d)
  }

  # Handle covariates
  if (!is.null(B)) {
    p <- ncol(B) / d
    AB <- cbind(B, A)
  } else {
    p <- 0
    AB <- A
  }

  # Initial precision matrix Q0 (with kappa_init, rho = I)
  rho_init <- diag(d)
  rhoxQ <- Q_rhoxQ(rho_init, kappa_init, priors$field.G)$rhoxQ

  if (p > 0) {
    Q0 <- Matrix::bdiag(
      Matrix::Diagonal(d * p, 1 / sigma_beta),
      rhoxQ
    )
  } else {
    Q0 <- rhoxQ
  }

  # Initialize x from data at observed locations (much better than zeros)
  y_safe <- no_zero_one(y, eps = 1e-4)
  x_init <- rep(0, N * d)
  # A is (n*d) x (N*d) block-diagonal; extract observation-to-grid mapping
  # from the first block (rows 1:n, cols 1:N)
  A_single <- A[1:n, 1:N, drop = FALSE]
  alr_y <- alr(y_safe)
  for (k in seq_len(d)) {
    offset <- (k - 1) * N
    x_k <- numeric(N)
    counts <- numeric(N)
    for (i in seq_len(n)) {
      j <- which(A_single[i, ] != 0)
      if (length(j) == 1) {
        x_k[j] <- x_k[j] + alr_y[i, k]
        counts[j] <- counts[j] + 1
      }
    }
    has_obs <- counts > 0
    x_k[has_obs] <- x_k[has_obs] / counts[has_obs]
    x_init[offset + seq_len(N)] <- x_k
  }

  # Initial state
  if (p > 0) {
    mu0 <- c(rep(0, p * d), x_init)
  } else {
    mu0 <- x_init
  }
  state0 <- c(mu0, alpha_init)

  if (verbose) cat("Optimizing [beta, x, alpha]...\n")

  # Objective function (negative log-posterior)
  obj_fn <- function(state) {
    x_vec <- state[1:(length(state) - 1)]
    alpha <- state[length(state)]

    if (alpha <= 0) return(1e10)

    result <- neg_log_lik(x_vec, alpha, AB, y, priors$alpha.a, priors$alpha.b,
                          Q0, d, w)
    return(result$L)
  }

  # Gradient function
  grad_fn <- function(state) {
    x_vec <- state[1:(length(state) - 1)]
    alpha <- state[length(state)]

    if (alpha <= 0) return(rep(0, length(state)))

    result <- neg_log_lik(x_vec, alpha, AB, y, priors$alpha.a, priors$alpha.b,
                          Q0, d, w)
    return(result$dL)
  }

  # Optimize using L-BFGS-B
  opt_result <- tryCatch({
    optim(
      par = state0,
      fn = obj_fn,
      gr = grad_fn,
      method = "L-BFGS-B",
      lower = c(rep(-Inf, length(state0) - 1), 0.1),
      upper = c(rep(Inf, length(state0) - 1), 1000),
      control = list(maxit = 500, factr = 1e7)
    )
  }, error = function(e) {
    if (verbose) cat("  Optimization failed, using defaults\n")
    list(par = state0, convergence = 1)
  })

  if (verbose) {
    if (opt_result$convergence == 0) {
      cat("  Converged.\n")
    } else {
      cat("  Warning: optimization may not have converged.\n")
    }
  }

  # Extract optimized values
  state_opt <- opt_result$par
  alpha_opt <- state_opt[length(state_opt)]

  if (p > 0) {
    beta_opt <- matrix(state_opt[1:(p * d)], nrow = p, ncol = d, byrow = FALSE)
    x_opt <- matrix(state_opt[(p * d + 1):(length(state_opt) - 1)],
                    nrow = N, ncol = d)
  } else {
    beta_opt <- NULL
    x_opt <- matrix(state_opt[1:(length(state_opt) - 1)], nrow = N, ncol = d)
  }

  # Optimize kappa given x
  if (verbose) cat("Optimizing kappa...\n")

  kappa_obj <- function(kappa) {
    if (kappa <= 0) return(1e10)
    log_posterior_kappa_x(x_opt, kappa, priors)
  }

  kappa_result <- tryCatch({
    optimize(kappa_obj, interval = c(0.01, 5), tol = 0.01)
  }, error = function(e) {
    list(minimum = kappa_init)
  })

  kappa_opt <- kappa_result$minimum
  if (verbose) cat("  kappa =", round(kappa_opt, 4), "\n")

  # Sample rho from conditional posterior
  if (verbose) cat("Sampling rho from conditional posterior...\n")
  rho_opt <- sample_invwishart(x_opt, kappa_opt, priors)

  if (verbose) {
    cat("\nOptimized initial values:\n")
    cat("  alpha =", round(alpha_opt, 2), "\n")
    cat("  kappa =", round(kappa_opt, 4), "\n")
    if (p > 0) {
      cat("  beta =\n")
      print(round(beta_opt, 3))
    }
  }

  # Return theta0
  theta0 <- list(
    alpha = alpha_opt,
    kappa = kappa_opt,
    rho = rho_opt,
    x = x_opt
  )

  if (p > 0) {
    theta0$beta <- beta_opt
  }

  return(theta0)
}
