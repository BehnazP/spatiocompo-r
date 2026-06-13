#' MCMC Sampling for Spatial Compositional Data
#'
#' Main function for Bayesian inference on spatial compositional data using
#' MCMC. Combines MALA for the latent field and alpha with Gibbs sampling
#' for the spatial parameters kappa and rho.
#'
#' @param theta0 List of initial values:
#'   \itemize{
#'     \item alpha: Dirichlet precision parameter (positive scalar)
#'     \item kappa: Spatial range parameter (positive scalar, or 0 for no spatial)
#'     \item rho: Covariance matrix (d x d positive definite)
#'     \item x: Latent field (N x d matrix)
#'     \item beta: Regression coefficients (p x d matrix, optional)
#'   }
#' @param priors List of prior hyperparameters:
#'   \itemize{
#'     \item alpha.a, alpha.b: Gamma(a, b) prior for alpha
#'     \item kappa.range, kappa.sigma: Gamma prior for kappa
#'     \item rho.df, rho.Sigma: Inverse-Wishart prior for rho
#'     \item field.G: Base precision matrix from [create_Q()]
#'     \item field.lambda: Eigenvalues from [comp_lambda()]
#'     \item logstep.MALA: Initial log step size for MALA
#'     \item logstep.Gibbs: Initial log step size for Gibbs
#'   }
#' @param y Observed compositional data (n x D matrix, rows sum to 1).
#' @param w Observation weights (scalar or vector of length n). Default is 1.
#' @param A Location matrix (n*d x N*d sparse matrix) connecting observations
#'   to grid cells. Default is identity.
#' @param B Covariate matrix (n x p) for regression effects. Default is NULL.
#' @param iter Number of MCMC iterations.
#' @param thin Thinning interval for x samples when iter > 1e5. Default is 50.
#' @param verbose Logical; if TRUE, show progress bar. Default is TRUE.
#'
#' @return A list with components:
#'   \item{alpha}{Vector of alpha samples (length iter)}
#'   \item{kappa}{Vector of kappa samples (length iter)}
#'   \item{rho}{Matrix of rho samples (d.d x iter)}
#'   \item{x}{Matrix of x samples (N.d x iter/thin)}
#'   \item{beta}{Matrix of beta samples (p.d x iter), if B provided}
#'   \item{MALA}{List with acceptance rate (acc), counts, and step sizes}
#'   \item{Gibbs}{List with acceptance rate (acc), counts, and step sizes}
#'
#' @details
#' The model assumes:
#' \deqn{y_i | z_i, \alpha \sim Dirichlet(\alpha \cdot z_i)}
#' \deqn{z_i = invalr(A x_i + B \beta)}
#' \deqn{x | \kappa, \rho \sim GMRF(0, \rho \otimes Q(\kappa))}
#'
#' The sampling algorithm alternates between:
#' \enumerate{
#'   \item MALA step for (x, beta, alpha) given (kappa, rho)
#'   \item Gibbs step for (kappa, rho) given x
#' }
#'
#' For long runs (iter > 1e5), the x samples are thinned to reduce memory.
#'
#' @references
#' Pirzamanbein, B., Lindstrom, J., Poska, A., Gaillard, M.-J. (2018).
#' Modelling spatial compositional data: Reconstructions of past land
#' cover and uncertainties. \emph{Spatial Statistics} 24: 14--31.
#' \doi{10.1016/j.spasta.2018.03.005}
#'
#' @seealso [create_priors()], [simulate_compo_data()], [plot_trace()]
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' sim <- simulate_compo_data(n_side = 10, kappa_true = 0.5,
#'                            alpha_true = 50, seed = 123)
#'
#' # Set up priors
#' priors <- create_priors(sim$G, d = 2)
#'
#' # Initial values
#' theta0 <- list(
#'   alpha = 30,
#'   kappa = 0.3,
#'   rho = diag(2),
#'   x = matrix(0, nrow = 100, ncol = 2)
#' )
#'
#' # Run MCMC
#' fit <- mcmc_sampling(theta0, priors, sim$y, iter = 1000)
#'
#' # Check results
#' plot_trace(fit)
#' }
#'
#' @export
mcmc_sampling <- function(theta0, priors, y, w = 1, A = NULL, B = NULL,
                          iter = 1000, thin = 50, verbose = TRUE) {

  # Dimensions
  if (is.vector(y)) y <- matrix(y, nrow = 1)
  n <- nrow(y)
  D <- ncol(y)
  d <- D - 1

  N <- nrow(theta0$x)

  # Default A is identity (observations at grid points)
  if (is.null(A)) {
    A <- Matrix::Diagonal(N * d)
  }

  # Handle covariates
  if (!is.null(B)) {
    p <- ncol(B) / d
    AB <- cbind(B, A)
    old <- c(as.vector(theta0$beta), as.vector(theta0$x))
  } else {
    p <- 0
    AB <- A
    old <- as.vector(theta0$x)
  }

  # Ensure w is a vector
  if (length(w) == 1) w <- rep(w, n)

  # Initialize storage
  MCMC <- list()
  MCMC$alpha <- numeric(iter)
  MCMC$alpha[1] <- theta0$alpha

  MCMC$kappa <- numeric(iter)
  MCMC$kappa[1] <- theta0$kappa

  MCMC$rho <- matrix(0, nrow = d * d, ncol = iter)
  MCMC$rho[, 1] <- as.vector(theta0$rho)

  # Thin x samples for long runs
  if (iter > 1e5) {
    x_iter <- iter / thin
  } else {
    x_iter <- iter
    thin <- 1
  }
  MCMC$x <- matrix(0, nrow = N * d, ncol = x_iter)
  MCMC$x[, 1] <- as.vector(theta0$x)

  if (p > 0) {
    MCMC$beta <- matrix(0, nrow = p * d, ncol = iter)
    MCMC$beta[, 1] <- as.vector(theta0$beta)

    # Horseshoe hyperparameters
    sigma_beta <- rep(1, p * d)
    MCMC$l2 <- matrix(0, nrow = p, ncol = iter)
    MCMC$l2[, 1] <- 1
    MCMC$tau2 <- numeric(iter)
    MCMC$tau2[1] <- 1
  }

  # Step sizes
  MCMC$Gibbs.logstep <- numeric(iter)
  MCMC$Gibbs.logstep[1] <- priors$logstep.Gibbs

  MCMC$MALA.logstep <- numeric(iter)
  MCMC$MALA.logstep[1] <- priors$logstep.MALA

  # Acceptance counts
  MCMC$Gibbs.count <- numeric(iter)
  MCMC$MALA.count <- numeric(iter)

  # Progress tracking
  if (verbose) {
    mcmc_start_time <- proc.time()[3]
    is_terminal <- interactive() || isatty(stdout())
    # In terminal: update ~100 times with \r (overwrite line)
    # In file/pipe: update ~50 times with \n (new lines)
    progress_interval <- if (is_terminal) max(1, floor(iter / 100))
                         else max(1, floor(iter / 50))

    format_duration <- function(secs) {
      if (secs < 60) return(sprintf("%.0fs", secs))
      if (secs < 3600) return(sprintf("%.0fm %02.0fs", secs %/% 60, secs %% 60))
      sprintf("%.0fh %02.0fm", secs %/% 3600, (secs %% 3600) %/% 60)
    }

    show_progress <- function(i, iter, start_time) {
      elapsed <- proc.time()[3] - start_time
      done <- i - 1
      per_iter <- elapsed / done
      remaining <- per_iter * (iter - i)
      pct <- done / (iter - 1) * 100
      bar_width <- 30
      filled <- round(bar_width * pct / 100)
      bar <- paste0("[", strrep("=", filled), strrep(" ", bar_width - filled), "]")
      msg <- sprintf("MCMC %s %5.1f%% | iter %d/%d | %s elapsed | %s/iter | ~%s remaining",
                     bar, pct, i, iter,
                     format_duration(elapsed),
                     sprintf("%.1fs", per_iter),
                     format_duration(remaining))
      if (is_terminal) {
        cat("\r", msg, sep = "")
      } else {
        cat(msg, "\n", sep = "")
      }
      flush.console()
      flush(stdout())
    }
  }

  # Main MCMC loop
  for (i in 2:iter) {
    if (verbose && (i == 2 || i == 10 || i %% progress_interval == 0 || i == iter)) {
      show_progress(i, iter, mcmc_start_time)
    }

    rho <- matrix(MCMC$rho[, i - 1], nrow = d, ncol = d)

    # Build precision matrix
    rhoxQ_result <- Q_rhoxQ(rho, MCMC$kappa[i - 1], priors$field.G)
    rhoxQ <- rhoxQ_result$rhoxQ

    # Add beta regularization if covariates present
    if (p > 0) {
      Q0 <- Matrix::bdiag(
        Matrix::Diagonal(p * d, 1 / sigma_beta),
        rhoxQ
      )
    } else {
      Q0 <- rhoxQ
    }

    # MALA step for alpha, beta, x
    step_result <- mala_step(
      old, MCMC$alpha[i - 1], y, w, Q0, AB, priors, d,
      MCMC$MALA.logstep[i - 1], i
    )

    new <- step_result$sample
    MCMC$MALA.count[i] <- step_result$count
    MCMC$MALA.logstep[i] <- step_result$logstep

    MCMC$alpha[i] <- new[length(new)]

    if (p > 0) {
      MCMC$beta[, i] <- new[1:(d * p)]
    }

    # Store x (with thinning)
    if (iter > 1e5) {
      if (i %% thin == 0) {
        MCMC$x[, i / thin] <- new[(d * p + 1):((N + p) * d)]
      }
    } else {
      MCMC$x[, i] <- new[(d * p + 1):((N + p) * d)]
    }

    old <- new[1:(length(new) - 1)]

    # Extract x for Gibbs step
    tmp <- new[(d * p + 1):((N + p) * d)]
    x_mat <- matrix(tmp, nrow = N, ncol = d)

    # Gibbs step for kappa and rho
    if (theta0$kappa != 0) {
      gibbs_result <- gibbs_kappa_rho(
        x_mat, MCMC$kappa[i - 1], rho, priors,
        MCMC$Gibbs.logstep[i - 1], i
      )

      MCMC$kappa[i] <- gibbs_result$sample$kappa
      MCMC$rho[, i] <- as.vector(gibbs_result$sample$rho)
      MCMC$Gibbs.count[i] <- gibbs_result$count
      MCMC$Gibbs.logstep[i] <- gibbs_result$logstep
    } else {
      # No spatial structure: sample rho directly
      MCMC$rho[, i] <- as.vector(sample_invwishart(x_mat, 0, priors))
      MCMC$kappa[i] <- 0
      MCMC$Gibbs.logstep[i] <- 0
    }

    # Sample horseshoe hyperparameters for beta
    if (p > 0) {
      # Auxiliary variables
      temp <- 1 + 1 / MCMC$l2[, i - 1]
      phi <- 1 / rgamma(p, 1, temp)
      temp <- 1 + 1 / MCMC$tau2[i - 1]
      xi <- 1 / rgamma(1, 1, temp)

      # Sample lambda
      beta_mat <- matrix(MCMC$beta[, i], nrow = d, ncol = p)
      tempa <- (1 + d) / 2
      tempb <- 1 / phi + colSums(beta_mat^2) / (2 * MCMC$tau2[i - 1])
      MCMC$l2[, i] <- 1 / rgamma(p, tempa, tempb)

      # Sample tau
      tempa <- (1 + p * d) / 2
      tempb <- 1 / xi + sum(colSums(beta_mat^2) / MCMC$l2[, i]) / 2
      MCMC$tau2[i] <- 1 / rgamma(1, tempa, tempb)

      sigma_beta <- rep(sqrt(MCMC$tau2[i] * MCMC$l2[, i]), d)
    }
  }

  if (verbose) {
    show_progress(iter, iter, mcmc_start_time)
    total_time <- proc.time()[3] - mcmc_start_time
    done_msg <- sprintf("MCMC complete: %d iterations in %s (%.2fs/iter)",
                        iter, format_duration(total_time),
                        total_time / (iter - 1))
    if (is_terminal) cat("\n")
    cat(done_msg, "\n")
    flush.console()
    flush(stdout())
  }

  # Compute acceptance rates
  MCMC$MALA.acc <- sum(MCMC$MALA.count) / iter
  MCMC$Gibbs.acc <- sum(MCMC$Gibbs.count) / iter

  # Organize output
  MCMC$MALA <- list(
    acc = MCMC$MALA.acc,
    count = MCMC$MALA.count,
    logstep = MCMC$MALA.logstep
  )
  MCMC$Gibbs <- list(
    acc = MCMC$Gibbs.acc,
    count = MCMC$Gibbs.count,
    logstep = MCMC$Gibbs.logstep
  )

  # Clean up
  MCMC$MALA.acc <- NULL
  MCMC$MALA.count <- NULL
  MCMC$MALA.logstep <- NULL
  MCMC$Gibbs.acc <- NULL
  MCMC$Gibbs.count <- NULL
  MCMC$Gibbs.logstep <- NULL

  return(MCMC)
}
