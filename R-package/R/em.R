#' EM (Block Coordinate Ascent) Estimation for Spatial Compositional Data
#'
#' Finds MAP point estimates for the spatial compositional model using
#' block coordinate ascent. Much faster than MCMC for obtaining point
#' estimates, reusing existing gradient and likelihood functions.
#'
#' @param theta0 List of initial values (same format as [mcmc_sampling()]):
#'   \itemize{
#'     \item alpha: Dirichlet precision parameter (positive scalar)
#'     \item kappa: Spatial range parameter (positive scalar, or 0 for no spatial)
#'     \item rho: Covariance matrix (d x d positive definite)
#'     \item x: Latent field (N x d matrix)
#'     \item beta: Regression coefficients (p x d matrix, optional)
#'   }
#' @param priors List of prior hyperparameters from [create_priors()].
#' @param y Observed compositional data (n x D matrix, rows sum to 1).
#' @param w Observation weights (scalar or vector of length n). Default is 1.
#' @param A Location matrix (n*d x N*d sparse matrix). Default is identity.
#' @param B Covariate matrix (n x p) for regression effects. Default is NULL.
#' @param max_iter Maximum number of EM iterations. Default is 100.
#' @param tol Convergence tolerance on relative change in neg log-posterior.
#'   Default is 1e-6.
#' @param verbose Logical; if TRUE, show progress. Default is TRUE.
#' @param lbfgs_maxit Maximum iterations for L-BFGS-B in each EM step.
#'   Default is 200.
#' @param compute_se Logical; if TRUE, compute standard errors at convergence.
#'   Default is TRUE.
#'
#' @return A list with components:
#'   \item{alpha}{Estimated Dirichlet precision}
#'   \item{kappa}{Estimated spatial range}
#'   \item{rho}{Estimated d x d covariance matrix}
#'   \item{x}{Estimated latent field (N x d matrix)}
#'   \item{beta}{Estimated regression coefficients (p x d matrix, if B provided)}
#'   \item{convergence}{List with neg_log_post trace, converged flag, n_iter}
#'   \item{se}{List of standard errors (if compute_se=TRUE)}
#'
#' @details
#' Each EM iteration cycles through:
#' \enumerate{
#'   \item Optimize (x, beta, alpha) jointly via L-BFGS-B
#'   \item Optimize kappa via 1D bounded optimization
#'   \item Set rho = Inverse-Wishart mode (closed form)
#'   \item Update horseshoe hyperparameters (tau2, l2) at conditional modes
#' }
#'
#' Convergence is checked via relative change in the negative log-posterior.
#'
#' @seealso [mcmc_sampling()], [find_initial_values()]
#'
#' @examples
#' \dontrun{
#' sim <- simulate_compo_data(n_side = 10, kappa_true = 0.5,
#'                            alpha_true = 50, seed = 123)
#' priors <- create_priors(sim$G, d = 2)
#' theta0 <- find_initial_values(sim$y, priors)
#' fit <- em_estimation(theta0, priors, sim$y)
#' }
#'
#' @export
em_estimation <- function(theta0, priors, y, w = 1, A = NULL, B = NULL,
                          max_iter = 100, tol = 1e-6, verbose = TRUE,
                          lbfgs_maxit = 200, compute_se = TRUE) {

  # --- Dimensions (mirrors mcmc_sampling lines 92-114) ---
  if (is.vector(y)) y <- matrix(y, nrow = 1)
  n <- nrow(y)
  D <- ncol(y)
  d <- D - 1

  N <- nrow(theta0$x)

  if (is.null(A)) {
    A <- Matrix::Diagonal(N * d)
  }

  if (!is.null(B)) {
    p <- ncol(B) / d
    AB <- cbind(B, A)
    state_xb <- c(as.vector(theta0$beta), as.vector(theta0$x))
  } else {
    p <- 0
    AB <- A
    state_xb <- as.vector(theta0$x)
  }

  if (length(w) == 1) w <- rep(w, n)

  # --- Initialize state ---
  alpha <- theta0$alpha
  kappa <- theta0$kappa
  rho <- theta0$rho
  if (is.vector(rho)) rho <- matrix(rho, nrow = d, ncol = d)

  # Horseshoe hyperparameters
  if (p > 0) {
    tau2 <- 1
    l2 <- rep(1, p)
    sigma_beta <- rep(1, p * d)
  }

  # --- Convergence tracking ---
  neg_log_post_trace <- numeric(max_iter)
  converged <- FALSE

  if (verbose) {
    em_start_time <- proc.time()[3]
    cat("EM estimation started\n")
  }

  # --- Main EM loop ---
  for (iter in seq_len(max_iter)) {

    # Build precision matrix Q0
    rhoxQ_result <- Q_rhoxQ(rho, kappa, priors$field.G)
    rhoxQ <- rhoxQ_result$rhoxQ

    if (p > 0) {
      Q0 <- Matrix::bdiag(
        Matrix::Diagonal(p * d, 1 / sigma_beta),
        rhoxQ
      )
    } else {
      Q0 <- rhoxQ
    }

    # === Step 1: Optimize (x, beta, alpha) via L-BFGS-B ===
    state <- c(state_xb, alpha)

    # Objective: returns both value and gradient (cache to avoid double computation)
    cached_result <- NULL
    cached_state <- NULL

    obj_fn <- function(s) {
      x_vec <- s[1:(length(s) - 1)]
      a <- s[length(s)]
      if (a <= 0) return(1e10)
      result <- neg_log_lik(x_vec, a, AB, y, priors$alpha.a, priors$alpha.b,
                            Q0, d, w)
      # Cache for gradient
      cached_result <<- result
      cached_state <<- s
      return(result$L)
    }

    grad_fn <- function(s) {
      x_vec <- s[1:(length(s) - 1)]
      a <- s[length(s)]
      if (a <= 0) return(rep(0, length(s)))

      # Use cached result if state matches
      if (!is.null(cached_state) && identical(s, cached_state)) {
        return(cached_result$dL)
      }

      result <- neg_log_lik(x_vec, a, AB, y, priors$alpha.a, priors$alpha.b,
                            Q0, d, w)
      cached_result <<- result
      cached_state <<- s
      return(result$dL)
    }

    opt_result <- tryCatch({
      optim(
        par = state,
        fn = obj_fn,
        gr = grad_fn,
        method = "L-BFGS-B",
        lower = c(rep(-Inf, length(state) - 1), 0.1),
        upper = c(rep(Inf, length(state) - 1), 1000),
        control = list(maxit = lbfgs_maxit, factr = 1e7)
      )
    }, error = function(e) {
      if (verbose) cat("  L-BFGS-B failed at iter", iter, ":", e$message, "\n")
      list(par = state, value = 1e10, convergence = 1)
    })

    state_opt <- opt_result$par
    alpha <- state_opt[length(state_opt)]
    state_xb <- state_opt[1:(length(state_opt) - 1)]

    # Extract x matrix for subsequent steps
    if (p > 0) {
      x_vec <- state_xb[(p * d + 1):((N + p) * d)]
    } else {
      x_vec <- state_xb
    }
    x_mat <- matrix(x_vec, nrow = N, ncol = d)

    # === Step 2: Optimize kappa ===
    if (theta0$kappa != 0) {
      kappa_obj <- function(k) {
        if (k <= 0) return(1e10)
        log_posterior_kappa_x(x_mat, k, priors)
      }

      kappa_result <- tryCatch({
        optimize(kappa_obj, interval = c(0.01, 5), tol = 0.01)
      }, error = function(e) {
        list(minimum = kappa)
      })
      kappa <- kappa_result$minimum
    }

    # === Step 3: Set rho = IW mode (closed form) ===
    Q_marg <- Q_rhoxQ(NULL, kappa, priors$field.G)$Q
    post_scale <- priors$rho.Sigma + as.matrix(t(x_mat) %*% Q_marg %*% x_mat)
    # IW mode: (Sigma_post) / (df_post + d + 1)
    post_df <- priors$rho.df + N
    rho <- post_scale / (post_df + d + 1)

    # === Step 4: Update horseshoe hyperparameters at conditional modes ===
    if (p > 0) {
      beta_vec <- state_xb[1:(p * d)]
      beta_mat <- matrix(beta_vec, nrow = d, ncol = p)

      # l2 mode: InvGamma((1 + d)/2, 1/phi + colSums(beta^2)/(2*tau2))
      # Mode of InvGamma(a, b) = b / (a + 1)
      # phi mode: InvGamma(1, 1 + 1/l2) -> mode = (1 + 1/l2) / 2
      phi <- (1 + 1 / l2) / 2
      a_l2 <- (1 + d) / 2
      b_l2 <- 1 / phi + colSums(beta_mat^2) / (2 * tau2)
      l2 <- b_l2 / (a_l2 + 1)

      # tau2 mode: InvGamma((1 + p*d)/2, 1/xi + sum(colSums(beta^2)/l2)/2)
      # xi mode: InvGamma(1, 1 + 1/tau2) -> mode = (1 + 1/tau2) / 2
      xi <- (1 + 1 / tau2) / 2
      a_tau2 <- (1 + p * d) / 2
      b_tau2 <- 1 / xi + sum(colSums(beta_mat^2) / l2) / 2
      tau2 <- b_tau2 / (a_tau2 + 1)

      sigma_beta <- rep(sqrt(tau2 * l2), d)
    }

    # === Compute neg log-posterior for convergence check ===
    # Rebuild Q0 with updated rho and kappa for accurate evaluation
    rhoxQ_new <- Q_rhoxQ(rho, kappa, priors$field.G)$rhoxQ
    if (p > 0) {
      Q0_new <- Matrix::bdiag(
        Matrix::Diagonal(p * d, 1 / sigma_beta),
        rhoxQ_new
      )
    } else {
      Q0_new <- rhoxQ_new
    }

    nll_result <- neg_log_lik(state_xb, alpha, AB, y,
                               priors$alpha.a, priors$alpha.b,
                               Q0_new, d, w)
    # Add kappa prior and rho prior contributions
    logkappa <- log(kappa)
    kappa_prior <- (priors$kappa.range - 1) * logkappa - priors$kappa.sigma * kappa

    logdetQ <- log_det_Q(kappa, priors$field.lambda)
    IxtQx <- priors$rho.Sigma + as.matrix(t(x_mat) %*% Q_marg %*% x_mat)
    R <- chol(IxtQx)
    logdet_IxtQx <- 2 * sum(log(diag(R)))
    rho_kappa_term <- 0.5 * d * logdetQ - 0.5 * (N + priors$rho.df) * logdet_IxtQx

    neg_log_post <- nll_result$L - kappa_prior - rho_kappa_term
    neg_log_post_trace[iter] <- neg_log_post

    # === Convergence check ===
    if (iter > 1) {
      rel_change <- abs(neg_log_post_trace[iter] - neg_log_post_trace[iter - 1]) /
        (abs(neg_log_post_trace[iter - 1]) + 1e-10)

      if (verbose && (iter <= 5 || iter %% 10 == 0)) {
        elapsed <- proc.time()[3] - em_start_time
        cat(sprintf("  EM iter %3d | neg_log_post = %.2f | rel_change = %.2e | %.1fs\n",
                    iter, neg_log_post, rel_change, elapsed))
      }

      if (rel_change < tol) {
        converged <- TRUE
        if (verbose) {
          elapsed <- proc.time()[3] - em_start_time
          cat(sprintf("  Converged at iter %d (rel_change = %.2e, %.1fs)\n",
                      iter, rel_change, elapsed))
        }
        break
      }
    } else {
      if (verbose) {
        cat(sprintf("  EM iter %3d | neg_log_post = %.2f\n", iter, neg_log_post))
      }
    }
  }

  if (!converged && verbose) {
    elapsed <- proc.time()[3] - em_start_time
    cat(sprintf("  EM did not converge after %d iterations (%.1fs)\n",
                max_iter, elapsed))
  }

  # --- Build result ---
  result <- list(
    alpha = alpha,
    kappa = kappa,
    rho = rho,
    x = x_mat,
    convergence = list(
      neg_log_post = neg_log_post_trace[1:min(iter, max_iter)],
      converged = converged,
      n_iter = min(iter, max_iter)
    )
  )

  if (p > 0) {
    result$beta <- matrix(state_xb[1:(p * d)], nrow = p, ncol = d, byrow = FALSE)
  }

  # --- Standard errors ---
  if (compute_se) {
    if (verbose) cat("Computing standard errors...\n")
    result$se <- compute_em_se(result, priors, y, w, AB, Q0_new, d, p, N)
    if (verbose) cat("  Done.\n")
  }

  if (verbose) {
    total_time <- proc.time()[3] - em_start_time
    cat(sprintf("EM complete: %d iterations in %.1fs\n",
                result$convergence$n_iter, total_time))
  }

  return(result)
}


#' Compute Standard Errors at EM Convergence
#'
#' Computes standard errors via numerical second derivatives at the
#' converged MAP estimates.
#'
#' @param fit Converged EM result list.
#' @param priors Prior specification.
#' @param y Observed data.
#' @param w Observation weights.
#' @param AB Combined covariate + observation matrix.
#' @param Q0 Current precision matrix.
#' @param d Number of ALR components.
#' @param p Number of covariate groups.
#' @param N Number of grid cells.
#'
#' @return A list with standard errors for alpha, kappa, and beta.
#'
#' @keywords internal
compute_em_se <- function(fit, priors, y, w, AB, Q0, d, p, N) {
  se <- list()

  # SE for alpha: 1D numerical second derivative
  alpha_hat <- fit$alpha
  state_xb <- if (p > 0) {
    c(as.vector(fit$beta), as.vector(fit$x))
  } else {
    as.vector(fit$x)
  }

  h_alpha <- max(1e-4, abs(alpha_hat) * 1e-4)
  f_plus <- neg_log_lik(state_xb, alpha_hat + h_alpha, AB, y,
                         priors$alpha.a, priors$alpha.b, Q0, d, w)$L
  f_minus <- neg_log_lik(state_xb, alpha_hat - h_alpha, AB, y,
                          priors$alpha.a, priors$alpha.b, Q0, d, w)$L
  f_center <- neg_log_lik(state_xb, alpha_hat, AB, y,
                           priors$alpha.a, priors$alpha.b, Q0, d, w)$L
  d2f_alpha <- (f_plus - 2 * f_center + f_minus) / h_alpha^2
  se$alpha <- if (d2f_alpha > 0) 1 / sqrt(d2f_alpha) else NA

  # SE for kappa: 1D numerical second derivative on marginal posterior
  kappa_hat <- fit$kappa
  x_mat <- fit$x
  if (kappa_hat > 0) {
    h_kappa <- max(1e-4, kappa_hat * 1e-4)
    f_plus_k <- log_posterior_kappa_x(x_mat, kappa_hat + h_kappa, priors)
    f_minus_k <- log_posterior_kappa_x(x_mat, kappa_hat - h_kappa, priors)
    f_center_k <- log_posterior_kappa_x(x_mat, kappa_hat, priors)
    d2f_kappa <- (f_plus_k - 2 * f_center_k + f_minus_k) / h_kappa^2
    se$kappa <- if (d2f_kappa > 0) 1 / sqrt(d2f_kappa) else NA
  } else {
    se$kappa <- NA
  }

  # SE for beta: numerical Hessian of the p*d block
  if (p > 0) {
    beta_vec <- as.vector(fit$beta)
    n_beta <- p * d
    se_beta <- numeric(n_beta)

    for (j in seq_len(n_beta)) {
      h_j <- max(1e-4, abs(beta_vec[j]) * 1e-4)
      state_plus <- state_xb
      state_minus <- state_xb
      state_plus[j] <- state_xb[j] + h_j
      state_minus[j] <- state_xb[j] - h_j

      f_plus_b <- neg_log_lik(state_plus, alpha_hat, AB, y,
                               priors$alpha.a, priors$alpha.b, Q0, d, w)$L
      f_minus_b <- neg_log_lik(state_minus, alpha_hat, AB, y,
                                priors$alpha.a, priors$alpha.b, Q0, d, w)$L
      d2f_j <- (f_plus_b - 2 * f_center + f_minus_b) / h_j^2
      se_beta[j] <- if (d2f_j > 0) 1 / sqrt(d2f_j) else NA
    }

    se$beta <- matrix(se_beta, nrow = p, ncol = d, byrow = FALSE)
  }

  return(se)
}
