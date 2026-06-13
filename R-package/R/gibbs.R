#' Gibbs Sampler for Spatial Parameters
#'
#' Performs a Metropolis-within-Gibbs step to sample the spatial range
#' parameter kappa and covariance matrix rho given the latent field x.
#'
#' @param x Matrix (N x d) of the current latent field values.
#' @param kappa Current value of the spatial range parameter.
#' @param rho Current d x d covariance matrix.
#' @param priors List containing prior hyperparameters:
#'   \itemize{
#'     \item kappa.range, kappa.sigma: Gamma prior shape and rate for kappa
#'     \item rho.df, rho.Sigma: Inverse-Wishart prior degrees of freedom and scale
#'     \item field.G: Base precision matrix from create_Q()
#'     \item field.lambda: Eigenvalues for log-determinant computation
#'   }
#' @param logstep0 Current log step size for adaptation.
#' @param iter Current MCMC iteration number.
#'
#' @return A list with components:
#'   \item{sample}{List with $kappa and $rho containing proposed values}
#'   \item{count}{1 if proposal accepted, 0 if rejected}
#'   \item{logstep}{Updated log step size}
#'
#' @details
#' The algorithm uses the fact that conditioned on x and kappa, the posterior
#' of rho is an inverse-Wishart distribution (conjugate prior). Thus:
#' 1. Sample kappa using a log-normal random walk proposal
#' 2. If accepted, sample rho from its conditional inverse-Wishart
#' 3. If rejected, keep the old values
#'
#' The step size is adapted to target an acceptance rate of 0.44.
#'
#' @seealso [sample_invwishart()], [log_posterior_kappa_x()], [mcmc_sampling()]
#'
#' @export
gibbs_kappa_rho <- function(x, kappa, rho, priors, logstep0, iter) {

  step <- exp(logstep0 / 2)

  if (kappa <= 0) {
    acc_prob <- 0
  } else {
    old <- kappa

    # Log-normal random walk proposal for kappa
    prop <- log(old) + step * rnorm(1)
    new <- exp(prop)

    # Log-posterior at old and new kappa
    logpxy <- -log_posterior_kappa_x(x, old, priors)
    logpyx <- -log_posterior_kappa_x(x, new, priors)

    # Jacobian for log-transform: proposal is symmetric in log-space
    # but need to account for transformation back
    acc_prob <- exp(min(0, logpyx - logpxy + log(new) - log(old)))
  }

  U <- runif(1)
  if (U < acc_prob) {
    sample <- list(
      kappa = new,
      rho = sample_invwishart(x, new, priors)
    )
    count <- 1
  } else {
    sample <- list(
      kappa = kappa,
      rho = rho
    )
    count <- 0
  }

  # Adapt step size (target acceptance rate 0.44)
  logstep <- logstep0 + (iter + 1)^(-0.5) * (acc_prob - 0.44)

  return(list(sample = sample, count = count, logstep = logstep))
}


#' Sample from Inverse-Wishart Conditional Posterior
#'
#' Samples the covariance matrix rho from its conditional inverse-Wishart
#' posterior given the latent field x and spatial range kappa.
#'
#' @param x Matrix (N x d) of latent field values.
#' @param kappa Spatial range parameter.
#' @param priors List containing prior hyperparameters:
#'   \itemize{
#'     \item field.G: Base precision matrix
#'     \item rho.df: Prior degrees of freedom
#'     \item rho.Sigma: Prior scale matrix
#'   }
#'
#' @return A d x d positive definite covariance matrix.
#'
#' @details
#' With the conjugate inverse-Wishart prior rho ~ IW(df, Sigma), the
#' conditional posterior given x and kappa is:
#' \deqn{rho | x, kappa \sim IW(df + N, Sigma + x^T Q x)}
#' where Q is the marginal precision matrix for one component.
#'
#' @seealso [gibbs_kappa_rho()]
#'
#' @export
sample_invwishart <- function(x, kappa, priors) {
  N <- nrow(x)

  # Get marginal Q
  Q <- Q_rhoxQ(NULL, kappa, priors$field.G)$Q

  # Posterior parameters
  post_df <- priors$rho.df + N
  post_scale <- priors$rho.Sigma + t(x) %*% Q %*% x

  # Sample from inverse-Wishart
  rho <- rinvwishart(nu = post_df, S = as.matrix(post_scale))

  return(rho)
}


#' Marginal Log-Posterior of Kappa Given X
#'
#' Computes the negative log-posterior of the spatial range parameter kappa
#' marginalized over the covariance matrix rho.
#'
#' @param x Matrix (N x d) of latent field values.
#' @param kappa Spatial range parameter (positive).
#' @param priors List containing prior hyperparameters:
#'   \itemize{
#'     \item kappa.range, kappa.sigma: Gamma prior shape and rate
#'     \item rho.df, rho.Sigma: Inverse-Wishart prior parameters
#'     \item field.G: Base precision matrix
#'     \item field.lambda: Eigenvalues for log-determinant
#'   }
#'
#' @return Negative log-posterior value (for minimization).
#'
#' @details
#' The marginal posterior of kappa given x is obtained by integrating out rho:
#' \deqn{p(kappa | x) \propto |Q|^{d/2} p(kappa) / |Sigma + x^T Q x|^{(N+df)/2}}
#'
#' The log-determinant of Q is computed efficiently using precomputed
#' eigenvalues via [log_det_Q()].
#'
#' @seealso [gibbs_kappa_rho()], [log_det_Q()]
#'
#' @export
log_posterior_kappa_x <- function(x, kappa, priors) {
  N <- nrow(x)
  d <- ncol(x)

  # Marginal precision Q
  Q <- Q_rhoxQ(NULL, kappa, priors$field.G)$Q

  # Log-determinant of Q using eigenvalues
  logdetQ <- log_det_Q(kappa, priors$field.lambda)

  # Posterior scale for rho
  IxtQx <- priors$rho.Sigma + as.matrix(t(x) %*% Q %*% x)

  # Log-determinant via Cholesky
  R <- chol(IxtQx)
  logdet <- 2 * sum(log(diag(R)))

  # Gamma prior on kappa: shape = range, rate = sigma
  logkappa <- log(kappa)
  pkappa <- (priors$kappa.range - 1) * logkappa - priors$kappa.sigma * kappa

  # Negative log-posterior
  logp <- -(0.5 * d * logdetQ - 0.5 * (N + priors$rho.df) * logdet + pkappa)

  return(logp)
}


#' Gradient of Log-Posterior of Kappa (w.r.t. log kappa)
#'
#' Computes the gradient of the marginal log-posterior of kappa with respect
#' to log(kappa), for use in gradient-based kappa samplers.
#'
#' @param x Matrix (N x d) of latent field values.
#' @param kappa Spatial range parameter (positive).
#' @param priors List containing prior hyperparameters.
#'
#' @return Gradient of log p(kappa | x) with respect to log(kappa).
#'
#' @keywords internal
grad_log_posterior_kappa <- function(x, kappa, priors) {
  N <- nrow(x)
  d <- ncol(x)
  lambda <- priors$field.lambda
  G <- priors$field.G

  # d/d(kappa) log|Q| = sum 2*(2*kappa^3 + kappa*lambda_i) / (kappa^4 + 2*kappa^2*lambda_i + lambda_i^2)
  xi <- kappa^4 + 2 * kappa^2 * lambda + lambda^2
  dlogdetQ_dkappa <- sum(2 * (2 * kappa^3 + kappa * lambda) / xi)  # * kappa for chain rule below

  # d/d(kappa) x^T Q(kappa) x
  # Q = kappa^4 I + 2*kappa^2 G + G^T G
  # dQ/dkappa = 4*kappa^3 I + 4*kappa G
  dQ_x <- 4 * kappa^3 * x + 4 * kappa * as.matrix(G %*% x)  # N x d
  xtdQx <- as.matrix(t(x) %*% dQ_x)  # d x d

  # Posterior scale
  Q <- Q_rhoxQ(NULL, kappa, priors$field.G)$Q
  IxtQx <- priors$rho.Sigma + as.matrix(t(x) %*% Q %*% x)
  IxtQx_inv <- solve(IxtQx)

  # d/d(kappa) log|Sigma + x^T Q x| = tr(IxtQx^{-1} * x^T dQ/dk x)
  dlogdet_dkappa <- sum(IxtQx_inv * xtdQx)  # tr(A*B) = sum(A*B)

  # Gamma prior: d/d(kappa) [(a-1)*log(kappa) - b*kappa]
  dpkappa_dkappa <- (priors$kappa.range - 1) / kappa - priors$kappa.sigma

  # Gradient of log-posterior w.r.t. kappa
  grad_kappa <- 0.5 * d * dlogdetQ_dkappa -
    0.5 * (N + priors$rho.df) * dlogdet_dkappa +
    dpkappa_dkappa

  # Chain rule: d/d(log kappa) = kappa * d/d(kappa)
  grad_kappa * kappa
}

