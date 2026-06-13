#' Metropolis-Adjusted Langevin Algorithm Step
#'
#' Performs one step of the Metropolis-Adjusted Langevin Algorithm (MALA)
#' to propose new values for the latent field x and precision parameter alpha.
#'
#' @param x Current value of the ALR-transformed latent field (vector).
#' @param alpha Current value of the Dirichlet precision parameter.
#' @param y Observed compositional data (n x D matrix).
#' @param w Observation weights (scalar or vector).
#' @param Q GMRF precision matrix for the latent field.
#' @param A Location matrix connecting observations to grid cells.
#' @param priors List containing prior hyperparameters:
#'   \itemize{
#'     \item alpha.a, alpha.b: Gamma prior parameters for alpha
#'   }
#' @param d Number of ALR components (D-1).
#' @param logstep0 Current log step size for adaptation.
#' @param iter Current MCMC iteration number.
#'
#' @return A list with components:
#'   \item{sample}{Proposed values c(x, alpha)}
#'   \item{count}{1 if proposal accepted, 0 if rejected}
#'   \item{logstep}{Updated log step size}
#'
#' @details
#' MALA uses gradient information to make proposals that are more likely
#' to be in high-probability regions. The algorithm:
#' 1. Computes the Fisher information matrix
#' 2. Proposes new values using a Langevin diffusion step
#' 3. Accepts/rejects using Metropolis-Hastings ratio
#'
#' The step size is adapted to target an acceptance rate of 0.56.
#'
#' @references
#' Girolami, M. and Calderhead, B. (2011). Riemann manifold Langevin and
#' Hamiltonian Monte Carlo methods. J. R. Statist. Soc. B, 73(2):123-214.
#'
#' Pirzamanbein, B., Lindstrom, J., Poska, A., Gaillard, M.-J. (2018).
#' Modelling spatial compositional data: Reconstructions of past land
#' cover and uncertainties. \emph{Spatial Statistics} 24: 14--31.
#' \doi{10.1016/j.spasta.2018.03.005}
#'
#' @seealso [mcmc_sampling()], [fisher_info()]
#'
#' @export
mala_step <- function(x, alpha, y, w, Q, A, priors, d, logstep0, iter) {

  step <- exp(logstep0 / 2)
  old <- c(x, alpha)

  # Compute Fisher information and proposal mean
  fi_result <- fisher_info(x, alpha, A, y, priors$alpha.a, priors$alpha.b,
                           Q, step, d, w)
  RFI0 <- fi_result$RFI
  m0 <- fi_result$m
  p0 <- fi_result$p

  # Generate proposal
  epsilon <- rnorm(nrow(RFI0))
  new_perm <- m0 + backsolve(RFI0, step * epsilon)

  # Unpermute
  new <- numeric(length(new_perm))
  new[p0] <- new_perm
  alphanew <- new[length(new)]

  # Check validity

if (alphanew <= 0) {
    acc_prob <- 0
  } else {
    # Log-posterior at old and new points
    nll_old <- neg_log_lik(old[1:(length(old) - 1)], old[length(old)],
                           A, y, priors$alpha.a, priors$alpha.b, Q, d, w)
    logpxy <- -nll_old$L

    nll_new <- neg_log_lik(new[1:(length(new) - 1)], new[length(new)],
                           A, y, priors$alpha.a, priors$alpha.b, Q, d, w)
    logpyx <- -nll_new$L

    if (is.finite(logpyx)) {
      # Compute Fisher info at new point
      fi_new <- fisher_info(new[1:(length(new) - 1)], new[length(new)],
                            A, y, priors$alpha.a, priors$alpha.b, Q, step, d, w)
      RFI <- fi_new$RFI
      m <- fi_new$m
      p <- fi_new$p

      # Permute old for comparison
      old_perm <- old[p]

      # Proposal density q(old|new)
      # log q = log|RFI| - 0.5 * ||RFI/step * (old - m)||^2
      logqxy <- sum(log(diag(RFI))) -
        sum(((RFI / step) %*% (old_perm - m))^2) / 2

      # Proposal density q(new|old)
      # = log|RFI0| - 0.5 * ||epsilon||^2
      logqyx <- sum(log(diag(RFI0))) - sum(epsilon^2) / 2

      # Acceptance probability
      # MH ratio: p(new) * q(old|new) / (p(old) * q(new|old))
      log_acc <- min(0, logpyx - logpxy + logqxy - logqyx)
      acc_prob <- exp(log_acc)
    } else {
      acc_prob <- 0
    }
  }

  # Accept or reject
  U <- runif(1)
  if (U < acc_prob) {
    sample <- new
    count <- 1
  } else {
    sample <- c(x, alpha)
    count <- 0
  }

  # Adapt step size (target acceptance rate 0.56)
  logstep <- logstep0 + (iter + 1)^(-0.5) * (acc_prob - 0.56)

  return(list(sample = sample, count = count, logstep = logstep))
}
