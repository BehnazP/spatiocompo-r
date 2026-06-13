#' Summarize MCMC Results
#'
#' Computes posterior summaries from MCMC output.
#'
#' @param mcmc Output from [mcmc_sampling()].
#' @param burnin Number of initial samples to discard. Default is iter/4.
#' @param probs Quantiles to compute. Default is c(0.025, 0.5, 0.975).
#'
#' @return A list with summaries for each parameter:
#'   \item{alpha}{Named vector with mean, sd, and quantiles}
#'   \item{kappa}{Named vector with mean, sd, and quantiles}
#'   \item{rho}{Matrix of summaries for each element of rho}
#'   \item{x}{Matrix of posterior means for the latent field}
#'   \item{acceptance}{List with MALA and Gibbs acceptance rates}
#'
#' @examples
#' \dontrun{
#' fit <- mcmc_sampling(theta0, priors, y, iter = 5000)
#' summary_mcmc(fit)
#' }
#'
#' @export
summary_mcmc <- function(mcmc, burnin = NULL, probs = c(0.025, 0.5, 0.975)) {

  iter <- length(mcmc$alpha)

  if (is.null(burnin)) {
    burnin <- floor(iter / 4)
  }

  post_idx <- (burnin + 1):iter

  # Helper function for summaries
  summarize <- function(x) {
    c(mean = mean(x),
      sd = sd(x),
      quantile(x, probs = probs))
  }

  # Alpha
  alpha_summary <- summarize(mcmc$alpha[post_idx])

  # Kappa
  kappa_summary <- summarize(mcmc$kappa[post_idx])

  # Rho
  d <- sqrt(nrow(mcmc$rho))
  rho_samples <- mcmc$rho[, post_idx]
  rho_summary <- t(apply(rho_samples, 1, summarize))
  rownames(rho_summary) <- paste0("rho[", rep(1:d, d), ",", rep(1:d, each = d), "]")

  # X (posterior means)
  x_samples <- mcmc$x
  if (ncol(mcmc$x) > length(post_idx)) {
    # Thinned samples
    x_mean <- rowMeans(x_samples)
  } else {
    # Compute which x samples correspond to post_idx
    # This depends on thinning
    thin_idx <- ceiling(post_idx / 50)
    thin_idx <- thin_idx[thin_idx <= ncol(x_samples)]
    x_mean <- rowMeans(x_samples[, unique(thin_idx), drop = FALSE])
  }

  result <- list(
    alpha = alpha_summary,
    kappa = kappa_summary,
    rho = rho_summary,
    x_mean = x_mean,
    acceptance = list(
      MALA = mcmc$MALA$acc,
      Gibbs = mcmc$Gibbs$acc
    ),
    n_iter = iter,
    n_post = length(post_idx)
  )

  result
}


#' Effective Sample Size
#'
#' Computes the effective sample size for MCMC samples, accounting for
#' autocorrelation.
#'
#' @param x Vector of MCMC samples.
#'
#' @return Effective sample size (positive number <= length(x)).
#'
#' @details
#' ESS is computed using the formula:
#' \deqn{ESS = n / (1 + 2 \sum_{k=1}^{K} \rho_k)}
#' where rho_k is the autocorrelation at lag k, summed until the
#' autocorrelations become negligible.
#'
#' Low ESS (< 100-400) suggests poor mixing and the need for longer chains
#' or improved sampling.
#'
#' @examples
#' # Uncorrelated samples
#' x <- rnorm(1000)
#' ess(x)  # Close to 1000
#'
#' # Correlated samples
#' x_corr <- filter(rnorm(1000), 0.9, method = "recursive")
#' ess(x_corr)  # Much less than 1000
#'
#' @export
ess <- function(x) {
  n <- length(x)
  if (n < 10) return(n)

  # Compute autocorrelations
  acf_vals <- acf(x, lag.max = min(n - 1, 500), plot = FALSE)$acf[-1]

  # Find where to truncate (first negative or near-zero)
  truncate_at <- which(acf_vals < 0.05)[1]
  if (is.na(truncate_at)) truncate_at <- length(acf_vals)

  # ESS formula
  rho_sum <- sum(acf_vals[1:truncate_at])
  ess_val <- n / (1 + 2 * rho_sum)

  return(max(1, ess_val))
}


#' Compute Posterior Predictive Samples
#'
#' Generates posterior predictive samples of compositional observations.
#'
#' @param mcmc Output from [mcmc_sampling()].
#' @param n_samples Number of posterior predictive samples. Default is 100.
#' @param burnin Burnin to discard. Default is iter/4.
#' @param A Location matrix (if different from identity).
#'
#' @return An array of dimension (N x D x n_samples) containing posterior
#'   predictive compositions.
#'
#' @details
#' For each posterior sample of (x, alpha), generates a new observation
#' from the Dirichlet likelihood. Useful for model checking by comparing
#' predictive distribution to observed data.
#'
#' @examples
#' \dontrun{
#' fit <- mcmc_sampling(theta0, priors, y, iter = 5000)
#' y_pred <- posterior_predictive(fit, n_samples = 50)
#'
#' # Compare observed vs predicted
#' hist(y[, 1], main = "Component 1", freq = FALSE)
#' lines(density(y_pred[, 1, ]), col = "red")
#' }
#'
#' @export
posterior_predictive <- function(mcmc, n_samples = 100, burnin = NULL, A = NULL) {

  iter <- length(mcmc$alpha)

  if (is.null(burnin)) {
    burnin <- floor(iter / 4)
  }

  post_idx <- (burnin + 1):iter
  sample_idx <- sample(post_idx, n_samples, replace = TRUE)

  # Dimensions
  N_d <- nrow(mcmc$x)
  d <- sqrt(nrow(mcmc$rho))
  D <- d + 1
  N <- N_d / d

  if (is.null(A)) {
    A <- Matrix::Diagonal(N * d)
  }

  n <- nrow(A) / d

  # Storage
  y_pred <- array(0, dim = c(n, D, n_samples))

  for (s in seq_len(n_samples)) {
    idx <- sample_idx[s]

    # Get x sample (handle thinning)
    x_idx <- ceiling(idx / 50)
    x_idx <- min(x_idx, ncol(mcmc$x))
    x <- mcmc$x[, x_idx]

    alpha <- mcmc$alpha[idx]

    # Transform to compositions
    Ax <- matrix(as.vector(A %*% x), nrow = n, ncol = d)
    z <- invalr(Ax)

    # Sample from Dirichlet likelihood
    for (i in seq_len(n)) {
      concentration <- alpha * z[i, ]
      y_pred[i, , s] <- rdirichlet(1, concentration)
    }
  }

  return(y_pred)
}


#' Check MCMC Convergence
#'
#' Performs basic convergence diagnostics on MCMC output.
#'
#' @param mcmc Output from [mcmc_sampling()].
#' @param split_chains If TRUE, split chain in half for Gelman-Rubin-like
#'   diagnostic. Default is TRUE.
#'
#' @return A data frame with columns:
#'   \item{parameter}{Parameter name}
#'   \item{ess}{Effective sample size}
#'   \item{rhat}{Potential scale reduction (should be < 1.1)}
#'
#' @details
#' Diagnostics include:
#' \itemize{
#'   \item ESS: Effective sample size (should be > 400)
#'   \item R-hat: Potential scale reduction comparing chain halves
#'     (should be < 1.1)
#' }
#'
#' @examples
#' \dontrun{
#' fit <- mcmc_sampling(theta0, priors, y, iter = 5000)
#' check_convergence(fit)
#' }
#'
#' @export
check_convergence <- function(mcmc, split_chains = TRUE) {

  iter <- length(mcmc$alpha)
  burnin <- floor(iter / 4)
  post_idx <- (burnin + 1):iter

  results <- data.frame(
    parameter = character(),
    ess = numeric(),
    rhat = numeric(),
    stringsAsFactors = FALSE
  )

  # Helper for split-chain R-hat
  compute_rhat <- function(x) {
    n <- length(x)
    mid <- n %/% 2
    chain1 <- x[1:mid]
    chain2 <- x[(mid + 1):n]

    # Between-chain variance
    m1 <- mean(chain1)
    m2 <- mean(chain2)
    B <- mid * var(c(m1, m2))

    # Within-chain variance
    W <- (var(chain1) + var(chain2)) / 2

    # R-hat
    sqrt((W * (mid - 1) / mid + B / mid) / W)
  }

  # Alpha
  alpha_post <- mcmc$alpha[post_idx]
  results <- rbind(results, data.frame(
    parameter = "alpha",
    ess = ess(alpha_post),
    rhat = if (split_chains) compute_rhat(alpha_post) else NA
  ))

  # Kappa
  kappa_post <- mcmc$kappa[post_idx]
  results <- rbind(results, data.frame(
    parameter = "kappa",
    ess = ess(kappa_post),
    rhat = if (split_chains) compute_rhat(kappa_post) else NA
  ))

  # Rho elements
  d <- sqrt(nrow(mcmc$rho))
  for (i in 1:d) {
    for (j in 1:d) {
      idx <- (i - 1) * d + j
      rho_ij <- mcmc$rho[idx, post_idx]
      results <- rbind(results, data.frame(
        parameter = paste0("rho[", i, ",", j, "]"),
        ess = ess(rho_ij),
        rhat = if (split_chains) compute_rhat(rho_ij) else NA
      ))
    }
  }

  return(results)
}
