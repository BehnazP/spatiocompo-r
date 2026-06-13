#' Random Dirichlet Distribution
#'
#' Generates random samples from the Dirichlet distribution.
#'
#' @param n Number of samples to generate.
#' @param alpha Vector of concentration parameters (all positive).
#'
#' @return A matrix with n rows and length(alpha) columns, where each row
#'   is a sample from Dirichlet(alpha).
#'
#' @details
#' The Dirichlet distribution is sampled by:
#' 1. Generate x_k ~ Gamma(alpha_k, 1) for each k
#' 2. Normalize: y_k = x_k / sum(x)
#'
#' @examples
#' # Sample from Dirichlet(2, 3, 5)
#' samples <- spatiocompo:::rdirichlet(100, c(2, 3, 5))
#' colMeans(samples)  # Should be close to c(0.2, 0.3, 0.5)
#'
#' @keywords internal
rdirichlet <- function(n, alpha) {
  k <- length(alpha)
  x <- matrix(rgamma(n * k, rep(alpha, each = n), 1), nrow = n, ncol = k)
  x / rowSums(x)
}


#' Random Inverse Wishart Distribution
#'
#' Generates a random sample from the inverse Wishart distribution.
#'
#' @param nu Degrees of freedom (must be > dim - 1).
#' @param S Scale matrix (positive definite).
#'
#' @return A positive definite matrix drawn from IW(nu, S).
#'   The expected value is S / (nu - p - 1) where p is the dimension.
#'
#' @details
#' Uses the Bartlett decomposition: sample \eqn{W \sim Wishart(\nu, S^{-1})}, return \eqn{W^{-1}}.
#'
#' @examples
#' # Sample from IW(10, I_2)
#' # E[X] = I / (10 - 2 - 1) = I / 7
#' set.seed(123)
#' samples <- replicate(1000, spatiocompo:::rinvwishart(10, diag(2)), simplify = FALSE)
#' mean_diag <- mean(sapply(samples, function(x) x[1,1]))
#' # Should be close to 1/7 = 0.143
#'
#' @keywords internal
rinvwishart <- function(nu, S) {
  p <- nrow(S)

  if (nu <= p - 1) {
    stop("Degrees of freedom must be > p - 1")
  }

  # Upper Cholesky of S: L_S' L_S = S
  L_S <- chol(S)

  # Lower triangular T for Bartlett decomposition
  # T_ii ~ sqrt(chi^2(nu - i + 1)), T_ij ~ N(0,1) for i > j
  T <- matrix(0, p, p)
  T[lower.tri(T)] <- rnorm(p * (p - 1) / 2)
  diag(T) <- sqrt(rchisq(p, df = nu - (1:p) + 1))

  # For W ~ Wishart(nu, S^{-1}): W = L_V T T' L_V' where L_V = chol(S^{-1})
  # X = W^{-1} = B' B where B = T^{-1} L_S^T = forwardsolve(T, t(L_S))
  B <- forwardsolve(T, t(L_S))

  # Return X = B' B
  crossprod(B)
}
