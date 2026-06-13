#' Gradient of Dirichlet Log-Likelihood
#'
#' Computes the gradient of the Dirichlet log-likelihood with respect to
#' the ALR-transformed latent field.
#'
#' @param z A matrix of compositional parameters (n x D) in (0,1).
#' @param alpha Positive precision parameter of the Dirichlet distribution.
#' @param w Observation weights, either scalar (uniform) or vector of length n.
#' @param y Observed compositional data (n x D) in (0,1).
#'
#' @return A list with components:
#'   \item{g}{Vector of length n*d containing the gradient with respect to
#'     the ALR-transformed field}
#'   \item{dg}{Matrix (n x D) of partial derivatives with respect to z}
#'
#' @details
#' The Dirichlet log-likelihood is:
#' \deqn{\log p(y|z,\alpha) = \sum_i [\log\Gamma(\alpha w_i) -
#'   \sum_k \log\Gamma(\alpha w_i z_{ik}) + \sum_k (\alpha w_i z_{ik} - 1)\log y_{ik}]}
#'
#' The gradient with respect to z_k is:
#' \deqn{\frac{\partial}{\partial z_k} = -\alpha w \psi(\alpha w z_k) + \alpha w \log y_k}
#'
#' The gradient with respect to the ALR-transformed field x is obtained by
#' chain rule using the Jacobian from [dalr()].
#'
#' @seealso [hessian_dirichlet()] for second derivatives, [neg_log_lik()] for
#'   the full posterior.
#'
#' @examples
#' z <- matrix(c(0.3, 0.5, 0.2), nrow = 1)
#' y <- matrix(c(0.35, 0.45, 0.2), nrow = 1)
#' grad <- gradient_dirichlet(z, alpha = 10, w = 1, y)
#'
#' @export
gradient_dirichlet <- function(z, alpha, w, y) {
  if (is.vector(z)) z <- matrix(z, nrow = 1)
  if (is.vector(y)) y <- matrix(y, nrow = 1)

  n <- nrow(y)
  D <- ncol(z)
  d <- D - 1

  # Ensure w is a vector
  if (length(w) == 1) w <- rep(w, n)

  # Jacobian of inverse ALR
  dz <- dalr(z)

  # Weighted compositions
  z_w <- z * w

  # Partial derivatives w.r.t. z
  dg <- matrix(0, nrow = n, ncol = D)
  for (l in seq_len(D)) {
    dg[, l] <- -alpha * w * digamma(alpha * z_w[, l]) + alpha * w * log(y[, l])
  }

  # Chain rule: gradient w.r.t. ALR-transformed x
  # g = sum over D of (dg * dz)
  # Replicate dg d times vertically (like MATLAB's repmat(dg, [d, 1]))
  dg_expanded <- dg[rep(1:n, d), , drop = FALSE]
  g <- rowSums(dg_expanded * dz)

  return(list(g = g, dg = dg))
}


#' Hessian of Dirichlet Log-Likelihood
#'
#' Computes the Hessian matrix (second derivatives) of the Dirichlet
#' log-likelihood with respect to the ALR-transformed latent field.
#'
#' @param z A matrix of compositional parameters (n x D) in (0,1).
#' @param alpha Positive precision parameter of the Dirichlet distribution.
#' @param w Observation weights, either scalar or vector of length n.
#' @param y Observed compositional data (n x D) in (0,1).
#'
#' @return A list with components:
#'   \item{H}{Sparse Hessian matrix of dimension (n*d) x (n*d)}
#'   \item{d2gii}{Matrix (n x d) of diagonal second derivatives}
#'   \item{d2gij}{Matrix (n x d) of off-diagonal second derivatives}
#'
#' @details
#' The Hessian is computed using the chain rule:
#' \deqn{H = \frac{\partial^2 \log p}{\partial x \partial x^T} =
#'   \frac{\partial z}{\partial x}^T \frac{\partial^2 \log p}{\partial z^2}
#'   \frac{\partial z}{\partial x} + ...}
#'
#' The resulting Hessian is sparse with a band structure determined by
#' the number of components.
#'
#' @seealso [gradient_dirichlet()] for first derivatives, [fisher_info()]
#'   for the expected Fisher information.
#'
#' @examples
#' z <- matrix(c(0.3, 0.5, 0.2), nrow = 1)
#' y <- matrix(c(0.35, 0.45, 0.2), nrow = 1)
#' H <- hessian_dirichlet(z, alpha = 10, w = 1, y)
#'
#' @export
hessian_dirichlet <- function(z, alpha, w, y) {
  if (is.vector(z)) z <- matrix(z, nrow = 1)
  if (is.vector(y)) y <- matrix(y, nrow = 1)

  n <- nrow(z)
  D <- ncol(z)
  d <- D - 1

  if (length(w) == 1) w <- rep(w, n)

  dz <- dalr(z)
  z_w <- z * w

  # Second derivatives of Dirichlet w.r.t. z
  d2g <- matrix(0, nrow = n, ncol = D)
  for (k in seq_len(D)) {
    d2g[, k] <- -alpha^2 * w^2 * trigamma(alpha * z_w[, k])
  }

  # Expand for chain rule computation (replicate rows d times)
  bigd2g <- d2g[rep(1:n, d), , drop = FALSE]

  # First term: d2g * dz^2
  d1 <- bigd2g * dz^2
  d2 <- rowSums(d1)

  # Get first derivatives for second term
  grad_result <- gradient_dirichlet(z, alpha, w, y)
  dg <- grad_result$dg

  bigdg <- dg[rep(1:n, d), , drop = FALSE]

  # Second term involving first derivatives
  Z1 <- numeric(n * d)
  for (k in seq_len(d)) {
    idx <- ((k - 1) * n + 1):(k * n)
    Z1[idx] <- 1 - 2 * z[, k]
  }

  C <- matrix(rep(Z1, D), nrow = n * d, ncol = D)
  dzC <- dz * C
  d3 <- rowSums(bigdg * dzC)

  # Diagonal elements of Hessian
  d4 <- d2 + d3
  d2gii <- matrix(d4, nrow = n, ncol = d)

  # Off-diagonal elements (cross terms between components)
  # For d=2 case (D=3)
  if (d == 2) {
    h0 <- dz[1:n, , drop = FALSE] * dz[(n + 1):(2 * n), , drop = FALSE]
    h1 <- h0[rep(1:n, d), , drop = FALSE]

    d5 <- rowSums(h1 * bigd2g)

    # Cross-derivative term
    Z0 <- cbind(
      -z[, 1] * z[, 2] * Z1[1:n],
      -z[, 1] * z[, 2] * Z1[(n + 1):(2 * n)],
      2 * z[, 1] * z[, 2] * z[, 3]
    )
    Z2 <- Z0[rep(1:n, d), , drop = FALSE]

    d6 <- rowSums(Z2 * bigdg)
    d7 <- d5 + d6
    d2gij <- matrix(d7, nrow = n, ncol = d)

    # Build sparse Hessian matrix
    # Structure: block tridiagonal for d=2
    G <- cbind(
      d2gij[, 1], d2gii[, 1], rep(0, n),
      rep(0, n), d2gii[, 2], d2gij[, 1]
    )
    G <- matrix(G, nrow = 2 * n, ncol = 3, byrow = FALSE)

    H <- Matrix::bandSparse(
      n = n * d,
      k = c(-n, 0, n),
      diagonals = list(G[1:n, 1], c(G[1:n, 2], G[(n + 1):(2 * n), 2]),
                       G[(n + 1):(2 * n), 3]),
      symmetric = FALSE
    )
  } else {
    # General case - simplified diagonal approximation
    d2gij <- matrix(0, nrow = n, ncol = d)
    H <- Matrix::Diagonal(x = d4)
  }

  return(list(H = H, d2gii = d2gii, d2gij = d2gij))
}


#' Negative Log-Posterior and Gradient
#'
#' Computes the negative log-posterior (for minimization) and its gradient
#' for the spatial compositional model.
#'
#' @param x Vector of ALR-transformed latent field values, length N*d.
#' @param alpha Positive precision parameter.
#' @param A Location matrix connecting observations to grid cells.
#' @param y Observed compositional data (n x D).
#' @param b,c Hyperparameters for the Gamma prior on alpha.
#' @param Q GMRF precision matrix.
#' @param d Number of ALR components (D-1).
#' @param w Observation weights.
#'
#' @return A list with components:
#'   \item{L}{Negative log-posterior value}
#'   \item{dL}{Gradient vector of length N*d + 1 (including alpha)}
#'
#' @details
#' The posterior combines:
#' \itemize{
#'   \item Dirichlet likelihood: \eqn{y_i \sim Dir(\alpha \cdot z_i)}
#'   \item GMRF prior: \eqn{x \sim N(0, Q^{-1})}
#'   \item Gamma prior on alpha: \eqn{\alpha \sim Gamma(b, c)}
#' }
#'
#' @seealso [mcmc_sampling()] which uses this function.
#'
#' @examples
#' # Small example
#' G <- create_Q(c(5, 5), alpha = 1)
#' Qfull <- Q_rhoxQ(NULL, kappa = 0.5, G)$Q
#' n <- 25
#' x <- rnorm(n * 2)  # d = 2
#' y <- matrix(runif(n * 3), n, 3)
#' y <- y / rowSums(y)
#' A <- Matrix::Diagonal(n * 2)
#' result <- neg_log_lik(x, alpha = 10, A, y, b = 1, c = 0.01,
#'                       Q = Matrix::bdiag(Qfull, Qfull), d = 2, w = 1)
#'
#' @export
neg_log_lik <- function(x, alpha, A, y, b, c, Q, d, w) {
  if (length(w) == 1) w <- rep(w, nrow(y))

  n <- nrow(y)

  # Transform to compositional space
  x1 <- matrix(as.vector(A %*% x), nrow = n, ncol = d)
  Az <- invalr(x1)
  Az_w <- Az * w

  # Standard Dirichlet log-likelihood
  logpy <- lgamma(alpha * w) -
    rowSums(lgamma(alpha * Az_w)) +
    rowSums((alpha * Az_w - 1) * log(y))
  sumlogpy <- sum(logpy)

  # Gradient via Dirichlet
  grad_result <- gradient_dirichlet(Az, alpha, w, y)

  # Gradient w.r.t. alpha
  dlogy <- w * digamma(alpha * w) -
    rowSums(Az_w * digamma(alpha * Az_w)) +
    rowSums(Az_w * log(y))

  g <- grad_result$g

  # GMRF prior
  logpx <- -0.5 * as.numeric(t(x) %*% Q %*% x)

  # Gamma prior on alpha
  logpg <- (b - 1) * log(alpha) - c * alpha

  # Negative log-posterior
  L <- -(sumlogpy + logpx + logpg)

  # Gradient w.r.t. x
  dLx <- as.vector(Matrix::t(A) %*% g) - as.vector(Q %*% x)

  # Gradient w.r.t. alpha
  sumdlogy <- sum(dlogy)
  dLalpha <- sumdlogy + (b - 1) / alpha - c

  # Negative gradient
  dL <- -c(dLx, dLalpha)

  return(list(L = L, dL = dL))
}


#' Fisher Information Matrix
#'
#' Computes the Fisher information matrix for the MALA proposal, including
#' the Cholesky factor and mean of the Langevin proposal.
#'
#' @param x Vector of ALR-transformed latent field values.
#' @param alpha Positive precision parameter.
#' @param A Location matrix.
#' @param y Observed compositional data.
#' @param b,c Hyperparameters for alpha prior.
#' @param Q GMRF precision matrix.
#' @param step Step size for MALA.
#' @param d Number of ALR components.
#' @param w Observation weights.
#'
#' @return A list with components:
#'   \item{RFI}{Upper triangular Cholesky factor of the Fisher information}
#'   \item{m}{Mean of the MALA proposal}
#'   \item{p}{Permutation vector from AMD ordering}
#'
#' @details
#' The Fisher information combines:
#' \itemize{
#'   \item Negative Hessian of the Dirichlet log-likelihood
#'   \item GMRF precision matrix Q
#'   \item Second derivative w.r.t. alpha from the gamma prior
#' }
#'
#' AMD (approximate minimum degree) ordering is used to reduce fill-in
#' during the Cholesky factorization.
#'
#' @seealso [mala_step()] which uses this function.
#'
#' @export
fisher_info <- function(x, alpha, A, y, b, c, Q, step, d, w) {
  if (length(w) == 1) w <- rep(w, nrow(y))

  n <- nrow(y)

  # Transform to compositional space
  Ax <- matrix(as.vector(A %*% x), nrow = n, ncol = d)
  Az <- invalr(Ax)
  Az_w <- Az * w

  # Hessian of likelihood
  hess_result <- hessian_dirichlet(Az, alpha, w, y)
  H <- hess_result$H

  # Fisher information for x (negative Hessian + prior precision)
  X2 <- -Matrix::t(A) %*% H %*% A + Q

  # Fisher information for alpha and cross term
  dzsj <- dalr(Az)
  D <- ncol(Az)

  # Standard Dirichlet version
  alpha2 <- sum(-w^2 * trigamma(alpha * w) +
                  rowSums(Az_w^2 * trigamma(alpha * Az_w))) +
    (b - 1) / alpha^2

  # Expand Az_w to (n*d) x D by repeating rows
  Az_w_expanded <- matrix(0, nrow = n * d, ncol = D)
  for (k in seq_len(d)) {
    Az_w_expanded[((k - 1) * n + 1):(k * n), ] <- Az_w
  }

  # Expand w to (n*d)
  w_expanded <- rep(w, d)

  # Xalpha = sum of (alpha * w * Az_w * dz * trigamma(alpha * Az_w))
  Xalpha_terms <- alpha * w_expanded * Az_w_expanded
  Xalpha_psi <- trigamma(alpha * Az_w_expanded)
  Xalpha <- rowSums(Xalpha_terms * dzsj * Xalpha_psi)

  # Build full Fisher information matrix
  FI <- rbind(
    cbind(X2, Matrix::t(A) %*% Xalpha),
    c(as.vector(Matrix::t(Xalpha) %*% A), alpha2)
  )

  # Convert to dense matrix for Cholesky
  FI <- as.matrix(FI)

  # Use identity ordering (like MATLAB amd, but simpler for small problems)
  p <- seq_len(nrow(FI))
  FI_perm <- FI

  # Cholesky factorization with regularization if needed (like MATLAB)
  RFI <- tryCatch(
    chol(FI_perm),
    error = function(e) NULL
  )

  if (is.null(RFI)) {
    # Add small diagonal for numerical stability (same as MATLAB: 1e-3)
    RFI <- tryCatch(
      chol(FI_perm + 1e-3 * diag(nrow(FI_perm))),
      error = function(e) NULL
    )
  }

  if (is.null(RFI)) {
    # Fall back to diagonal approximation if still failing
    RFI <- diag(sqrt(pmax(diag(FI_perm), 1e-6)))
  }

  # Compute negative log-posterior gradient
  nll_result <- neg_log_lik(x, alpha, A, y, b, c, Q, d, w)
  dl <- -nll_result$dL  # Gradient of log-posterior
  dl <- dl[p]

  # MALA proposal mean
  old <- c(x, alpha)[p]
  m <- old + 0.5 * step^2 * backsolve(RFI, forwardsolve(t(RFI), dl))

  return(list(RFI = RFI, m = m, p = p))
}
