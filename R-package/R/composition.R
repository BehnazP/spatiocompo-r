#' Additive Log-Ratio Transformation
#'
#' Transforms compositional data from the simplex to real space using the
#' additive log-ratio (ALR) transformation. The last component is used as
#' the reference.
#'
#' @param compo A matrix of compositional data with n rows (observations) and
#'   D columns (components). Each row must sum to 1, and all values must be
#'   in the interval (0, 1).
#'
#' @return A matrix of transformed values with n rows and d = D-1 columns,
#'   where values are in the interval (-Inf, Inf).
#'
#' @details
#' The ALR transformation for component k is defined as:
#' \deqn{y_k = \log(x_k / x_D)}
#' where \eqn{x_D} is the reference (last) component.
#'
#' @seealso [invalr()] for the inverse transformation, [dalr()] for the Jacobian.
#'
#' @examples
#' # Create compositional data
#' compo <- matrix(c(0.3, 0.5, 0.2,
#'                   0.4, 0.4, 0.2), nrow = 2, byrow = TRUE)
#' # Transform to real space
#' y <- alr(compo)
#'
#' # Verify inverse
#' all.equal(compo, invalr(y))
#'
#' @export
alr <- function(compo) {
  if (is.vector(compo)) {
    compo <- matrix(compo, nrow = 1)
  }

  n <- nrow(compo)
  D <- ncol(compo)
  d <- D - 1

  # Check validity
  if (any(compo <= 0 | compo >= 1)) {
    warning("Compositional data should be strictly between 0 and 1")
  }

  transform <- matrix(0, nrow = n, ncol = d)
  for (k in seq_len(d)) {
    transform[, k] <- log(compo[, k] / compo[, D])
  }

  return(transform)
}


#' Inverse Additive Log-Ratio Transformation
#'
#' Transforms data from real space back to the simplex using the inverse
#' additive log-ratio (ALR) transformation.
#'
#' @param x A matrix of ALR-transformed data with n rows (observations) and
#'   d columns. Values can be any real number.
#'
#' @return A matrix of compositional data with n rows and D = d+1 columns.
#'   Each row sums to 1, and all values are in (0, 1).
#'
#' @details
#' The inverse ALR transformation is:
#' \deqn{x_k = \exp(y_k) / (1 + \sum_{j=1}^{d} \exp(y_j))} for k = 1, ..., d
#' \deqn{x_D = 1 / (1 + \sum_{j=1}^{d} \exp(y_j))}
#'
#' @seealso [alr()] for the forward transformation.
#'
#' @examples
#' # ALR-transformed data
#' y <- matrix(c(0.5, -0.3,
#'               0.2, 0.1), nrow = 2, byrow = TRUE)
#' # Transform back to compositions
#' compo <- invalr(y)
#' rowSums(compo)  # Should all be 1
#'
#' @export
invalr <- function(x) {
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }

  n <- nrow(x)
  d <- ncol(x)
  D <- d + 1

  transform <- matrix(0, nrow = n, ncol = D)

  # Compute sum of exp(x) for numerical stability
  sum_exp <- rowSums(exp(x))
  denom <- 1 + sum_exp

  for (k in seq_len(d)) {
    transform[, k] <- exp(x[, k]) / denom
  }
  transform[, D] <- 1 / denom

  return(transform)
}


#' Jacobian of the Inverse ALR Transformation
#'
#' Computes the Jacobian matrix of the inverse ALR transformation,
#' i.e., the partial derivatives dz/dx where z is compositional and
#' x is ALR-transformed.
#'
#' @param z A matrix of compositional data with n rows (observations) and
#'   D columns (components).
#'
#' @return A matrix of size (n*d) x D containing the Jacobian. The matrix
#'   is structured in blocks where rows (k-1)*n+1 to k*n correspond to
#'   the derivatives of component k with respect to all D components.
#'
#' @details
#' The derivatives are:
#' \itemize{
#'   \item \eqn{dz_k/dx_j = z_k(1 - z_k)} if k = j
#'   \item \eqn{dz_k/dx_j = -z_k z_j} if k != j
#' }
#'
#' @seealso [alr()], [invalr()]
#'
#' @examples
#' compo <- matrix(c(0.3, 0.5, 0.2), nrow = 1)
#' J <- dalr(compo)
#'
#' @export
dalr <- function(z) {
  if (is.vector(z)) {
    z <- matrix(z, nrow = 1)
  }

  n <- nrow(z)
  D <- ncol(z)
  d <- D - 1

  dhkj <- matrix(0, nrow = n * d, ncol = D)

  for (k in seq_len(d)) {
    for (j in seq_len(D)) {
      row_idx <- ((k - 1) * n + 1):(k * n)
      if (k == j) {
        dhkj[row_idx, j] <- z[, k] * (1 - z[, j])
      } else {
        dhkj[row_idx, j] <- -z[, k] * z[, j]
      }
    }
  }

  return(dhkj)
}


#' Aitchison Distance Between Compositions
#'
#' Computes the Aitchison distance between two compositional vectors or
#' between rows of two compositional matrices.
#'
#' @param x,y Compositional vectors or matrices. If matrices, they must
#'   have the same dimensions, and the distance is computed row-wise.
#'
#' @return A numeric value (if vectors) or vector (if matrices) of
#'   Aitchison distances.
#'
#' @details
#' The Aitchison distance is defined as the Euclidean distance between
#' CLR (centered log-ratio) transformed compositions:
#' \deqn{d_A(x, y) = \sqrt{\sum_{k=1}^{D} (\log(x_k/g(x)) - \log(y_k/g(y)))^2}}
#' where g(x) is the geometric mean of x.
#'
#' @examples
#' x <- c(0.3, 0.5, 0.2)
#' y <- c(0.4, 0.4, 0.2)
#' compo_dist(x, y)
#'
#' @export
compo_dist <- function(x, y) {
  if (is.vector(x)) x <- matrix(x, nrow = 1)
  if (is.vector(y)) y <- matrix(y, nrow = 1)

  if (!all(dim(x) == dim(y))) {
    stop("x and y must have the same dimensions")
  }

  # CLR transformation
  clr_x <- log(x) - rowMeans(log(x))
  clr_y <- log(y) - rowMeans(log(y))

  # Euclidean distance in CLR space
  sqrt(rowSums((clr_x - clr_y)^2))
}


#' Handle Boundary Values in Compositional Data
#'
#' Replaces exact zeros and ones in compositional data with small positive
#' values to ensure the data is strictly in the open simplex (0, 1).
#'
#' @param compo A matrix or vector of compositional data.
#' @param eps Small value to use for replacement. Default is 1e-6.
#'
#' @return Compositional data with zeros replaced by eps and ones replaced
#'   by 1-eps, then renormalized to sum to 1.
#'
#' @details
#' Exact zeros and ones cause problems for log-ratio transformations.
#' This function applies simple replacement followed by renormalization.
#' For more sophisticated approaches, consider Bayesian imputation or
#' the multiplicative replacement strategy.
#'
#' @examples
#' compo <- matrix(c(0, 0.5, 0.5,
#'                   1, 0, 0), nrow = 2, byrow = TRUE)
#' no_zero_one(compo)
#'
#' @export
no_zero_one <- function(compo, eps = 1e-6) {
  if (is.vector(compo)) {
    compo <- matrix(compo, nrow = 1)
    was_vector <- TRUE
  } else {
    was_vector <- FALSE
  }

  # Replace zeros and ones
  compo[compo <= 0] <- eps
  compo[compo >= 1] <- 1 - eps

  # Renormalize rows to sum to 1
  compo <- compo / rowSums(compo)

  if (was_vector) {
    compo <- as.vector(compo)
  }

  return(compo)
}
