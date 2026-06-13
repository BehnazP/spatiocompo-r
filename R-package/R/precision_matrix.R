#' Create GMRF Precision Matrix for Regular Grid
#'
#' Creates a precision matrix Q for a Gaussian Markov Random Field on a
#' regular grid with Neumann (zero-flux) or torus (periodic) boundary conditions.
#'
#' @param sz Integer vector specifying grid dimensions, e.g., c(n1, n2) for 2D.
#' @param alpha Order of the Laplace operator: 1 for first-order, 2 for
#'   second-order (Matern with nu=1). Default is 2.
#' @param a Positive weight vector for anisotropy in each dimension.
#'   Scalar values are repeated for all dimensions. Default is 1.
#' @param neumann Logical; if TRUE (default), use Neumann boundary conditions,
#'   otherwise use torus (periodic) boundaries.
#'
#' @return A sparse precision matrix Q of dimension prod(sz) x prod(sz),
#'   stored as a Matrix::dgCMatrix.
#'
#' @details
#' The precision matrix corresponds to the discretization of:
#' \deqn{(-\Delta)^{\alpha/2} u = \epsilon}
#' on a regular grid. For alpha=2, this gives a Matern-like field with
#' smoothness nu=1.
#'
#' The anisotropy weights allow different correlation lengths in each
#' dimension: Q = (sum_i a_i * W_i)^(alpha/2) where W_i is the
#' second-difference matrix in dimension i.
#'
#' @seealso [Q_rhoxQ()] for the full GMRF precision with spatial range.
#'
#' @examples
#' # 2D grid 10x10 with second-order Laplacian
#' Q <- create_Q(c(10, 10), alpha = 2)
#' dim(Q)  # 100 x 100
#'
#' # 1D grid
#' Q1d <- create_Q(20, alpha = 1)
#'
#' @export
create_Q <- function(sz, alpha = 2, a = 1, neumann = TRUE) {
  # Input validation
  if (!alpha %in% c(1, 2)) {
    stop("alpha should be 1 or 2")
  }

  sz <- as.integer(sz)

  # Remove singleton dimensions
  valid <- sz > 1
  if (!all(valid)) {
    sz <- sz[valid]
    if (length(a) > 1) {
      a <- a[valid]
    }
  }

  ndim <- length(sz)

  # Expand a to match dimensions
  if (length(a) == 1) {
    a <- rep(a, ndim)
  }
  if (length(a) != ndim) {
    stop("a must be scalar or have length equal to number of dimensions")
  }
  if (any(a < 0)) {
    stop("a must be non-negative")
  }

  # Create W matrix for each dimension
  W_list <- lapply(sz, function(n) create_W(n, neumann))

  # Use Kronecker products to expand to full grid
  if (ndim == 1) {
    W_expanded <- W_list
  } else {
    W_expanded <- vector("list", ndim)

    # First dimension
    W_expanded[[1]] <- Matrix::kronecker(
      Matrix::Diagonal(prod(sz[2:ndim])),
      W_list[[1]]
    )

    # Middle dimensions
    if (ndim > 2) {
      for (i in 2:(ndim - 1)) {
        W_expanded[[i]] <- Matrix::kronecker(
          Matrix::Diagonal(prod(sz[(i + 1):ndim])),
          Matrix::kronecker(W_list[[i]], Matrix::Diagonal(prod(sz[1:(i - 1)])))
        )
      }
    }

    # Last dimension
    W_expanded[[ndim]] <- Matrix::kronecker(
      W_list[[ndim]],
      Matrix::Diagonal(prod(sz[1:(ndim - 1)]))
    )
  }

  # Accumulate into Q matrix with anisotropy weights
  Q <- a[1] * W_expanded[[1]]
  if (ndim > 1) {
    for (i in 2:ndim) {
      Q <- Q + a[i] * W_expanded[[i]]
    }
  }

  # Apply Laplacian order
  if (alpha == 2) {
    Q <- Q %*% Q
  }

  return(Q)
}


#' Create Second-Difference Matrix for 1D
#'
#' Internal helper function to create the second-difference (Laplacian)
#' matrix for a 1D grid.
#'
#' @param n Grid size.
#' @param neumann If TRUE, use Neumann boundaries; otherwise torus.
#'
#' @return Sparse matrix of dimension n x n.
#'
#' @keywords internal
create_W <- function(n, neumann = TRUE) {
  # Main diagonal and off-diagonals
  diag_vals <- rep(2, n)
  off_vals <- rep(-1, n - 1)

  # Create symmetric tridiagonal matrix
  # Use symmetric = FALSE and specify all diagonals explicitly
  W <- Matrix::bandSparse(
    n = n,
    k = c(-1, 0, 1),
    diagonals = list(off_vals, diag_vals, off_vals),
    symmetric = FALSE
  )
  # Ensure symmetry
  W <- Matrix::forceSymmetric(W)

  if (neumann) {
    # Neumann: modify corner values
    W[1, 1] <- 1
    W[n, n] <- 1
  } else {
    # Torus: add wrap-around connections
    W[1, n] <- -1
    W[n, 1] <- -1
  }

  return(W)
}


#' GMRF Precision Matrix with Spatial Covariance
#'
#' Computes the full GMRF precision matrix including the spatial range
#' parameter kappa and optionally the cross-component covariance structure.
#'
#' @param rho A d x d positive definite covariance matrix for the d-1
#'   transformed compositional components. Set to NULL to return only
#'   the marginal Q matrix.
#' @param kappa Positive spatial range parameter. Larger values give
#'   shorter correlation range (smoother fields).
#' @param G The base precision matrix from [create_Q()], typically with
#'   alpha = 1.
#'
#' @return A list with components:
#'   \item{Q}{The marginal precision matrix for one component, dimension N x N}
#'   \item{rhoxQ}{The full precision matrix kron(inv(rho), Q) of dimension
#'     Nd x Nd, only returned if rho is provided}
#'
#' @details
#' The precision matrix is computed as:
#' \deqn{Q = \kappa^4 I + 2\kappa^2 G + G^T G}
#' where G is the first-order Laplacian. This corresponds to a Matern field
#' with smoothness nu = 1.
#'
#' When rho is provided, the full precision for the multivariate field is
#' kron(inv(rho), Q), which models cross-component correlations.
#'
#' @seealso [create_Q()] for creating G, [log_det_Q()] for efficient
#'   log-determinant computation.
#'
#' @examples
#' # Create base matrix G
#' G <- create_Q(c(10, 10), alpha = 1)
#'
#' # Marginal precision with kappa = 0.5
#' result <- Q_rhoxQ(NULL, kappa = 0.5, G)
#' Q <- result$Q
#'
#' # With 2x2 covariance matrix
#' rho <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
#' result <- Q_rhoxQ(rho, kappa = 0.5, G)
#'
#' @export
Q_rhoxQ <- function(rho, kappa, G) {
  N <- nrow(G)

  # Compute Q = kappa^4 * I + 2*kappa^2 * G + G'G
  Q <- kappa^4 * Matrix::Diagonal(N) + 2 * kappa^2 * G + Matrix::t(G) %*% G

  result <- list(Q = Q)

  if (!is.null(rho)) {
    # Compute kron(inv(rho), Q)
    rho_inv <- solve(rho)
    result$rhoxQ <- Matrix::kronecker(rho_inv, Q)
  }

  return(result)
}


#' Eigenvalues for Efficient Log-Determinant
#'
#' Computes the eigenvalues of the base Laplacian matrix G, which can be
#' reused to efficiently compute log|Q| for different values of kappa.
#'
#' @param G The base precision matrix from [create_Q()] with alpha = 1.
#'
#' @return A numeric vector of eigenvalues.
#'
#' @details
#' For the precision matrix Q = kappa^4*I + 2*kappa^2*G + G'G, if G has
#' eigenvalues lambda_i, then Q has eigenvalues:
#' \deqn{\xi_i = \kappa^4 + 2\kappa^2 \lambda_i + \lambda_i^2}
#' and log|Q| = sum(log(xi_i)).
#'
#' Computing the eigenvalues once and reusing them is much faster than
#' repeatedly computing log-determinants via Cholesky factorization.
#'
#' @seealso [log_det_Q()] which uses these eigenvalues.
#'
#' @examples
#' G <- create_Q(c(10, 10), alpha = 1)
#' lambda <- comp_lambda(G)
#' log_det_Q(kappa = 0.5, lambda)
#'
#' @export
comp_lambda <- function(G) {
  # For symmetric matrices, use eigen
  # Note: G is sparse, so we convert to dense for eigen
  # For very large matrices, consider using RSpectra package

  if (nrow(G) > 5000) {
    warning("Large matrix - eigenvalue computation may be slow")
  }

  lambda <- eigen(as.matrix(G), symmetric = TRUE, only.values = TRUE)$values

  return(lambda)
}


#' Log-Determinant of GMRF Precision Matrix
#'
#' Efficiently computes the log-determinant of the GMRF precision matrix Q
#' using precomputed eigenvalues.
#'
#' @param kappa Positive spatial range parameter.
#' @param lambda Vector of eigenvalues from [comp_lambda()].
#'
#' @return The log-determinant log|Q|.
#'
#' @details
#' Uses the formula:
#' \deqn{\log|Q| = \sum_i \log(\kappa^4 + 2\kappa^2 \lambda_i + \lambda_i^2)}
#'
#' This is O(N) once eigenvalues are computed, compared to \eqn{O(N^{1.5})} for
#' sparse Cholesky.
#'
#' @seealso [comp_lambda()] for computing eigenvalues.
#'
#' @examples
#' G <- create_Q(c(10, 10), alpha = 1)
#' lambda <- comp_lambda(G)
#' log_det_Q(kappa = 0.5, lambda)
#' log_det_Q(kappa = 1.0, lambda)
#'
#' @export
log_det_Q <- function(kappa, lambda) {
  xi <- kappa^4 + 2 * kappa^2 * lambda + lambda^2
  sum(log(xi))
}
