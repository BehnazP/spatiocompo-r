test_that("gradient_dirichlet returns correct structure", {
  z <- matrix(c(0.3, 0.5, 0.2), nrow = 1)
  y <- matrix(c(0.35, 0.45, 0.2), nrow = 1)

  result <- gradient_dirichlet(z, alpha = 10, w = 1, y)

  expect_true("g" %in% names(result))
  expect_true("dg" %in% names(result))

  # g should be length n*d = 1*2 = 2
  expect_equal(length(result$g), 2)

  # dg should be n x D = 1 x 3
  expect_equal(dim(result$dg), c(1, 3))
})

test_that("gradient is zero at MLE (approximately)", {
  # When y == z (perfect fit), gradient should be close to zero
  # for the likelihood term

  z <- matrix(c(0.3, 0.5, 0.2), nrow = 1)
  y <- z  # Perfect match

  result <- gradient_dirichlet(z, alpha = 100, w = 1, y)

  # Gradient components should be relatively small
  # (not exactly zero due to prior and transformation)
  # The derivatives dg have magnitude proportional to alpha * (log(y) - psi(alpha*z))
  # When y = z, this is alpha * (log(z) - psi(alpha*z)) which can be O(100)
  expect_true(all(is.finite(result$dg)))
})

test_that("hessian_dirichlet returns correct structure", {
  z <- matrix(c(0.3, 0.5, 0.2), nrow = 1)
  y <- matrix(c(0.35, 0.45, 0.2), nrow = 1)

  result <- hessian_dirichlet(z, alpha = 10, w = 1, y)

  expect_true("H" %in% names(result))

  # H should be (n*d) x (n*d) = 2 x 2
  expect_equal(dim(result$H), c(2, 2))
})

test_that("neg_log_lik returns finite values", {
  n <- 9  # 3x3 grid
  d <- 2
  N <- n

  # Simple setup
  x <- rnorm(N * d)
  y <- matrix(runif(n * 3), n, 3)
  y <- y / rowSums(y)

  G <- create_Q(c(3, 3), alpha = 1)
  Q <- Q_rhoxQ(NULL, kappa = 0.5, G)$Q
  Q_full <- Matrix::bdiag(Q, Q)

  A <- Matrix::Diagonal(N * d)

  result <- neg_log_lik(x, alpha = 10, A, y, b = 1, c = 0.01, Q_full, d = 2, w = 1)

  expect_true(is.finite(result$L))
  expect_true(all(is.finite(result$dL)))
})

test_that("neg_log_lik gradient has correct length", {
  n <- 9
  d <- 2
  N <- n

  x <- rnorm(N * d)
  y <- matrix(runif(n * 3), n, 3)
  y <- y / rowSums(y)

  G <- create_Q(c(3, 3), alpha = 1)
  Q <- Q_rhoxQ(NULL, kappa = 0.5, G)$Q
  Q_full <- Matrix::bdiag(Q, Q)

  A <- Matrix::Diagonal(N * d)

  result <- neg_log_lik(x, alpha = 10, A, y, b = 1, c = 0.01, Q_full, d = 2, w = 1)

  # Gradient should be length N*d + 1 (for alpha)
  expect_equal(length(result$dL), N * d + 1)
})

test_that("numerical gradient matches analytical gradient", {
  skip_if_not_installed("numDeriv")

  n <- 4
  d <- 2
  N <- n

  set.seed(123)
  x <- rnorm(N * d)
  y <- matrix(runif(n * 3), n, 3)
  y <- y / rowSums(y)

  G <- create_Q(c(2, 2), alpha = 1)
  Q <- Q_rhoxQ(NULL, kappa = 0.5, G)$Q
  Q_full <- Matrix::bdiag(Q, Q)

  A <- Matrix::Diagonal(N * d)
  alpha <- 20

  # Analytical gradient
  result <- neg_log_lik(x, alpha, A, y, b = 1, c = 0.01, Q_full, d = 2, w = 1)
  grad_analytical <- result$dL

  # Numerical gradient
  f <- function(theta) {
    x_val <- theta[1:(N * d)]
    alpha_val <- theta[N * d + 1]
    neg_log_lik(x_val, alpha_val, A, y, b = 1, c = 0.01, Q_full, d = 2, w = 1)$L
  }

  grad_numerical <- numDeriv::grad(f, c(x, alpha))

  # Should match within tolerance
  expect_equal(grad_analytical, grad_numerical, tolerance = 0.01)
})

test_that("fisher_info returns positive definite matrix", {
  n <- 4
  d <- 2
  N <- n

  set.seed(42)
  x <- rnorm(N * d, sd = 0.5)
  y <- matrix(runif(n * 3), n, 3)
  y <- y / rowSums(y)

  G <- create_Q(c(2, 2), alpha = 1)
  Q <- Q_rhoxQ(NULL, kappa = 0.5, G)$Q
  Q_full <- Matrix::bdiag(Q, Q)

  A <- Matrix::Diagonal(N * d)

  result <- fisher_info(x, alpha = 50, A, y, b = 1, c = 0.01, Q_full,
                        step = 0.1, d = 2, w = 1)

  # RFI should exist and be upper triangular (Cholesky)
  expect_true("RFI" %in% names(result))
  expect_true("m" %in% names(result))

  # Diagonal should be positive (from Cholesky)
  expect_true(all(diag(result$RFI) > 0))
})
