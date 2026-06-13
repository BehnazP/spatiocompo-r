test_that("create_Q produces symmetric matrix", {
  Q <- create_Q(c(5, 5), alpha = 2)

  # Check symmetry
  expect_true(Matrix::isSymmetric(Q))
})

test_that("create_Q has correct dimensions", {
  Q1 <- create_Q(c(10, 10), alpha = 2)
  expect_equal(dim(Q1), c(100, 100))

  Q2 <- create_Q(c(5, 5, 5), alpha = 1)
  expect_equal(dim(Q2), c(125, 125))

  Q3 <- create_Q(20, alpha = 1)  # 1D
  expect_equal(dim(Q3), c(20, 20))
})

test_that("create_Q with alpha=1 is different from alpha=2", {
  Q1 <- create_Q(c(5, 5), alpha = 1)
  Q2 <- create_Q(c(5, 5), alpha = 2)

  expect_false(isTRUE(all.equal(as.matrix(Q1), as.matrix(Q2))))
})

test_that("create_Q respects Neumann vs torus boundaries", {
  Q_neumann <- create_Q(c(5, 5), alpha = 1, neumann = TRUE)
  Q_torus <- create_Q(c(5, 5), alpha = 1, neumann = FALSE)

  # They should be different
  expect_false(isTRUE(all.equal(as.matrix(Q_neumann), as.matrix(Q_torus))))

  # Neumann has zero row sums at boundaries
  # Torus is doubly stochastic (row sums constant)
})

test_that("create_Q with anisotropy weights works", {
  Q_iso <- create_Q(c(5, 5), alpha = 1, a = c(1, 1))
  Q_aniso <- create_Q(c(5, 5), alpha = 1, a = c(1, 2))

  expect_false(isTRUE(all.equal(as.matrix(Q_iso), as.matrix(Q_aniso))))
})

test_that("Q_rhoxQ returns correct structure", {
  G <- create_Q(c(5, 5), alpha = 1)
  rho <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  result <- Q_rhoxQ(rho, kappa = 0.5, G)

  # Q should be 25 x 25

  expect_equal(dim(result$Q), c(25, 25))

  # rhoxQ should be 50 x 50 (25*2)
  expect_equal(dim(result$rhoxQ), c(50, 50))
})

test_that("Q_rhoxQ with NULL rho returns only Q", {
  G <- create_Q(c(5, 5), alpha = 1)
  result <- Q_rhoxQ(NULL, kappa = 0.5, G)

  expect_true("Q" %in% names(result))
  expect_false("rhoxQ" %in% names(result))
})

test_that("Q_rhoxQ produces positive definite matrix", {
  G <- create_Q(c(5, 5), alpha = 1)
  result <- Q_rhoxQ(NULL, kappa = 0.5, G)

  # All eigenvalues should be positive
  eig <- eigen(as.matrix(result$Q), only.values = TRUE)$values
  expect_true(all(eig > 0))
})

test_that("comp_lambda returns correct number of eigenvalues", {
  G <- create_Q(c(5, 5), alpha = 1)
  lambda <- comp_lambda(G)

  expect_equal(length(lambda), 25)
})

test_that("log_det_Q is consistent with direct computation", {
  G <- create_Q(c(5, 5), alpha = 1)
  lambda <- comp_lambda(G)
  kappa <- 0.5

  # Via eigenvalues
  logdet1 <- log_det_Q(kappa, lambda)

  # Direct computation via Cholesky
  Q <- Q_rhoxQ(NULL, kappa, G)$Q
  R <- chol(Q + 1e-10 * Matrix::Diagonal(25))  # Small regularization
  logdet2 <- 2 * sum(log(Matrix::diag(R)))

  expect_equal(logdet1, logdet2, tolerance = 0.01)
})

test_that("log_det_Q increases with kappa", {
  G <- create_Q(c(5, 5), alpha = 1)
  lambda <- comp_lambda(G)

  ld1 <- log_det_Q(0.1, lambda)
  ld2 <- log_det_Q(0.5, lambda)
  ld3 <- log_det_Q(1.0, lambda)

  expect_true(ld1 < ld2)
  expect_true(ld2 < ld3)
})
