test_that("mala_step returns correct structure", {
  set.seed(123)
  n <- 9
  d <- 2
  N <- n

  x <- rnorm(N * d, sd = 0.5)
  y <- matrix(runif(n * 3), n, 3)
  y <- y / rowSums(y)

  G <- create_Q(c(3, 3), alpha = 1)
  rho <- diag(2)
  Q <- Q_rhoxQ(rho, kappa = 0.5, G)$rhoxQ

  A <- Matrix::Diagonal(N * d)

  priors <- list(alpha.a = 1, alpha.b = 0.01)

  result <- mala_step(x, alpha = 50, y, w = 1, Q, A, priors, d, logstep0 = 0, iter = 10)

  expect_true("sample" %in% names(result))
  expect_true("count" %in% names(result))
  expect_true("logstep" %in% names(result))

  # Sample should have length N*d + 1
  expect_equal(length(result$sample), N * d + 1)

  # Count should be 0 or 1
  expect_true(result$count %in% c(0, 1))
})

test_that("gibbs_kappa_rho returns correct structure", {
  set.seed(123)
  N <- 25
  d <- 2

  x <- matrix(rnorm(N * d, sd = 0.5), nrow = N, ncol = d)
  kappa <- 0.5
  rho <- diag(d)

  G <- create_Q(c(5, 5), alpha = 1)
  lambda <- comp_lambda(G)

  priors <- list(
    kappa.range = 2,
    kappa.sigma = 1,
    rho.df = 3,
    rho.Sigma = diag(d),
    field.G = G,
    field.lambda = lambda
  )

  result <- gibbs_kappa_rho(x, kappa, rho, priors, logstep0 = 0, iter = 10)

  expect_true("sample" %in% names(result))
  expect_true(result$sample$kappa > 0)
  expect_equal(dim(result$sample$rho), c(d, d))
  expect_true(result$count %in% c(0, 1))
})

test_that("sample_invwishart returns positive definite matrix", {
  set.seed(123)
  N <- 25
  d <- 2

  x <- matrix(rnorm(N * d), nrow = N, ncol = d)
  G <- create_Q(c(5, 5), alpha = 1)

  priors <- list(
    field.G = G,
    rho.df = 3,
    rho.Sigma = diag(d)
  )

  rho <- sample_invwishart(x, kappa = 0.5, priors)

  expect_equal(dim(rho), c(d, d))

  # Check positive definite
  eig <- eigen(rho, only.values = TRUE)$values
  expect_true(all(eig > 0))

  # Check symmetric
  expect_equal(rho, t(rho), tolerance = 1e-10)
})

test_that("log_posterior_kappa_x returns finite values", {
  set.seed(123)
  N <- 25
  d <- 2

  x <- matrix(rnorm(N * d, sd = 0.5), nrow = N, ncol = d)
  G <- create_Q(c(5, 5), alpha = 1)
  lambda <- comp_lambda(G)

  priors <- list(
    kappa.range = 2,
    kappa.sigma = 1,
    rho.df = 3,
    rho.Sigma = diag(d),
    field.G = G,
    field.lambda = lambda
  )

  logp <- log_posterior_kappa_x(x, kappa = 0.5, priors)

  expect_true(is.finite(logp))
})

test_that("mcmc_sampling runs without error", {
  skip_on_cran()  # Skip on CRAN due to time

  set.seed(123)

  # Small simulation
  sim <- simulate_compo_data(n_side = 5, kappa_true = 0.5,
                             alpha_true = 50, seed = 123)

  priors <- create_priors(sim$G, d = 2)

  theta0 <- list(
    alpha = 30,
    kappa = 0.3,
    rho = diag(2),
    x = matrix(0, nrow = 25, ncol = 2)
  )

  # Short run
  result <- mcmc_sampling(theta0, priors, sim$y, iter = 20, verbose = FALSE)

  expect_true("alpha" %in% names(result))
  expect_true("kappa" %in% names(result))
  expect_true("rho" %in% names(result))
  expect_true("x" %in% names(result))

  expect_equal(length(result$alpha), 20)
  expect_equal(length(result$kappa), 20)
})

test_that("mcmc_sampling produces reasonable acceptance rates", {
  skip_on_cran()

  set.seed(456)

  sim <- simulate_compo_data(n_side = 5, kappa_true = 0.5,
                             alpha_true = 50, seed = 456)

  priors <- create_priors(sim$G, d = 2)

  theta0 <- find_initial_values(sim$y, priors, verbose = FALSE)

  result <- mcmc_sampling(theta0, priors, sim$y, iter = 50, verbose = FALSE)

  # Acceptance rates should be between 0 and 1
  expect_true(result$MALA$acc >= 0 && result$MALA$acc <= 1)
  expect_true(result$Gibbs$acc >= 0 && result$Gibbs$acc <= 1)
})
