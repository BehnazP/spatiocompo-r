test_that("simulate_compo_data returns correct structure", {
  sim <- simulate_compo_data(n_side = 10, seed = 123)

  expect_true("y" %in% names(sim))
  expect_true("x_true" %in% names(sim))
  expect_true("z_true" %in% names(sim))
  expect_true("G" %in% names(sim))
  expect_true("params" %in% names(sim))
  expect_true("grid" %in% names(sim))
})

test_that("simulate_compo_data produces valid dimensions", {
  n_side <- 15
  D <- 3
  sim <- simulate_compo_data(n_side = n_side, D = D, seed = 123)

  N <- n_side^2
  d <- D - 1

  expect_equal(dim(sim$y), c(N, D))
  expect_equal(dim(sim$x_true), c(N, d))
  expect_equal(dim(sim$z_true), c(N, D))
  expect_equal(dim(sim$G), c(N, N))
})

test_that("simulate_compo_data produces valid compositions", {
  sim <- simulate_compo_data(n_side = 10, seed = 123)

  # y should be in (0, 1)
  expect_true(all(sim$y > 0 & sim$y < 1))

  # Rows should sum to 1
  expect_equal(rowSums(sim$y), rep(1, 100), tolerance = 1e-10)

  # Same for z_true
  expect_true(all(sim$z_true > 0 & sim$z_true < 1))
  expect_equal(rowSums(sim$z_true), rep(1, 100), tolerance = 1e-10)
})

test_that("simulate_compo_data is reproducible with seed", {
  sim1 <- simulate_compo_data(n_side = 10, seed = 42)
  sim2 <- simulate_compo_data(n_side = 10, seed = 42)

  expect_equal(sim1$y, sim2$y)
  expect_equal(sim1$x_true, sim2$x_true)
})

test_that("simulate_compo_data varies with seed", {
  sim1 <- simulate_compo_data(n_side = 10, seed = 1)
  sim2 <- simulate_compo_data(n_side = 10, seed = 2)

  expect_false(isTRUE(all.equal(sim1$y, sim2$y)))
})

test_that("simulate_compo_data stores true parameters", {
  kappa_true <- 0.7
  alpha_true <- 100

  sim <- simulate_compo_data(n_side = 10, kappa_true = kappa_true,
                             alpha_true = alpha_true, seed = 123)

  expect_equal(sim$params$kappa, kappa_true)
  expect_equal(sim$params$alpha, alpha_true)
})

test_that("simulate_landcover returns correct structure", {
  sim <- simulate_landcover(n_side = 10, seed = 123)

  expect_true("categories" %in% names(sim))
  expect_equal(sim$categories, c("Forest", "Open", "Other"))
  expect_equal(ncol(sim$y), 3)
})

test_that("simulate_landcover with gradient pattern works", {
  sim <- simulate_landcover(n_side = 20, pattern = "gradient", seed = 123)

  expect_true("trend" %in% names(sim))

  # With gradient, northern cells should have more forest (z[,1])
  coords <- sim$grid$coords
  north_idx <- which(coords[, 2] > 15)
  south_idx <- which(coords[, 2] < 5)

  # On average, forest proportion should be higher in north
  # (this is probabilistic, so we use a weak test)
  expect_true(is.finite(mean(sim$z_true[north_idx, 1])))
  expect_true(is.finite(mean(sim$z_true[south_idx, 1])))
})

test_that("simulate_landcover with different patterns", {
  sim_grad <- simulate_landcover(n_side = 15, pattern = "gradient", seed = 123)
  sim_patch <- simulate_landcover(n_side = 15, pattern = "patchy", seed = 123)
  sim_rand <- simulate_landcover(n_side = 15, pattern = "random", seed = 123)

  # All should produce valid data
  expect_true(all(sim_grad$y > 0))
  expect_true(all(sim_patch$y > 0))
  expect_true(all(sim_rand$y > 0))

  # Patterns should differ
  expect_false(isTRUE(all.equal(sim_grad$trend, sim_rand$trend)))
})

test_that("simulate_compo_data works with custom rho", {
  rho_custom <- matrix(c(2, 0.5, 0.5, 1), 2, 2)

  sim <- simulate_compo_data(n_side = 10, rho_true = rho_custom, seed = 123)

  expect_equal(sim$params$rho, rho_custom)
})
