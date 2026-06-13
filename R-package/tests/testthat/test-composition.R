test_that("alr transformation works correctly", {
  # Simple case
  compo <- matrix(c(0.3, 0.5, 0.2), nrow = 1)
  result <- alr(compo)

  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), 1)

  # Verify: alr[1] = log(0.3/0.2), alr[2] = log(0.5/0.2)
  expect_equal(result[1, 1], log(0.3 / 0.2))
  expect_equal(result[1, 2], log(0.5 / 0.2))
})

test_that("invalr is inverse of alr", {
  # Generate random compositions
  set.seed(42)
  compo <- matrix(runif(15), nrow = 5, ncol = 3)
  compo <- compo / rowSums(compo)

  # Transform and back
  transformed <- alr(compo)
  recovered <- invalr(transformed)

  expect_equal(compo, recovered, tolerance = 1e-10)
})

test_that("invalr produces valid compositions", {
  # Test with extreme values
  x <- matrix(c(-5, 5, 0, 0, -10, 10), nrow = 3, ncol = 2)
  result <- invalr(x)

  # All values should be in (0, 1)
  expect_true(all(result > 0))
  expect_true(all(result < 1))

  # Rows should sum to 1
  expect_equal(rowSums(result), rep(1, 3), tolerance = 1e-10)
})

test_that("dalr has correct dimensions", {
  compo <- matrix(c(0.3, 0.5, 0.2,
                    0.4, 0.4, 0.2), nrow = 2, byrow = TRUE)

  J <- dalr(compo)

  # Should be (n*d) x D = (2*2) x 3 = 4 x 3
  expect_equal(dim(J), c(4, 3))
})

test_that("dalr satisfies derivative identities", {
  # For compositional data z with d = D-1
  # dz_k/dx_j = z_k(1 - z_k) if k = j
  # dz_k/dx_j = -z_k * z_j if k != j

  z <- matrix(c(0.3, 0.5, 0.2), nrow = 1)
  J <- dalr(z)

  # Check specific values
  # J[1,1] = dz_1/dx_1 = z_1(1-z_1) = 0.3 * 0.7 = 0.21
  expect_equal(J[1, 1], 0.3 * 0.7, tolerance = 1e-10)

  # J[1,2] = dz_1/dx_2 = -z_1 * z_2 = -0.3 * 0.5 = -0.15
  expect_equal(J[1, 2], -0.3 * 0.5, tolerance = 1e-10)
})

test_that("compo_dist is symmetric and zero for identical compositions", {
  x <- c(0.3, 0.5, 0.2)
  y <- c(0.4, 0.4, 0.2)

  # Distance to self should be 0
  expect_equal(compo_dist(x, x), 0, tolerance = 1e-10)

  # Symmetry
  expect_equal(compo_dist(x, y), compo_dist(y, x), tolerance = 1e-10)
})

test_that("no_zero_one handles boundary values", {
  compo <- matrix(c(0, 0.5, 0.5,
                    1, 0, 0,
                    0.3, 0.3, 0.4), nrow = 3, byrow = TRUE)

  result <- no_zero_one(compo)

  # No zeros or ones
  expect_true(all(result > 0))
  expect_true(all(result < 1))

  # Still sums to 1
  expect_equal(rowSums(result), rep(1, 3), tolerance = 1e-10)
})

test_that("alr handles vector input", {
  compo <- c(0.3, 0.5, 0.2)
  result <- alr(compo)

  expect_equal(dim(result), c(1, 2))
})

test_that("invalr handles vector input", {
  x <- c(0.5, -0.5)
  result <- invalr(x)

  expect_equal(dim(result), c(1, 3))
  expect_equal(sum(result), 1, tolerance = 1e-10)
})
