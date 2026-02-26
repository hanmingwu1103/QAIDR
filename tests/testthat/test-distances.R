test_that("idist_euclidean returns correct shape", {
  C <- matrix(rnorm(12), 4, 3)
  R <- matrix(runif(12, 0.1, 0.5), 4, 3)
  D <- idist_euclidean(C, R)

  expect_true(is.matrix(D))
  expect_equal(dim(D), c(4, 4))
  expect_equal(diag(D), rep(0, 4))
  # Symmetric

  expect_equal(D, t(D))
})

test_that("idist_hausdorff returns correct shape", {
  C <- matrix(rnorm(12), 4, 3)
  R <- matrix(runif(12, 0.1, 0.5), 4, 3)
  D <- idist_hausdorff(C, R)

  expect_equal(dim(D), c(4, 4))
  expect_equal(D, t(D))
  expect_true(all(D >= 0))
})

test_that("idist_ichino_yaguchi returns correct shape", {
  C <- matrix(rnorm(12), 4, 3)
  R <- matrix(runif(12, 0.1, 0.5), 4, 3)
  D <- idist_ichino_yaguchi(C, R)

  expect_equal(dim(D), c(4, 4))
  expect_equal(D, t(D))
  expect_true(all(D >= 0))
})

test_that("idist_wasserstein returns correct shape", {
  C <- matrix(rnorm(12), 4, 3)
  R <- matrix(runif(12, 0.1, 0.5), 4, 3)
  D <- idist_wasserstein(C, R)

  expect_equal(dim(D), c(4, 4))
  expect_equal(D, t(D))
  expect_true(all(D >= 0))
})

test_that("idist_wasserstein with zero radii equals Euclidean", {
  C <- matrix(c(0, 1, 0, 0, 0, 1), nrow = 3, ncol = 2)
  R <- matrix(0, nrow = 3, ncol = 2)

  D_wass <- idist_wasserstein(C, R)
  D_euc <- as.matrix(dist(C))

  expect_equal(D_wass, D_euc, tolerance = 1e-10)
})

test_that("idist dispatcher works for all metrics", {
  C <- matrix(rnorm(12), 4, 3)
  R <- matrix(runif(12, 0.1, 0.5), 4, 3)

  for (met in c("Int-Euclidean", "Hausdorff", "Ichino-Yaguchi", "Wasserstein")) {
    D <- idist(C, R, met)
    expect_equal(dim(D), c(4, 4))
  }
})

test_that("idist works with interval_data objects", {
  C <- matrix(rnorm(12), 4, 3)
  R <- matrix(runif(12, 0.1, 0.5), 4, 3)
  x <- interval_data(C, R)

  D <- idist(x, metric = "Wasserstein")
  expect_equal(dim(D), c(4, 4))
})

test_that("idist errors on unknown metric", {
  C <- matrix(1:4, 2, 2)
  R <- matrix(0.1, 2, 2)
  expect_error(idist(C, R, "FooBar"), "Unknown metric")
})

test_that("distance functions are non-negative", {
  set.seed(123)
  C <- matrix(rnorm(30), 10, 3)
  R <- matrix(runif(30, 0, 1), 10, 3)

  expect_true(all(idist_euclidean(C, R) >= 0))
  expect_true(all(idist_hausdorff(C, R) >= 0))
  expect_true(all(idist_ichino_yaguchi(C, R) >= 0))
  expect_true(all(idist_wasserstein(C, R) >= 0))
})
