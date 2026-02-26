test_that("coranking_indices returns named vector of correct length", {
  set.seed(42)
  X <- matrix(rnorm(50), 10, 5)
  Y <- matrix(rnorm(20), 10, 2)
  Dh <- as.matrix(dist(X))
  Dl <- as.matrix(dist(Y))

  result <- coranking_indices(Dh, Dl, K = 3)

  expect_length(result, 6)
  expect_named(result, c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC"))
})

test_that("quality indices are in [0, 1]", {
  set.seed(42)
  X <- matrix(rnorm(50), 10, 5)
  Y <- matrix(rnorm(20), 10, 2)
  Dh <- as.matrix(dist(X))
  Dl <- as.matrix(dist(Y))

  result <- coranking_indices(Dh, Dl, K = 3)

  expect_true(result["Q_TC"] >= 0 && result["Q_TC"] <= 1)
  expect_true(result["Q_RE"] >= 0 && result["Q_RE"] <= 1)
  expect_true(result["Q_LC"] >= 0 && result["Q_LC"] <= 1)
})

test_that("behavior indices are in [-1, 1]", {
  set.seed(42)
  X <- matrix(rnorm(50), 10, 5)
  Y <- matrix(rnorm(20), 10, 2)
  Dh <- as.matrix(dist(X))
  Dl <- as.matrix(dist(Y))

  result <- coranking_indices(Dh, Dl, K = 3)

  expect_true(result["B_TC"] >= -1 && result["B_TC"] <= 1)
  expect_true(result["B_RE"] >= -1 && result["B_RE"] <= 1)
  expect_true(result["B_LC"] >= -1 && result["B_LC"] <= 1)
})

test_that("perfect embedding gives Q = 1", {
  set.seed(42)
  X <- matrix(rnorm(30), 10, 3)
  Dh <- as.matrix(dist(X))

  # Same distances => perfect preservation

  result <- coranking_indices(Dh, Dh, K = 3)

  expect_equal(unname(result["Q_TC"]), 1, tolerance = 1e-10)
  expect_equal(unname(result["Q_RE"]), 1, tolerance = 1e-10)
  expect_equal(unname(result["Q_LC"]), 1, tolerance = 1e-10)
})

test_that("coranking_indices rejects invalid K", {
  Dh <- matrix(0, 5, 5)
  Dl <- matrix(0, 5, 5)
  expect_error(coranking_indices(Dh, Dl, K = 0))
  expect_error(coranking_indices(Dh, Dl, K = 4))
})

test_that("coranking_indices works with K = 1", {
  set.seed(42)
  X <- matrix(rnorm(50), 10, 5)
  Y <- matrix(rnorm(20), 10, 2)
  Dh <- as.matrix(dist(X))
  Dl <- as.matrix(dist(Y))

  result <- coranking_indices(Dh, Dl, K = 1)
  expect_length(result, 6)
})
