test_that("assess_quality returns qaidr_assessment object", {
  set.seed(42)
  C <- matrix(rnorm(50), 10, 5)
  R <- matrix(runif(50, 0, 0.5), 10, 5)
  x <- interval_data(C, R)

  # Create mock projections
  proj <- structure(
    list(
      "mock" = list(
        C = matrix(rnorm(20), 10, 2),
        R = matrix(0, 10, 2),
        type = "Point"
      )
    ),
    class = "idr_projections"
  )

  result <- assess_quality(x, proj, K = 3,
                           metrics = "Wasserstein", perm_test = FALSE)

  expect_s3_class(result, "qaidr_assessment")
  expect_true(is.data.frame(result$results))
  expect_equal(nrow(result$results), 1)
  expect_null(result$pvalues)
})

test_that("assess_quality with perm_test returns pvalues", {
  set.seed(42)
  C <- matrix(rnorm(50), 10, 5)
  R <- matrix(runif(50, 0, 0.5), 10, 5)
  x <- interval_data(C, R)

  proj <- structure(
    list(
      "mock" = list(
        C = matrix(rnorm(20), 10, 2),
        R = matrix(0, 10, 2),
        type = "Point"
      )
    ),
    class = "idr_projections"
  )

  result <- assess_quality(x, proj, K = 3,
                           metrics = "Wasserstein",
                           perm_test = TRUE, n_perm = 19)

  expect_s3_class(result, "qaidr_assessment")
  expect_true(is.data.frame(result$pvalues))
  expect_equal(nrow(result$pvalues), 1)
})

test_that("perm_test returns correct structure", {
  set.seed(42)
  Dh <- as.matrix(dist(matrix(rnorm(50), 10, 5)))
  Dl <- as.matrix(dist(matrix(rnorm(20), 10, 2)))

  pt <- perm_test(Dh, Dl, K = 3, n_perm = 19)

  expect_type(pt, "list")
  expect_named(pt, c("vals", "pQ", "pB"))
  expect_length(pt$vals, 6)
  expect_length(pt$pQ, 3)
  expect_length(pt$pB, 3)
  # P-values in [0, 1]
  expect_true(all(pt$pQ >= 0 & pt$pQ <= 1))
  expect_true(all(pt$pB >= 0 & pt$pB <= 1))
})

test_that("k_profiles returns correct data frame", {
  set.seed(42)
  C <- matrix(rnorm(50), 10, 5)
  R <- matrix(runif(50, 0, 0.5), 10, 5)
  x <- interval_data(C, R)

  proj <- structure(
    list(
      "mock" = list(
        C = matrix(rnorm(20), 10, 2),
        R = matrix(0, 10, 2),
        type = "Point"
      )
    ),
    class = "idr_projections"
  )

  profiles <- k_profiles(x, proj, K_max = 3, metrics = "Wasserstein")

  expect_true(is.data.frame(profiles))
  expect_equal(nrow(profiles), 3)  # K_max = 3 * 1 method * 1 metric
  expect_true(all(c("Method", "Metric", "K", "Q_TC", "B_TC") %in% names(profiles)))
})

test_that("print.qaidr_assessment works", {
  result <- structure(
    list(
      results = data.frame(
        IDR = "test", Metric = "Wasserstein",
        Q_TC = 0.5, B_TC = 0.1,
        Q_RE = 0.6, B_RE = -0.1,
        Q_LC = 0.7, B_LC = 0.05
      ),
      pvalues = NULL,
      K = 5
    ),
    class = "qaidr_assessment"
  )

  expect_output(print(result), "QAIDR Assessment")
})
