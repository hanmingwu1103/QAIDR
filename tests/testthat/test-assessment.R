# Tests for assessment, permutation inference, and multiplicity adjustment
# (Phase B, Stage 2).

mk_interval <- function(n = 12, p = 4, seed = 42) {
  set.seed(seed)
  interval_data(matrix(rnorm(n * p), n, p),
                matrix(runif(n * p, 0.05, 0.5), n, p))
}

mk_proj <- function(n = 12, d = 2, seed = 43, type = "Point") {
  set.seed(seed)
  structure(list(mock = list(C = matrix(rnorm(n * d), n, d),
                             R = matrix(if (type == "Point") 0 else runif(n * d, 0.05, 0.3), n, d),
                             type = type)),
            class = "idr_projections")
}

test_that("assess_quality returns results incl. the center-only baseline row", {
  x <- mk_interval()
  res <- assess_quality(x, mk_proj(), K = 3, metrics = "Wasserstein",
                        baseline = TRUE)
  expect_s3_class(res, "qaidr_assessment")
  expect_setequal(res$results$Metric, c("Wasserstein", "Centers-Euclidean"))
  expect_null(res$pvalues)
})

test_that("assess_quality baseline can be disabled", {
  x <- mk_interval()
  res <- assess_quality(x, mk_proj(), K = 3, metrics = "Wasserstein",
                        baseline = FALSE)
  expect_equal(res$results$Metric, "Wasserstein")
})

test_that("permutation p-values respect the 1/(m+1) floor and agree across entry points", {
  x <- mk_interval(n = 14)
  proj <- mk_proj(n = 14)
  m <- 39
  res <- assess_quality(x, proj, K = 3, metrics = "Wasserstein",
                        baseline = FALSE, perm_test = TRUE, n_perm = m,
                        seed = 202)
  p_all <- unlist(res$pvalues[, c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")])
  expect_true(all(p_all >= 1 / (m + 1) - 1e-12))
  expect_true(all(p_all <= 1))
  padj <- unlist(res$pvalues_adj[, c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")])
  # single-step min-P adjustment can only increase p-values
  expect_true(all(padj >= p_all - 1e-12))

  # perm_test and assess_quality share the engine: same seeded stream, same p
  Dh <- idist(x$centers, x$radii, "Wasserstein")
  Dl <- as.matrix(dist(proj$mock$C))
  pt <- perm_test(Dh, Dl, K = 3, n_perm = m, seed = 202)
  expect_equal(unname(pt$pQ["TC"]), res$pvalues$Q_TC[1], tolerance = 1e-12)
  expect_equal(unname(pt$pB["LC"]), res$pvalues$B_LC[1], tolerance = 1e-12)
})

test_that("identity projection yields the minimal attainable p-value", {
  # Perfect correspondence: observed Q indices are at the maximum, so
  # p must be exactly 1/(m+1) for the quality indices.
  set.seed(77)
  n <- 15
  C <- matrix(rnorm(n * 2), n, 2)
  x <- interval_data(cbind(C, C), matrix(0.0, n, 4) + 1e-6)
  proj <- structure(list(ident = list(C = C, R = matrix(1e-6, n, 2),
                                      type = "Interval")),
                    class = "idr_projections")
  m <- 19
  res <- assess_quality(x, proj, K = 3, metrics = "Wasserstein",
                        baseline = FALSE, perm_test = TRUE, n_perm = m,
                        seed = 5)
  expect_equal(res$pvalues$Q_LC[1], 1 / (m + 1), tolerance = 1e-12)
})

test_that("lambda and nu are honored end-to-end", {
  x <- mk_interval(n = 10)
  D1 <- idist(x$centers, x$radii, "Int-Euclidean", lambda = 0)
  D2 <- idist(x$centers, x$radii, "Int-Euclidean", lambda = 1)
  expect_false(isTRUE(all.equal(D1, D2)))
  D3 <- idist(x$centers, x$radii, "Ichino-Yaguchi", nu = 0)
  D4 <- idist(x$centers, x$radii, "Ichino-Yaguchi", nu = 0.5)
  expect_false(isTRUE(all.equal(D3, D4)))
  expect_error(idist(x$centers, x$radii, "Int-Euclidean", lambda = 2))
  # canonical IY range is [0, 0.5]; values above it are rejected
  expect_error(idist(x$centers, x$radii, "Ichino-Yaguchi", nu = 0.75))
  # boundary values work
  expect_silent(idist(x$centers, x$radii, "Int-Euclidean", lambda = 0))
  expect_silent(idist(x$centers, x$radii, "Ichino-Yaguchi", nu = 0.5))
})

test_that("k_profiles evaluates all K on one co-ranking matrix", {
  x <- mk_interval()
  profiles <- k_profiles(x, mk_proj(), K_max = 4, metrics = "Wasserstein",
                         baseline = FALSE)
  expect_true(is.data.frame(profiles))
  expect_equal(nrow(profiles), 4)
  expect_true(all(c("Method", "Metric", "K", "Q_TC", "B_LC") %in% names(profiles)))
  # K = N-2 upper bound enforced
  expect_error(k_profiles(x, mk_proj(), K_max = 11, metrics = "Wasserstein"),
               "K")
})

test_that("print.qaidr_assessment works with and without p-values", {
  x <- mk_interval()
  res <- assess_quality(x, mk_proj(), K = 3, metrics = "Wasserstein",
                        baseline = FALSE)
  expect_output(print(res), "QAIDR Assessment")
  res2 <- assess_quality(x, mk_proj(), K = 3, metrics = "Wasserstein",
                         baseline = FALSE, perm_test = TRUE, n_perm = 19,
                         seed = 1)
  expect_output(print(res2), "min-P")
})
