# Tests for the adjudicated co-ranking implementation (Phase B, Stage 2).
# Definitions per Lee & Verleysen (2009), Neurocomputing 72:1431-1443,
# as adjudicated in Octopus_Workflow/reports/REVISION_WORKLOG.md (Stage 1).

tiefree_data <- function(n, p, seed) {
  set.seed(seed)
  matrix(rnorm(n * p), n, p)
}

# Independent first-principles recomputation of all six indices directly from
# rank pairs (no co-ranking matrix), used as an oracle against the package.
oracle_indices <- function(Dh, Dl, K) {
  N <- nrow(Dh)
  rk <- function(D) {
    R <- matrix(0L, N, N)
    for (i in seq_len(N)) R[i, -i] <- rank(D[i, -i], ties.method = "first")
    R
  }
  Rh <- rk(Dh); Rl <- rk(Dl)
  G <- if (K < N / 2) N * K * (2 * N - 3 * K - 1) else N * (N - K) * (N - K - 1)
  H <- N * sum(abs(N - 2 * seq_len(K) + 1) / seq_len(K))
  intr <- extr <- Wn <- Wv <- ul <- un <- ux <- 0
  for (i in seq_len(N)) for (j in seq_len(N)) {
    if (i == j) next
    k <- Rh[i, j]; l <- Rl[i, j]
    if (k > K && l <= K) intr <- intr + (k - K)
    if (k <= K && l > K) extr <- extr + (l - K)
    if (l <= K) Wn <- Wn + abs(k - l) / l
    if (k <= K) Wv <- Wv + abs(k - l) / k
    if (k <= K && l <= K) {
      ul <- ul + 1
      if (l < k) un <- un + 1
      if (l > k) ux <- ux + 1
    }
  }
  M_T <- 1 - 2 * intr / G; M_C <- 1 - 2 * extr / G
  c(Q_TC = (M_T + M_C) / 2, B_TC = M_C - M_T,
    Q_RE = 1 - (Wn / H + Wv / H) / 2, B_RE = Wn / H - Wv / H,
    Q_LC = ul / (N * K), B_LC = ux / (N * K) - un / (N * K))
}

test_that("coranking_indices returns the six named indices", {
  X <- tiefree_data(10, 5, 42); Y <- tiefree_data(10, 2, 43)
  result <- coranking_indices(as.matrix(dist(X)), as.matrix(dist(Y)), K = 3)
  expect_length(result, 6)
  expect_named(result, c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC"))
})

test_that("identity embedding: all Q exactly 1, all B exactly 0", {
  X <- tiefree_data(20, 4, 1)
  D <- as.matrix(dist(X))
  for (K in c(1, 5, 10, 18)) {
    v <- coranking_indices(D, D, K)
    expect_equal(unname(v[c("Q_TC", "Q_RE", "Q_LC")]), c(1, 1, 1))
    expect_equal(unname(v[c("B_TC", "B_RE", "B_LC")]), c(0, 0, 0))
  }
})

test_that("package indices equal the first-principles oracle", {
  X <- tiefree_data(25, 5, 2)
  Y <- tiefree_data(25, 2, 3)
  Dh <- as.matrix(dist(X)); Dl <- as.matrix(dist(Y))
  for (K in c(2, 7, 12, 20)) {
    expect_equal(unname(coranking_indices(Dh, Dl, K)[1:6]),
                 unname(oracle_indices(Dh, Dl, K)),
                 tolerance = 1e-12)
  }
})

test_that("co-ranking mass and margin invariants hold", {
  X <- tiefree_data(15, 3, 4); Y <- tiefree_data(15, 2, 5)
  Rh <- rank_matrix(as.matrix(dist(X)))
  Rl <- rank_matrix(as.matrix(dist(Y)))
  Q <- coranking_matrix(Rh, Rl)
  N <- 15
  expect_equal(sum(Q), N * (N - 1))
  expect_equal(unname(rowSums(Q)), rep(N, N - 1))
  expect_equal(unname(colSums(Q)), rep(N, N - 1))
})

test_that("self-exclusion: duplicate observations do not lose mass", {
  # The pre-revision implementation ranked the diagonal and dropped pairs for
  # duplicates, so a perfect embedding could score Q_LC < 1.
  X <- tiefree_data(12, 3, 6)
  X[2, ] <- X[1, ]                    # exact duplicate pair
  D <- as.matrix(dist(X))
  # duplicates create off-diagonal zero-distance ties -> primary policy errors
  expect_error(rank_matrix(D), "tie")
  # seeded fallback keeps the full mass invariant and perfect self-agreement
  set.seed(99)
  R1 <- rank_matrix(D, ties = "random")
  Q <- coranking_matrix(R1, R1)
  expect_equal(sum(Q), 12 * 11)
  expect_equal(unname(indices_from_coranking(Q, 4)[c("Q_TC", "Q_RE", "Q_LC")]),
               c(1, 1, 1))
})

test_that("tie policy: error is default, random fallback is seeded", {
  D <- matrix(c(0, 1, 1, 2,
                1, 0, 2, 1,
                1, 2, 0, 1,
                2, 1, 1, 0), 4, 4)
  expect_error(rank_matrix(D), "tie")
  set.seed(7); Ra <- rank_matrix(D, ties = "random")
  set.seed(7); Rb <- rank_matrix(D, ties = "random")
  expect_identical(Ra, Rb)
})

test_that("degenerate intervals reduce every dissimilarity to point Euclidean", {
  X <- tiefree_data(18, 4, 8)
  R0 <- matrix(0, 18, 4)
  for (met in c("Int-Euclidean", "Hausdorff", "Ichino-Yaguchi", "Wasserstein")) {
    Dm <- idist(X, R0, met)
    expect_equal(unname(Dm), unname(as.matrix(dist(X))), tolerance = 1e-10,
                 label = met)
  }
})

test_that("null expectations match adjudicated formulas (D.2-1)", {
  skip_on_cran()
  set.seed(11)
  N <- 40; K <- 5; m <- 4000
  Rh <- rank_matrix(as.matrix(dist(tiefree_data(N, 3, 12))))
  Rl <- rank_matrix(as.matrix(dist(tiefree_data(N, 2, 13))))
  nulls <- perm_null_indices(Rh, Rl, K, m = m)
  se <- apply(nulls, 2, sd) / sqrt(m)
  expect_lt(abs(mean(nulls[, "Q_LC"]) - K / (N - 1)), 3 * se["Q_LC"])
  expect_lt(abs(mean(nulls[, "B_LC"])), 3 * se["B_LC"])
  expect_lt(abs(mean(nulls[, "B_TC"])), 3 * se["B_TC"])
  expect_lt(abs(mean(nulls[, "B_RE"])), 3 * se["B_RE"])
})

test_that("B_LC sign convention: positive for intrusive, negative for extrusive", {
  # Strong canonical intrusion: five well-separated clusters superimposed in
  # the embedding (cross-cluster pairs become severe false neighbors).
  # B_LC (= Lee-Verleysen's B_NX) is a mild-error tendency and is validated
  # here on strong distortions; near zero it is intentionally noisy.
  set.seed(1)
  K <- 8; n_per <- 20; n_cl <- 5
  offs <- matrix(rnorm(n_cl * 2, sd = 30), n_cl, 2)
  X <- do.call(rbind, lapply(seq_len(n_cl), function(g)
    sweep(matrix(rnorm(n_per * 2), n_per, 2), 2, offs[g, ], "+")))
  Y <- X - offs[rep(seq_len(n_cl), each = n_per), ]
  Y <- Y + matrix(rnorm(nrow(Y) * 2, sd = 1e-6), nrow(Y))
  v1 <- coranking_indices(as.matrix(dist(X)), as.matrix(dist(Y)), K)
  expect_gt(v1[["B_TC"]], 0); expect_gt(v1[["B_LC"]], 0.05)

  # Extrusive: tear a circle open to a line.
  set.seed(3)
  th <- sort(runif(100, 0, 2 * pi))
  v2 <- coranking_indices(as.matrix(dist(cbind(cos(th), sin(th)))),
                          as.matrix(dist(th)), K)
  expect_lt(v2[["B_TC"]], 0); expect_lt(v2[["B_LC"]], 0)
})

test_that("B_LC respects its tight range (K-1)/K without clamping", {
  X <- tiefree_data(16, 3, 31); Y <- tiefree_data(16, 2, 32)
  for (K in c(2, 5, 10)) {
    v <- coranking_indices(as.matrix(dist(X)), as.matrix(dist(Y)), K)
    expect_lte(abs(v[["B_LC"]]), (K - 1) / K + 1e-12)
  }
})

test_that("fast permutation engine equals the reference implementation", {
  X <- tiefree_data(30, 4, 51); Y <- tiefree_data(30, 2, 52)
  Rh <- rank_matrix(as.matrix(dist(X)))
  Rl <- rank_matrix(as.matrix(dist(Y)))
  for (K in c(3, 10)) {
    set.seed(777)
    fast <- perm_null_indices(Rh, Rl, K, m = 25)
    set.seed(777)
    slow <- perm_null_indices_reference(Rh, Rl, K, m = 25)
    expect_equal(fast, slow, tolerance = 1e-12)
  }
})

test_that("min-P families share joint permutation draws across metrics", {
  # Explicit perms are honored deterministically
  X <- tiefree_data(20, 3, 61); Y <- tiefree_data(20, 2, 62)
  Rh <- rank_matrix(as.matrix(dist(X)))
  Rl <- rank_matrix(as.matrix(dist(Y)))
  perms <- t(vapply(1:15, function(s) sample.int(20), integer(20)))
  a <- perm_null_indices(Rh, Rl, K = 4, perms = perms)
  b <- perm_null_indices(Rh, Rl, K = 4, perms = perms)
  expect_identical(a, b)

  # In assess_quality, two metrics of the same family must be evaluated
  # under the same draws: with an interval geometry whose four dissimilarities
  # induce the same rank system (radii constant), the null draws must agree
  # exactly across metric slices.
  set.seed(63)
  C <- matrix(rnorm(18 * 3), 18, 3)
  x <- interval_data(C, matrix(0.1, 18, 3))
  proj <- structure(list(mock = list(C = C[, 1:2], R = matrix(0.1, 18, 2),
                                     type = "Interval")),
                    class = "idr_projections")
  res <- assess_quality(x, proj, K = 4, metrics = c("Hausdorff", "Wasserstein"),
                        baseline = FALSE, perm_test = TRUE, n_perm = 25,
                        seed = 99)
  nd <- res$null_draws$mock
  expect_equal(nd[, "Hausdorff", ], nd[, "Wasserstein", ], tolerance = 1e-12)
})

test_that("invalid inputs are rejected", {
  X <- tiefree_data(10, 3, 41)
  D <- as.matrix(dist(X))
  expect_error(rank_matrix(D[, -1]), "square")
  Dn <- D; Dn[1, 2] <- Dn[2, 1] <- -1
  expect_error(rank_matrix(Dn), "negative")
  Df <- D; Df[1, 2] <- Df[2, 1] <- NA
  expect_error(rank_matrix(Df), "non-finite")
  Dd <- D; diag(Dd)[1] <- 0.5
  expect_error(rank_matrix(Dd), "diagonal")
  Ds <- D; Ds[1, 2] <- Ds[1, 2] + 1
  expect_error(rank_matrix(Ds), "symmetric")
  Q <- coranking_matrix(rank_matrix(D), rank_matrix(D))
  expect_error(indices_from_coranking(Q, 9), "K")
})

# ---- Prompt E (2026-07-17): tolerance-block randomization regression set ----
# Tie blocks detected within the relative tolerance must be randomized as
# whole blocks; base R's random tie method randomizes only exact equals.

sym_D4 <- function(d12, d13, d14, d23 = 3, d24 = 5, d34 = 7) {
  D <- matrix(0, 4, 4)
  D[upper.tri(D)] <- c(d12, d13, d23, d14, d24, d34)
  D + t(D)
}

test_that("error policy catches a default-tolerance near tie", {
  D <- sym_D4(1, 1 + 5e-13, 9)
  expect_error(rank_matrix(D), "tie")
})

test_that("near-tie randomization reaches both orders and is seeded", {
  D <- sym_D4(1, 1 + 5e-13, 9)
  ord <- vapply(1:32, function(seed) {
    set.seed(seed)
    R <- rank_matrix(D, ties = "random")
    sign(R[1, 2] - R[1, 3])
  }, numeric(1))
  expect_setequal(unique(ord), c(-1, 1))
  set.seed(17); a <- rank_matrix(D, ties = "random")
  set.seed(17); b <- rank_matrix(D, ties = "random")
  expect_identical(a, b)
})

test_that("a tolerance-connected triple is randomized as one block", {
  # Endpoints differ by 1.5e-12 (> tol), but both adjacent gaps are links.
  D <- sym_D4(1, 1 + 0.75e-12, 1 + 1.50e-12)
  perms <- vapply(1:128, function(seed) {
    set.seed(seed)
    paste(rank_matrix(D, ties = "random")[1, 2:4], collapse = "")
  }, character(1))
  expect_setequal(unique(perms),
                  c("123", "132", "213", "231", "312", "321"))
  # The far endpoints can reverse, proving the chain was not split.
  expect_true(any(vapply(1:128, function(seed) {
    set.seed(seed)
    R <- rank_matrix(D, ties = "random")
    R[1, 2] > R[1, 4]
  }, logical(1))))
})

test_that("exact ties still randomize", {
  D <- sym_D4(1, 1, 9)
  ord <- vapply(1:32, function(seed) {
    set.seed(seed)
    R <- rank_matrix(D, ties = "random")
    sign(R[1, 2] - R[1, 3])
  }, numeric(1))
  expect_setequal(unique(ord), c(-1, 1))
})

test_that("tie-free ranks are bit-for-bit pre-patch and use no RNG", {
  set.seed(901)
  D <- as.matrix(dist(matrix(rnorm(40), 10, 4)))
  old <- matrix(0L, 10, 10)
  for (i in seq_len(10)) {
    old[i, -i] <- as.integer(rank(D[i, -i], ties.method = "first"))
  }
  attr(old, "n_ties") <- 0L

  set.seed(2718); before <- .Random.seed
  got <- rank_matrix(D, ties = "random")
  after <- .Random.seed
  expect_identical(got, old)
  expect_identical(after, before)
})

test_that("every row is a valid self-excluded rank permutation under ties", {
  D <- sym_D4(1, 1 + 5e-13, 9)
  set.seed(29)
  R <- rank_matrix(D, ties = "random")
  expect_identical(diag(R), rep(0L, 4))
  for (i in seq_len(4)) expect_identical(sort(R[i, -i]), 1:3)
})

test_that("the IY tangency is exceptional under both tie policies", {
  # Prompt D/E regression: d_IY(O,A)^2 - d_IY(O,B)^2 = (1 - 5 nu)^2,
  # double root at nu = 0.2 (isolated tie without sign reversal).
  centers <- rbind(c(0, 0), c(6, sqrt(45) / 2),
                   c(3 * sqrt(80) / 4, 0), c(50, 70))
  radii <- rbind(c(0, 0), c(0, sqrt(45) / 2),
                 c(sqrt(80) / 4, 0), c(0, 0))
  D <- idist_ichino_yaguchi(centers, radii, nu = 0.2)
  expect_error(rank_matrix(D, ties = "error"), "tie")
  ord <- vapply(1:32, function(seed) {
    set.seed(seed)
    R <- rank_matrix(D, ties = "random")
    sign(R[1, 2] - R[1, 3])
  }, numeric(1))
  expect_setequal(unique(ord), c(-1, 1))
})

test_that("n = 3 is supported end-to-end as Theorem thm:nullcal permits", {
  X <- matrix(c(0, 0, 1, 0, 0.3, 0.9), 3, 2, byrow = TRUE)
  D <- as.matrix(dist(X))
  R <- rank_matrix(D)
  expect_identical(diag(R), rep(0L, 3))
  for (i in 1:3) expect_identical(sort(R[i, -i]), 1:2)
  v <- coranking_indices(D, D, K = 1)   # K = 1 = n - 2
  expect_equal(unname(v[c("Q_TC", "Q_RE", "Q_LC")]), c(1, 1, 1))
  expect_equal(unname(v[c("B_TC", "B_RE", "B_LC")]), c(0, 0, 0))
})
