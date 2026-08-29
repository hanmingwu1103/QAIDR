# Oracle tests for the symmetric randomization min-P adjustment (Prompt H,
# contract 1.1/1.3). The literal triple-loop oracle below is written directly
# from the frozen algorithm; the numeric fixtures in minp_oracle_vectors.rds
# were produced by an independent implementation (OpenAI Codex seat) and are
# consumed here as exact rational expectations (integer numerators over B).

vec <- readRDS(test_path("minp_oracle_vectors.rds"))

# Literal contract oracle: p[b,h] = (1/B) sum_r 1{T[r,h] >= T[b,h]},
# M[b] = min_h p[b,h], padj[h] = (1/B) sum_b 1{M[b] <= p[1,h]} (row 1 observed).
oracle_minp <- function(T) {
  B <- nrow(T); H <- ncol(T)
  p <- matrix(NA_real_, B, H)
  for (b in seq_len(B)) for (h in seq_len(H)) {
    cnt <- 0L
    for (r in seq_len(B)) if (T[r, h] >= T[b, h]) cnt <- cnt + 1L
    p[b, h] <- cnt / B
  }
  M <- apply(p, 1, min)
  padj <- vapply(seq_len(H), function(h) sum(M <= p[1, h]) / B, numeric(1))
  list(p = p, M = M, padj = padj)
}

# Family-pooled literal oracle on production-shaped inputs: for each index
# family, pool {Q_F, B_F} across all G metric rows into one B x 2G array
# (|.| applied to behavior columns) and run the contract oracle.
oracle_p_minP <- function(obs_mat, null_arr) {
  m <- dim(null_arr)[1]
  fams <- list(c("Q_TC", "B_TC"), c("Q_RE", "B_RE"), c("Q_LC", "B_LC"))
  padj <- obs_mat; padj[] <- NA_real_
  for (fam in fams) {
    G <- nrow(obs_mat)
    T <- matrix(NA_real_, m + 1, G * 2)
    for (g in seq_len(G)) for (j in 1:2) {
      v <- c(obs_mat[g, fam[j]], null_arr[, g, fam[j]])
      if (startsWith(fam[j], "B")) v <- abs(v)
      T[, (g - 1L) * 2L + j] <- v
    }
    res <- oracle_minp(T)
    for (g in seq_len(G)) for (j in 1:2)
      padj[g, fam[j]] <- res$padj[(g - 1L) * 2L + j]
  }
  padj
}

test_that("literal oracle reproduces the independent random-array vectors exactly", {
  for (id in names(vec$random$T)) {
    got <- oracle_minp(vec$random$T[[id]])
    exp <- vec$random$expected[[id]]
    expect_identical(got$p * exp$B, exp$p_num, info = id)
    expect_identical(got$M * exp$B, exp$M_num, info = id)
    expect_identical(got$padj * exp$B, exp$padj_num, info = id)
  }
})

test_that("literal oracle reproduces the edge-case vectors (m=1, duplicates, all-equal, family sizes 1 and 4)", {
  for (id in names(vec$edges$T)) {
    got <- oracle_minp(vec$edges$T[[id]])
    exp <- vec$edges$expected[[id]]
    expect_identical(got$padj * exp$B, exp$padj_num, info = id)
  }
})

test_that("one-hypothesis reduction: family of size 1 returns the marginal symmetric p", {
  T <- vec$edges$T$family1_with_ties
  got <- oracle_minp(T)
  expect_identical(got$padj, got$p[1, 1])
})

test_that("complete-orbit enumeration controls weak FWER at every attainable alpha", {
  for (id in names(vec$complete_orbits)) {
    orb <- vec$complete_orbits[[id]]
    T <- orb$T; B <- nrow(T)
    min_padj <- numeric(B)
    for (s in seq_len(B)) {
      ord <- c(s, seq_len(B)[-s])
      min_padj[s] <- min(oracle_minp(T[ord, , drop = FALSE])$padj)
    }
    expect_equal(min_padj, orb$observed$min_padj, tolerance = 1e-12, info = id)
    for (k in seq_len(B)) {
      fwer <- mean(min_padj <= k / B)
      expect_lte(fwer, k / B + 1e-12)
      expect_equal(fwer, orb$fwer$fwer[orb$fwer$alpha_num == k],
                   tolerance = 1e-12, info = paste(id, "alpha", k))
    }
  }
})

test_that("Prompt G m=19 counterexample: old formula fails, symmetric formula controls", {
  ce <- vec$counterexample
  expect_gt(ce$old_formula_fwer, ce$alpha)          # 0.0883 > 0.05
  expect_lte(ce$new_formula_fwer, ce$alpha)         # 0.0108 <= 0.05
  expect_true(ce$new_controls)
  # Exhaustive in-test reproduction of the small labelled grid (q=3, B=4,
  # H=2 independent coordinates): enumerate all q^(B*H) arrays.
  q <- ce$small_grid_q; B <- ce$small_grid_B; H <- 2L
  grid <- as.matrix(expand.grid(rep(list(seq_len(q)), B * H)))
  hit <- 0L
  for (i in seq_len(nrow(grid))) {
    T <- matrix(grid[i, ], B, H)
    if (min(oracle_minp(T)$padj) <= 1 / B) hit <- hit + 1L
  }
  fwer <- hit / nrow(grid)
  expect_equal(fwer, ce$small_grid_enumerated_fwer, tolerance = 1e-12)
  expect_equal(fwer, ce$small_grid_formula_fwer, tolerance = 1e-12)
  expect_lte(fwer, 1 / B)
})

mk_qb <- function(m, G, seed, tied = TRUE) {
  set.seed(seed)
  idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  draw <- function(k) if (tied) sample(-3:5, k, replace = TRUE) / 4 else rnorm(k)
  obs <- matrix(draw(G * 6), G, 6, dimnames = list(NULL, idx))
  nul <- array(draw(m * G * 6), c(m, G, 6), dimnames = list(NULL, NULL, idx))
  list(obs = obs, nul = nul)
}

test_that("production .p_minP agrees exactly with the family-pooled literal oracle", {
  for (seed in c(11, 12, 13)) {
    d <- mk_qb(m = 7, G = 3, seed = seed, tied = TRUE)
    expect_equal(QAIDR:::.p_minP(d$obs, d$nul), oracle_p_minP(d$obs, d$nul),
                 tolerance = 1e-15)
  }
  d <- mk_qb(m = 19, G = 4, seed = 14, tied = FALSE)
  expect_equal(QAIDR:::.p_minP(d$obs, d$nul), oracle_p_minP(d$obs, d$nul),
               tolerance = 1e-15)
})

test_that("production .p_minP matches the independent production-adapter vector", {
  ad <- vec$production_adapter
  idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  G <- nrow(ad$obs_mat)
  obs <- matrix(NA_real_, G, 6, dimnames = list(NULL, idx))
  nul <- array(NA_real_, c(dim(ad$null_arr)[1], G, 6),
               dimnames = list(NULL, NULL, idx))
  # Adapter carries one two-column family; replicate it into all three
  # families (computed independently, so expectations are unchanged).
  for (f in c("TC", "RE", "LC")) {
    obs[, paste0("Q_", f)] <- ad$obs_mat[, 1]
    obs[, paste0("B_", f)] <- ad$obs_mat[, 2]
    nul[, , paste0("Q_", f)] <- ad$null_arr[, , 1]
    nul[, , paste0("B_", f)] <- ad$null_arr[, , 2]
  }
  padj <- QAIDR:::.p_minP(obs, nul)
  want <- ad$expected$padj             # labels g1_Q, g1_B, g2_Q, g2_B
  expect_equal(unname(padj[1, "Q_TC"]), unname(want[1]), tolerance = 1e-15)
  expect_equal(unname(padj[1, "B_TC"]), unname(want[2]), tolerance = 1e-15)
  expect_equal(unname(padj[2, "Q_TC"]), unname(want[3]), tolerance = 1e-15)
  expect_equal(unname(padj[2, "B_TC"]), unname(want[4]), tolerance = 1e-15)
})

test_that("observed symmetric p equals the plus-one marginal p and padj >= praw", {
  for (seed in c(21, 22)) {
    d <- mk_qb(m = 9, G = 3, seed = seed, tied = TRUE)
    praw <- QAIDR:::.p_unadjusted(d$obs, d$nul)
    padj <- QAIDR:::.p_minP(d$obs, d$nul)
    # p[0,h] = (1 + #{null >= obs})/(m+1) is exactly the symmetric row-0 p.
    for (g in seq_len(3)) for (col in colnames(d$obs)) {
      v <- c(d$obs[g, col], d$nul[, g, col])
      if (startsWith(col, "B")) v <- abs(v)
      expect_equal(sum(v >= v[1]) / length(v), unname(praw[g, col]),
                   tolerance = 1e-15)
    }
    expect_true(all(padj >= praw - 1e-15))
  }
})

test_that("quality uses the upper tail and behavior the absolute two-sided tail", {
  m <- 9; idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  obs <- matrix(0, 1, 6, dimnames = list(NULL, idx))
  nul <- array(0, c(m, 1, 6), dimnames = list(NULL, NULL, idx))
  nul[, 1, ] <- seq(0.1, 0.9, length.out = m)   # all null values in (0,1)
  obs[1, "Q_TC"] <- 1                            # top of upper tail
  obs[1, "B_TC"] <- -1                           # extreme by |.| only
  obs[1, c("Q_RE", "Q_LC")] <- -1                # bottom of upper tail
  obs[1, c("B_RE", "B_LC")] <- 0                 # center: not extreme
  padj <- QAIDR:::.p_minP(obs, nul)
  expect_equal(unname(padj[1, "Q_TC"]), 1 / (m + 1), tolerance = 1e-15)
  expect_equal(unname(padj[1, "B_TC"]), 1 / (m + 1), tolerance = 1e-15)
  expect_equal(unname(padj[1, "Q_RE"]), 1, tolerance = 1e-15)
  expect_equal(unname(padj[1, "B_RE"]), 1, tolerance = 1e-15)
})

test_that("metric-row relabeling equivariance", {
  d <- mk_qb(m = 7, G = 4, seed = 31, tied = TRUE)
  perm <- c(3, 1, 4, 2)
  p1 <- QAIDR:::.p_minP(d$obs, d$nul)
  p2 <- QAIDR:::.p_minP(d$obs[perm, , drop = FALSE],
                        d$nul[, perm, , drop = FALSE])
  expect_equal(p2, p1[perm, , drop = FALSE], tolerance = 1e-15)
})

test_that("deliberate-failure self-test: an exclusive-comparison variant is caught", {
  broken_minp <- function(T) {
    B <- nrow(T); H <- ncol(T)
    p <- matrix(NA_real_, B, H)
    for (b in seq_len(B)) for (h in seq_len(H))
      p[b, h] <- sum(T[, h] > T[b, h]) / B      # strict: WRONG
    M <- apply(p, 1, min)
    vapply(seq_len(H), function(h) sum(M <= p[1, h]) / B, numeric(1))
  }
  mismatch <- FALSE
  for (id in names(vec$random$T)) {
    exp <- vec$random$expected[[id]]
    if (!isTRUE(all.equal(broken_minp(vec$random$T[[id]]) * exp$B,
                          exp$padj_num, tolerance = 1e-12)))
      mismatch <- TRUE
  }
  expect_true(mismatch)
})
