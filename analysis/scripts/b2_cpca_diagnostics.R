# B2 C-PCA mechanism diagnostics (preregistered A-D4, frozen 2026-07-16).
# Scenario I replications 1..25 re-derived from original seeds. C-PCA only.
# W recomputed with the exact symbolicDA expression W = eigen(cor(centers))
# $vectors[,1:2]; implementation identity C_out = C_in W, R_out = R_in |W|
# checked by max_abs_error (subspace/sign handled by comparing against the
# actually fitted projection via the wrapper output, using |.| invariants).
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.R"))
save_sessioninfo("b2_cpca_diagnostics")

N_SUB <- 25L
pseudoF <- function(v, g) {
  gm <- tapply(v, g, mean); n_g <- tapply(v, g, length); mu <- mean(v)
  b <- sum(n_g * (gm[levels(g)] - mu)^2) / (nlevels(g) - 1)
  w <- sum((v - gm[g])^2) / (length(v) - nlevels(g))
  b / w
}

rows <- list()
for (r in seq_len(N_SUB)) {
  seed_r <- SEED_BASE + 1000L + r
  x <- gen_scenario1(seed_r); xs <- standardize(x)
  g <- x$labels
  set.seed(seed_r + 500000L)
  rm <- run_methods_timed(xs, methods = "C-PCA")
  pr <- rm$projections[["C-PCA"]]

  ## exact recomputation of W per the symbolicDA source
  W <- eigen(cor(xs$centers))$vectors[, 1:2]
  C_pred <- xs$centers %*% W
  R_pred <- xs$radii %*% abs(W)
  ## subspace/sign handling: align predicted to fitted via per-column sign
  sgn <- sign(colSums(C_pred * pr$C)); sgn[sgn == 0] <- 1
  err_C <- max(abs(sweep(C_pred, 2, sgn, "*") - pr$C))
  err_R <- max(abs(R_pred - pr$R))          # radii are sign-invariant

  ## T1: association of |w_st| with the group-informative radius profile
  ## (per-variable mean radius correlates with |w| columns?)
  rad_profile <- colMeans(xs$radii)
  t1 <- max(abs(cor(abs(W), rad_profile)))
  ## T2: paired separation, projected radii vs projected centers
  f_rad <- pseudoF(rowMeans(pr$R), g)
  f_cen <- pseudoF(rowMeans(pr$C), g)
  ## T3: top-2 eigenvalue share of cor(centers) + isotropic-null reference
  ev <- eigen(cor(xs$centers), only.values = TRUE)$values
  share <- sum(ev[1:2]) / sum(ev)
  set.seed(SEED_BASE + 7900L + r)
  evn <- eigen(cor(matrix(rnorm(300 * 5), 300, 5)), only.values = TRUE)$values
  share_null <- sum(evn[1:2]) / sum(evn)

  rows[[r]] <- data.frame(rep = r, err_C = err_C, err_R = err_R, t1_maxcor = t1,
                          F_radii = f_rad, F_centers = f_cen,
                          F_diff = f_rad - f_cen,
                          eig_share = share, eig_share_null = share_null)
  if (r %% 5 == 0) cat("rep", r, "/", N_SUB, "\n")
}
d <- do.call(rbind, rows)
crit <- list(
  identity_max_err = max(d$err_C, d$err_R),
  t1_interval = quantile(d$t1_maxcor, c(0.025, 0.975)),
  F_diff_q025 = quantile(d$F_diff, 0.025),
  criterion_T2_pass = unname(quantile(d$F_diff, 0.025) > 0),
  eig_share_range = range(d$eig_share),
  eig_share_null_range = range(d$eig_share_null)
)
save_result(d, "b2_cpca_diagnostics",
            extra = c(list(prereg_sha256 = "2c25e3b2d39ee824043fdf5261e299f13773c5a7734799f7f35ba6257f2df67a"),
                      lapply(crit, unname)))
cat(sprintf("identity max err: C %.2e R %.2e | T1 max|cor| 95%% int [%.3f, %.3f] | F_diff q025 %.1f (pass=%s)\n",
            max(d$err_C), max(d$err_R), crit$t1_interval[1], crit$t1_interval[2],
            crit$F_diff_q025, crit$criterion_T2_pass))
cat(sprintf("eig share: data [%.3f, %.3f] vs isotropic null [%.3f, %.3f]\n",
            crit$eig_share_range[1], crit$eig_share_range[2],
            crit$eig_share_null_range[1], crit$eig_share_null_range[2]))
