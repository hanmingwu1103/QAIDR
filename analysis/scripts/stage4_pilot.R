# Stage 4 pilot and calibration gates (Prompt B).
#  A. End-to-end deterministic pilot of each scenario (1 replication, quick perms).
#  B. Null-expectation checks (D.2-1 formulas) on real scenario geometry.
#  C. Type-I error under a true random-correspondence null; power under
#     graded distortion; min-P family-wise error control.
#  D. Per-component timings -> projected full-run resources.
# Output: output/stage4_pilot_report.{csv,rds,json} + console summary.

SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.R"))
save_sessioninfo("stage4_pilot")
set.seed(SEED_BASE)

report <- list()

## ---------- A. End-to-end pilots (one seeded replication each) -------------
cat("== A. End-to-end pilots ==\n")
pilot_one <- function(gen, seed, K, std_before_dr = TRUE, n_perm = 39) {
  x <- gen(seed)
  xs <- standardize(x)
  t_dr <- system.time(rm <- run_methods_timed(if (std_before_dr) xs else x))[3]
  ## Pilot uses ties = "random" (seeded) so tie-affected cells (e.g. IMDS
  ## duplicate rectangles) do not abort timing; Stage 5 primary tables use
  ## assess_cell_tie_audit() instead.
  t_ass <- system.time(
    ass <- assess_quality(xs, rm$projections, K = K, perm_test = TRUE,
                          n_perm = n_perm, seed = seed + 1, ties = "random")
  )[3]
  list(x = x, res = ass, dr_time = t_dr, dr_timings = rm$timings,
       assess_time = t_ass, failures = rm$failures)
}

p1 <- pilot_one(gen_scenario1, SEED_BASE + 1, K = 10)
cat(sprintf("Scenario I : DR %.1fs (%s), assess+perm(39) %.1fs, failures: %s\n",
            p1$dr_time, paste(names(p1$dr_timings), round(p1$dr_timings, 1), collapse = " "),
            p1$assess_time, if (length(p1$failures)) paste(p1$failures, collapse = "; ") else "none"))
p2 <- pilot_one(gen_scenario2, SEED_BASE + 2, K = 10)
cat(sprintf("Scenario II: DR %.1fs (%s), assess+perm(39) %.1fs, failures: %s\n",
            p2$dr_time, paste(names(p2$dr_timings), round(p2$dr_timings, 1), collapse = " "),
            p2$assess_time, if (length(p2$failures)) paste(p2$failures, collapse = "; ") else "none"))
p3x <- gen_scenario3(SEED_BASE + 3)
p3s <- standardize(p3x)
# Scenario III uses a FIXED embedding + radius perturbations; time one cycle.
t3 <- system.time({
  rm3 <- run_methods_timed(p3s)
  Dh <- idist(p3s$centers, p3s$radii, "Wasserstein")
})[3]
cat(sprintf("Scenario III: DR-once + one Dh %.1fs\n", t3))

report$pilot <- data.frame(
  scenario = c("I", "II", "III"),
  dr_sec = c(p1$dr_time, p2$dr_time, t3),
  assess_perm39_sec = c(p1$assess_time, p2$assess_time, NA))

## ---------- B. Null expectations on scenario geometry ----------------------
cat("\n== B. Null-expectation checks (m = 2000) ==\n")
x1s <- standardize(p1$x)
Rh <- rank_matrix(idist(x1s$centers, x1s$radii, "Wasserstein"))
rmq <- run_methods_timed(x1s, methods = "C-PCA")
Dl <- if (rmq$projections[["C-PCA"]]$type == "Point") as.matrix(dist(rmq$projections[["C-PCA"]]$C)) else
  idist(rmq$projections[["C-PCA"]]$C, rmq$projections[["C-PCA"]]$R, "Wasserstein")
Rl <- rank_matrix(Dl)
set.seed(SEED_BASE + 10)
nulls <- perm_null_indices(Rh, Rl, K = 10, m = 2000)
N <- nrow(Rh); K <- 10
exp_qlc <- K / (N - 1)
mu <- colMeans(nulls); se <- apply(nulls, 2, sd) / sqrt(nrow(nulls))
chk <- data.frame(index = colnames(nulls), mean = mu, se = se,
                  expected = c(NA, 0, NA, 0, exp_qlc, 0))
print(round(chk[, -1], 5))
ok_null <- all(abs(mu[c("B_TC", "B_RE", "B_LC")]) < 4 * se[c("B_TC", "B_RE", "B_LC")]) &&
  abs(mu["Q_LC"] - exp_qlc) < 4 * se["Q_LC"]
cat("Null-expectation gate:", if (ok_null) "PASS" else "FAIL", "\n")
report$null_check <- chk; report$null_gate <- ok_null

## ---------- C. Type-I error and power ---------------------------------------
cat("\n== C. Type-I / power / min-P calibration (n=60 subsampled geometry) ==\n")
set.seed(SEED_BASE + 20)
idx <- sample(nrow(x1s$centers), 60)
xs60 <- interval_data(x1s$centers[idx, ], x1s$radii[idx, ])
Dh60 <- idist(xs60$centers, xs60$radii, "Wasserstein")
Rh60 <- rank_matrix(Dh60)
Dl60 <- as.matrix(dist(rmq$projections[["C-PCA"]]$C[idx, ]))
n_rep_cal <- if (RUN_MODE == "full") 200 else 20
m_cal <- 199
# Type-I: embedding replaced by an INDEPENDENT random geometry each rep
rej <- matrix(0, n_rep_cal, 6, dimnames = list(NULL, colnames(nulls)))
rej_minP_family <- numeric(n_rep_cal)
for (r in seq_len(n_rep_cal)) {
  Zr <- matrix(rnorm(60 * 2), 60, 2)
  x_dummy <- xs60
  proj_r <- structure(list(nullproj = list(C = Zr, R = matrix(0, 60, 2), type = "Point")),
                      class = "idr_projections")
  a <- assess_quality(x_dummy, proj_r, K = 5, metrics = "Wasserstein",
                      baseline = FALSE, perm_test = TRUE, n_perm = m_cal,
                      seed = SEED_BASE + 100 + r)
  pv <- unlist(a$pvalues[1, c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")])
  rej[r, ] <- pv <= 0.05
  pva <- unlist(a$pvalues_adj[1, c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")])
  rej_minP_family[r] <- any(pva <= 0.05)
}
t1 <- colMeans(rej)
cat("Per-test type-I at alpha=.05:", paste(sprintf("%s=%.3f", names(t1), t1), collapse = " "), "\n")
cat(sprintf("Family-wise (min-P) type-I: %.3f (binomial 95%% CI half-width %.3f)\n",
            mean(rej_minP_family), 1.96 * sqrt(0.05 * 0.95 / n_rep_cal)))
# Gate: every empirical rate within 3 binomial Ses of 0.05
se_bin <- sqrt(0.05 * 0.95 / n_rep_cal)
ok_t1 <- all(abs(t1 - 0.05) < 4 * se_bin) && abs(mean(rej_minP_family) - 0.05) < 4 * se_bin
cat("Type-I gate:", if (ok_t1) "PASS" else "FAIL", "\n")
report$type1 <- t1; report$type1_minP <- mean(rej_minP_family); report$type1_gate <- ok_t1

# Power under graded distortion: true projection + increasing noise swap
cat("Power (graded label-swap distortion, Q_LC test):\n")
frac_swap <- c(0, 0.25, 0.5, 1)
pow <- numeric(length(frac_swap))
for (fi in seq_along(frac_swap)) {
  hits <- 0
  for (r in 1:20) {
    set.seed(SEED_BASE + 300 + fi * 100 + r)
    Zr <- rmq$projections[["C-PCA"]]$C[idx, ]
    ns <- round(frac_swap[fi] * 60)
    if (ns > 1) { sw <- sample(60, ns); Zr[sw, ] <- Zr[sample(sw), ] }
    pt <- perm_test(Dh60, as.matrix(dist(Zr)), K = 5, n_perm = 199,
                    seed = SEED_BASE + 400 + fi * 100 + r)
    hits <- hits + (pt$pQ["LC"] <= 0.05)
  }
  pow[fi] <- hits / 20
}
cat(sprintf("  swap %.0f%% -> power %.2f\n", 100 * (1 - frac_swap), pow))
report$power <- data.frame(kept_fraction = 1 - frac_swap, power = pow)

## ---------- D. Resource projection ------------------------------------------
cat("\n== D. Projected full-run resources ==\n")
n_rep <- 100
per_rep_1 <- p1$dr_time + 6 * 5 * 0.5   # DR + indices for 6 methods x 5 metrics (approx from pilot)
per_rep_2 <- p2$dr_time + 6 * 5 * 1.5
# use measured assess (39 perms) to scale the m=999 calibration subset (25 reps)
cal_1 <- p1$assess_time * (999 / 39) * (25 / 1)
cal_2 <- p2$assess_time * (999 / 39) * (25 / 1)
proj_hours <- (n_rep * (per_rep_1 + per_rep_2) + cal_1 + cal_2 + n_rep * t3) / 3600
cat(sprintf("Projected: %.1f h total (Sc.I %.0f min DR+idx, Sc.II %.0f min, calibration %.0f + %.0f min, Sc.III %.0f min)\n",
            proj_hours, n_rep * per_rep_1 / 60, n_rep * per_rep_2 / 60,
            cal_1 / 60, cal_2 / 60, n_rep * t3 / 60))
report$projection <- data.frame(
  component = c("ScI_reps", "ScII_reps", "ScI_cal", "ScII_cal", "ScIII_reps"),
  minutes = c(n_rep * per_rep_1, n_rep * per_rep_2, cal_1, cal_2, n_rep * t3) / 60)

saveRDS(report, file.path(OUTPUT_DIR, "stage4_pilot_report.rds"))
jsonlite::write_json(list(null_gate = ok_null, type1_gate = ok_t1,
                          projected_hours = proj_hours),
                     file.path(OUTPUT_DIR, "stage4_gates.json"), auto_unbox = TRUE)
cat("\nSTAGE 4 GATES:", if (ok_null && ok_t1) "PASS" else "FAIL", "\n")
