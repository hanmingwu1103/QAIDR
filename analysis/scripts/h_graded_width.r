# Prompt H empirical addition B (H_PREREGISTRATION.md sec. 3 + A1.5/A1.6 +
# A2, frozen before results): graded width-signal Scenario-I extension.
# Levels (radii supports per group of 100):
#   L1 weak     U(0.1,0.2) / U(0.3,0.4) / U(0.5,0.7)
#   L2 moderate U(0.1,0.2) / U(0.5,0.7) / U(1.5,1.9)
#   L3 reference = production Scenario I: U(0.1,0.2) / U(1.0,1.2) / U(5,6)
# L3 reuses production replications 1..25 (data seeds SEED_BASE+1000+r, DR
# stream seed_r+500000) and is reconstructed FIRST in this cold process,
# before any L1/L2 fit (A1.5); recomputed pre-existing cells must agree with
# stored scenario1_per_rep within 1e-5 (A2 numerical-noise tolerance).
# Evaluations on the SAME embeddings: four named dissimilarities, center-only
# Euclidean, and one (center,radius)-Euclidean evaluation "CR-Euclidean"
# (rank-equivalent to endpoint Euclidean up to the constant sqrt(2);
# Wasserstein is (c, r/sqrt(3))-Euclidean). No additional radius weight (H.5).
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.r"))
save_sessioninfo("h_graded_width")

stopifnot(utils::packageVersion("QAIDR") >= "0.3.0")
K <- 10L
N_REP <- 25L
SEED_H <- 20260719L
L3_TOL <- 1e-5
pcix <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
t_start <- Sys.time()

gen_graded <- function(seed, lo, mid, hi) {
  ## Scenario I center geometry with identical RNG consumption; only the
  ## radii supports differ (same runif call count and order).
  set.seed(seed)
  N <- 300L; P <- 5L
  centers <- matrix(rnorm(N * P), N, P)
  radii <- rbind(
    matrix(runif(100 * P, lo[1], lo[2]), 100, P),
    matrix(runif(100 * P, mid[1], mid[2]), 100, P),
    matrix(runif(100 * P, hi[1], hi[2]), 100, P)
  )
  colnames(centers) <- colnames(radii) <- paste0("V", seq_len(P))
  interval_data(centers, radii,
                labels = factor(rep(c("small", "medium", "large"), each = 100)))
}
LEVELS <- list(
  L1 = list(lo = c(0.1, 0.2), mid = c(0.3, 0.4), hi = c(0.5, 0.7),
            seed_off = 8000L),
  L2 = list(lo = c(0.1, 0.2), mid = c(0.5, 0.7), hi = c(1.5, 1.9),
            seed_off = 8100L))

cr_dist <- function(C, R) as.matrix(stats::dist(cbind(C, R)))

eval_with_cr <- function(xs, projections, tie_seed) {
  cells <- evaluate_cells(xs, projections, K = K, tie_seed = tie_seed)
  Dh_cr <- cr_dist(xs$centers, xs$radii)
  extra <- lapply(names(projections), function(m) {
    pr <- projections[[m]]
    Dl <- if (pr$type == "Point") as.matrix(stats::dist(pr$C))
          else cr_dist(pr$C, pr$R)
    cell <- assess_cell_tie_audit(Dh_cr, Dl, K, seed = tie_seed)
    data.frame(Method = m, Metric = "CR-Euclidean", t(cell$vals),
               tie_affected = cell$tie_affected, tie_spread = cell$max_spread,
               stringsAsFactors = FALSE)
  })
  rbind(cells, do.call(rbind, extra))
}

per_stored <- readRDS(file.path(OUTPUT_DIR, "scenario1_per_rep.rds"))
all_rows <- list(); aa <- 0L
l3_dmax <- 0

## ---- L3 FIRST (production reconstruction, A1.5) ---------------------------
for (r in seq_len(N_REP)) {
  seed_r <- SEED_BASE + 1000L + r
  x <- gen_scenario1(seed_r)
  xs <- standardize(x)
  set.seed(seed_r + 500000L)
  rm <- run_methods_timed(xs)
  cells <- eval_with_cr(xs, rm$projections, tie_seed = seed_r)
  ## Agreement gate vs stored production per-rep values (pre-existing cells).
  o <- per_stored[per_stored$rep == r, ]
  o <- o[order(o$Method, o$Metric), ]
  n <- cells[cells$Metric != "CR-Euclidean", ]
  n <- n[order(n$Method, n$Metric), ]
  stopifnot(nrow(o) == nrow(n), all(o$Method == n$Method),
            all(o$Metric == n$Metric))
  dmax <- max(abs(as.matrix(o[pcix]) - as.matrix(n[pcix])))
  l3_dmax <- max(l3_dmax, dmax)
  if (dmax > L3_TOL) stop(sprintf(
    "L3 rep %d disagrees with scenario1_per_rep: max|delta| = %g > %g",
    r, dmax, L3_TOL))
  cells$level <- "L3"; cells$rep <- r
  aa <- aa + 1L; all_rows[[aa]] <- cells
  cat(sprintf("L3 rep %d/%d ok (max|delta|=%.2e) elapsed %.1f min\n", r,
              N_REP, dmax, as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}

## ---- L1, L2 ---------------------------------------------------------------
for (lv in names(LEVELS)) {
  cfg <- LEVELS[[lv]]
  for (r in seq_len(N_REP)) {
    seed_r <- SEED_H + cfg$seed_off + r
    x <- gen_graded(seed_r, cfg$lo, cfg$mid, cfg$hi)
    xs <- standardize(x)
    set.seed(seed_r + 500000L)
    rm <- run_methods_timed(xs)
    cells <- eval_with_cr(xs, rm$projections, tie_seed = seed_r)
    cells$level <- lv; cells$rep <- r
    aa <- aa + 1L; all_rows[[aa]] <- cells
    cat(sprintf("%s rep %d/%d elapsed %.1f min\n", lv, r, N_REP,
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  }
}

per <- do.call(rbind, all_rows)
save_result(per, "h_graded_width_per_rep",
            extra = list(K = K, n_rep = N_REP,
                         levels = "L1 weak / L2 moderate / L3 production",
                         l3_agreement_max_delta = l3_dmax, l3_tol = L3_TOL,
                         prereg = "H_PREREGISTRATION.md sec 3 (post-A2)"))

## ---- Primary estimand: paired Qc_interval[Wasserstein] - Qc_center --------
n_obj <- 300L
EQ <- local({
  n <- n_obj; Kk <- K
  GK <- if (Kk < n / 2) n * Kk * (2 * n - 3 * Kk - 1) else n * (n - Kk) * (n - Kk - 1)
  HK <- n * sum(abs(n - 2 * seq_len(Kk) + 1) / seq_len(Kk))
  j <- seq_len(Kk)
  EWn <- n / ((n - 1) * HK) * sum((1 / j) * (j * (j - 1) / 2 + (n - 1 - j) * (n - j) / 2))
  c(Q_TC = 1 - n * Kk * (n - 1 - Kk) * (n - Kk) / (GK * (n - 1)),
    Q_RE = 1 - EWn,
    Q_LC = Kk / (n - 1))
})
qc <- function(q, mu0) (q - mu0) / (1 - mu0)

contrast_rows <- list(); cc <- 0L
for (lv in unique(per$level)) for (m in unique(per$Method)) {
  base <- per[per$level == lv & per$Method == m & per$Metric == "Centers-Euclidean", ]
  base <- base[order(base$rep), ]
  for (met in setdiff(unique(per$Metric), "Centers-Euclidean")) {
    int <- per[per$level == lv & per$Method == m & per$Metric == met, ]
    int <- int[order(int$rep), ]
    stopifnot(nrow(int) == nrow(base), all(int$rep == base$rep))
    for (fam in c("Q_TC", "Q_RE", "Q_LC")) {
      d <- qc(int[[fam]], EQ[fam]) - qc(base[[fam]], EQ[fam])
      cc <- cc + 1L
      contrast_rows[[cc]] <- data.frame(
        level = lv, Method = m, Metric = met, index = fam,
        primary = (met == "Wasserstein" && fam == "Q_TC"),
        mean_qc_diff = mean(d), mcse = stats::sd(d) / sqrt(length(d)),
        n_rep = length(d), stringsAsFactors = FALSE)
    }
  }
}
summ <- do.call(rbind, contrast_rows)
save_result(summ, "h_graded_width_summary",
            extra = list(K = K, n = n_obj, null_means = as.list(EQ),
                         qc_def = "(Q - mu0)/(1 - mu0), mu0 closed-form null mean (thm:nullcal)",
                         primary = "Wasserstein Q_TC chance-adjusted paired diff vs center baseline",
                         prereg = "H_PREREGISTRATION.md sec 3 (post-A2)"))
rt <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
save_result(data.frame(component = "h_graded_width_total", seconds = rt),
            "h_graded_width_runtime")
cat("h_graded_width done:", format(Sys.time()), "\n")
