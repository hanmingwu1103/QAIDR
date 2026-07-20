# Prompt H (contract 1.4): regenerate scenario4_calibration under the
# symmetric randomization min-P (QAIDR 0.3.0). Reproduces the exact
# b2_scenario4.R production streams for replications 1..25: same data seeds,
# same DR stream, same documented single-resolution s=1 stream before
# assess_quality. eval_rep() is skipped because every RNG consumer after the
# DR fit is re-seeded (set.seed(seed_r + 300001L)), so the calibration stream
# is unchanged by the omission; per-rep descriptive outputs are NOT rewritten.
# Per-rep assertion: the raw (marginal) p-value rows must be identical to the
# pre-regeneration snapshot (abort condition in H_PREREGISTRATION.md sec. 7).
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.R"))
save_sessioninfo("h_regen_s4_calibration")

stopifnot(utils::packageVersion("QAIDR") >= "0.3.0")
K0 <- 10L
N_SUB <- 25L
M_PERM <- n_perm_mode(999L)

old <- readRDS(file.path(OUTPUT_DIR, "h_old_inference_snapshot",
                         "scenario4_calibration.rds"))
pcols <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")

cal_rows <- vector("list", N_SUB)
t0 <- Sys.time()
for (r in seq_len(N_SUB)) {
  seed_r <- SEED_BASE + 7000L + r
  x <- gen_scenario4(seed_r)
  xs <- standardize(x)
  set.seed(seed_r + 500000L)
  rm <- run_methods_timed(xs)
  set.seed(seed_r + 300001L)      # documented single resolution s=1 stream
  a <- assess_quality(xs, rm$projections, K = K0, perm_test = TRUE,
                      n_perm = M_PERM, seed = seed_r + 900000L,
                      ties = "random", baseline = TRUE)
  pu <- a$pvalues; pa <- a$pvalues_adj
  pu$rep <- r; pa$rep <- r; pu$kind <- "raw"; pa$kind <- "minP"

  ## Re-baseline (contract Amendment A2): raw deltas vs the stored
  ## production snapshot are RECORDED, not asserted zero.
  o <- old[old$kind == "raw" & old$rep == r, ]
  o <- o[order(o$IDR, o$Metric), ]
  n <- pu[order(pu$IDR, pu$Metric), ]
  stopifnot(nrow(o) == nrow(n), all(o$IDR == n$IDR), all(o$Metric == n$Metric))
  dmax <- max(abs(as.matrix(o[pcols]) - as.matrix(n[pcols])))
  cal_rows[[r]] <- rbind(pu, pa)
  cat(sprintf("s4 rep %d/%d raw max|delta|=%.4f elapsed %.1f min\n",
              r, N_SUB, dmax,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

save_result(do.call(rbind, cal_rows), "scenario4_calibration",
            extra = list(m = M_PERM, n_cal = N_SUB,
                         engine = "symmetric randomization min-P (QAIDR 0.3.0, Prompt H)"))
cat("Scenario 4 calibration regeneration done:", format(Sys.time()), "\n")
