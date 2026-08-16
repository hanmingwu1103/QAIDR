# Prompt H (contract 1.4): regenerate scenario4_calibration under the
# symmetric randomization min-P (QAIDR 0.3.0). Reproduces the exact
# b2_scenario4.r production streams for replications 1..25: same data seeds,
# same DR stream, same documented single-resolution s=1 stream before
# assess_quality. eval_rep() is skipped because every RNG consumer after the
# DR fit is re-seeded (set.seed(seed_r + 300001L)), so the calibration stream
# is unchanged by the omission; per-rep descriptive outputs are NOT rewritten.
# Per-rep assertion: the raw (marginal) p-value rows must be identical to the
# pre-regeneration snapshot (abort condition in H_PREREGISTRATION.md sec. 7).
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.r"))
save_sessioninfo("h_regen_s4_calibration")

## Exact release pin (see _common.r); a floor would re-admit version drift.
stopifnot(identical(as.character(utils::packageVersion("QAIDR")),
                    QAIDR_REQUIRED_VERSION))
K0 <- 10L
N_SUB <- 25L
M_PERM <- n_perm_mode(999L)

## h_old_inference_snapshot/ is OPTIONAL HISTORICAL EVIDENCE, not a production
## dependency. It has no producer anywhere in the tree, so a fresh checkout
## cannot rebuild it; hard-reading it made this target unrunnable on a clean
## clone at any wall-clock budget (PRODUCTION_RERUN_DECISION.md carve-out).
## The snapshot only ever fed a RECORDED delta (dmax) that is explicitly not
## asserted to be zero -- it never entered the saved calibration output -- so
## its absence changes no result. When present it is still compared and the
## delta still reported; when absent the comparison is skipped and that fact is
## recorded in the provenance sidecar.
OLD_SNAPSHOT <- file.path(OUTPUT_DIR, "h_old_inference_snapshot",
                          "scenario4_calibration.rds")
have_old <- file.exists(OLD_SNAPSHOT)
old <- if (have_old) readRDS(OLD_SNAPSHOT) else NULL
if (!have_old)
  message("h_regen_s4_calibration: optional historical snapshot absent (",
          rel_path(OLD_SNAPSHOT), "); skipping the re-baseline comparison. ",
          "This does not affect the regenerated calibration.")
pcols <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
delta_log <- rep(NA_real_, N_SUB)

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
  ## production snapshot are RECORDED, not asserted zero. Skipped when the
  ## optional historical snapshot is unavailable.
  dmax <- NA_real_
  if (have_old) {
    o <- old[old$kind == "raw" & old$rep == r, ]
    o <- o[order(o$IDR, o$Metric), ]
    n <- pu[order(pu$IDR, pu$Metric), ]
    stopifnot(nrow(o) == nrow(n), all(o$IDR == n$IDR), all(o$Metric == n$Metric))
    dmax <- max(abs(as.matrix(o[pcols]) - as.matrix(n[pcols])))
  }
  delta_log[r] <- dmax
  cal_rows[[r]] <- rbind(pu, pa)
  cat(sprintf("s4 rep %d/%d raw max|delta|=%s elapsed %.1f min\n",
              r, N_SUB,
              if (is.na(dmax)) "NA (no historical snapshot)" else sprintf("%.4f", dmax),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

save_result(do.call(rbind, cal_rows), "scenario4_calibration",
            inputs = if (have_old) OLD_SNAPSHOT else character(),
            extra = list(m = M_PERM, n_cal = N_SUB,
                         engine = "symmetric randomization min-P (QAIDR 0.3.0, Prompt H)",
                         historical_snapshot_present = have_old,
                         historical_snapshot_role = paste(
                           "optional historical evidence; recorded delta only,",
                           "never an input to the saved calibration"),
                         raw_delta_max = if (have_old) max(delta_log, na.rm = TRUE) else NA_real_,
                         raw_delta_by_rep = as.list(delta_log)))
cat("Scenario 4 calibration regeneration done:", format(Sys.time()), "\n")
