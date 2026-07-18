# Regenerate Scenario I/II permutation calibration under the corrected
# joint-draw (Westfall-Young) engine. Re-derives the SAME per-replicate data
# and embeddings from the original seeds (identical RNG offsets to
# scenario1.R / scenario2.R), reruns assess_quality (which now shares one set
# of permutation draws across all metrics per method), and overwrites
# scenario{1,2}_calibration outputs. Point estimates (per_rep/summary) are
# untouched: they never depended on the permutation engine.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.R"))
save_sessioninfo("regen_calibration")

stopifnot(utils::packageVersion("QAIDR") >= "0.2.0")
N_CAL <- 25L
M_PERM <- n_perm_mode(999L)

regen <- function(gen, seed_off, dr_seed_off, cal_seed_off, K, name, raw_dr = FALSE) {
  cal_rows <- vector("list", N_CAL)
  for (r in seq_len(N_CAL)) {
    seed_r <- SEED_BASE + seed_off + r
    x <- gen(seed_r)
    xs <- standardize(x)
    set.seed(seed_r + dr_seed_off)
    rm <- run_methods_timed(if (raw_dr) x else xs)
    a <- assess_quality(xs, rm$projections, K = K, perm_test = TRUE,
                        n_perm = M_PERM, seed = seed_r + cal_seed_off,
                        ties = "random", baseline = TRUE)
    pu <- a$pvalues; pa <- a$pvalues_adj
    pu$rep <- r; pa$rep <- r; pu$kind <- "raw"; pa$kind <- "minP"
    cal_rows[[r]] <- rbind(pu, pa)
    cat(name, "rep", r, "/", N_CAL, "\n")
  }
  save_result(do.call(rbind, cal_rows), paste0(name, "_calibration"),
              extra = list(m = M_PERM, n_cal = N_CAL,
                           engine = "joint-draws shared across metrics (Westfall-Young corrected)"))
}

regen(gen_scenario1, 1000L, 500000L, 900000L, K = 10L, name = "scenario1")
regen(gen_scenario2, 2000L, 500000L, 900000L, K = 10L, name = "scenario2", raw_dr = TRUE)
cat("Calibration regeneration done:", format(Sys.time()), "\n")
