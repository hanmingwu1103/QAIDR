# Table I (Scenario I): width-only structure. Monte Carlo replication with
# center-only baseline and m=999 min-P calibration on the first 25 replicates.
# Full mode: 100 replications. Quick mode: 3 replications, 19 perms.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.r"))
save_sessioninfo("scenario1")

K <- 10L
N_REP <- n_rep_mode(100L)
N_CAL <- min(25L, N_REP)
M_PERM <- n_perm_mode(999L)

all_rows <- list()
cal_rows <- list()
t_start <- Sys.time()
for (r in seq_len(N_REP)) {
  seed_r <- SEED_BASE + 1000L + r
  x <- gen_scenario1(seed_r)
  xs <- standardize(x)
  set.seed(seed_r + 500000L)          # DR stochasticity (UMAP) seeded per rep
  rm <- run_methods_timed(xs)
  if (length(rm$failures)) message("rep ", r, " failures: ",
                                   paste(rm$failures, collapse = "; "))
  cells <- evaluate_cells(xs, rm$projections, K = K, tie_seed = seed_r)
  cells$rep <- r
  all_rows[[r]] <- cells

  if (r <= N_CAL) {
    a <- assess_quality(xs, rm$projections, K = K, perm_test = TRUE,
                        n_perm = M_PERM, seed = seed_r + 900000L,
                        ties = "random", baseline = TRUE)
    pu <- a$pvalues; pa <- a$pvalues_adj
    pu$rep <- r; pa$rep <- r; pu$kind <- "raw"; pa$kind <- "minP"
    cal_rows[[r]] <- rbind(pu, pa)
  }
  if (r %% 10 == 0) cat(sprintf("rep %d/%d elapsed %.1f min\n", r, N_REP,
                                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}

per_rep <- do.call(rbind, all_rows)
idx_cols <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
agg <- aggregate(per_rep[idx_cols], per_rep[c("Method", "Metric")],
                 function(v) c(mean = mean(v), sd = sd(v), mcse = sd(v) / sqrt(length(v))))
summary_df <- do.call(data.frame, agg)
tie_flags <- aggregate(per_rep["tie_affected"], per_rep[c("Method", "Metric")], any)

save_result(per_rep, "scenario1_per_rep",
            extra = list(K = K, n_rep = N_REP, seed_base = SEED_BASE))
save_result(summary_df, "scenario1_summary",
            extra = list(K = K, n_rep = N_REP))
save_result(tie_flags, "scenario1_tie_flags")
if (length(cal_rows)) {
  save_result(do.call(rbind, cal_rows), "scenario1_calibration",
              extra = list(m = M_PERM, n_cal = N_CAL))
}
cat("Scenario 1 done:", format(Sys.time()), "\n")
