# Table II (Scenario II): Swiss roll, variance-dominated view.
# Primary protocol (prespecified to match the manuscript's stated mechanism):
#   DR on RAW data, evaluation on STANDARDIZED geometry.
# Factorial (Stage 5 item 6) on the first N_FACT replicates:
#   {raw, std} DR input x {raw, std} evaluation geometry.
# Full mode: 100 replications, m=999 calibration on 25. Quick: 3 reps, 19 perms.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.r"))
save_sessioninfo("scenario2")

K <- 10L
N_REP <- n_rep_mode(100L)
N_CAL <- min(25L, N_REP)
N_FACT <- min(25L, N_REP)
M_PERM <- n_perm_mode(999L)

all_rows <- list(); fact_rows <- list(); cal_rows <- list()
t_start <- Sys.time()
for (r in seq_len(N_REP)) {
  seed_r <- SEED_BASE + 2000L + r
  x_raw <- gen_scenario2(seed_r)
  x_std <- standardize(x_raw)

  set.seed(seed_r + 500000L)
  rm_raw <- run_methods_timed(x_raw)      # primary: DR on raw data
  if (length(rm_raw$failures)) message("rep ", r, " raw-DR failures: ",
                                       paste(rm_raw$failures, collapse = "; "))
  cells <- evaluate_cells(x_std, rm_raw$projections, K = K, tie_seed = seed_r)
  cells$rep <- r; cells$dr_input <- "raw"; cells$eval_geom <- "std"
  all_rows[[r]] <- cells

  if (r <= N_FACT) {
    set.seed(seed_r + 600000L)
    rm_std <- run_methods_timed(x_std)
    combos <- list(
      list(rm = rm_raw, ev = x_raw, dr = "raw", ge = "raw"),
      list(rm = rm_std, ev = x_raw, dr = "std", ge = "raw"),
      list(rm = rm_std, ev = x_std, dr = "std", ge = "std")
    )
    fr <- lapply(combos, function(cb) {
      cc <- evaluate_cells(cb$ev, cb$rm$projections, K = K, tie_seed = seed_r)
      cc$rep <- r; cc$dr_input <- cb$dr; cc$eval_geom <- cb$ge
      cc
    })
    fact_rows[[r]] <- do.call(rbind, fr)
  }

  if (r <= N_CAL) {
    a <- assess_quality(x_std, rm_raw$projections, K = K, perm_test = TRUE,
                        n_perm = M_PERM, seed = seed_r + 900000L,
                        ties = "random", baseline = TRUE)
    pu <- a$pvalues; pa <- a$pvalues_adj
    pu$rep <- r; pa$rep <- r; pu$kind <- "raw"; pa$kind <- "minP"
    cal_rows[[r]] <- rbind(pu, pa)
  }
  if (r %% 5 == 0) cat(sprintf("rep %d/%d elapsed %.1f min\n", r, N_REP,
                               as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}

per_rep <- do.call(rbind, all_rows)
idx_cols <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
agg <- aggregate(per_rep[idx_cols], per_rep[c("Method", "Metric")],
                 function(v) c(mean = mean(v), sd = sd(v), mcse = sd(v) / sqrt(length(v))))
save_result(per_rep, "scenario2_per_rep", extra = list(K = K, n_rep = N_REP))
save_result(do.call(data.frame, agg), "scenario2_summary", extra = list(K = K, n_rep = N_REP))
## Prompt E fix (2026-07-17): tables_to_tex.r expects scenario2_tie_flags.csv
## for the dagger markers (same any()-over-replications rule as Scenario I);
## this aggregation was missing, so Table II silently lost its daggers.
tie_flags <- aggregate(per_rep["tie_affected"], per_rep[c("Method", "Metric")], any)
save_result(tie_flags, "scenario2_tie_flags", extra = list(rule = "any over replications"))
if (length(fact_rows)) {
  fact <- rbind(cbind(per_rep[per_rep$rep <= N_FACT, ]), do.call(rbind, fact_rows))
  fagg <- aggregate(fact[idx_cols], fact[c("Method", "Metric", "dr_input", "eval_geom")], mean)
  save_result(fact, "scenario2_factorial_per_rep", extra = list(n_fact = N_FACT))
  save_result(fagg, "scenario2_factorial_summary", extra = list(n_fact = N_FACT))
}
if (length(cal_rows)) save_result(do.call(rbind, cal_rows), "scenario2_calibration",
                                  extra = list(m = M_PERM, n_cal = N_CAL))
cat("Scenario 2 done:", format(Sys.time()), "\n")
