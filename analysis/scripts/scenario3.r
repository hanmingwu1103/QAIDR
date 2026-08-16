# Table III (Scenario III): index stability under dense overlapping noise.
# Corrected protocol: R_DATA independently seeded datasets (full: 100), each
# with fixed embeddings and N_PERT radius perturbations; report mean absolute
# deviations and sign-flip rates with between-dataset variability.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.r"))
save_sessioninfo("scenario3")

K <- 10L
R_DATA <- n_rep_mode(100L)
N_PERT <- if (RUN_MODE == "full") 100L else 5L
EPS <- 0.25

idx_cols <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
out_rows <- list(); oo <- 0L
t_start <- Sys.time()

for (d in seq_len(R_DATA)) {
  seed_d <- SEED_BASE + 3000L + d
  x <- gen_scenario3(seed_d)
  xs <- standardize(x)
  set.seed(seed_d + 500000L)
  rm <- run_methods_timed(xs)
  if (length(rm$failures)) message("dataset ", d, " failures: ",
                                   paste(rm$failures, collapse = "; "))

  ## Baseline (unperturbed) indices per method x metric; embeddings FIXED.
  base <- evaluate_cells(xs, rm$projections, K = K, tie_seed = seed_d)

  ## Low-space distance matrices are fixed; precompute their ranks per cell.
  for (pert in seq_len(N_PERT)) {
    set.seed(seed_d + 700000L + pert)
    delta <- matrix(runif(length(x$radii), -EPS, EPS), nrow(x$radii))
    x_p <- interval_data(x$centers, x$radii * (1 + delta))
    xs_p <- standardize(x_p)
    cells_p <- evaluate_cells(xs_p, rm$projections, K = K, tie_seed = seed_d + pert)
    mg <- merge(cells_p, base, by = c("Method", "Metric"),
                suffixes = c("_pert", "_base"))
    oo <- oo + 1L
    out_rows[[oo]] <- data.frame(
      dataset = d, pert = pert,
      Method = mg$Method, Metric = mg$Metric,
      dQ_TC = abs(mg$Q_TC_pert - mg$Q_TC_base),
      dB_TC = abs(mg$B_TC_pert - mg$B_TC_base),
      dQ_RE = abs(mg$Q_RE_pert - mg$Q_RE_base),
      dB_RE = abs(mg$B_RE_pert - mg$B_RE_base),
      dQ_LC = abs(mg$Q_LC_pert - mg$Q_LC_base),
      dB_LC = abs(mg$B_LC_pert - mg$B_LC_base),
      flip_TC = sign(mg$B_TC_pert) != sign(mg$B_TC_base) & mg$B_TC_base != 0,
      flip_RE = sign(mg$B_RE_pert) != sign(mg$B_RE_base) & mg$B_RE_base != 0,
      flip_LC = sign(mg$B_LC_pert) != sign(mg$B_LC_base) & mg$B_LC_base != 0)
  }
  if (d %% 10 == 0) cat(sprintf("dataset %d/%d elapsed %.1f min\n", d, R_DATA,
                                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}

per <- do.call(rbind, out_rows)
dcols <- grep("^d|^flip", names(per), value = TRUE)
summ <- aggregate(per[dcols], per[c("Method", "Metric")],
                  function(v) c(mean = mean(v), sd = sd(v)))
save_result(per, "scenario3_per_perturbation",
            extra = list(K = K, R_data = R_DATA, n_pert = N_PERT, eps = EPS))
save_result(do.call(data.frame, summ), "scenario3_summary",
            extra = list(K = K, R_data = R_DATA, n_pert = N_PERT, eps = EPS))
cat("Scenario 3 done:", format(Sys.time()), "\n")
