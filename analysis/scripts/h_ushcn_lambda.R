# Prompt H empirical addition A (H_PREREGISTRATION.md sec. 2 + A1.4, frozen
# before results): USHCN Interval-Euclidean lambda mechanism audit.
# Cold R process: frozen 645-station snapshot, standardize, then exactly one
# seeded production embedding call (set.seed(SEED_BASE + 6000L) immediately
# before run_methods_timed) -- verified to reproduce the production raw
# p-values byte-identically in the Stage 3 regeneration.
# Grid: lambda in seq(0, 1, by = 0.05), 21 values, none omitted; every cell
# = all six methods; center baseline recorded once as lambda-invariant.
# S_lambda decomposes exactly as lambda*d_max + (1-lambda)*d_min per entry,
# so the lambda grid uses the two component matrices (bitwise identical to
# calling idist_euclidean(lambda) directly).
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.R"))
save_sessioninfo("h_ushcn_lambda")

stopifnot(utils::packageVersion("QAIDR") >= "0.3.0",
          exists("interval_data_from_dataSDA",
                 envir = asNamespace("QAIDR"), inherits = FALSE))
K <- 10L
LAMBDA_GRID <- seq(0, 1, by = 0.05)
TIE_SEED <- SEED_BASE + 6000L
t_start <- Sys.time()

snap <- readRDS(file.path(DATA_PROC, "ghcn_seasonal_intervals.rds"))
x <- interval_data(snap$centers, snap$radii)
xs <- standardize(x)
N <- nrow(snap$centers)
cat("USHCN snapshot:", N, "stations, year", snap$target_year, "\n")

set.seed(SEED_BASE + 6000L)
rm <- run_methods_timed(xs)
if (length(rm$failures)) stop("DR failures: ", paste(rm$failures, collapse = "; "))

## Component matrices: d_min (lambda = 0) and d_max (lambda = 1) per space.
comp <- list()
comp[["high"]] <- list(
  D0 = idist_euclidean(xs$centers, xs$radii, lambda = 0),
  D1 = idist_euclidean(xs$centers, xs$radii, lambda = 1))
for (m in names(rm$projections)) {
  pr <- rm$projections[[m]]
  comp[[m]] <- if (pr$type == "Point") {
    Dp <- as.matrix(stats::dist(pr$C)); list(D0 = Dp, D1 = Dp)
  } else {
    list(D0 = idist_euclidean(pr$C, pr$R, lambda = 0),
         D1 = idist_euclidean(pr$C, pr$R, lambda = 1))
  }
}

## Fixed aggregate d^-/d^+ diagnostics over ALL unordered i<j pairs (A1.4).
comp_rows <- lapply(names(comp), function(sp) {
  ut <- upper.tri(comp[[sp]]$D0)
  dm <- comp[[sp]]$D0[ut]; dp <- comp[[sp]]$D1[ut]
  elig <- dp > 0
  data.frame(space = sp,
             n_pairs = length(dm), n_eligible = sum(elig),
             mean_dminus = mean(dm), sd_dminus = stats::sd(dm),
             mean_dplus = mean(dp), sd_dplus = stats::sd(dp),
             mean_ratio = mean(dm[elig] / dp[elig]),
             sd_ratio = stats::sd(dm[elig] / dp[elig]),
             frac_dminus_zero = mean(dm == 0),
             stringsAsFactors = FALSE)
})
save_result(do.call(rbind, comp_rows), "h_ushcn_lambda_components",
            extra = list(n = N, prereg = "H_PREREGISTRATION.md sec 2 (post-A2)",
                         pair_def = "unordered off-diagonal i<j"))

## Lambda profile: all six methods x 21 lambda values, tie-audit protocol.
rows <- list(); ii <- 0L
for (lam in LAMBDA_GRID) {
  Dh <- lam * comp[["high"]]$D1 + (1 - lam) * comp[["high"]]$D0
  for (m in names(rm$projections)) {
    Dl <- lam * comp[[m]]$D1 + (1 - lam) * comp[[m]]$D0
    cell <- assess_cell_tie_audit(Dh, Dl, K, seed = TIE_SEED)
    ii <- ii + 1L
    rows[[ii]] <- data.frame(Method = m, lambda = lam, t(cell$vals),
                             tie_affected = cell$tie_affected,
                             tie_spread = cell$max_spread,
                             n_res = cell$n_res, stringsAsFactors = FALSE)
  }
  cat(sprintf("lambda %.2f done, elapsed %.1f min\n", lam,
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}
## Center-only baseline, lambda-invariant by construction, recorded once.
Dh_c <- as.matrix(stats::dist(xs$centers))
for (m in names(rm$projections)) {
  Dl_c <- as.matrix(stats::dist(rm$projections[[m]]$C))
  cell <- assess_cell_tie_audit(Dh_c, Dl_c, K, seed = TIE_SEED)
  ii <- ii + 1L
  rows[[ii]] <- data.frame(Method = m, lambda = NA_real_, t(cell$vals),
                           tie_affected = cell$tie_affected,
                           tie_spread = cell$max_spread,
                           n_res = cell$n_res, stringsAsFactors = FALSE)
}
prof <- do.call(rbind, rows)
save_result(prof, "h_ushcn_lambda_profile",
            extra = list(K = K, n = N, lambda_grid = LAMBDA_GRID,
                         tie_seed = TIE_SEED, n_res = 50L,
                         prereg = "H_PREREGISTRATION.md sec 2 (post-A2)",
                         baseline_row = "lambda = NA, center-only Euclidean"))
rt <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
save_result(data.frame(component = "h_ushcn_lambda_total", seconds = rt),
            "h_ushcn_lambda_runtime")
cat("h_ushcn_lambda done:", format(Sys.time()), "\n")
