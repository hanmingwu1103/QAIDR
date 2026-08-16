# Scalability benchmark (approved G.3 evidence): n in {200, 500, 1000, 2000},
# runtime and peak memory for distance construction, ranking, co-ranking
# indices, and a m=99 permutation block; DR methods timed where feasible.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.r"))
save_sessioninfo("benchmark")

NS <- if (RUN_MODE == "full") c(200L, 500L, 1000L, 2000L) else c(200L, 500L)
rows <- list(); ii <- 0L
for (n in NS) {
  set.seed(SEED_BASE + 5000L + n)
  centers <- matrix(rnorm(n * 5), n, 5)
  radii <- matrix(runif(n * 5, 0.1, 1), n, 5)
  x <- interval_data(centers, radii)

  t_dist <- system.time(Dh <- idist(centers, radii, "Wasserstein"))[3]
  t_rank <- system.time(Rh <- rank_matrix(Dh))[3]
  Zl <- centers[, 1:2] + matrix(rnorm(n * 2, sd = 0.1), n, 2)
  Dl <- as.matrix(dist(Zl)); Rl <- rank_matrix(Dl)
  t_idx <- system.time(v <- indices_from_coranking(coranking_matrix(Rh, Rl), K = 10))[3]
  t_perm <- system.time(perm_null_indices(Rh, Rl, K = 10, m = 99))[3]
  mem_mb <- 3 * 8 * n^2 / 1e6   # Dh + two rank matrices resident (approx)

  ## DR method timings (each guarded; skip gracefully on failure/timeouts)
  dr_times <- rep(NA_real_, length(METHODS)); names(dr_times) <- METHODS
  if (n <= 2000) {
    xs <- standardize(x)
    for (m in METHODS) {
      tm <- tryCatch(system.time(run_idr(xs, methods = m, verbose = FALSE))[3],
                     error = function(e) NA_real_)
      dr_times[m] <- tm
    }
  }
  ii <- ii + 1L
  rows[[ii]] <- data.frame(n = n, t_dist = t_dist, t_rank = t_rank,
                           t_indices = t_idx, t_perm99 = t_perm,
                           approx_matrix_mem_MB = mem_mb, t(dr_times))
  cat(sprintf("n=%d: dist %.2fs rank %.2fs idx %.3fs perm99 %.1fs\n",
              n, t_dist, t_rank, t_idx, t_perm))
}
bench <- do.call(rbind, rows)
save_result(bench, "scalability_benchmark", extra = list(sizes = NS))
cat("Benchmark done:", format(Sys.time()), "\n")
