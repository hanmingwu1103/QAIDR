# Corrected assessment-only scalability benchmark (Stage 3 / Stage 4 item 16).
#
# This is the SAME benchmark design as benchmark.r -- same sizes, same seeds,
# same generated data, same operations. Only the MEASUREMENT is corrected. It is
# not a new benchmark and not a new simulation.
#
# Defects in benchmark.r that this script fixes:
#   1. Dl (low-space distances) and Rl (low-space ranks) were constructed
#      OUTSIDE any system.time() call, so the reported total omitted two of the
#      components the manuscript claimed were timed.
#   2. A single timing pass, so no dispersion was available.
#   3. approx_matrix_mem_MB = 3*8*n^2/1e6 is a size CALCULATION from the matrix
#      dimensions, not a measurement. It was reported as though measured.
#   4. Only m=99 was run, while production inference uses m=999.
#
# Embedding (dimensionality reduction) time is deliberately EXCLUDED here: this
# is an assessment-only benchmark. The DR timings remain in the original
# scalability_benchmark.csv and are reported separately in the manuscript; they
# are dominated by the bundled implementations (SPCA ~11,845 s at n=2000) and
# would otherwise swamp the assessment cost this benchmark is measuring.

SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
save_sessioninfo("benchmark_corrected")

NS <- if (RUN_MODE == "full") c(200L, 500L, 1000L, 2000L) else c(200L, 500L)
N_REP <- as.integer(Sys.getenv("QAIDR_BENCH_REPS", "3"))
stopifnot(N_REP >= 1L)

## Peak memory below is R's gc max-used high-water mark: a measurement of this
## process's allocation peak, not the OS-level peak resident set size. Obtaining
## RSS requires an external observer, so it is reported as NA rather than
## approximated by a size calculation.

timeit <- function(expr) as.numeric(system.time(force(expr))[3])

rows <- list(); ii <- 0L
for (n in NS) {
  ## identical data generation to benchmark.r (same seed, same construction)
  set.seed(SEED_BASE + 5000L + n)
  centers <- matrix(rnorm(n * 5), n, 5)
  radii   <- matrix(runif(n * 5, 0.1, 1), n, 5)
  Zl_fixed <- centers[, 1:2] + matrix(rnorm(n * 2, sd = 0.1), n, 2)

  td <- tr <- tdl <- trl <- ti <- tp99 <- tp999 <- numeric(N_REP)

  for (r in seq_len(N_REP)) {
    ## every component the manuscript claims is timed IS timed
    td[r]  <- timeit(Dh <- idist(centers, radii, "Wasserstein"))
    tr[r]  <- timeit(Rh <- rank_matrix(Dh))
    tdl[r] <- timeit(Dl <- as.matrix(dist(Zl_fixed)))   # was untimed
    trl[r] <- timeit(Rl <- rank_matrix(Dl))             # was untimed
    ti[r]  <- timeit(v  <- indices_from_coranking(coranking_matrix(Rh, Rl), K = 10))
    tp99[r]  <- timeit(perm_null_indices(Rh, Rl, K = 10, m = 99L))
    tp999[r] <- timeit(perm_null_indices(Rh, Rl, K = 10, m = 999L))
  }

  ## measured peak R heap over one full assessment pass at this n
  invisible(gc(reset = TRUE, full = TRUE))
  Dh <- idist(centers, radii, "Wasserstein"); Rh <- rank_matrix(Dh)
  Dl <- as.matrix(dist(Zl_fixed));            Rl <- rank_matrix(Dl)
  v  <- indices_from_coranking(coranking_matrix(Rh, Rl), K = 10)
  invisible(perm_null_indices(Rh, Rl, K = 10, m = 99L))
  gpk <- gc(full = TRUE)
  measured_peak_mb <- sum(gpk[, 6])

  ## the figure benchmark.r reported, retained for comparison and labelled
  calculated_mem_mb <- 3 * 8 * n^2 / 1e6

  med <- function(x) as.numeric(stats::median(x))
  total99  <- med(td) + med(tr) + med(tdl) + med(trl) + med(ti) + med(tp99)
  total999 <- med(td) + med(tr) + med(tdl) + med(trl) + med(ti) + med(tp999)

  ii <- ii + 1L
  rows[[ii]] <- data.frame(
    n = n, reps = N_REP,
    t_dist_hi = med(td),  t_dist_hi_min = min(td),  t_dist_hi_max = max(td),
    t_rank_hi = med(tr),  t_rank_hi_min = min(tr),  t_rank_hi_max = max(tr),
    t_dist_lo = med(tdl), t_dist_lo_min = min(tdl), t_dist_lo_max = max(tdl),
    t_rank_lo = med(trl), t_rank_lo_min = min(trl), t_rank_lo_max = max(trl),
    t_indices = med(ti),  t_indices_min = min(ti),  t_indices_max = max(ti),
    t_perm99  = med(tp99),  t_perm99_min  = min(tp99),  t_perm99_max  = max(tp99),
    t_perm999 = med(tp999), t_perm999_min = min(tp999), t_perm999_max = max(tp999),
    total_assessment_m99  = total99,
    total_assessment_m999 = total999,
    measured_peak_r_heap_MB = measured_peak_mb,
    calculated_matrix_mem_MB = calculated_mem_mb,
    os_peak_rss_MB = NA_real_,
    embedding_time_included = FALSE
  )
  cat(sprintf(paste0("n=%d (%d reps): Dh %.2fs Rh %.2fs Dl %.2fs Rl %.2fs idx %.3fs ",
                     "perm99 %.2fs perm999 %.2fs | total(m=99) %.2fs total(m=999) %.2fs ",
                     "| measured peak heap %.0f MB (calculated %.0f MB)\n"),
              n, N_REP, med(td), med(tr), med(tdl), med(trl), med(ti),
              med(tp99), med(tp999), total99, total999,
              measured_peak_mb, calculated_mem_mb))
}

bench <- do.call(rbind, rows)
save_result(bench, "scalability_benchmark_corrected",
            extra = list(
              sizes = NS, reps = N_REP,
              supersedes = "scalability_benchmark",
              assessment_only = TRUE,
              embedding_excluded = paste(
                "DR/embedding time is excluded by design; see scalability_benchmark.csv",
                "for the stored per-method embedding timings"),
              memory_note = paste(
                "measured_peak_r_heap_MB is R's gc max-used high-water mark, a",
                "measurement; calculated_matrix_mem_MB is the 3*8*n^2 size",
                "calculation previously reported as if measured; os_peak_rss_MB",
                "is NA because OS-level RSS requires an external observer"),
              m_note = "m=99 is the benchmark setting; m=999 is production inference"))
cat("Corrected benchmark done:", format(Sys.time()), "\n")
