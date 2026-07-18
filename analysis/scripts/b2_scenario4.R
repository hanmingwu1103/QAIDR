# B2 Scenario IV (preregistered A-D2, frozen 2026-07-16): heterogeneous/
# adversarial stress scenario. 100 replications; K grid {5,10,20,40} evaluated
# for replications 1..25 in the same pass; min-P calibration (m=999, joint
# draws) on replications 1..25 under one documented tie resolution (s=1).
# Pilot mode: set B2_PILOT=<n> to run n replications and report ONLY wall
# time, memory floor, and tie counts (quality-blind feasibility audit).
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.R"))

PILOT <- as.integer(Sys.getenv("B2_PILOT", "0"))
if (PILOT == 0) save_sessioninfo("b2_scenario4")

K0 <- 10L
K_GRID <- c(5L, 10L, 20L, 40L)      # {0.5,1,2,4} x K0 per prereg A-D1
N_REP <- if (PILOT > 0) PILOT else n_rep_mode(100L)
N_SUB <- 25L
N_RES <- 50L
M_PERM <- n_perm_mode(999L)

count_row_ties <- function(D, tol = 1e-12) {
  n <- nrow(D); tot <- 0L
  for (i in seq_len(n)) {
    ds <- sort(D[i, -i])
    tot <- tot + sum(diff(ds) <= tol * pmax(ds[-1], .Machine$double.xmin))
  }
  tot
}

## Tie-aware multi-resolution evaluation of all cells at all K in K_GRID.
## Shares the high-space rank resolution across methods within each metric
## and resolution. Returns list(cells = long df over resolutions-mean,
## ties = per-cell tie info).
eval_rep <- function(x_eval, projections, kmax_needed, tie_seed) {
  mets <- c(METRICS, "Centers-Euclidean")
  Dh <- list(); Dl <- list()
  for (met in mets) {
    Dh[[met]] <- if (met == "Centers-Euclidean") as.matrix(stats::dist(x_eval$centers))
                 else idist(x_eval$centers, x_eval$radii, met)
  }
  for (m in names(projections)) {
    pr <- projections[[m]]
    for (met in mets) {
      Dl[[paste(m, met)]] <- if (met == "Centers-Euclidean" || pr$type == "Point")
        as.matrix(stats::dist(pr$C)) else idist(pr$C, pr$R, met)
    }
  }
  hi_ties <- vapply(mets, function(met) count_row_ties(Dh[[met]]), integer(1))
  acc <- list(); aa <- 0L
  for (met in mets) {
    for (s in seq_len(N_RES)) {
      set.seed(tie_seed + s)
      Rh <- rank_matrix(Dh[[met]], ties = "random")
      for (m in names(projections)) {
        Rl <- rank_matrix(Dl[[paste(m, met)]], ties = "random")
        Qm <- coranking_matrix(Rh, Rl)
        for (K in K_GRID[K_GRID <= nrow(Rh) - 2]) {
          aa <- aa + 1L
          acc[[aa]] <- data.frame(Method = m, Metric = met, K = K, res = s,
                                  t(indices_from_coranking(Qm, K)))
        }
      }
    }
  }
  long <- do.call(rbind, acc)
  idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  cells <- aggregate(long[idx], long[c("Method", "Metric", "K")], mean)
  spread <- aggregate(long[idx], long[c("Method", "Metric", "K")],
                      function(v) diff(range(v)))
  names(spread)[4:9] <- paste0(idx, "_spread")
  lo_ties <- vapply(names(Dl), function(k) count_row_ties(Dl[[k]]), integer(1))
  list(cells = merge(cells, spread), hi_ties = hi_ties, lo_ties = lo_ties)
}

all_rows <- list(); cal_rows <- list(); recs <- list()
t0 <- Sys.time()
for (r in seq_len(N_REP)) {
  seed_r <- SEED_BASE + 7000L + r
  x <- gen_scenario4(seed_r)
  xs <- standardize(x)
  set.seed(seed_r + 500000L)
  rm <- run_methods_timed(xs)
  need_grid <- (r <= N_SUB) || (PILOT > 0)
  t_eval <- system.time(
    ev <- eval_rep(xs, rm$projections,
                   kmax_needed = if (need_grid) max(K_GRID) else K0,
                   tie_seed = seed_r + 300000L)
  )[3]
  cells <- ev$cells
  cells$rep <- r
  all_rows[[r]] <- if (need_grid) cells else cells[cells$K == K0, ]
  recs[[r]] <- list(contaminated = attr(x, "contaminated"),
                    degenerate = attr(x, "degenerate"),
                    donors = attr(x, "donors"), recipients = attr(x, "recipients"),
                    failures = rm$failures, hi_ties = ev$hi_ties,
                    lo_ties = ev$lo_ties, dr_sec = sum(rm$timings),
                    eval_sec = t_eval)
  if (PILOT > 0) {
    cat(sprintf("PILOT rep %d: DR %.0fs eval %.0fs; hi-ties per metric: %s; lo-tie cells>0: %d/30\n",
                r, sum(rm$timings), t_eval,
                paste(ev$hi_ties, collapse = "/"), sum(ev$lo_ties > 0)))
    next
  }
  if (r <= N_SUB) {
    set.seed(seed_r + 300001L)      # documented single resolution s=1 stream
    a <- assess_quality(xs, rm$projections, K = K0, perm_test = TRUE,
                        n_perm = M_PERM, seed = seed_r + 900000L,
                        ties = "random", baseline = TRUE)
    pu <- a$pvalues; pa <- a$pvalues_adj
    pu$rep <- r; pa$rep <- r; pu$kind <- "raw"; pa$kind <- "minP"
    cal_rows[[r]] <- rbind(pu, pa)
  }
  if (r %% 5 == 0) cat(sprintf("rep %d/%d elapsed %.1f min\n", r, N_REP,
                               as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

if (PILOT > 0) {
  cat(sprintf("PILOT projection: %.1f h for 100 reps (grid reps 1-25 heavier)\n",
              (25 * mean(vapply(recs, function(z) z$dr_sec + z$eval_sec, 0)) +
               75 * mean(vapply(recs, function(z) z$dr_sec + z$eval_sec * 0.3, 0))) / 3600))
  quit(status = 0)
}

per <- do.call(rbind, all_rows)
idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
prim <- per[per$K == K0, ]
summ <- aggregate(prim[idx], prim[c("Method", "Metric")],
                  function(v) c(mean = mean(v), sd = sd(v),
                                mcse = sd(v) / sqrt(length(v)),
                                q025 = unname(quantile(v, 0.025)),
                                q975 = unname(quantile(v, 0.975))))
save_result(per, "scenario4_per_rep",
            extra = list(K0 = K0, K_grid = K_GRID, n_rep = N_REP, n_res = N_RES,
                         prereg_sha256 = "2c25e3b2d39ee824043fdf5261e299f13773c5a7734799f7f35ba6257f2df67a"))
save_result(do.call(data.frame, summ), "scenario4_summary",
            extra = list(K0 = K0, n_rep = N_REP))
## E2 contrast: interval-average minus baseline, per method per replication
base <- prim[prim$Metric == "Centers-Euclidean", c("Method", "rep", "Q_TC")]
intv <- aggregate(Q_TC ~ Method + rep, prim[prim$Metric != "Centers-Euclidean", ], mean)
ctr <- merge(intv, base, by = c("Method", "rep"), suffixes = c("_int", "_base"))
ctr$gain <- ctr$Q_TC_int - ctr$Q_TC_base
save_result(aggregate(gain ~ Method, ctr, function(v) c(mean = mean(v), sd = sd(v))) |>
              do.call(what = data.frame),
            "scenario4_baseline_contrast")
saveRDS(recs, file.path(OUTPUT_DIR, "scenario4_records.rds"))
if (length(cal_rows)) save_result(do.call(rbind, cal_rows), "scenario4_calibration",
                                  extra = list(m = M_PERM, n_cal = N_SUB))
cat("Scenario 4 done:", format(Sys.time()), "\n")
