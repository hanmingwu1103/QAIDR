# B2 K-profile evidence (preregistered A-D1, frozen 2026-07-16).
# Scenarios I-III: replications 1..25 re-derived from the original seeds
# (identical RNG offsets to scenario{1,2,3}.R); K grid = {0.5,1,2,4} x K0,
# K0 = 10; all 6 methods x (4 dissimilarities + baseline) x 6 indices.
# Scenario IV grid cells are produced by b2_scenario4.R in the same pass and
# are summarized here if present. Kendall-tau stability per prereg.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.R"))
save_sessioninfo("b2_k_profiles")

K0 <- 10L; K_GRID <- c(5L, 10L, 20L, 40L); N_SUB <- 25L
idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")

profile_rep <- function(x_eval, projections, tie_seed) {
  mets <- c(METRICS, "Centers-Euclidean")
  rows <- list(); ii <- 0L
  for (met in mets) {
    Dh <- if (met == "Centers-Euclidean") as.matrix(stats::dist(x_eval$centers))
          else idist(x_eval$centers, x_eval$radii, met)
    for (m in names(projections)) {
      pr <- projections[[m]]
      Dl <- if (met == "Centers-Euclidean" || pr$type == "Point")
        as.matrix(stats::dist(pr$C)) else idist(pr$C, pr$R, met)
      cell <- tryCatch({
        Rh <- rank_matrix(Dh); Rl <- rank_matrix(Dl)
        Qm <- coranking_matrix(Rh, Rl)
        lapply(K_GRID, function(K) indices_from_coranking(Qm, K))
      }, error = function(e) {   # structural ties -> audited multi-resolution mean
        res <- array(NA_real_, c(50L, length(K_GRID), 6L))
        for (s in 1:50) {
          set.seed(tie_seed + s)
          Rh <- rank_matrix(Dh, ties = "random"); Rl <- rank_matrix(Dl, ties = "random")
          Qm <- coranking_matrix(Rh, Rl)
          for (ki in seq_along(K_GRID)) res[s, ki, ] <- indices_from_coranking(Qm, K_GRID[ki])
        }
        lapply(seq_along(K_GRID), function(ki) setNames(colMeans(res[, ki, ]), idx))
      })
      for (ki in seq_along(K_GRID)) {
        ii <- ii + 1L
        rows[[ii]] <- data.frame(Method = m, Metric = met, K = K_GRID[ki],
                                 t(cell[[ki]]))
      }
    }
  }
  do.call(rbind, rows)
}

run_scen <- function(gen, seed_off, dr_off, raw_dr, label) {
  out <- list()
  for (r in seq_len(N_SUB)) {
    seed_r <- SEED_BASE + seed_off + r
    x <- gen(seed_r); xs <- standardize(x)
    set.seed(seed_r + dr_off)
    rm <- run_methods_timed(if (raw_dr) x else xs)
    pr <- profile_rep(xs, rm$projections, tie_seed = seed_r + 300000L)
    pr$rep <- r; out[[r]] <- pr
    if (r %% 5 == 0) cat(label, "rep", r, "/", N_SUB, "\n")
  }
  do.call(rbind, out)
}

res <- list(
  I   = run_scen(gen_scenario1, 1000L, 500000L, FALSE, "Sc I"),
  II  = run_scen(gen_scenario2, 2000L, 500000L, TRUE,  "Sc II"),
  III = run_scen(gen_scenario3, 3000L, 500000L, FALSE, "Sc III")
)
long <- do.call(rbind, Map(function(d, s) { d$scenario <- s; d }, res, names(res)))
## include Scenario IV grid rows if b2_scenario4.R has produced them
s4 <- file.path(OUTPUT_DIR, "scenario4_per_rep.csv")
if (file.exists(s4)) {
  d4 <- read.csv(s4); d4 <- d4[d4$K %in% K_GRID & d4$rep <= N_SUB,
                               c("Method", "Metric", "K", idx, "rep")]
  d4$scenario <- "IV"; long <- rbind(long, d4)
}
summ <- aggregate(long[idx], long[c("scenario", "Method", "Metric", "K")],
                  function(v) c(mean = mean(v), sd = sd(v),
                                mcse = sd(v) / sqrt(length(v)),
                                q025 = unname(quantile(v, 0.025)),
                                q975 = unname(quantile(v, 0.975))))

## Kendall tau of method ordering at K vs K0, per replication/metric/index
taus <- list(); tt <- 0L
for (sc in unique(long$scenario)) for (met in unique(long$Metric))
  for (jj in idx) for (K in setdiff(K_GRID, K0)) {
    for (r in seq_len(N_SUB)) {
      a <- long[long$scenario == sc & long$Metric == met & long$rep == r, ]
      v0 <- a[a$K == K0, ]; vK <- a[a$K == K, ]
      if (!nrow(v0) || !nrow(vK)) next
      mm <- intersect(v0$Method, vK$Method)
      if (length(mm) < 3) next
      tau <- suppressWarnings(cor(v0[[jj]][match(mm, v0$Method)],
                                  vK[[jj]][match(mm, vK$Method)],
                                  method = "kendall"))
      tt <- tt + 1L
      taus[[tt]] <- data.frame(scenario = sc, Metric = met, index = jj,
                               K = K, rep = r, tau = tau)
    }
  }
tau_df <- do.call(rbind, taus)
tau_summ <- aggregate(tau ~ scenario + Metric + index + K, tau_df,
                      function(v) c(mean = mean(v, na.rm = TRUE),
                                    q025 = unname(quantile(v, 0.025, na.rm = TRUE)),
                                    q975 = unname(quantile(v, 0.975, na.rm = TRUE))))

save_result(long, "b2_k_profiles_per_rep",
            extra = list(K_grid = K_GRID, K0 = K0, n_sub = N_SUB,
                         prereg_sha256 = "2c25e3b2d39ee824043fdf5261e299f13773c5a7734799f7f35ba6257f2df67a"))
save_result(do.call(data.frame, summ), "b2_k_profiles_summary")
save_result(do.call(data.frame, tau_summ), "b2_k_stability_kendall")
cat("B2 K profiles done:", format(Sys.time()), "\n")
