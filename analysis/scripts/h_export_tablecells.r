# Prompt H Stage 5: export per-cell table data for the manuscript tables.
# For simulation tables: rejection fractions across the 25 calibration
# replicates (adjusted p <= 0.05) under the corrected symmetric min-P, plus
# Clopper-Pearson 95% intervals (archived), plus chance-adjusted Qc of the
# displayed mean for the three Q columns. For Face/USHCN: stars from the
# regenerated adjusted p-values, Qc of the displayed values, and (Face) the
# re-baselined index values themselves.
# Usage: Rscript h_export_tablecells.r <which>  (s1|s2|s4|face|ushcn)
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
which_tab <- commandArgs(trailingOnly = TRUE)[1]

pcix <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
qcols <- c("Q_TC", "Q_RE", "Q_LC")

mu0 <- function(n, K) {
  GK <- if (K < n / 2) n * K * (2 * n - 3 * K - 1) else n * (n - K) * (n - K - 1)
  HK <- n * sum(abs(n - 2 * seq_len(K) + 1) / seq_len(K))
  j <- seq_len(K)
  EWn <- n / ((n - 1) * HK) * sum((1 / j) * (j * (j - 1) / 2 + (n - 1 - j) * (n - j) / 2))
  c(Q_TC = 1 - n * K * (n - 1 - K) * (n - K) / (GK * (n - 1)),
    Q_RE = 1 - EWn, Q_LC = K / (n - 1))
}

frac_tab <- function(cal_file, n, K) {
  cal <- readRDS(file.path(OUTPUT_DIR, cal_file))
  pa <- cal[cal$kind == "minP", ]
  agg <- aggregate(pa[pcix], pa[c("IDR", "Metric")],
                   function(p) mean(p <= 0.05))
  cnt <- aggregate(pa[pcix], pa[c("IDR", "Metric")],
                   function(p) sum(p <= 0.05))
  nrep <- length(unique(pa$rep))
  cp <- function(x) c(stats::binom.test(x, nrep)$conf.int)
  rows <- list(); k <- 0L
  for (i in seq_len(nrow(agg))) for (col in pcix) {
    k <- k + 1L
    ci <- cp(cnt[i, col])
    rows[[k]] <- data.frame(IDR = agg$IDR[i], Metric = agg$Metric[i],
                            col = col, frac = agg[i, col],
                            cp_lo = ci[1], cp_hi = ci[2],
                            qc = NA_real_, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  m0 <- mu0(n, K)
  out
}

if (which_tab %in% c("s1", "s2", "s4")) {
  cfg <- list(s1 = list(f = "scenario1_calibration.rds", n = 300L, K = 10L,
                        sum = "scenario1_summary.rds", name = "sim1"),
              s2 = list(f = "scenario2_calibration.rds", n = 800L, K = 10L,
                        sum = "scenario2_summary.rds", name = "sim2"),
              s4 = list(f = "scenario4_calibration.rds", n = 300L, K = 10L,
                        sum = "scenario4_summary.rds", name = "sim4"))[[which_tab]]
  out <- frac_tab(cfg$f, cfg$n, cfg$K)
  ## Qc of the displayed (100-rep) mean for Q columns
  m0 <- mu0(cfg$n, cfg$K)
  summ <- readRDS(file.path(OUTPUT_DIR, cfg$sum))
  mean_cols <- grep("\\.mean$", names(summ), value = TRUE)
  for (i in seq_len(nrow(out))) {
    if (out$col[i] %in% qcols) {
      r <- summ[summ$Method == out$IDR[i] & summ$Metric == out$Metric[i], ]
      mc <- paste0(out$col[i], ".mean")
      if (!mc %in% names(summ)) mc <- out$col[i]
      v <- if (is.data.frame(r[[mc]])) r[[mc]][[1]] else r[[mc]]
      out$qc[i] <- (v - m0[out$col[i]]) / (1 - m0[out$col[i]])
    }
  }
  utils::write.csv(out, file.path(OUTPUT_DIR, paste0("h_cells_", cfg$name, ".csv")),
                   row.names = FALSE)
  cat("wrote h_cells_", cfg$name, ".csv  rows:", nrow(out), "\n", sep = "")
} else if (which_tab == "face") {
  res <- readRDS(file.path(OUTPUT_DIR, "face_assessment_results.rds"))
  adj <- readRDS(file.path(OUTPUT_DIR, "face_pvalues_minP.rds"))
  old <- readRDS(file.path(OUTPUT_DIR, "h_old_inference_snapshot",
                           "face_assessment_results.rds"))
  m0 <- mu0(27L, 5L)
  rows <- list(); k <- 0L
  for (i in seq_len(nrow(res))) for (col in pcix) {
    k <- k + 1L
    a <- adj[adj$IDR == res$IDR[i] & adj$Metric == res$Metric[i], col]
    o <- old[old$IDR == res$IDR[i] & old$Metric == res$Metric[i], col]
    rows[[k]] <- data.frame(
      IDR = res$IDR[i], Metric = res$Metric[i], col = col,
      value = res[i, col], star = a <= 0.05,
      qc = if (col %in% qcols) (res[i, col] - m0[col]) / (1 - m0[col]) else NA_real_,
      changed = round(res[i, col], 3) != round(o, 3),
      stringsAsFactors = FALSE)
  }
  utils::write.csv(do.call(rbind, rows),
                   file.path(OUTPUT_DIR, "h_cells_face.csv"), row.names = FALSE)
  cat("wrote h_cells_face.csv\n")
} else if (which_tab == "ushcn") {
  res <- readRDS(file.path(OUTPUT_DIR, "midscale_assessment_results.rds"))
  idx <- readRDS(file.path(OUTPUT_DIR, "midscale_indices.rds"))
  adj <- readRDS(file.path(OUTPUT_DIR, "midscale_pvalues_minP.rds"))
  m0 <- mu0(645L, 10L)
  rows <- list(); k <- 0L
  for (i in seq_len(nrow(idx))) for (col in pcix) {
    k <- k + 1L
    a <- adj[adj$IDR == idx$Method[i] & adj$Metric == idx$Metric[i], col]
    rows[[k]] <- data.frame(
      IDR = idx$Method[i], Metric = idx$Metric[i], col = col,
      value = idx[i, col], star = a <= 0.05,
      dagger = idx$tie_affected[i],
      qc = if (col %in% qcols) (idx[i, col] - m0[col]) / (1 - m0[col]) else NA_real_,
      stringsAsFactors = FALSE)
  }
  utils::write.csv(do.call(rbind, rows),
                   file.path(OUTPUT_DIR, "h_cells_ushcn.csv"), row.names = FALSE)
  cat("wrote h_cells_ushcn.csv\n")
}
