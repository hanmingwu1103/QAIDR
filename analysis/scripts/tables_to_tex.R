# Render manuscript tables from pipeline CSVs into LaTeX fragments
# (analysis/output/tex/). No number is ever hand-transcribed: fragments are
# spliced verbatim into the manuscript. All generated content is wrapped for
# red highlighting at the row level.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
TEX_DIR <- file.path(OUTPUT_DIR, "tex")
dir.create(TEX_DIR, showWarnings = FALSE)

METHOD_ORDER <- c("C-PCA", "V-PCA", "MR-PCA", "SPCA", "IMDS", "Int-UMAP")
METRIC_ORDER <- c("Int-Euclidean", "Hausdorff", "Ichino-Yaguchi", "Wasserstein",
                  "Centers-Euclidean")
METRIC_LABEL <- c("Int-Euclidean" = "Int-Euclidean", "Hausdorff" = "Hausdorff",
                  "Ichino-Yaguchi" = "Ichino-Yaguchi", "Wasserstein" = "Wasserstein",
                  "Centers-Euclidean" = "Centers (baseline)")

## Compact SD notation: parenthetical value is the SD in units of 1e-3,
## rounded to an integer (stated in every caption). Keeps six mean(SD)
## columns within the text width without boxing (sn-jnl tabulars cannot be
## wrapped in resizebox/adjustbox).
fmt_ms <- function(m, s) sprintf("%.3f{\\scriptsize(%d)}", m, round(s * 1000))

## Significance stars from calibration: star iff min-P adjusted p <= 0.05 in
## >= 90% of calibration replicates for that (method, metric, index).
star_map <- function(cal_file) {
  if (!file.exists(cal_file)) return(NULL)
  cal <- read.csv(cal_file)
  cal <- cal[cal$kind == "minP", ]
  idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  ag <- aggregate(cal[idx], cal[c("IDR", "Metric")],
                  function(p) mean(p <= 0.05) >= 0.9)
  ag
}

render_mc_table <- function(summary_csv, cal_csv, tie_csv, out_name) {
  d <- read.csv(summary_csv)
  st <- star_map(cal_csv)
  tf <- if (!is.null(tie_csv) && file.exists(tie_csv)) read.csv(tie_csv) else NULL
  idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  lines <- character(0)
  for (m in METHOD_ORDER) {
    first <- TRUE
    for (g in METRIC_ORDER) {
      r <- d[d$Method == m & d$Metric == g, ]
      if (!nrow(r)) next
      cells <- vapply(idx, function(j) {
        v <- fmt_ms(r[[paste0(j, ".mean")]], r[[paste0(j, ".sd")]])
        starred <- FALSE
        if (!is.null(st)) {
          sr <- st[st$IDR == m & st$Metric == g, ]
          if (nrow(sr) && isTRUE(sr[[j]])) starred <- TRUE
        }
        if (starred) v <- paste0(v, "$^{*}$")
        v
      }, character(1))
      dag <- ""
      if (!is.null(tf)) {
        tr <- tf[tf$Method == m & tf$Metric == g, ]
        if (nrow(tr) && isTRUE(tr$tie_affected)) dag <- "$^{\\dagger}$"
      }
      lead <- if (first) sprintf("\\textbf{%s}", m) else ""
      cells_red <- sprintf("\\textcolor{red}{%s}", cells)
      lines <- c(lines, sprintf("\t\t\t%s & \\textcolor{red}{%s%s} & %s \\\\",
                                lead, METRIC_LABEL[g], dag,
                                paste(cells_red, collapse = " & ")))
      first <- FALSE
    }
    lines <- c(lines, "\t\t\t\\midrule")
  }
  lines <- head(lines, -1)
  writeLines(lines, file.path(TEX_DIR, out_name))
  cat("wrote", out_name, ":", length(lines), "lines\n")
}

## Table I
if (file.exists(file.path(OUTPUT_DIR, "scenario1_summary.csv"))) {
  render_mc_table(file.path(OUTPUT_DIR, "scenario1_summary.csv"),
                  file.path(OUTPUT_DIR, "scenario1_calibration.csv"),
                  file.path(OUTPUT_DIR, "scenario1_tie_flags.csv"),
                  "table1_body.tex")
}
## Table II
if (file.exists(file.path(OUTPUT_DIR, "scenario2_summary.csv"))) {
  render_mc_table(file.path(OUTPUT_DIR, "scenario2_summary.csv"),
                  file.path(OUTPUT_DIR, "scenario2_calibration.csv"),
                  file.path(OUTPUT_DIR, "scenario2_tie_flags.csv"),
                  "table2_body.tex")
}
## Table III (different columns: MAE + flip rates)
if (file.exists(file.path(OUTPUT_DIR, "scenario3_summary.csv"))) {
  d <- read.csv(file.path(OUTPUT_DIR, "scenario3_summary.csv"))
  cols <- c("dQ_TC", "dB_TC", "flip_TC", "dQ_RE", "dB_RE", "flip_RE", "dQ_LC", "dB_LC")
  lines <- character(0)
  for (m in METHOD_ORDER) {
    first <- TRUE
    for (g in METRIC_ORDER[1:4]) {
      r <- d[d$Method == m & d$Metric == g, ]
      if (!nrow(r)) next
      cells <- vapply(cols, function(j) {
        mu <- r[[paste0(j, ".mean")]]
        sdv <- r[[paste0(j, ".sd")]]
        if (startsWith(j, "flip")) sprintf("%.2f", mu) else fmt_ms(mu, sdv)
      }, character(1))
      lead <- if (first) sprintf("\\textbf{%s}", m) else ""
      cells_red <- sprintf("\\textcolor{red}{%s}", cells)
      lines <- c(lines, sprintf("\t\t\t%s & \\textcolor{red}{%s} & %s \\\\",
                                lead, METRIC_LABEL[g], paste(cells_red, collapse = " & ")))
      first <- FALSE
    }
    lines <- c(lines, "\t\t\t\\midrule")
  }
  writeLines(head(lines, -1), file.path(TEX_DIR, "table3_body.tex"))
  cat("wrote table3_body.tex\n")
}
## Table IV (Face) + mid-scale: single-run indices with min-P stars
render_single_table <- function(idx_csv, minp_csv, out_name) {
  d <- read.csv(idx_csv); pm <- read.csv(minp_csv)
  idxn <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  lines <- character(0)
  for (m in METHOD_ORDER) {
    first <- TRUE
    for (g in METRIC_ORDER) {
      r <- d[d$Method == m & d$Metric == g, ]
      if (!nrow(r)) next
      p <- pm[pm$IDR == m & pm$Metric == g, ]
      cells <- vapply(idxn, function(j) {
        v <- sprintf("%.3f", r[[j]])
        if (nrow(p) && p[[j]] <= 0.05) v <- paste0(v, "$^{*}$")
        v
      }, character(1))
      dag <- if (isTRUE(r$tie_affected)) "$^{\\dagger}$" else ""
      lead <- if (first) sprintf("\\textbf{%s}", m) else ""
      cells_red <- sprintf("\\textcolor{red}{%s}", cells)
      lines <- c(lines, sprintf("\t\t\t%s & \\textcolor{red}{%s%s} & %s \\\\",
                                lead, METRIC_LABEL[g], dag, paste(cells_red, collapse = " & ")))
      first <- FALSE
    }
    lines <- c(lines, "\t\t\t\\midrule")
  }
  writeLines(head(lines, -1), file.path(TEX_DIR, out_name))
  cat("wrote", out_name, "\n")
}
if (file.exists(file.path(OUTPUT_DIR, "face_indices.csv"))) {
  render_single_table(file.path(OUTPUT_DIR, "face_indices.csv"),
                      file.path(OUTPUT_DIR, "face_pvalues_minP.csv"),
                      "table4_body.tex")
}
if (file.exists(file.path(OUTPUT_DIR, "midscale_indices.csv"))) {
  render_single_table(file.path(OUTPUT_DIR, "midscale_indices.csv"),
                      file.path(OUTPUT_DIR, "midscale_pvalues_minP.csv"),
                      "midscale_body.tex")
}
## Benchmark table
if (file.exists(file.path(OUTPUT_DIR, "scalability_benchmark.csv"))) {
  b <- read.csv(file.path(OUTPUT_DIR, "scalability_benchmark.csv"))
  lines <- vapply(seq_len(nrow(b)), function(i) {
    vals <- c(sprintf("%d", b$n[i]), sprintf("%.2f", b$t_dist[i]),
              sprintf("%.2f", b$t_rank[i]),
              sprintf("%.2f", b$t_indices[i] + b$t_perm99[i] * 999 / 99),
              sprintf("%.1f", b$approx_matrix_mem_MB[i]),
              sprintf("%.0f", sum(b[i, intersect(METHOD_ORDER, names(b))], na.rm = TRUE)))
    paste0("\t\t\t", paste(sprintf("\\textcolor{red}{%s}", vals), collapse = " & "), " \\\\")
  }, character(1))
  writeLines(lines, file.path(TEX_DIR, "benchmark_body.tex"))
  cat("wrote benchmark_body.tex\n")
}
## Scenario IV (B2)
if (file.exists(file.path(OUTPUT_DIR, "scenario4_summary.csv"))) {
  render_mc_table(file.path(OUTPUT_DIR, "scenario4_summary.csv"),
                  file.path(OUTPUT_DIR, "scenario4_calibration.csv"),
                  NULL, "table_sc4_body.tex")
}
## B2 K-stability (compact: per scenario x K, worst and median Kendall tau
## of the method ordering vs K0, across all metrics x indices)
if (file.exists(file.path(OUTPUT_DIR, "b2_k_stability_kendall.csv"))) {
  kt <- read.csv(file.path(OUTPUT_DIR, "b2_k_stability_kendall.csv"))
  agg <- aggregate(tau.mean ~ scenario + K, kt,
                   function(v) c(min = min(v), med = median(v)))
  agg <- do.call(data.frame, agg)
  lines <- vapply(seq_len(nrow(agg)), function(i) {
    vals <- sprintf("\\textcolor{red}{%s}",
                    c(agg$scenario[i], agg$K[i],
                      sprintf("%.2f", agg$tau.mean.min[i]),
                      sprintf("%.2f", agg$tau.mean.med[i])))
    paste0("\t\t\t", paste(vals, collapse = " & "), " \\\\")
  }, character(1))
  writeLines(lines, file.path(TEX_DIR, "table_kstab_body.tex"))
  cat("wrote table_kstab_body.tex\n")
}
## B2 sensitivity (Face exact rows + grid rows)
if (file.exists(file.path(OUTPUT_DIR, "b2_sensitivity.csv"))) {
  sd_ <- read.csv(file.path(OUTPUT_DIR, "b2_sensitivity.csv"))
  rngc <- grep("^range_", names(sd_))
  sd_$max_range <- apply(sd_[rngc], 1, max)
  lines <- vapply(seq_len(nrow(sd_)), function(i) {
    vals <- sprintf("\\textcolor{red}{%s}",
                    c(sd_$arm[i], if (sd_$par[i] == "lambda") "$\\lambda$" else "$\\nu$",
                      sd_$Method[i],
                      ifelse(is.na(sd_$n_breakpoints[i]), "grid",
                             format(sd_$n_breakpoints[i], big.mark = "{,}")),
                      format(sd_$n_distinct_rankings[i], big.mark = "{,}"),
                      sprintf("%.3f", sd_$max_range[i])))
    paste0("\t\t\t", paste(vals, collapse = " & "), " \\\\")
  }, character(1))
  writeLines(lines, file.path(TEX_DIR, "table_sens_body.tex"))
  cat("wrote table_sens_body.tex\n")
}
cat("tables_to_tex done\n")
