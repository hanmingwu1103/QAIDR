# Prompt H delta gate (contract 1.4): assert that regeneration changed ONLY
# adjusted p-values; report the old-vs-new adjusted-p delta and star flips.
# Usage: Rscript h_delta_check.r <output_dir>
args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args)) args[1] else "."
SNAP <- file.path(OUT, "h_old_inference_snapshot")
pc <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")

same <- function(f) {
  o <- readRDS(file.path(SNAP, f)); n <- readRDS(file.path(OUT, f))
  isTRUE(all.equal(o, n, tolerance = 1e-15, check.attributes = FALSE))
}
adj_delta <- function(f, label) {
  o <- readRDS(file.path(SNAP, f)); n <- readRDS(file.path(OUT, f))
  om <- as.matrix(o[pc]); nm <- as.matrix(n[pc])
  d <- abs(om - nm)
  cat(sprintf("%s: max|delta|=%.6f cells_changed=%d/%d\n", label, max(d),
              sum(d > 1e-12), length(d)))
  os <- om <= 0.05; ns <- nm <= 0.05
  ch <- which(os != ns, arr.ind = TRUE)
  cat(sprintf("%s: star_flips=%d\n", label, nrow(ch)))
  if (nrow(ch)) for (i in seq_len(nrow(ch)))
    cat(sprintf("  flip %s | %s | %s : %.4f -> %.4f\n",
                o$IDR[ch[i, 1]], o$Metric[ch[i, 1]], pc[ch[i, 2]],
                om[ch[i, 1], ch[i, 2]], nm[ch[i, 1], ch[i, 2]]))
  invisible(NULL)
}
cal_check <- function(f, label) {
  o <- readRDS(file.path(SNAP, f)); n <- readRDS(file.path(OUT, f))
  or <- o[o$kind == "raw", ]; nr <- n[n$kind == "raw", ]
  or <- or[order(or$rep, or$IDR, or$Metric), ]
  nr <- nr[order(nr$rep, nr$IDR, nr$Metric), ]
  raw_ok <- nrow(or) == nrow(nr) &&
    isTRUE(all.equal(as.matrix(or[pc]), as.matrix(nr[pc]),
                     tolerance = 1e-12, check.attributes = FALSE))
  cat(sprintf("%s: raw_identical=%s\n", label, raw_ok))
  oa <- o[o$kind == "minP", ]; na <- n[n$kind == "minP", ]
  oa <- oa[order(oa$rep, oa$IDR, oa$Metric), ]
  na <- na[order(na$rep, na$IDR, na$Metric), ]
  om <- as.matrix(oa[pc]); nm <- as.matrix(na[pc])
  d <- abs(om - nm)
  rej_o <- colMeans(matrix(apply(om <= 0.05, 1, any), ncol = 1))
  cat(sprintf("%s: adj max|delta|=%.6f cells_changed=%d/%d any-rej reps old=%d new=%d of %d\n",
              label, max(d), sum(d > 1e-12), length(d),
              sum(tapply(apply(om <= 0.05, 1, any), oa$rep, any)),
              sum(tapply(apply(nm <= 0.05, 1, any), na$rep, any)),
              length(unique(oa$rep))))
  invisible(NULL)
}

if (file.exists(file.path(OUT, "face_pvalues_raw.rds"))) {
  cat("== Face ==\n")
  cat("raw identical:", same("face_pvalues_raw.rds"), "\n")
  cat("results identical:", same("face_assessment_results.rds"), "\n")
  cat("indices identical:", same("face_indices.rds"), "\n")
  cat("k_profiles identical:", same("face_k_profiles.rds"), "\n")
  adj_delta("face_pvalues_minP.rds", "face_adj")
}
if (file.exists(file.path(OUT, "midscale_pvalues_raw.rds"))) {
  cat("== USHCN midscale ==\n")
  cat("raw identical:", same("midscale_pvalues_raw.rds"), "\n")
  cat("results identical:", same("midscale_assessment_results.rds"), "\n")
  cat("indices identical:", same("midscale_indices.rds"), "\n")
  cat("k_profiles identical:", same("midscale_k_profiles.rds"), "\n")
  adj_delta("midscale_pvalues_minP.rds", "midscale_adj")
}
for (s in c("scenario1", "scenario2", "scenario4")) {
  f <- paste0(s, "_calibration.rds")
  if (file.exists(file.path(OUT, f))) {
    cat("==", s, "calibration ==\n")
    cal_check(f, s)
  }
}
