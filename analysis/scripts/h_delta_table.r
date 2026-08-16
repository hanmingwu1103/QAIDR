# Prompt H (contract 1.4 + A2): complete old-vs-new inference delta table.
# For Face/USHCN: per-cell old vs new adjusted p and star; raw deltas
# recorded (A2 re-baseline). For the scenario calibrations: per-cell old vs
# new rejection fraction across the 25 replicates (adjusted p <= 0.05), and
# max |delta| of raw and adjusted p over replicates.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
save_sessioninfo("h_delta_table")

SNAP <- file.path(OUTPUT_DIR, "h_old_inference_snapshot")
pcix <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
rows <- list(); rr <- 0L
add <- function(...) { rr <<- rr + 1L; rows[[rr]] <<- data.frame(..., stringsAsFactors = FALSE) }

single <- function(ds, raw_f, adj_f) {
  or <- readRDS(file.path(SNAP, raw_f)); nr <- readRDS(file.path(OUTPUT_DIR, raw_f))
  oa <- readRDS(file.path(SNAP, adj_f)); na <- readRDS(file.path(OUTPUT_DIR, adj_f))
  for (i in seq_len(nrow(oa))) for (col in pcix) {
    o <- oa[i, col]; n <- na[na$IDR == oa$IDR[i] & na$Metric == oa$Metric[i], col]
    orv <- or[i, col]; nrv <- nr[nr$IDR == oa$IDR[i] & nr$Metric == oa$Metric[i], col]
    add(dataset = ds, IDR = oa$IDR[i], Metric = oa$Metric[i], test = col,
        old_raw = orv, new_raw = nrv, raw_delta = nrv - orv,
        old_adj = o, new_adj = n, adj_delta = n - o,
        old_marker = o <= 0.05, new_marker = n <= 0.05,
        marker_flip = (o <= 0.05) != (n <= 0.05))
  }
}
calib <- function(ds, f) {
  o <- readRDS(file.path(SNAP, f)); n <- readRDS(file.path(OUTPUT_DIR, f))
  om <- o[o$kind == "minP", ]; nm <- n[n$kind == "minP", ]
  orw <- o[o$kind == "raw", ]; nrw <- n[n$kind == "raw", ]
  key <- unique(om[c("IDR", "Metric")])
  for (i in seq_len(nrow(key))) for (col in pcix) {
    so <- om[om$IDR == key$IDR[i] & om$Metric == key$Metric[i], ]
    sn <- nm[nm$IDR == key$IDR[i] & nm$Metric == key$Metric[i], ]
    so <- so[order(so$rep), ]; sn <- sn[order(sn$rep), ]
    ro <- orw[orw$IDR == key$IDR[i] & orw$Metric == key$Metric[i], ]
    rn <- nrw[nrw$IDR == key$IDR[i] & nrw$Metric == key$Metric[i], ]
    ro <- ro[order(ro$rep), ]; rn <- rn[order(rn$rep), ]
    fo <- mean(so[[col]] <= 0.05); fn <- mean(sn[[col]] <= 0.05)
    add(dataset = ds, IDR = key$IDR[i], Metric = key$Metric[i], test = col,
        old_raw = NA_real_, new_raw = NA_real_,
        raw_delta = max(abs(rn[[col]] - ro[[col]])),
        old_adj = fo, new_adj = fn, adj_delta = max(abs(sn[[col]] - so[[col]])),
        old_marker = fo >= 0.9, new_marker = fn >= 0.9,
        marker_flip = (fo >= 0.9) != (fn >= 0.9))
  }
}

single("Face", "face_pvalues_raw.rds", "face_pvalues_minP.rds")
single("USHCN", "midscale_pvalues_raw.rds", "midscale_pvalues_minP.rds")
calib("Scenario1", "scenario1_calibration.rds")
calib("Scenario2", "scenario2_calibration.rds")
if (file.exists(file.path(OUTPUT_DIR, "h_s4_done.flag")) ||
    !identical(readRDS(file.path(SNAP, "scenario4_calibration.rds")),
               readRDS(file.path(OUTPUT_DIR, "scenario4_calibration.rds"))))
  calib("Scenario4", "scenario4_calibration.rds")

del <- do.call(rbind, rows)
save_result(del, "h_delta_table",
            extra = list(note = paste(
              "A2 re-baseline: for Face/USHCN, old = stored production values,",
              "new = committed-pipeline QAIDR 0.3.0 regeneration; raw deltas",
              "reflect the pre-existing pipeline re-baselining, adjusted deltas",
              "additionally reflect the symmetric min-P correction. For",
              "calibrations, old/new_adj columns hold the rejection fraction",
              "across 25 replicates and old/new_marker the legacy 90% star rule",
              "for comparison; the manuscript now displays fractions directly."),
              prereg = "H_PREREGISTRATION.md sec 1.4 (post-A2)"))
cat("rows:", nrow(del), "\n")
cat("Face marker flips:", sum(del$marker_flip[del$dataset == "Face"]), "\n")
cat("USHCN marker flips:", sum(del$marker_flip[del$dataset == "USHCN"]), "\n")
cat("Scenario fraction-marker flips:",
    sum(del$marker_flip[grepl("Scenario", del$dataset)]), "\n")
