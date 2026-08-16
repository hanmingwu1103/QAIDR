pc <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
SNAP <- "analysis/output/h_old_inference_snapshot"
for (f in c("face_indices.rds", "face_pvalues_raw.rds")) {
  o <- readRDS(file.path(SNAP, f))
  n <- readRDS(file.path("analysis/output", f))
  key <- if ("Method" %in% names(o)) "Method" else "IDR"
  d <- abs(as.matrix(o[pc]) - as.matrix(n[pc]))
  agg <- aggregate(rowSums(d > 0), list(m = o[[key]]), sum)
  cat("--", f, "cells>0 by method:\n"); print(agg)
}
oe <- readRDS(file.path(SNAP, "..", "face_effective_params.rds"))
cat("stored seed:", oe$seed, "\n")
cat("dr timings stored:", paste(names(oe$dr_timings),
    round(oe$dr_timings, 2), collapse = " "), "\n")
