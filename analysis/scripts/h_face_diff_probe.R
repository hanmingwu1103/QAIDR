pc <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
SNAP <- "analysis/output/h_old_inference_snapshot"
for (f in c("midscale_indices.rds", "midscale_pvalues_raw.rds",
            "midscale_assessment_results.rds", "midscale_k_profiles.rds")) {
  o <- readRDS(file.path(SNAP, f))
  n <- readRDS(file.path("analysis/output", f))
  num <- intersect(pc, names(o))
  d <- abs(as.matrix(o[num]) - as.matrix(n[num]))
  cat(f, "max|d|=", format(max(d), digits = 6), "cells>0:", sum(d > 0),
      "of", length(d), "\n")
  if (max(d) > 0) {
    w <- which(d == max(d), arr.ind = TRUE)[1, ]
    lab <- if ("IDR" %in% names(o)) o$IDR[w[1]] else o$IDR[w[1]]
    cat("  worst:", lab, o$Metric[w[1]],
        if ("K" %in% names(o)) paste0("K=", o$K[w[1]]) else "",
        num[w[2]], as.matrix(o[num])[w[1], w[2]], "vs",
        as.matrix(n[num])[w[1], w[2]], "\n")
  }
  ae <- all.equal(o, n, check.attributes = FALSE)
  if (!isTRUE(ae)) print(utils::head(ae, 4))
}
