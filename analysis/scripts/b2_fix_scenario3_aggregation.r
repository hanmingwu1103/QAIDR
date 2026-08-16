# B2 repair (Codex audit finding): scenario3_summary.csv aggregated SDs over
# all 10,000 dataset x perturbation rows while the manuscript caption claims
# between-dataset SDs. Recompute correctly from the untouched
# scenario3_per_perturbation.csv: first average within dataset (dataset-level
# means), then take mean/SD/MCSE across the 100 dataset means. The defective
# summary is preserved as *_pooled_deprecated.csv for audit.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
save_sessioninfo("b2_fix_scenario3")

per <- read.csv(file.path(OUTPUT_DIR, "scenario3_per_perturbation.csv"))
old <- file.path(OUTPUT_DIR, "scenario3_summary.csv")
if (file.exists(old)) {
  file.copy(old, file.path(OUTPUT_DIR, "scenario3_summary_pooled_deprecated.csv"),
            overwrite = TRUE)
}
dcols <- grep("^d|^flip", names(per), value = TRUE)
## dataset-level means first
ds <- aggregate(per[dcols], per[c("Method", "Metric", "dataset")], mean)
## then across-dataset summaries
summ <- aggregate(ds[dcols], ds[c("Method", "Metric")],
                  function(v) c(mean = mean(v), sd = sd(v),
                                mcse = sd(v) / sqrt(length(v))))
out <- do.call(data.frame, summ)
save_result(out, "scenario3_summary",
            inputs = file.path(OUTPUT_DIR, "scenario3_per_perturbation.csv"),
            extra = list(aggregation = "between-dataset SD of dataset means (B2 repair)",
                         n_datasets = length(unique(per$dataset)),
                         prereg = "B2_PREREGISTRATION.md A-Repairs"))
cat("Rebuilt scenario3_summary.csv with between-dataset aggregation.\n")
## sanity echo for the two Codex-cited cells
r <- out[out$Method == "C-PCA" & out$Metric == "Hausdorff", ]
cat(sprintf("C-PCA/Hausdorff: dQ_TC mean=%.4f sd=%.6f ; dQ_RE mean=%.4f sd=%.6f\n",
            r$dQ_TC.mean, r$dQ_TC.sd, r$dQ_RE.mean, r$dQ_RE.sd))
