# Table IV + Figures 1-5 (Face data): regenerated with corrected indices,
# center-only baseline, m=999 permutations, min-P adjustment, K profiles,
# effective hyperparameters and seeds recorded.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.R"))
save_sessioninfo("face")

K <- 5L
M_PERM <- n_perm_mode(999L)
if (!requireNamespace("dataSDA", quietly = TRUE)) {
  stop("Package 'dataSDA' is required for the Face analysis.")
}
utils::data("face.iGAP", package = "dataSDA", envir = environment())
face_labels <- sub("[[:digit:]]+$", "", rownames(face.iGAP))
face_data <- interval_data_from_dataSDA(face.iGAP, labels = face_labels)
x <- standardize(face_data)

set.seed(SEED_BASE + 4000L)
rm <- run_methods_timed(x)
if (length(rm$failures)) message("Face DR failures: ", paste(rm$failures, collapse = "; "))

## Primary table (tie-audited cells) + inference
cells <- evaluate_cells(x, rm$projections, K = K, tie_seed = SEED_BASE + 4000L)
a <- assess_quality(x, rm$projections, K = K, perm_test = TRUE, n_perm = M_PERM,
                    seed = SEED_BASE + 4001L, ties = "random", baseline = TRUE)
save_result(cells, "face_indices", extra = list(K = K))
save_result(a$results, "face_assessment_results", extra = list(K = K, m = M_PERM))
save_result(a$pvalues, "face_pvalues_raw", extra = list(m = M_PERM))
save_result(a$pvalues_adj, "face_pvalues_minP", extra = list(m = M_PERM))
saveRDS(a, file.path(OUTPUT_DIR, "face_assessment_full.rds"))

## K profiles (1..25) incl. baseline
prof <- k_profiles(x, rm$projections, K_max = 25, ties = "random", baseline = TRUE)
save_result(prof, "face_k_profiles")

## Effective hyperparameters (incl. what RSDA actually honored for Int-UMAP)
eff <- lapply(rm$projections, function(p) attr(p, "effective_umap_config"))
saveRDS(list(effective_umap = eff, dr_timings = rm$timings,
             seed = SEED_BASE + 4000L),
        file.path(OUTPUT_DIR, "face_effective_params.rds"))

## Figures: K-profile plots per metric + 2-D projections.
## A null device is opened so that grob construction (ggplotGrob inside
## arrangeGrob) cannot open the default device and emit Rplots.pdf.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  grDevices::pdf(NULL)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  for (met in c(METRICS, "Centers-Euclidean")) {
    sub <- prof[prof$Metric == met, ]
    pdf_file <- file.path(OUTPUT_DIR, paste0("RealData_QBK_", gsub("[^A-Za-z-]", "", met), ".pdf"))
    tryCatch({
      p <- plot_k_profiles(sub)
      if (is.list(p) && !grid::is.grob(p) && length(p) >= 1) p <- p[[1]]
      ggplot2::ggsave(pdf_file, p, width = 12, height = 7)
    }, error = function(e) message("profile plot failed for ", met, ": ",
                                   conditionMessage(e)))
  }
  tryCatch({
    p2 <- plot_projections(rm$projections, labels = face_data$labels)
    ggplot2::ggsave(file.path(OUTPUT_DIR, "RealData_2D_Projections.pdf"), p2,
                    width = 12, height = 8)
  }, error = function(e) message("projection plot failed: ", conditionMessage(e)))
}
cat("Face analysis done:", format(Sys.time()), "\n")
