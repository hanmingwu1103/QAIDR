# Mid-scale application: USHCN seasonal temperature intervals (selected arm,
# author-approved n=645; selection frozen in REVISION_WORKLOG.md BEFORE any
# index computation). Full corrected analysis: 6 DR methods x (4 interval
# dissimilarities + center baseline), K profiles, m=999 permutation inference
# with min-P adjustment, runtime/memory notes.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.R"))
save_sessioninfo("midscale_ghcn")

snap <- readRDS(file.path(DATA_PROC, "ghcn_seasonal_intervals.rds"))
x <- interval_data(snap$centers, snap$radii)
xs <- standardize(x)
N <- nrow(snap$centers)
K <- 10L
M_PERM <- n_perm_mode(999L)
cat("USHCN snapshot:", N, "stations, year", snap$target_year, "\n")

set.seed(SEED_BASE + 6000L)
t_dr <- system.time(rm <- run_methods_timed(xs))[3]
cat("DR total:", round(t_dr, 1), "s; per-method:",
    paste(names(rm$timings), round(rm$timings, 1), collapse = " "), "\n")
if (length(rm$failures)) message("failures: ", paste(rm$failures, collapse = "; "))

t_cells <- system.time(
  cells <- evaluate_cells(xs, rm$projections, K = K, tie_seed = SEED_BASE + 6000L)
)[3]
t_inf <- system.time(
  a <- assess_quality(xs, rm$projections, K = K, perm_test = TRUE, n_perm = M_PERM,
                      seed = SEED_BASE + 6001L, ties = "random", baseline = TRUE)
)[3]
prof <- k_profiles(xs, rm$projections, K_max = 25, ties = "random", baseline = TRUE)

save_result(cells, "midscale_indices", extra = list(K = K, n = N, year = snap$target_year))
save_result(a$results, "midscale_assessment_results", extra = list(K = K, m = M_PERM))
save_result(a$pvalues, "midscale_pvalues_raw")
save_result(a$pvalues_adj, "midscale_pvalues_minP")
save_result(prof, "midscale_k_profiles")
save_result(data.frame(component = c("DR_all6", "indices_30cells", "inference_m999"),
                       seconds = c(t_dr, t_cells, t_inf)),
            "midscale_runtime", extra = list(n = N, p = ncol(snap$centers)))
saveRDS(list(effective_umap = lapply(rm$projections, function(p) attr(p, "effective_umap_config")),
             dr_timings = rm$timings, seed = SEED_BASE + 6000L),
        file.path(OUTPUT_DIR, "midscale_effective_params.rds"))
cat("Mid-scale analysis done:", format(Sys.time()), "\n")
