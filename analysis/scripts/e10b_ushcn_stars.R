# E10b follow-up (Prompt E 2026-07-17): USHCN inference re-run under the
# corrected tolerance-block tie randomization. The tie-impact audit found
# that all USHCN IMDS tie blocks are tolerance near-ties (not exact
# duplicates), so the production p-values/stars must be re-derived under the
# corrected policy and compared. Production call replicated exactly
# (seed SEED_BASE + 6001, m = 999, min-P adjustment).
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.R"))

K <- 10L
snap <- readRDS(file.path(DATA_PROC, "ghcn_seasonal_intervals.rds"))
x <- interval_data(snap$centers, snap$radii)
xs <- standardize(x)
set.seed(SEED_BASE + 6000L)
rm <- run_methods_timed(xs)
a <- assess_quality(xs, rm$projections, K = K, perm_test = TRUE, n_perm = 999L,
                    seed = SEED_BASE + 6001L, ties = "random", baseline = TRUE)
new_adj <- a$pvalues_adj
old_adj <- read.csv(file.path(OUTPUT_DIR, "midscale_pvalues_minP.csv"))
idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
mg <- merge(old_adj, new_adj, by = c("IDR", "Metric"),
            suffixes = c(".old", ".new"))
cat("aligned rows:", nrow(mg), "of", nrow(old_adj), "\n")
tot <- 0L
for (ix in idx) {
  po <- mg[[paste0(ix, ".old")]]; pn <- mg[[paste0(ix, ".new")]]
  fl <- which((po <= 0.05) != (pn <= 0.05))
  cat(sprintf("%-5s max |dp| = %.4g  star flips = %d\n",
              ix, max(abs(po - pn)), length(fl)))
  if (length(fl)) print(mg[fl, c("IDR", "Metric",
                                 paste0(ix, c(".old", ".new")))])
  tot <- tot + length(fl)
}
cat("TOTAL STAR FLIPS:", tot, "\n")
write.csv(mg, file.path(OUTPUT_DIR, "e10b_ushcn_stars_compare.csv"),
          row.names = FALSE)
cat("E10B_USHCN_STARS_DONE\n")
