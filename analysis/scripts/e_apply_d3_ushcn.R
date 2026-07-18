# Apply author decision D3(b) (Prompt E, 2026-07-17): adopt the
# corrected-policy USHCN inference. Re-runs the production inference call
# under the corrected tolerance-block tie randomization (production seeds)
# and REPLACES the midscale p-value artifacts. Index values, K-profiles,
# and runtime artifacts are NOT touched (their display values are proven
# invariant; runtimes are machine-time measurements quoted in the
# manuscript). Provenance recorded in the artifact metadata.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.R"))
save_sessioninfo("e_apply_d3_ushcn")

K <- 10L
snap <- readRDS(file.path(DATA_PROC, "ghcn_seasonal_intervals.rds"))
x <- interval_data(snap$centers, snap$radii)
xs <- standardize(x)
set.seed(SEED_BASE + 6000L)
rm <- run_methods_timed(xs)
a <- assess_quality(xs, rm$projections, K = K, perm_test = TRUE, n_perm = 999L,
                    seed = SEED_BASE + 6001L, ties = "random", baseline = TRUE)
prov <- list(K = K, m = 999L,
  provenance = paste("Prompt E D3(b) 2026-07-17: corrected tolerance-block",
                     "tie randomization (rank_matrix patch); production",
                     "call and seeds otherwise identical; supersedes the",
                     "2026-07-16 artifacts computed under the defective",
                     "near-tie randomization"))
save_result(a$pvalues, "midscale_pvalues_raw", extra = prov)
save_result(a$pvalues_adj, "midscale_pvalues_minP", extra = prov)
## sanity: exactly the one expected star change vs the archived comparison
cmp <- read.csv(file.path(OUTPUT_DIR, "e10b_ushcn_stars_compare.csv"))
idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
new <- a$pvalues_adj
flips <- 0L
for (ix in idx) {
  mg <- merge(cmp[, c("IDR", "Metric", paste0(ix, ".old"))], new[, c("IDR", "Metric", ix)],
              by = c("IDR", "Metric"))
  flips <- flips + sum((mg[[paste0(ix, ".old")]] <= 0.05) != (mg[[ix]] <= 0.05))
}
cat("star changes vs published:", flips, "(expected 1)\n")
stopifnot(flips == 1L)
cat("D3_APPLIED\n")
