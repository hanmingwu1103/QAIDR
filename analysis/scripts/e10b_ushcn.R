# E10b tie-impact audit, USHCN flagged cells (Prompt E 2026-07-17).
# Four flagged cells: IMDS x {Int-Euclidean, Hausdorff, Ichino-Yaguchi,
# Centers-Euclidean}. Single dataset (no replications).
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.R"))
source(file.path(SCRIPT_DIR, "e10b_lib.R"))

K <- 10L
snap <- readRDS(file.path(DATA_PROC, "ghcn_seasonal_intervals.rds"))
x <- interval_data(snap$centers, snap$radii)
xs <- standardize(x)
set.seed(SEED_BASE + 6000L)
rm <- run_methods_timed(xs)
tie_seed <- SEED_BASE + 6000L

mi <- read.csv(file.path(OUTPUT_DIR, "midscale_indices.csv"))
fl <- mi[mi$tie_affected == TRUE, ]
cat("USHCN flagged cells:", nrow(fl), "\n")
geom_h <- list(C = xs$centers, R = xs$radii)

blocks_all <- list(); vals_all <- list(); ii <- 0L
for (q in seq_len(nrow(fl))) {
  m <- fl$Method[q]; met <- fl$Metric[q]
  pr <- rm$projections[[m]]
  Dh <- dh_for(xs, met); Dl <- dl_for(pr, met)
  a <- audit_cell("USHCN", NA_integer_, m, met, Dh, Dl, geom_h,
                  list(C = pr$C, R = pr$R), K, tie_seed = tie_seed)
  stored <- fl[q, c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")]
  a$vals$replication_fidelity <-
    max(abs(unlist(stored) - unlist(a$vals[paste0("pre_", names(stored))])))
  a$vals$display_changed <- any(
    sprintf("%.3f", unlist(a$vals[paste0("pre_", names(stored))])) !=
    sprintf("%.3f", unlist(a$vals[paste0("post_", names(stored))])))
  ii <- ii + 1L
  blocks_all[[ii]] <- a$blocks; vals_all[[ii]] <- a$vals
  cat(m, met, ": blocks", if (is.null(a$blocks)) 0 else nrow(a$blocks),
      "fidelity", a$vals$replication_fidelity,
      "max delta", a$vals$max_abs_mean_delta, "\n")
}
blocks <- do.call(rbind, blocks_all)
vals <- do.call(rbind, vals_all)
write.csv(blocks, file.path(OUTPUT_DIR, "e10b_ushcn_blocks.csv"), row.names = FALSE)
write.csv(vals, file.path(OUTPUT_DIR, "e10b_ushcn_values.csv"), row.names = FALSE)
cat("\n=== USHCN block classification ===\n")
print(table(blocks$classification))
cat("cells with display-level change:", sum(vals$display_changed), "\n")
cat("E10B_USHCN_VERDICT:",
    if (all(blocks$classification == "exact_duplicate") &&
        sum(vals$display_changed) == 0) "PASS" else "REVIEW", "\n")
