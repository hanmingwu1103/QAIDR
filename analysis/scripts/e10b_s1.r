# E10b tie-impact audit, Scenario I flagged cells (Prompt E 2026-07-17).
# IMDS is deterministic and RNG-independent (verified: identical output
# across RNG states; consumes no RNG), so IMDS projections are regenerated
# solo. The single Int-UMAP cell (rep 17, Hausdorff) is regenerated with the
# production stream: C-PCA warm (rep > 1), set.seed(seed_r + 500000), full
# METHODS sequence. Replication fidelity is gated against the stored
# per-replication values in every audited cell.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.r"))
source(file.path(SCRIPT_DIR, "e10b_lib.r"))

K <- 10L
pr_all <- read.csv(file.path(OUTPUT_DIR, "scenario1_per_rep.csv"))
fl <- pr_all[pr_all$tie_affected == TRUE, ]
cat("Scenario I flagged rows:", nrow(fl), "over reps:",
    length(unique(fl$rep)), "\n")

## warm up the one-time C-PCA RNG initialization (verified: first call only)
xw <- standardize(gen_scenario1(SEED_BASE + 999999L))
invisible(run_idr(xw, methods = "C-PCA", verbose = FALSE))

blocks_all <- list(); vals_all <- list(); ii <- 0L
t0 <- Sys.time()
for (r in sort(unique(fl$rep))) {
  seed_r <- SEED_BASE + 1000L + r
  xs <- standardize(gen_scenario1(seed_r))
  flr <- fl[fl$rep == r, ]
  geom_h <- list(C = xs$centers, R = xs$radii)
  projs <- list()
  if (any(flr$Method == "IMDS"))
    projs[["IMDS"]] <- run_idr(xs, methods = "IMDS", verbose = FALSE)[["IMDS"]]
  if (any(flr$Method == "Int-UMAP")) {
    set.seed(seed_r + 500000L)
    projs[["Int-UMAP"]] <- run_methods_timed(xs)$projections[["Int-UMAP"]]
  }
  for (q in seq_len(nrow(flr))) {
    m <- flr$Method[q]; met <- flr$Metric[q]
    pr <- projs[[m]]
    Dh <- dh_for(xs, met); Dl <- dl_for(pr, met)
    a <- audit_cell("ScenarioI", r, m, met, Dh, Dl, geom_h,
                    list(C = pr$C, R = pr$R), K, tie_seed = seed_r)
    stored <- flr[q, c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")]
    a$vals$replication_fidelity <-
      max(abs(unlist(stored) -
              unlist(a$vals[paste0("pre_", names(stored))])))
    ii <- ii + 1L
    blocks_all[[ii]] <- a$blocks; vals_all[[ii]] <- a$vals
  }
  if (r %% 10 == 0) cat(sprintf("rep %d done, %.1f min elapsed\n", r,
    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
blocks <- do.call(rbind, blocks_all)
vals <- do.call(rbind, vals_all)
write.csv(blocks, file.path(OUTPUT_DIR, "e10b_s1_blocks.csv"), row.names = FALSE)
write.csv(vals, file.path(OUTPUT_DIR, "e10b_s1_values.csv"), row.names = FALSE)

cat("\n=== Scenario I classification summary ===\n")
print(table(blocks$classification))
cat("max replication-fidelity error:", max(vals$replication_fidelity), "\n")
cat("max |pre-post| 50-resolution mean delta:", max(vals$max_abs_mean_delta), "\n")

## Table-level impact: substitute post-patch means and re-aggregate.
sm_idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
imp <- list(); pp <- 0L
for (cell in unique(paste(fl$Method, fl$Metric, sep = "||"))) {
  mm <- strsplit(cell, "||", fixed = TRUE)[[1]]
  rows <- pr_all$Method == mm[1] & pr_all$Metric == mm[2]
  newpr <- pr_all[rows, ]
  for (q in seq_len(nrow(vals))) {
    if (vals$Method[q] == mm[1] && vals$Metric[q] == mm[2]) {
      tr <- newpr$rep == vals$replication[q]
      newpr[tr, sm_idx] <- vals[q, paste0("post_", sm_idx)]
    }
  }
  for (ix in sm_idx) {
    old_m <- mean(pr_all[rows, ix]); new_m <- mean(newpr[, ix])
    old_s <- sd(pr_all[rows, ix]);   new_s <- sd(newpr[, ix])
    pp <- pp + 1L
    imp[[pp]] <- data.frame(Method = mm[1], Metric = mm[2], index = ix,
      old_mean_fmt = sprintf("%.3f", old_m), new_mean_fmt = sprintf("%.3f", new_m),
      old_sd_fmt = sprintf("%.0f", old_s * 1000), new_sd_fmt = sprintf("%.0f", new_s * 1000),
      mean_delta = new_m - old_m, sd_delta = new_s - old_s)
  }
}
impact <- do.call(rbind, imp)
impact$display_changed <- impact$old_mean_fmt != impact$new_mean_fmt |
                          impact$old_sd_fmt != impact$new_sd_fmt
write.csv(impact, file.path(OUTPUT_DIR, "e10b_s1_impact.csv"), row.names = FALSE)
cat("cells with display-level change:", sum(impact$display_changed), "\n")
cat("max |mean delta|:", max(abs(impact$mean_delta)), "\n")
cat("max replication-fidelity:", max(vals$replication_fidelity), "\n")
cat("E10B_S1_VERDICT:",
    if (all(blocks$classification == "exact_duplicate") &&
        sum(impact$display_changed) == 0 &&
        max(vals$replication_fidelity) < 1e-9) "PASS" else "REVIEW", "\n")
