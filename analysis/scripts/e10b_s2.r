# E10b tie-impact audit, Scenario II flagged cells (Prompt E 2026-07-17).
# All flagged methods (IMDS, C-PCA, SPCA) are deterministic and
# RNG-independent, so projections are regenerated solo per replication (DR on
# RAW data; evaluation on STANDARDIZED geometry, as in production).
# Fidelity is gated against stored per-replication values.
# Env REP_CHUNK="a:b" restricts to flagged reps in [a, b] and suffixes
# outputs "_a_b" (for parallel chunking); the merge step is e10b_s2_merge.r.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.r"))
source(file.path(SCRIPT_DIR, "e10b_lib.r"))

K <- 10L
pr_all <- read.csv(file.path(OUTPUT_DIR, "scenario2_per_rep.csv"))
pr_all <- pr_all[pr_all$dr_input == "raw" & pr_all$eval_geom == "std", ]
fl <- pr_all[pr_all$tie_affected == TRUE, ]
chunk <- Sys.getenv("REP_CHUNK", "")
sfx <- ""
if (nzchar(chunk)) {
  ab <- as.integer(strsplit(chunk, ":")[[1]])
  fl <- fl[fl$rep >= ab[1] & fl$rep <= ab[2], ]
  sfx <- paste0("_", ab[1], "_", ab[2])
}
cat("Scenario II flagged rows in scope:", nrow(fl), "reps:",
    length(unique(fl$rep)), "\n")

blocks_all <- list(); vals_all <- list(); ii <- 0L
t0 <- Sys.time()
for (r in sort(unique(fl$rep))) {
  seed_r <- SEED_BASE + 2000L + r
  x_raw <- gen_scenario2(seed_r)
  xs <- standardize(x_raw)              # evaluation geometry
  flr <- fl[fl$rep == r, ]
  geom_h <- list(C = xs$centers, R = xs$radii)
  projs <- list()
  for (m in unique(flr$Method))
    projs[[m]] <- run_idr(x_raw, methods = m, verbose = FALSE)[[m]]
  for (q in seq_len(nrow(flr))) {
    m <- flr$Method[q]; met <- flr$Metric[q]
    pr <- projs[[m]]
    Dh <- dh_for(xs, met); Dl <- dl_for(pr, met)
    a <- audit_cell("ScenarioII", r, m, met, Dh, Dl, geom_h,
                    list(C = pr$C, R = pr$R), K, tie_seed = seed_r)
    stored <- flr[q, c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")]
    a$vals$replication_fidelity <-
      max(abs(unlist(stored) -
              unlist(a$vals[paste0("pre_", names(stored))])))
    ii <- ii + 1L
    blocks_all[[ii]] <- a$blocks; vals_all[[ii]] <- a$vals
  }
  cat(sprintf("rep %d done, %.1f min elapsed\n", r,
    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
blocks <- do.call(rbind, blocks_all)
vals <- do.call(rbind, vals_all)
write.csv(blocks, file.path(OUTPUT_DIR, paste0("e10b_s2_blocks", sfx, ".csv")),
          row.names = FALSE)
write.csv(vals, file.path(OUTPUT_DIR, paste0("e10b_s2_values", sfx, ".csv")),
          row.names = FALSE)
cat("\n=== chunk classification summary ===\n")
print(table(blocks$classification))
cat("max replication-fidelity error:", max(vals$replication_fidelity), "\n")
cat("max |pre-post| mean delta:", max(vals$max_abs_mean_delta), "\n")
cat("E10B_S2_CHUNK_DONE", sfx, "\n")
