# Merge the four e10b_s2 chunk outputs and compute the Scenario II
# table-level impact (Prompt E 2026-07-17).
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))

sfx <- c("_1_25", "_26_50", "_51_75", "_76_100")
blocks <- do.call(rbind, lapply(sfx, function(s)
  read.csv(file.path(OUTPUT_DIR, paste0("e10b_s2_blocks", s, ".csv")))))
vals <- do.call(rbind, lapply(sfx, function(s)
  read.csv(file.path(OUTPUT_DIR, paste0("e10b_s2_values", s, ".csv")))))
write.csv(blocks, file.path(OUTPUT_DIR, "e10b_s2_blocks.csv"), row.names = FALSE)
write.csv(vals, file.path(OUTPUT_DIR, "e10b_s2_values.csv"), row.names = FALSE)
cat("=== Scenario II classification summary ===\n")
print(table(blocks$classification))
cat("blocks by method:\n"); print(table(blocks$Method, blocks$classification))
cat("max replication-fidelity error:", max(vals$replication_fidelity), "\n")
cat("max |pre-post| mean delta:", max(vals$max_abs_mean_delta), "\n")

pr_all <- read.csv(file.path(OUTPUT_DIR, "scenario2_per_rep.csv"))
pr_all <- pr_all[pr_all$dr_input == "raw" & pr_all$eval_geom == "std", ]
fl <- pr_all[pr_all$tie_affected == TRUE, ]
idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
imp <- list(); pp <- 0L
for (cell in unique(paste(fl$Method, fl$Metric, sep = "||"))) {
  mm <- strsplit(cell, "||", fixed = TRUE)[[1]]
  rows <- pr_all$Method == mm[1] & pr_all$Metric == mm[2]
  newpr <- pr_all[rows, ]
  for (q in seq_len(nrow(vals))) {
    if (vals$Method[q] == mm[1] && vals$Metric[q] == mm[2]) {
      tr <- newpr$rep == vals$replication[q]
      newpr[tr, idx] <- vals[q, paste0("post_", idx)]
    }
  }
  for (ix in idx) {
    old_m <- mean(pr_all[rows, ix]); new_m <- mean(newpr[, ix])
    old_s <- sd(pr_all[rows, ix]);   new_s <- sd(newpr[, ix])
    pp <- pp + 1L
    imp[[pp]] <- data.frame(Method = mm[1], Metric = mm[2], index = ix,
      old_mean_fmt = sprintf("%.3f", old_m), new_mean_fmt = sprintf("%.3f", new_m),
      old_sd_fmt = sprintf("%.0f", old_s * 1000), new_sd_fmt = sprintf("%.0f", new_s * 1000),
      mean_delta = new_m - old_m)
  }
}
impact <- do.call(rbind, imp)
impact$display_changed <- impact$old_mean_fmt != impact$new_mean_fmt |
                          impact$old_sd_fmt != impact$new_sd_fmt
write.csv(impact, file.path(OUTPUT_DIR, "e10b_s2_impact.csv"), row.names = FALSE)
cat("cells with display-level change:", sum(impact$display_changed), "\n")
cat("max |mean delta|:", max(abs(impact$mean_delta)), "\n")
cat("E10B_S2_VERDICT:",
    if (all(blocks$classification %in%
            c("exact_duplicate", "tolerance_near_tie")) &&
        sum(impact$display_changed) == 0 &&
        max(vals$replication_fidelity) < 1e-9)
      "VALUES-PASS (classification per report)" else "REVIEW", "\n")
