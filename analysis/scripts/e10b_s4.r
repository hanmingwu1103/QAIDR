# E10b tie-impact audit, Scenario IV (Prompt E 2026-07-17).
# Part 1 (all 100 reps, records-based): verify that every recorded high-space
# tie link is exactly accounted for by the built-in duplicate pairs
# (hi_ties == (n-2) * n_duplicate_pairs for every metric and replication),
# so no tolerance near-tie existed in any Scenario IV high space; summarize
# lo_ties. Part 2 (subsample recomputation): regenerate reps with production
# seeds, classify every block in both spaces, and compare pre/post-patch
# 50-resolution means under the production eval_rep stream structure (K=10).
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.r"))
source(file.path(SCRIPT_DIR, "e10b_lib.r"))

K0 <- 10L; N_RES <- 50L
rec <- readRDS(file.path(OUTPUT_DIR, "scenario4_records.rds"))

## Part 1: records identity over all 100 reps
ok_hi <- vapply(seq_along(rec), function(r) {
  e <- rec[[r]]
  n_dup <- length(e$donors)
  all(e$hi_ties == (300L - 2L) * n_dup)
}, logical(1))
lo_pos <- vapply(seq_along(rec), function(r) sum(rec[[r]]$lo_ties > 0), integer(1))
lo_max <- vapply(seq_along(rec), function(r) max(rec[[r]]$lo_ties), numeric(1))
cat("Part 1: hi_ties == 298 * n_dup(intended) for all metrics in", sum(ok_hi),
    "of", length(rec), "reps\n")
cat("lo_ties: cells > 0 per rep, range:", min(lo_pos), "-", max(lo_pos),
    "; max lo_ties value:", max(lo_max), "\n")

## Reps that miss the intended-pair identity: contamination can redraw one
## member of a duplicate pair (fewer links), or a far-away anchor can
## compress relative gaps into genuine tolerance near-ties (extra links).
## Non-conforming reps get full block diagnosis + value-impact treatment in
## Part 2 below.
part1_ok <- ok_hi
extra_reps <- integer(0)
for (r in which(!ok_hi)) {
  e <- rec[[r]]
  x <- gen_scenario4(SEED_BASE + 7000L + r)
  eff <- sum(vapply(seq_along(e$donors), function(k) {
    d <- e$donors[k]; rcp <- e$recipients[k]
    identical(x$centers[d, ], x$centers[rcp, ]) &&
      identical(x$radii[d, ], x$radii[rcp, ])
  }, logical(1)))
  ok_eff <- all(e$hi_ties == (300L - 2L) * eff)
  cat(sprintf("rep %d: intended pairs %d, effective %d, hi_ties %s -> %s\n",
              r, length(e$donors), eff, paste(unique(e$hi_ties), collapse = "/"),
              if (ok_eff) "MATCHES effective identity"
              else "EXTRA LINKS - full diagnosis in Part 2"))
  part1_ok[r] <- ok_eff
  if (!ok_eff) extra_reps <- c(extra_reps, r)
}

## Part 2: subsample recomputation with production seeds; non-conforming
## reps are always included and get table-level impact substitution.
SUB <- sort(unique(c(1L, 26L, 51L, 76L, extra_reps)))
mets <- c(METRICS, "Centers-Euclidean")
blocks_all <- list(); vals_all <- list(); bb <- 0L
for (r in SUB) {
  seed_r <- SEED_BASE + 7000L + r
  x <- gen_scenario4(seed_r)
  xs <- standardize(x)
  set.seed(seed_r + 500000L)
  rm <- run_methods_timed(xs)
  tie_seed <- seed_r + 300000L
  geom_h <- list(C = xs$centers, R = xs$radii)

  Dh <- lapply(mets, function(met) dh_for(xs, met)); names(Dh) <- mets
  Dl <- list()
  for (m in names(rm$projections)) for (met in mets)
    Dl[[paste(m, met)]] <- dl_for(rm$projections[[m]], met)

  ## block classification, both spaces, all cells
  for (met in mets) {
    b <- extract_blocks(Dh[[met]], geom_h)
    if (!is.null(b)) { bb <- bb + 1L
      blocks_all[[bb]] <- cbind(scenario = "ScenarioIV", replication = r,
        Method = "(high space)", Metric = met, space = "high", b) }
    for (m in names(rm$projections)) {
      pr <- rm$projections[[m]]
      bl <- extract_blocks(Dl[[paste(m, met)]], list(C = pr$C, R = pr$R))
      if (!is.null(bl)) { bb <- bb + 1L
        blocks_all[[bb]] <- cbind(scenario = "ScenarioIV", replication = r,
          Method = m, Metric = met, space = "low", bl) }
    }
  }

  ## pre/post value comparison under the production stream structure
  for (routine in c("pre", "post")) {
    rkfun <- if (routine == "pre") rank_matrix_prepatch else rank_matrix
    acc <- list(); aa <- 0L
    for (met in mets) {
      for (s in seq_len(N_RES)) {
        set.seed(tie_seed + s)
        Rh <- rkfun(Dh[[met]], ties = "random")
        for (m in names(rm$projections)) {
          Rl <- rkfun(Dl[[paste(m, met)]], ties = "random")
          aa <- aa + 1L
          acc[[aa]] <- data.frame(Method = m, Metric = met, res = s,
            t(indices_from_coranking(coranking_matrix(Rh, Rl), K0)))
        }
      }
    }
    long <- do.call(rbind, acc)
    idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
    cells <- aggregate(long[idx], long[c("Method", "Metric")], mean)
    names(cells)[3:8] <- paste0(routine, "_", idx)
    if (routine == "pre") pre_cells <- cells else post_cells <- cells
  }
  mg <- merge(pre_cells, post_cells, by = c("Method", "Metric"))
  mg$replication <- r
  idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  mg$max_abs_mean_delta <- apply(
    abs(mg[paste0("pre_", idx)] - mg[paste0("post_", idx)]), 1, max)
  vals_all[[length(vals_all) + 1L]] <- mg
  cat("rep", r, "done; max |pre-post| mean delta:",
      max(mg$max_abs_mean_delta), "\n")
}
blocks <- do.call(rbind, blocks_all)
vals <- do.call(rbind, vals_all)
write.csv(blocks, file.path(OUTPUT_DIR, "e10b_s4_blocks.csv"), row.names = FALSE)
write.csv(vals, file.path(OUTPUT_DIR, "e10b_s4_values.csv"), row.names = FALSE)
cat("\n=== Scenario IV block classification (audited reps) ===\n")
print(table(blocks$classification, blocks$replication))

## Table-level impact for reps carrying non-exact blocks: substitute
## post-patch means into the per-rep records and re-aggregate the published
## mean(SD) cells (K = 10 primary rows).
idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
pr_all <- read.csv(file.path(OUTPUT_DIR, "scenario4_per_rep.csv"))
pr_all <- pr_all[pr_all$K == 10, ]
imp <- list(); pp <- 0L
for (cell in unique(paste(vals$Method, vals$Metric, sep = "||"))) {
  mm <- strsplit(cell, "||", fixed = TRUE)[[1]]
  rows <- pr_all$Method == mm[1] & pr_all$Metric == mm[2]
  newpr <- pr_all[rows, ]
  fid <- numeric(0)
  for (q in seq_len(nrow(vals))) {
    if (vals$Method[q] == mm[1] && vals$Metric[q] == mm[2] &&
        vals$replication[q] %in% extra_reps) {
      tr <- newpr$rep == vals$replication[q]
      fid <- c(fid, max(abs(unlist(newpr[tr, idx]) -
                            unlist(vals[q, paste0("pre_", idx)]))))
      newpr[tr, idx] <- vals[q, paste0("post_", idx)]
    }
  }
  if (!length(fid)) next
  for (ix in idx) {
    old_m <- mean(pr_all[rows, ix]); new_m <- mean(newpr[, ix])
    old_s <- sd(pr_all[rows, ix]);   new_s <- sd(newpr[, ix])
    pp <- pp + 1L
    imp[[pp]] <- data.frame(Method = mm[1], Metric = mm[2], index = ix,
      old_mean_fmt = sprintf("%.3f", old_m), new_mean_fmt = sprintf("%.3f", new_m),
      old_sd_fmt = sprintf("%.0f", old_s * 1000), new_sd_fmt = sprintf("%.0f", new_s * 1000),
      mean_delta = new_m - old_m, replication_fidelity = max(fid))
  }
}
if (length(imp)) {
  impact <- do.call(rbind, imp)
  impact$display_changed <- impact$old_mean_fmt != impact$new_mean_fmt |
                            impact$old_sd_fmt != impact$new_sd_fmt
  write.csv(impact, file.path(OUTPUT_DIR, "e10b_s4_impact.csv"), row.names = FALSE)
  cat("impact rows:", nrow(impact),
      "; display-level changes:", sum(impact$display_changed),
      "; max |mean delta|:", max(abs(impact$mean_delta)),
      "; max replication fidelity error:", max(impact$replication_fidelity), "\n")
} else impact <- data.frame(display_changed = logical(0))

nonexact <- blocks[blocks$classification != "exact_duplicate", ]
cat("non-exact blocks in audited reps:", nrow(nonexact),
    "(reps:", paste(unique(nonexact$replication), collapse = ","), ")\n")
cat("calibration subset (reps 1-25) carries non-exact blocks:",
    any(nonexact$replication <= 25), "-> stars",
    if (any(nonexact$replication <= 25)) "MUST BE re-audited" else "unaffected",
    "\n")
cat("E10B_S4_VERDICT:",
    if (all(part1_ok | seq_along(part1_ok) %in% extra_reps) &&
        sum(impact$display_changed) == 0 &&
        !any(nonexact$replication <= 25)) "PASS-WITH-DISCLOSED-NEAR-TIES"
    else "REVIEW", "\n")
