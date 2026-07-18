# Apply author decision D4(b) (Prompt E, 2026-07-17): adopt the
# corrected-policy Scenario II calibration. Re-runs the production
# calibration (reps 1-25, m = 999, production seeds and call order) under
# the corrected tolerance-block tie randomization, saving BOTH raw and
# min-P adjusted p-values in the production format, and REPLACES
# scenario2_calibration (before-copies preserved by the caller).
# Asserts: the min-P part is identical to the audit run
# (e10b_s2_cal_new.csv, same seeds) and exactly one star flips vs the
# published calibration (IMDS x Wasserstein x B_RE, 22 -> 23 of 25).
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.R"))
save_sessioninfo("e_apply_d4_s2")

K <- 10L
cal_rows <- list()
t0 <- Sys.time()
for (r in 1:25) {
  seed_r <- SEED_BASE + 2000L + r
  x_raw <- gen_scenario2(seed_r)
  xs <- standardize(x_raw)
  set.seed(seed_r + 500000L)
  rm <- run_methods_timed(x_raw)        # production: DR on RAW data
  a <- assess_quality(xs, rm$projections, K = K, perm_test = TRUE,
                      n_perm = 999L, seed = seed_r + 900000L,
                      ties = "random", baseline = TRUE)
  pu <- a$pvalues; pa <- a$pvalues_adj
  pu$rep <- r; pa$rep <- r; pu$kind <- "raw"; pa$kind <- "minP"
  cal_rows[[r]] <- rbind(pu, pa)
  cat(sprintf("rep %d done, %.1f min elapsed\n", r,
      as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
new_cal <- do.call(rbind, cal_rows)

idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
## identity vs the audit run (same seeds -> must agree exactly)
aud <- read.csv(file.path(OUTPUT_DIR, "e10b_s2_cal_new.csv"))
mp <- new_cal[new_cal$kind == "minP", ]
mg <- merge(mp, aud, by = c("IDR", "Metric", "rep"), suffixes = c("", ".aud"))
stopifnot(nrow(mg) == nrow(aud),
          max(abs(as.matrix(mg[idx]) - as.matrix(mg[paste0(idx, ".aud")]))) < 1e-12)
cat("min-P identity vs audit run: OK\n")

## star delta vs published calibration: exactly one flip
old <- read.csv(file.path(OUTPUT_DIR, "scenario2_calibration.csv"))
oldm <- old[old$kind == "minP", ]
flips <- character(0)
for (ix in idx) {
  so <- aggregate((oldm[[ix]] <= 0.05), oldm[c("IDR", "Metric")], sum)
  sn <- aggregate((mp[[ix]] <= 0.05), mp[c("IDR", "Metric")], sum)
  cc <- merge(so, sn, by = c("IDR", "Metric"), suffixes = c(".o", ".n"))
  fl <- cc[(cc$x.o >= 23) != (cc$x.n >= 23), ]
  if (nrow(fl)) flips <- c(flips, paste(ix, fl$IDR, fl$Metric))
}
cat("star flips vs published:", length(flips), "->", paste(flips, collapse = "; "), "\n")
stopifnot(identical(flips, "B_RE IMDS Wasserstein"))

save_result(new_cal, "scenario2_calibration",
            extra = list(m = 999L, n_cal = 25L,
  provenance = paste("Prompt E D4(b) 2026-07-17: corrected tolerance-block",
                     "tie randomization (rank_matrix patch); production",
                     "call and seeds otherwise identical; supersedes the",
                     "2026-07-16 artifact computed under the defective",
                     "near-tie randomization")))
cat("D4_APPLIED\n")
