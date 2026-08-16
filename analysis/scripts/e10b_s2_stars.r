# E10b follow-up (Prompt E 2026-07-17): Scenario II calibration inference
# re-run under the corrected tie randomization. The block audit found
# thousands of low-space tolerance near-ties per replication, including in
# the calibration subset (reps 1-25), so the per-replication adjusted
# p-values and the >= 23/25 star rule must be compared old vs new.
# Production protocol replicated (scenario2.r): reps ascending from a cold
# session, per rep: xs = standardize(gen_scenario1(seed_r));
# set.seed(seed_r + 500000); run_methods_timed(xs);
# assess_quality(..., n_perm = 999, seed = seed_r + 900000, ties = "random").
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.r"))

K0 <- 10L
rows <- list()
t0 <- Sys.time()
for (r in 1:25) {
  seed_r <- SEED_BASE + 2000L + r
  x_raw <- gen_scenario2(seed_r)
  xs <- standardize(x_raw)
  set.seed(seed_r + 500000L)
  rm <- run_methods_timed(x_raw)   # production: DR on RAW data
  a <- assess_quality(xs, rm$projections, K = K0, perm_test = TRUE,
                      n_perm = 999L, seed = seed_r + 900000L,
                      ties = "random", baseline = TRUE)
  pa <- a$pvalues_adj; pa$rep <- r; pa$kind <- "minP"
  rows[[r]] <- pa
  cat(sprintf("rep %d done, %.1f min elapsed\n", r,
      as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
new_cal <- do.call(rbind, rows)
write.csv(new_cal, file.path(OUTPUT_DIR, "e10b_s2_cal_new.csv"), row.names = FALSE)

old_cal <- read.csv(file.path(OUTPUT_DIR, "scenario2_calibration.csv"))
old_cal <- old_cal[old_cal$kind == "minP", ]
idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
mg <- merge(old_cal, new_cal, by = c("IDR", "Metric", "rep"),
            suffixes = c(".old", ".new"))
cat("aligned calibration rows:", nrow(mg), "of", nrow(old_cal), "\n")
## star rule: adjusted p <= 0.05 in >= 23 of the 25 calibration reps
tot_flips <- 0L
for (ix in idx) {
  so <- aggregate((mg[[paste0(ix, ".old")]] <= 0.05),
                  mg[c("IDR", "Metric")], sum)
  sn <- aggregate((mg[[paste0(ix, ".new")]] <= 0.05),
                  mg[c("IDR", "Metric")], sum)
  cmp <- merge(so, sn, by = c("IDR", "Metric"), suffixes = c(".old", ".new"))
  cmp$star_old <- cmp$x.old >= 23; cmp$star_new <- cmp$x.new >= 23
  fl <- which(cmp$star_old != cmp$star_new)
  cat(sprintf("%-5s max |count delta| = %d  star flips = %d\n", ix,
              max(abs(cmp$x.old - cmp$x.new)), length(fl)))
  if (length(fl)) print(cmp[fl, ])
  tot_flips <- tot_flips + length(fl)
}
cat("TOTAL S2 STAR FLIPS:", tot_flips, "\n")
cat("E10B_S2_STARS_DONE\n")
