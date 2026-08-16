# Attribution probe: regenerate scenario4 calibration rep 1 raw p-values and
# compare to the stored production file, printing per-method mismatch counts
# and the package version/path used. Writes nothing to output/.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.r"))

cat("QAIDR", as.character(utils::packageVersion("QAIDR")), "at",
    find.package("QAIDR"), "\n")
K0 <- 10L; M_PERM <- 999L
r <- 1L
seed_r <- SEED_BASE + 7000L + r
x <- gen_scenario4(seed_r)
xs <- standardize(x)
set.seed(seed_r + 500000L)
rm <- run_methods_timed(xs)
set.seed(seed_r + 300001L)
a <- assess_quality(xs, rm$projections, K = K0, perm_test = TRUE,
                    n_perm = M_PERM, seed = seed_r + 900000L,
                    ties = "random", baseline = TRUE)
pu <- a$pvalues
old <- readRDS(file.path(OUTPUT_DIR, "h_old_inference_snapshot",
                         "scenario4_calibration.rds"))
o <- old[old$kind == "raw" & old$rep == r, ]
o <- o[order(o$IDR, o$Metric), ]
n <- pu[order(pu$IDR, pu$Metric), ]
pc <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
d <- abs(as.matrix(o[pc]) - as.matrix(n[pc]))
cat("rep1 raw cells differing:", sum(d > 1e-12), "of", length(d),
    "max|d|=", max(d), "\n")
agg <- aggregate(rowSums(d > 1e-12), list(Method = o$IDR), sum)
print(agg)
saveRDS(n, file.path(Sys.getenv("H_PROBE_OUT", tempdir()),
                     paste0("s4_rep1_raw_", gsub("[^0-9.]", "",
                            as.character(utils::packageVersion("QAIDR"))), ".rds")))
