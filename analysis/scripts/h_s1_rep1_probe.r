# Attribution probe: scenario 1 replication 1 under the current committed
# pipeline vs (a) stored scenario1_per_rep indices (B2-era production, basis
# of the graded-width L3 reuse assertion A1.5) and (b) stored
# scenario1_calibration raw p-values (Prompt E regeneration). Writes nothing.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.r"))

cat("QAIDR", as.character(utils::packageVersion("QAIDR")), "at",
    find.package("QAIDR"), "\n")
K <- 10L; M_PERM <- 999L
r <- 1L
seed_r <- SEED_BASE + 1000L + r
x <- gen_scenario1(seed_r)
xs <- standardize(x)
set.seed(seed_r + 500000L)
rm <- run_methods_timed(xs)
cells <- evaluate_cells(xs, rm$projections, K = K, tie_seed = seed_r)
a <- assess_quality(xs, rm$projections, K = K, perm_test = TRUE,
                    n_perm = M_PERM, seed = seed_r + 900000L,
                    ties = "random", baseline = TRUE)
pc <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")

per <- readRDS(file.path(OUTPUT_DIR, "scenario1_per_rep.rds"))
o <- per[per$rep == r, ]
o <- o[order(o$Method, o$Metric), ]
n <- cells[order(cells$Method, cells$Metric), ]
d <- abs(as.matrix(o[pc]) - as.matrix(n[pc]))
cat("rep1 per_rep index cells differing:", sum(d > 1e-12), "of", length(d),
    "max|d|=", max(d), "\n")
print(aggregate(rowSums(d > 1e-12), list(Method = o$Method), sum))

old <- readRDS(file.path(OUTPUT_DIR, "h_old_inference_snapshot",
                         "scenario1_calibration.rds"))
oc <- old[old$kind == "raw" & old$rep == r, ]
oc <- oc[order(oc$IDR, oc$Metric), ]
nc <- a$pvalues[order(a$pvalues$IDR, a$pvalues$Metric), ]
d2 <- abs(as.matrix(oc[pc]) - as.matrix(nc[pc]))
cat("rep1 calibration raw cells differing:", sum(d2 > 1e-12), "of",
    length(d2), "max|d|=", max(d2), "\n")
print(aggregate(rowSums(d2 > 1e-12), list(Method = oc$IDR), sum))
