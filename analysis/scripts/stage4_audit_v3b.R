# Prompt C supplementary multiplicity audit at the PRODUCTION m = 999.
# Motivated by the frozen v3 audit's finding that m = 199 min-P is
# anti-conservative (familywise ~0.067-0.070); this run determines whether the
# production configuration (m = 999, as used for every manuscript star) is
# calibrated. Fresh disjoint seed stream SEED_BASE + 9,700,000 + r (data),
# + 9,800,000 + r (permutations); 2,000 replications; same Wilson CI rule.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.R"))
save_sessioninfo("stage4_audit_v3b")

N_REP <- if (RUN_MODE == "full") 2000L else 50L
N_OBJ <- 60L; K <- 5L; M_CAL <- 999L; ALPHA <- 0.05
set.seed(SEED_BASE + 20L)
x1 <- standardize(gen_scenario1(SEED_BASE + 1L))
sub <- sample(nrow(x1$centers), N_OBJ)
xs <- interval_data(x1$centers[sub, ], x1$radii[sub, ])
fams <- list(TC = c("Q_TC", "B_TC"), RE = c("Q_RE", "B_RE"), LC = c("Q_LC", "B_LC"))
fam_rej <- matrix(0, N_REP, 3, dimnames = list(NULL, names(fams)))
t0 <- Sys.time()
for (r in seq_len(N_REP)) {
  set.seed(SEED_BASE + 9700000L + r)
  Zr <- matrix(rnorm(N_OBJ * 2), N_OBJ, 2)
  proj_r <- structure(list(nullproj = list(C = Zr, R = matrix(0, N_OBJ, 2),
                                           type = "Point")),
                      class = "idr_projections")
  a <- assess_quality(xs, proj_r, K = K, perm_test = TRUE, n_perm = M_CAL,
                      seed = SEED_BASE + 9800000L + r, ties = "error",
                      baseline = TRUE)
  pa <- a$pvalues_adj
  for (f in names(fams)) fam_rej[r, f] <- any(as.matrix(pa[fams[[f]]]) <= ALPHA)
  if (r %% 250 == 0) cat("rep", r, "/", N_REP, "elapsed",
                         round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
}
wilson <- function(x, n, z = 1.96) {
  p <- x / n; c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n),
                (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n))
}
res <- data.frame(family = names(fams), rejections = colSums(fam_rej),
                  rate = colMeans(fam_rej))
cis <- t(vapply(res$rejections, wilson, numeric(2), n = N_REP))
res$ci_lo <- cis[, 1]; res$ci_hi <- cis[, 2]
res$pass <- (res$ci_lo <= ALPHA & ALPHA <= res$ci_hi) | res$ci_hi < ALPHA
save_result(res, "stage4v3b_audit_m999",
            extra = list(n_rep = N_REP, m = M_CAL,
                         seed_stream = "SEED_BASE+9.7e6/9.8e6+r",
                         motivation = "m=199 audit anti-conservative; production m=999 audited"))
print(res, digits = 4)
cat("AUDIT m=999", if (all(res$pass)) "PASS" else "FAIL", "\n")
if (!all(res$pass)) quit(status = 1)
