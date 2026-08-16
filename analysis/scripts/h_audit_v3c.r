# Prompt H (contract 1.4): rerun of the frozen stage4_audit_v3b familywise
# calibration audit under the symmetric randomization min-P (QAIDR 0.3.0).
# Identical frozen design (same fixed 60-object Scenario-I subsample, K = 5,
# production m = 999, ties = "error", complete random-correspondence null via
# fresh point projections, 2,000 replications, Wilson 95% CI rule); only the
# replication seed streams are fresh and disjoint per the preregistration:
# base 20260719 with the v3b offsets (data +9,700,000 + r; perms
# +9,800,000 + r).
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.r"))
save_sessioninfo("h_audit_v3c")

stopifnot(utils::packageVersion("QAIDR") >= "0.3.0")
SEED_H <- 20260719L
N_REP <- if (RUN_MODE == "full") 2000L else 50L
N_OBJ <- 60L; K <- 5L; M_CAL <- 999L; ALPHA <- 0.05
set.seed(SEED_BASE + 20L)                       # frozen dataset, unchanged
x1 <- standardize(gen_scenario1(SEED_BASE + 1L))
sub <- sample(nrow(x1$centers), N_OBJ)
xs <- interval_data(x1$centers[sub, ], x1$radii[sub, ])
fams <- list(TC = c("Q_TC", "B_TC"), RE = c("Q_RE", "B_RE"), LC = c("Q_LC", "B_LC"))
fam_rej <- matrix(0, N_REP, 3, dimnames = list(NULL, names(fams)))
t0 <- Sys.time()
for (r in seq_len(N_REP)) {
  set.seed(SEED_H + 9700000L + r)               # fresh disjoint data stream
  Zr <- matrix(rnorm(N_OBJ * 2), N_OBJ, 2)
  proj_r <- structure(list(nullproj = list(C = Zr, R = matrix(0, N_OBJ, 2),
                                           type = "Point")),
                      class = "idr_projections")
  a <- assess_quality(xs, proj_r, K = K, perm_test = TRUE, n_perm = M_CAL,
                      seed = SEED_H + 9800000L + r, ties = "error",
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
res <- data.frame(family = names(fams),
                  rejections = colSums(fam_rej),
                  rate = colMeans(fam_rej))
cis <- t(vapply(res$rejections, wilson, numeric(2), n = N_REP))
res$ci_lo <- cis[, 1]; res$ci_hi <- cis[, 2]
res$pass <- (res$ci_lo <= ALPHA & ALPHA <= res$ci_hi) | res$ci_hi < ALPHA
save_result(res, "h_audit_v3c_familywise",
            extra = list(n_rep = N_REP, m = M_CAL, alpha = ALPHA,
                         engine = "symmetric randomization min-P (QAIDR 0.3.0)",
                         seed_stream = "20260719+9.7e6+r (data), +9.8e6+r (perms)",
                         prereg = "H_PREREGISTRATION.md sec 1.4 (post-A1)",
                         rule = "Wilson 95% CI contains alpha or lies below"))
print(res, digits = 4)
cat("h_audit_v3c done:", format(Sys.time()), "\n")
