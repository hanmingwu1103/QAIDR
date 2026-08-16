# Prompt C Phase 2.3 multiplicity audit (frozen plan, SHA-256 3d15e9fd...):
# 2,000 independent replications under the COMPLETE random-correspondence
# null with a fresh, disjoint seed stream; Wilson 95% CI acceptance rule.
# Writes stage4v3_audit_* ONLY; does not touch stage4_gates.json or the v2
# calibration CSV.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.r"))
save_sessioninfo("stage4_audit_v3")

N_REP <- if (RUN_MODE == "full") 2000L else 50L
N_OBJ <- 60L; K <- 5L; M_CAL <- 199L; ALPHA <- 0.05
set.seed(SEED_BASE + 20L)
x1 <- standardize(gen_scenario1(SEED_BASE + 1L))
sub <- sample(nrow(x1$centers), N_OBJ)
xs <- interval_data(x1$centers[sub, ], x1$radii[sub, ])

fams <- list(TC = c("Q_TC", "B_TC"), RE = c("Q_RE", "B_RE"), LC = c("Q_LC", "B_LC"))
idx <- unlist(fams, use.names = FALSE)
fam_rej <- matrix(0, N_REP, 3, dimnames = list(NULL, names(fams)))
per_rej <- NULL
t0 <- Sys.time()
for (r in seq_len(N_REP)) {
  set.seed(SEED_BASE + 9000000L + r)          # fresh, disjoint stream
  Zr <- matrix(rnorm(N_OBJ * 2), N_OBJ, 2)
  proj_r <- structure(list(nullproj = list(C = Zr, R = matrix(0, N_OBJ, 2),
                                           type = "Point")),
                      class = "idr_projections")
  a <- assess_quality(xs, proj_r, K = K, perm_test = TRUE, n_perm = M_CAL,
                      seed = SEED_BASE + 9500000L + r, ties = "error",
                      baseline = TRUE)
  pu <- a$pvalues; pa <- a$pvalues_adj
  rj <- as.matrix(pu[idx]) <= ALPHA
  if (is.null(per_rej)) per_rej <- rj * 0
  per_rej <- per_rej + rj
  for (f in names(fams)) fam_rej[r, f] <- any(as.matrix(pa[fams[[f]]]) <= ALPHA)
  if (r %% 250 == 0) cat("rep", r, "/", N_REP, "elapsed",
                         round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
                         "min\n")
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
per_rates <- per_rej / N_REP
save_result(res, "stage4v3_audit_familywise",
            extra = list(n_rep = N_REP, m = M_CAL, alpha = ALPHA,
                         seed_stream = "SEED_BASE+9e6+r (data), +9.5e6+r (perms)",
                         plan_sha256 = "3d15e9fd4da697765418c87c193fa09ba04cfa3d42f413d8931af680a50a9730",
                         rule = "Wilson 95% CI contains alpha or lies below"))
save_result(as.data.frame(per_rates), "stage4v3_audit_pertest",
            extra = list(n_rep = N_REP))
print(res, digits = 4)
cat("AUDIT", if (all(res$pass)) "PASS" else "FAIL", "\n")
if (!all(res$pass)) quit(status = 1)
