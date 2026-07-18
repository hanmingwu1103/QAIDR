# B2 Stage 4 gate v2: min-P familywise calibration under the CORRECTED
# joint-draw engine, exercising the full claimed family ({Q,B} x 4 interval
# dissimilarities + center baseline) per index family (TC/RE/LC), with a
# machine-readable evidence record (Codex schema). The 2026-07-15 record is
# preserved as deprecated audit history in stage4_gates_v1_deprecated.json.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.R"))
save_sessioninfo("stage4_gate2")

N_REP_CAL <- if (RUN_MODE == "full") 200L else 20L
N_OBJ <- 60L; K <- 5L; M_CAL <- 199L; ALPHA <- 0.05
set.seed(SEED_BASE + 20L)
x1 <- standardize(gen_scenario1(SEED_BASE + 1L))
sub <- sample(nrow(x1$centers), N_OBJ)
xs <- interval_data(x1$centers[sub, ], x1$radii[sub, ])

fams <- list(TC = c("Q_TC", "B_TC"), RE = c("Q_RE", "B_RE"), LC = c("Q_LC", "B_LC"))
per_test_rej <- NULL; fam_rej <- matrix(0, N_REP_CAL, 3,
                                        dimnames = list(NULL, names(fams)))
cal_rows <- list()
t0 <- Sys.time()
for (r in seq_len(N_REP_CAL)) {
  set.seed(SEED_BASE + 5000000L + r)
  Zr <- matrix(rnorm(N_OBJ * 2), N_OBJ, 2)          # independent null embedding
  proj_r <- structure(list(nullproj = list(C = Zr, R = matrix(0, N_OBJ, 2),
                                           type = "Point")),
                      class = "idr_projections")
  a <- assess_quality(xs, proj_r, K = K, perm_test = TRUE, n_perm = M_CAL,
                      seed = SEED_BASE + 6000000L + r, ties = "error",
                      baseline = TRUE)
  pu <- a$pvalues; pa <- a$pvalues_adj
  idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  rej_r <- as.matrix(pu[idx]) <= ALPHA
  rownames(rej_r) <- pu$Metric
  if (is.null(per_test_rej)) per_test_rej <- rej_r * 0
  per_test_rej <- per_test_rej + rej_r
  for (f in names(fams)) {
    fam_rej[r, f] <- any(as.matrix(pa[fams[[f]]]) <= ALPHA)
  }
  pu$rep <- r; pa$rep <- r; pu$kind <- "raw"; pa$kind <- "minP"
  cal_rows[[r]] <- rbind(pu, pa)
  if (r %% 50 == 0) cat("rep", r, "/", N_REP_CAL, "\n")
}
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cal <- do.call(rbind, cal_rows)
save_result(cal, "stage4v2_calibration_per_rep",
            extra = list(m = M_CAL, n_rep = N_REP_CAL, n = N_OBJ, K = K))
cal_hash <- digest::digest(file.path(OUTPUT_DIR, "stage4v2_calibration_per_rep.csv"),
                           file = TRUE, algo = "sha256")

se_bin <- sqrt(ALPHA * (1 - ALPHA) / N_REP_CAL)
fam_rates <- colMeans(fam_rej)
per_rates <- per_test_rej / N_REP_CAL
rule <- "every familywise rate within 4 binomial SEs of alpha; per-test rates <= alpha + 4 SE"
fam_pass <- all(abs(fam_rates - ALPHA) < 4 * se_bin)
per_pass <- all(per_rates <= ALPHA + 4 * se_bin)

## null-expectation check (re-run under current engine for the record)
Rh <- rank_matrix(idist(xs$centers, xs$radii, "Wasserstein"))
set.seed(SEED_BASE + 10L)
Zl <- matrix(rnorm(N_OBJ * 2), N_OBJ, 2)
Rl <- rank_matrix(as.matrix(dist(Zl)))
nulls <- perm_null_indices(Rh, Rl, K, m = 4000)
mu <- colMeans(nulls); se <- apply(nulls, 2, sd) / sqrt(nrow(nulls))
null_ok <- abs(mu["Q_LC"] - K / (N_OBJ - 1)) < 4 * se["Q_LC"] &&
  all(abs(mu[c("B_TC", "B_RE", "B_LC")]) < 4 * se[c("B_TC", "B_RE", "B_LC")])

pkg_src <- list.files(file.path(ANALYSIS_DIR, "..", "R"), full.names = TRUE)
gate <- list(
  schema = "stage4-gate-v2",
  created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
  engine = "joint-draw Westfall-Young min-P",
  assessment_R_sha256 = digest::digest(file.path(ANALYSIS_DIR, "..", "R", "assessment.R"),
                                       file = TRUE, algo = "sha256"),
  package_version = as.character(utils::packageVersion("QAIDR")),
  design = list(n_rep = N_REP_CAL, n = N_OBJ, K = K, m = M_CAL, alpha = ALPHA,
                metrics = c(METRICS, "Centers-Euclidean"), baseline = TRUE,
                families = fams,
                seed_rule = "data SEED_BASE+20/+1; null embeddings SEED_BASE+5e6+r; permutations SEED_BASE+6e6+r"),
  per_test_rates = as.list(as.data.frame(per_rates)),
  familywise_rates = as.list(fam_rates),
  familywise_ci_halfwidth_95 = 1.96 * se_bin,
  acceptance_rule = rule,
  gates = list(null_expectation = unname(null_ok),
               per_test_type1 = unname(per_pass),
               familywise_type1 = unname(fam_pass)),
  null_check = list(means = as.list(mu), mc_se = as.list(se),
                    expected_Q_LC = K / (N_OBJ - 1)),
  calibration_csv = "stage4v2_calibration_per_rep.csv",
  calibration_csv_sha256 = cal_hash,
  runtime_sec = elapsed,
  overall_pass = unname(null_ok && per_pass && fam_pass)
)
old <- file.path(OUTPUT_DIR, "stage4_gates.json")
if (file.exists(old)) file.copy(old, file.path(OUTPUT_DIR, "stage4_gates_v1_deprecated.json"),
                                overwrite = TRUE)
jsonlite::write_json(gate, file.path(OUTPUT_DIR, "stage4_gates.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = 6)
## Hash sidecar so the quick-mode runner can validate the gate file itself
## (not just parse it) before accepting a cached run.
jsonlite::write_json(
  list(file = "stage4_gates.json",
       sha256 = digest::digest(file.path(OUTPUT_DIR, "stage4_gates.json"),
                               file = TRUE, algo = "sha256"),
       provenance_status = "contemporaneous"),
  file.path(OUTPUT_DIR, "stage4_gates_sha256.json"),
  auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("Familywise rates: TC %.3f RE %.3f LC %.3f (CI half-width %.3f); overall_pass=%s\n",
            fam_rates["TC"], fam_rates["RE"], fam_rates["LC"], 1.96 * se_bin,
            gate$overall_pass))
if (!gate$overall_pass) quit(status = 1)
