# Master runner for the QA-IDR manuscript analysis pipeline (B2 revision).
#
# Usage:
#   QAIDR_RUN_MODE=quick Rscript analysis/make.R    # validation pass (_QUICK outputs)
#   QAIDR_RUN_MODE=full  Rscript analysis/make.R    # full manuscript rebuild
#   QAIDR_PROJECT_ROOT=<02-QA-IDR root>             # optional override
#
# Contract (B2 Part 2.1): dependency-ordered, fail-fast. Full mode recomputes
# every stage and ignores _QUICK artifacts. Quick mode may reuse an existing
# FULL artifact only when its recorded output hash still matches the file on
# disk (reported as CACHED_FULL_VALIDATED - a validated cached run, NOT a
# clean rebuild); otherwise it runs the stage in quick mode. Time expectations
# are printed before every stage (basis: measured under 4-way CPU contention
# during Phase B, or projected/unmeasured as stated).

SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
ANALYSIS_ROOT <- dirname(SCRIPT_FILE)
SCRIPTS <- file.path(ANALYSIS_ROOT, "scripts")
OUT <- file.path(ANALYSIS_ROOT, "output")
MODE <- Sys.getenv("QAIDR_RUN_MODE", "full")

# `out` is the primary artifact (post-run existence check); `arts` is the
# COMPLETE list of provenance-tracked artifacts the stage writes. Quick-mode
# cache validation checks every entry in `arts` (CSV and RDS hashes against
# the recorded sidecar), not just the primary, plus any gate JSON.
stages <- list(
  list(id = "stage4_gates",    script = "stage4_gate2.R",
       out = "stage4v2_calibration_per_rep",
       arts = "stage4v2_calibration_per_rep", gate_json = "stage4_gates.json",
       expect = "40-80 min (measured basis: 200 reps x m=199)", critical = TRUE),
  list(id = "ushcn_pilot",     script = "midscale_pilot_ghcn.R",
       out = "midscale_pilot_ghcn_rubric", arts = "midscale_pilot_ghcn_rubric",
       expect = "6-10 min + one-time 316 MB download if data-raw absent"),
  list(id = "sc2_warnings",    script = "scenario2_warning_diagnostic.R",
       out = "scenario2_warning_inventory", arts = "scenario2_warning_inventory",
       expect = "10-20 min (2 diagnostic replications)"),
  list(id = "scenario1",       script = "scenario1.R",  out = "scenario1_summary",
       arts = c("scenario1_summary", "scenario1_per_rep",
                "scenario1_calibration", "scenario1_tie_flags"),
       expect = "~4.2 h measured (100 reps + calibration)"),
  list(id = "scenario2",       script = "scenario2.R",  out = "scenario2_summary",
       arts = c("scenario2_summary", "scenario2_per_rep", "scenario2_calibration",
                "scenario2_factorial_per_rep", "scenario2_factorial_summary",
                "scenario2_tie_flags"),
       expect = "~20.5 h measured (100 reps + factorial + calibration)"),
  list(id = "scenario3",       script = "scenario3.R",  out = "scenario3_per_perturbation",
       arts = "scenario3_per_perturbation",
       expect = "~31.2 h measured (100 datasets x 100 perturbations)"),
  list(id = "scenario3_agg",   script = "b2_fix_scenario3_aggregation.R",
       out = "scenario3_summary", arts = "scenario3_summary", expect = "<1 min"),
  list(id = "face",            script = "face.R",       out = "face_indices",
       arts = c("face_indices", "face_assessment_results", "face_pvalues_raw",
                "face_pvalues_minP", "face_k_profiles"),
       expect = "5-10 min"),
  list(id = "benchmark",       script = "benchmark.R",  out = "scalability_benchmark",
       arts = "scalability_benchmark",
       expect = "~4.6 h measured (dominated by n=2000 DR methods)"),
  list(id = "ushcn",           script = "midscale_ghcn.R", out = "midscale_indices",
       arts = c("midscale_indices", "midscale_assessment_results",
                "midscale_pvalues_raw", "midscale_pvalues_minP",
                "midscale_k_profiles", "midscale_runtime"),
       expect = "~10 min measured"),
  list(id = "b2_scenario4",    script = "b2_scenario4.R", out = "scenario4_summary",
       arts = c("scenario4_summary", "scenario4_per_rep", "scenario4_calibration",
                "scenario4_baseline_contrast"),
       expect = "~8-10 h (pilot-projected; 100 reps, 50-resolution tie audit)"),
  list(id = "b2_k_profiles",   script = "b2_k_profiles.R", out = "b2_k_profiles_summary",
       arts = c("b2_k_profiles_summary", "b2_k_profiles_per_rep"),
       expect = "~4-6 h (75 re-derived embedding sets; measured ~4.3 h)"),
  list(id = "b2_k_merge",      script = "b2_k_profiles_merge.R", out = "b2_k_stability_kendall",
       arts = c("b2_k_stability_kendall", "b2_k_profiles_summary",
                "b2_k_profiles_per_rep"),
       expect = "<2 min (re-aggregation incl. Scenario IV)"),
  list(id = "b2_sensitivity",  script = "b2_sensitivity.R", out = "b2_sensitivity",
       arts = c("b2_sensitivity", "b2_sensitivity_exceptional_summary",
                "b2_sensitivity_endpoint_audit"),
       expect = "~1-2 h (Face exact enumeration + grids; unmeasured)"),
  list(id = "b2_cpca",         script = "b2_cpca_diagnostics.R", out = "b2_cpca_diagnostics",
       arts = "b2_cpca_diagnostics",
       expect = "~40 min (25 C-PCA refits; unmeasured)"),
  list(id = "regen_calibration", script = "regen_calibration.R", out = "scenario1_calibration",
       arts = c("scenario1_calibration", "scenario2_calibration"),
       expect = "~2.5 h (50 re-derived embedding sets + m=999 inference; QAIDR >= 0.3.0)"),
  list(id = "h_s4_calibration", script = "h_regen_s4_calibration.R", out = "scenario4_calibration",
       arts = "scenario4_calibration",
       expect = "~1 h (25 re-derived embedding sets + m=999 inference; QAIDR >= 0.3.0)"),
  list(id = "h_audit_v3c",     script = "h_audit_v3c.R", out = "h_audit_v3c_familywise",
       arts = "h_audit_v3c_familywise",
       expect = "~20 min (2,000 null replicates at m=999, symmetric min-P)"),
  list(id = "h_ushcn_lambda",  script = "h_ushcn_lambda.R", out = "h_ushcn_lambda_profile",
       arts = c("h_ushcn_lambda_profile", "h_ushcn_lambda_components",
                "h_ushcn_lambda_runtime"),
       expect = "~15 min measured (frozen 21-point lambda grid, tie-audit path)"),
  list(id = "h_graded_width",  script = "h_graded_width.R", out = "h_graded_width_summary",
       arts = c("h_graded_width_per_rep", "h_graded_width_summary",
                "h_graded_width_runtime"),
       expect = "~4-6 h (75 embedding sets: L3 reconstruction first, then L1/L2)"),
  list(id = "h_tie_policy",    script = "h_tie_policy_audit.R", out = "h_tie_policy_audit",
       arts = "h_tie_policy_audit",
       expect = "~5-6 h (A3-scoped flagged cells; 50 vs 500 vs coupled-priority)"),
  list(id = "tables",          script = "tables_to_tex.R", out = NULL,
       expect = "<1 min", produces_dir = "tex"),
  list(id = "verify",          script = "verify_manuscript_numbers.R", out = NULL,
       expect = "<2 min (includes deliberate-failure self-test)")
)

`%||%` <- function(x, y) if (is.null(x)) y else x
sha256_file <- function(f) digest::digest(f, file = TRUE, algo = "sha256")

# Validate ONE provenance-tracked artifact: sidecar present, CSV and RDS both
# present with hashes equal to the recorded ones. Returns "" if valid, else a
# reason string.
validate_artifact <- function(name) {
  mj <- file.path(OUT, paste0(name, "_meta.json"))
  if (!file.exists(mj)) return(paste0(name, ": no provenance sidecar"))
  m <- jsonlite::read_json(mj)
  for (ext in c("csv", "rds")) {
    rec <- m$output_hashes[[ext]]
    if (is.null(rec)) next
    f <- file.path(OUT, paste0(name, ".", ext))
    if (!file.exists(f)) return(paste0(name, ".", ext, ": missing"))
    if (!identical(rec, sha256_file(f)))
      return(paste0(name, ".", ext, ": hash mismatch vs recorded provenance"))
  }
  ""
}

# Validate ALL artifacts of a stage (primary AND secondary) plus its gate
# JSON. This is the quick-mode cache gate: any failure forces a live rerun.
validate_stage_cache <- function(st) {
  for (a in st$arts) {
    r <- validate_artifact(a)
    if (nzchar(r)) return(r)
  }
  if (!is.null(st$gate_json)) {
    gj <- file.path(OUT, st$gate_json)
    if (!file.exists(gj)) return(paste0(st$gate_json, ": missing"))
    sc <- file.path(OUT, sub("\\.json$", "_sha256.json", st$gate_json))
    if (!file.exists(sc)) return(paste0(basename(sc), ": missing hash sidecar"))
    rec <- tryCatch(jsonlite::read_json(sc)$sha256, error = function(e) NULL)
    if (is.null(rec) || !identical(rec, sha256_file(gj)))
      return(paste0(st$gate_json, ": hash mismatch vs recorded sidecar"))
    g <- tryCatch(jsonlite::read_json(gj), error = function(e) NULL)
    if (is.null(g) || !isTRUE(g$overall_pass))
      return(paste0(st$gate_json, ": overall_pass is not TRUE"))
  }
  ""
}

log_dir <- file.path(OUT, "logs"); dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
run_manifest <- list()
cat(sprintf("QAIDR pipeline | mode=%s | %d stages | root=%s\n", MODE, length(stages),
            ANALYSIS_ROOT))
cat("NOTE: sequential full rebuild is ~70+ hours of compute on the reference laptop;\n")
cat("      long stages can be run in parallel R sessions (see README_REPRODUCIBILITY.md).\n\n")

for (st in stages) {
  cat(sprintf("== %-15s %s | expect: %s\n", st$id, st$script, st$expect))
  cached <- FALSE
  if (MODE == "quick" && !is.null(st$out)) {
    reason <- validate_stage_cache(st)
    if (!nzchar(reason)) {
      cat(sprintf("   CACHED_FULL_VALIDATED (%d artifact(s)%s validated against recorded provenance)\n",
                  length(st$arts), if (is.null(st$gate_json)) "" else " + gate JSON"))
      run_manifest[[st$id]] <- list(status = "CACHED_FULL_VALIDATED",
                                    artifacts = file.path(OUT, paste0(st$arts, ".csv")))
      cached <- TRUE
    } else {
      cat("   cache REJECTED (", reason, ") -> running stage in quick mode\n", sep = "")
    }
  }
  if (!cached) {
    lf <- file.path(log_dir, sub("\\.R$", ".log", st$script))
    t0 <- Sys.time()
    status <- system2(file.path(R.home("bin"), "Rscript"),
                      args = shQuote(file.path(SCRIPTS, st$script)),
                      stdout = lf, stderr = lf)
    el <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
    cat(sprintf("   -> exit %d in %.1f min (log %s)\n", status, el, basename(lf)))
    run_manifest[[st$id]] <- list(status = if (status == 0) "RAN" else "FAILED",
                                  exit = status, minutes = el, log = lf)
    if (status != 0) {
      jsonlite::write_json(run_manifest, file.path(OUT, "run_manifest.json"),
                           auto_unbox = TRUE, pretty = TRUE)
      stop("Pipeline stopped at stage '", st$id, "'; see ", lf)
    }
    suffix <- if (MODE == "quick") "_QUICK" else ""
    if (!is.null(st$out)) {
      expect_csv <- file.path(OUT, paste0(st$out, suffix, ".csv"))
      if (!file.exists(expect_csv) || file.size(expect_csv) == 0)
        stop("Stage '", st$id, "' finished but expected output missing/empty: ", expect_csv)
    }
    if (!is.null(st$gate_json)) {
      g <- jsonlite::read_json(file.path(OUT, st$gate_json))
      if (!isTRUE(g$overall_pass)) stop("Gate failed in ", st$gate_json)
      cat("   gate: overall_pass = TRUE (engine: ", g$engine, ")\n", sep = "")
    }
  }
}
jsonlite::write_json(run_manifest, file.path(OUT, "run_manifest.json"),
                     auto_unbox = TRUE, pretty = TRUE)
cat("\nPipeline complete:", format(Sys.time()), "- see run_manifest.json\n")
