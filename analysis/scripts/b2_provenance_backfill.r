# B2 Part 2.4: provenance manifests + metadata backfill.
# Creates analysis/provenance/{package_source_manifest,seed_manifest}.json and
# upgrades every existing *_meta.json to named SHA-256 input-hash maps per the
# Codex scheme. Backfilled records are marked provenance_status="backfilled"
# and record the CURRENT (unchanged) output hashes; they do not claim the
# current script hash was the original producer hash (the original run's
# script may have since been edited for path portability). Numeric outputs
# are never modified. RDS-only and PDF outputs receive new sidecars.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))

## ---- package source manifest ----------------------------------------------
pkg_files <- c(file.path(PACKAGE_ROOT, "DESCRIPTION"),
               file.path(PACKAGE_ROOT, "NAMESPACE"),
               list.files(file.path(PACKAGE_ROOT, "R"), full.names = TRUE, pattern = "\\.R$"))
psm <- list(package = "QAIDR",
            version = as.character(utils::packageVersion("QAIDR")),
            created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
            files = setNames(lapply(pkg_files, sha256_file),
                             vapply(pkg_files, rel_path, character(1))))
jsonlite::write_json(psm, file.path(PROV_DIR, "package_source_manifest.json"),
                     auto_unbox = TRUE, pretty = TRUE)

## ---- external CRAN data manifest ------------------------------------------
external_data_manifest <- list()
if (requireNamespace("dataSDA", quietly = TRUE)) {
  e <- new.env(parent = emptyenv())
  utils::data("face.iGAP", package = "dataSDA", envir = e)
  utils::data("cars.int", package = "dataSDA", envir = e)
  external_data_manifest <- list(
    package = "dataSDA",
    version = as.character(utils::packageVersion("dataSDA")),
    objects = list(
      face.iGAP = digest::digest(e$face.iGAP, algo = "sha256", serialize = TRUE),
      cars.int = digest::digest(e$cars.int, algo = "sha256", serialize = TRUE)
    ),
    note = "SHA-256 digests use digest serialization; analysis uses face.iGAP."
  )
}
jsonlite::write_json(external_data_manifest,
                     file.path(PROV_DIR, "external_data_manifest.json"),
                     auto_unbox = TRUE, pretty = TRUE)

## ---- seed manifest ----------------------------------------------------------
seeds <- list(
  SEED_BASE = 20260715,
  rules = list(
    scenario1 = "data SEED_BASE+1000+r; DR +500000; tie audits +r or tie_seed=seed_r; calibration permutations +900000",
    scenario2 = "data SEED_BASE+2000+r; DR raw +500000, std +600000; calibration +900000",
    scenario3 = "data SEED_BASE+3000+d; DR +500000; perturbation +700000+pert",
    face = "DR SEED_BASE+4000; inference +4001",
    benchmark = "per size SEED_BASE+5000+n",
    ushcn = "DR SEED_BASE+6000; inference +6001",
    scenario4 = "data SEED_BASE+7000+r; DR +500000; tie resolutions +300000+s; calibration +900000",
    cpca_isotropic_null = "SEED_BASE+7900+r",
    stage4_gate_v2 = "data SEED_BASE+20/+1; null embeddings +5e6+r; permutations +6e6+r"
  ))
jsonlite::write_json(seeds, file.path(PROV_DIR, "seed_manifest.json"),
                     auto_unbox = TRUE, pretty = TRUE)

## ---- design manifests -------------------------------------------------------
for (d in list(
  list(id = "phaseB", file = file.path(PROJECT_ROOT, "Octopus_Workflow", "PROMPT_B_IMPLEMENTATION.md")),
  list(id = "B2", file = file.path(ANALYSIS_DIR, "B2_PREREGISTRATION.md")))) {
  jsonlite::write_json(list(design = d$id, path = rel_path(d$file),
                            sha256 = sha256_file(d$file),
                            recorded = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
                       file.path(PROV_DIR, paste0("design_", d$id, ".json")),
                       auto_unbox = TRUE, pretty = TRUE)
}

## ---- input families ---------------------------------------------------------
S <- function(...) file.path(ANALYSIS_DIR, "scripts", ...)
GH <- file.path(DATA_RAW, "ghcnd")
fam <- list(
  "^stage4_pilot"        = c(S("stage4_pilot.r"), S("scenarios.r")),
  "^stage4v2_"           = c(S("stage4_gate2.r"), S("scenarios.r")),
  "^scenario1_"          = c(S("scenario1.r"), S("scenarios.r")),
  "^scenario2_warning"   = c(S("scenario2_warning_diagnostic.r"), S("scenarios.r")),
  "^scenario2_"          = c(S("scenario2.r"), S("scenarios.r")),
  "^scenario3_per"       = c(S("scenario3.r"), S("scenarios.r")),
  "^scenario3_summary"   = c(S("b2_fix_scenario3_aggregation.r"),
                             file.path(OUTPUT_DIR, "scenario3_per_perturbation.csv")),
  "^scenario4_"          = c(S("b2_scenario4.r"), S("scenarios.r"),
                             file.path(ANALYSIS_DIR, "B2_PREREGISTRATION.md")),
  "^scalability_"        = c(S("benchmark.r"), S("scenarios.r")),
  "^face_"               = c(S("face.r"), S("scenarios.r"),
                             file.path(PROV_DIR, "external_data_manifest.json")),
  "^midscale_pilot"      = c(S("midscale_pilot_ghcn.r"),
                             file.path(GH, "ghcnd_hcn.tar.gz"),
                             file.path(GH, "ghcnd-stations.txt"),
                             file.path(GH, "ghcnd-version.txt"),
                             file.path(GH, "readme.txt")),
  "^midscale_"           = c(S("midscale_ghcn.r"), S("scenarios.r"),
                             file.path(DATA_PROC, "ghcn_seasonal_intervals.rds"),
                             file.path(GH, "ghcnd_hcn.tar.gz")),
  "^b2_k_"               = c(S("b2_k_profiles.r"), S("scenarios.r"),
                             file.path(ANALYSIS_DIR, "B2_PREREGISTRATION.md")),
  "^b2_sensitivity"      = c(S("b2_sensitivity.r"), S("scenarios.r"),
                             file.path(ANALYSIS_DIR, "B2_PREREGISTRATION.md"),
                             file.path(PROV_DIR, "external_data_manifest.json")),
  "^b2_cpca"             = c(S("b2_cpca_diagnostics.r"), S("scenarios.r"),
                             file.path(ANALYSIS_DIR, "B2_PREREGISTRATION.md"),
                             file.path(PROJECT_ROOT, "R-packages", "symbolicDA", "R", "PCA.SDA.r"))
)

metas <- list.files(OUTPUT_DIR, pattern = "_meta\\.json$", full.names = TRUE)
n_up <- 0L
for (mj in metas) {
  m <- jsonlite::read_json(mj)
  nm <- sub("_meta\\.json$", "", basename(mj))
  hit <- NULL
  for (pat in names(fam)) if (grepl(pat, nm)) { hit <- fam[[pat]]; break }
  if (is.null(hit)) next
  already_named <- is.list(m$input_hashes) && length(names(m$input_hashes) %||% character(0)) > 0
  ih <- if (already_named) m$input_hashes else list()
  for (f in hit) if (file.exists(f)) ih[[rel_path(f)]] <- sha256_file(f)
  ih[[rel_path(file.path(PROV_DIR, "package_source_manifest.json"))]] <-
    sha256_file(file.path(PROV_DIR, "package_source_manifest.json"))
  ih[[rel_path(file.path(PROV_DIR, "seed_manifest.json"))]] <-
    sha256_file(file.path(PROV_DIR, "seed_manifest.json"))
  m$input_hashes <- ih
  if (!already_named) m$provenance_status <- m$provenance_status %||% "backfilled"
  csv <- file.path(OUTPUT_DIR, paste0(nm, ".csv"))
  rds <- file.path(OUTPUT_DIR, paste0(nm, ".rds"))
  m$output_hashes <- list(csv = if (file.exists(csv)) sha256_file(csv) else NULL,
                          rds = if (file.exists(rds)) sha256_file(rds) else NULL)
  jsonlite::write_json(m, mj, auto_unbox = TRUE, pretty = TRUE)
  n_up <- n_up + 1L
}
cat("Upgraded", n_up, "metadata sidecars.\n")

## ---- sidecars for RDS-only / PDF artifacts ---------------------------------
extra_art <- c("face_assessment_full.rds", "face_effective_params.rds",
               "midscale_effective_params.rds", "stage4_pilot_report.rds",
               "scenario4_records.rds",
               "RealData_QBK_Int-Euclidean.pdf", "RealData_QBK_Hausdorff.pdf",
               "RealData_QBK_Ichino-Yaguchi.pdf", "RealData_QBK_Wasserstein.pdf",
               "RealData_QBK_Centers-Euclidean.pdf", "RealData_2D_Projections.pdf")
n_side <- 0L
for (a in extra_art) {
  f <- file.path(OUTPUT_DIR, a)
  if (!file.exists(f)) next
  fam_hit <- if (grepl("^midscale", a)) fam[["^midscale_"]]
             else if (grepl("^stage4", a)) fam[["^stage4_pilot"]]
             else if (grepl("^scenario4", a)) fam[["^scenario4_"]]
             else fam[["^face_"]]
  ih <- list(); for (ff in fam_hit) if (file.exists(ff)) ih[[rel_path(ff)]] <- sha256_file(ff)
  jsonlite::write_json(list(name = a, provenance_status = "backfilled",
                            created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
                            input_hashes = ih,
                            output_hashes = list(file = sha256_file(f)),
                            note = "sidecar added by b2_provenance_backfill.r"),
                       file.path(OUTPUT_DIR, paste0(a, ".prov.json")),
                       auto_unbox = TRUE, pretty = TRUE)
  n_side <- n_side + 1L
}
## GHCN processed snapshot sidecar
gs <- file.path(DATA_PROC, "ghcn_seasonal_intervals.rds")
if (file.exists(gs)) {
  ih <- list(); for (ff in fam[["^midscale_pilot"]]) if (file.exists(ff)) ih[[rel_path(ff)]] <- sha256_file(ff)
  jsonlite::write_json(list(name = "ghcn_seasonal_intervals.rds",
                            provenance_status = "backfilled",
                            input_hashes = ih,
                            output_hashes = list(file = sha256_file(gs))),
                       paste0(gs, ".prov.json"), auto_unbox = TRUE, pretty = TRUE)
  n_side <- n_side + 1L
}
cat("Wrote", n_side, "artifact sidecars.\n")
