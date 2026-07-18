# Shared infrastructure for the QA-IDR manuscript analysis pipeline (B2 rev).
#
# Path contract (B2 Part 2.2): no user-specific absolute paths. Entry scripts
# define SCRIPT_DIR (see header snippet) and source this file from it; the
# project root is derived by walking up from the scripts directory
# (scripts -> analysis -> QAIDR -> Rcode-QAIDR -> <project root>), and can be
# overridden with the environment variable QAIDR_PROJECT_ROOT. The sentinel
# Rcode-QAIDR/QAIDR/DESCRIPTION is validated. Library paths honor standard
# R_LIBS_USER; an optional QAIDR_R_LIB is prepended only when explicitly set.
#
# Provenance contract (B2 Part 2.4): save_result() writes named SHA-256
# input-hash maps, output hashes, producer, run mode, seeds where supplied,
# and provenance_status ("contemporaneous" here; "backfilled" is used only by
# the dedicated backfill utility).

`%||%` <- function(x, y) if (is.null(x)) y else x

.qaidr_gitbash_fix <- function(p) {
  if (.Platform$OS.type == "windows" && grepl("^/[A-Za-z]/", p))
    p <- paste0(toupper(substr(p, 2, 2)), ":", substr(p, 3, nchar(p)))
  p
}

if (!exists("SCRIPT_DIR") || !nzchar(SCRIPT_DIR %||% "")) {
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) {
    SCRIPT_DIR <- dirname(normalizePath(.qaidr_gitbash_fix(sub("^--file=", "", fa[1])),
                                        winslash = "/", mustWork = TRUE))
  } else stop("_common.R: define SCRIPT_DIR before source() or run via Rscript")
}

PROJECT_ROOT <- Sys.getenv("QAIDR_PROJECT_ROOT", "")
if (!nzchar(PROJECT_ROOT)) {
  PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", "..", "..", ".."),
                                winslash = "/", mustWork = TRUE)
} else {
  PROJECT_ROOT <- normalizePath(.qaidr_gitbash_fix(PROJECT_ROOT),
                                winslash = "/", mustWork = TRUE)
}
if (!file.exists(file.path(PROJECT_ROOT, "Rcode-QAIDR", "QAIDR", "DESCRIPTION")))
  stop("_common.R: sentinel Rcode-QAIDR/QAIDR/DESCRIPTION not found under ",
       PROJECT_ROOT, "; set QAIDR_PROJECT_ROOT to the 02-QA-IDR root")

PACKAGE_ROOT <- file.path(PROJECT_ROOT, "Rcode-QAIDR", "QAIDR")
ANALYSIS_DIR <- file.path(PACKAGE_ROOT, "analysis")
OUTPUT_DIR <- file.path(ANALYSIS_DIR, "output")
DATA_RAW <- file.path(ANALYSIS_DIR, "data-raw")
DATA_PROC <- file.path(ANALYSIS_DIR, "data-processed")
PROV_DIR <- file.path(ANALYSIS_DIR, "provenance")
MANUSCRIPT_TEX <- file.path(PROJECT_ROOT, "LaTeX", "QA-IntervalDR_20260830_SADM",
                            "LaTeX", "QA-IntervalDR.tex")
MANUSCRIPT_DIR <- dirname(MANUSCRIPT_TEX)
for (d in c(OUTPUT_DIR, DATA_RAW, DATA_PROC, PROV_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

qlib <- Sys.getenv("QAIDR_R_LIB", "")
if (nzchar(qlib) && dir.exists(qlib)) .libPaths(c(qlib, .libPaths()))

RUN_MODE <- Sys.getenv("QAIDR_RUN_MODE", "full")
stopifnot(RUN_MODE %in% c("full", "quick"))
SEED_BASE <- 20260715L

suppressMessages(library(QAIDR))
stopifnot(utils::packageVersion("QAIDR") >= "0.2.0")

METHODS <- c("C-PCA", "V-PCA", "MR-PCA", "SPCA", "IMDS", "Int-UMAP")
METRICS <- c("Int-Euclidean", "Hausdorff", "Ichino-Yaguchi", "Wasserstein")

sha256_file <- function(f) digest::digest(f, file = TRUE, algo = "sha256")

#' Relative path (to PROJECT_ROOT) for provenance records (fixed-prefix match).
rel_path <- function(f) {
  f <- tryCatch(normalizePath(f, winslash = "/", mustWork = FALSE), error = function(e) f)
  if (startsWith(f, paste0(PROJECT_ROOT, "/"))) substring(f, nchar(PROJECT_ROOT) + 2)
  else f
}

#' Save a result table (CSV + RDS) plus a provenance sidecar with NAMED
#' SHA-256 input hashes.
save_result <- function(obj, name, inputs = character(), extra = list(),
                        seeds = NULL) {
  stopifnot(is.character(name), length(name) == 1)
  if (RUN_MODE == "quick") name <- paste0(name, "_QUICK")
  t0 <- get0(".qaidr_t0", envir = globalenv(), ifnotfound = Sys.time())
  csv <- file.path(OUTPUT_DIR, paste0(name, ".csv"))
  rds <- file.path(OUTPUT_DIR, paste0(name, ".rds"))
  if (is.data.frame(obj)) utils::write.csv(obj, csv, row.names = FALSE)
  saveRDS(obj, rds)
  producer <- tryCatch(rel_path(get0("SCRIPT_FILE", ifnotfound =
    file.path(SCRIPT_DIR, "unknown.R"))), error = function(e) NA_character_)
  ih <- list()
  for (f in inputs) if (file.exists(f)) ih[[rel_path(f)]] <- sha256_file(f)
  cfile <- file.path(SCRIPT_DIR, "_common.R")
  ih[[rel_path(cfile)]] <- sha256_file(cfile)
  psm <- file.path(PROV_DIR, "package_source_manifest.json")
  if (file.exists(psm)) ih[[rel_path(psm)]] <- sha256_file(psm)
  meta <- list(
    name = name,
    producer = producer,
    provenance_status = "contemporaneous",
    run_mode = RUN_MODE,
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    max_gc_mb = sum(gc()[, 6]),
    r_version = R.version.string,
    packages = as.list(vapply(c("QAIDR", "dataSDA", "symbolicDA", "RSDA",
                                "umap", "coRanking"),
                              function(p) tryCatch(as.character(utils::packageVersion(p)),
                                                   error = function(e) NA_character_),
                              character(1))),
    input_hashes = ih,
    output_hashes = list(csv = if (file.exists(csv)) sha256_file(csv) else NULL,
                         rds = sha256_file(rds)),
    seeds = seeds,
    extra = extra
  )
  jsonlite::write_json(meta, file.path(OUTPUT_DIR, paste0(name, "_meta.json")),
                       auto_unbox = TRUE, pretty = TRUE)
  invisible(list(csv = csv, rds = rds))
}

save_sessioninfo <- function(name) {
  if (RUN_MODE == "quick") name <- paste0(name, "_QUICK")
  writeLines(c(utils::capture.output(sessionInfo()),
               "", paste("libPaths:", paste(.libPaths(), collapse = "; ")),
               paste("PROJECT_ROOT:", PROJECT_ROOT)),
             file.path(OUTPUT_DIR, paste0(name, "_sessionInfo.txt")))
}

n_perm_mode <- function(full = 999L, quick = 19L) if (RUN_MODE == "full") full else quick
n_rep_mode <- function(full = 100L, quick = 3L) if (RUN_MODE == "full") full else quick

assess_cell_tie_audit <- function(Dh, Dl, K, n_res = 50L, seed = 1L) {
  primary <- tryCatch(
    list(vals = coranking_indices(Dh, Dl, K), tie_affected = FALSE,
         max_spread = 0, n_res = 0L),
    error = function(e) NULL)
  if (!is.null(primary)) return(primary)
  vals <- matrix(NA_real_, n_res, 6)
  for (s in seq_len(n_res)) {
    set.seed(seed + s)
    vals[s, ] <- coranking_indices(Dh, Dl, K, ties = "random")
  }
  colnames(vals) <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  list(vals = colMeans(vals), tie_affected = TRUE,
       max_spread = max(apply(vals, 2, function(v) diff(range(v)))),
       n_res = n_res)
}

evaluate_cells <- function(x_eval, projections, K, lambda = 0.5, nu = 0.5,
                           metrics = METRICS, tie_seed = 1L) {
  rows <- list(); ii <- 0L
  Dh_cache <- list()
  for (met in c(metrics, "Centers-Euclidean")) {
    Dh_cache[[met]] <- if (met == "Centers-Euclidean")
      as.matrix(stats::dist(x_eval$centers))
    else idist(x_eval$centers, x_eval$radii, met, lambda = lambda, nu = nu)
  }
  for (m in names(projections)) {
    pr <- projections[[m]]
    for (met in c(metrics, "Centers-Euclidean")) {
      Dl <- if (met == "Centers-Euclidean" || pr$type == "Point")
        as.matrix(stats::dist(pr$C))
      else idist(pr$C, pr$R, met, lambda = lambda, nu = nu)
      cell <- assess_cell_tie_audit(Dh_cache[[met]], Dl, K, seed = tie_seed)
      ii <- ii + 1L
      rows[[ii]] <- data.frame(Method = m, Metric = met, t(cell$vals),
                               tie_affected = cell$tie_affected,
                               tie_spread = cell$max_spread,
                               stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

.qaidr_t0 <- Sys.time()
