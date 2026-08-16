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
  } else stop("_common.r: define SCRIPT_DIR before source() or run via Rscript")
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
  stop("_common.r: sentinel Rcode-QAIDR/QAIDR/DESCRIPTION not found under ",
       PROJECT_ROOT, "; set QAIDR_PROJECT_ROOT to the 02-QA-IDR root")

PACKAGE_ROOT <- file.path(PROJECT_ROOT, "Rcode-QAIDR", "QAIDR")
ANALYSIS_DIR <- file.path(PACKAGE_ROOT, "analysis")
OUTPUT_DIR <- file.path(ANALYSIS_DIR, "output")
DATA_RAW <- file.path(ANALYSIS_DIR, "data-raw")
DATA_PROC <- file.path(ANALYSIS_DIR, "data-processed")
PROV_DIR <- file.path(ANALYSIS_DIR, "provenance")
## Case-exact paths. On disk the working directory is lowercase `latex/` and the
## master is lowercase `qa-intervaldr.tex`. The previous spelling ("LaTeX",
## "QA-IntervalDR.tex") resolved only on case-insensitive filesystems and would
## fail on a case-sensitive host. assert_case_exact() re-checks each component
## against the real directory listing, so a future drift fails here rather than
## silently on someone else's machine.
assert_case_exact <- function(path, what = "path") {
  ## Walk the path AS WRITTEN. normalizePath() would silently canonicalize the
  ## spelling to whatever is on disk (Windows and macOS are case-insensitive),
  ## which is exactly the difference this function exists to detect.
  path <- gsub("\\\\", "/", path)
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  if (length(parts) < 2L) return(invisible(TRUE))
  cur <- parts[1]
  if (!nzchar(cur)) {
    cur <- "/"
  } else if (grepl("^[A-Za-z]:$", cur)) {
    ## On Windows "D:" means "the current directory on drive D", not its root;
    ## list.files("D:") would therefore enumerate the wrong directory.
    cur <- paste0(cur, "/")
  }
  for (p in parts[-1]) {
    entries <- tryCatch(list.files(cur, all.files = TRUE, no.. = TRUE),
                        error = function(e) character())
    if (!length(entries)) break
    if (!(p %in% entries)) {
      hit <- entries[tolower(entries) == tolower(p)]
      stop("_common.r: case-inexact ", what, ": component '", p, "' under '",
           cur, "' ",
           if (length(hit))
             paste0("is actually spelled '", hit[1],
                    "'. Fix the literal in the script; do not rely on a ",
                    "case-insensitive filesystem.")
           else "does not exist.")
    }
    cur <- file.path(cur, p)
  }
  invisible(TRUE)
}

MANUSCRIPT_TEX <- file.path(PROJECT_ROOT, "LaTeX", "QA-IntervalDR_20260830_SADM",
                            "latex", "qa-intervaldr.tex")
MANUSCRIPT_DIR <- dirname(MANUSCRIPT_TEX)
if (file.exists(MANUSCRIPT_TEX)) assert_case_exact(MANUSCRIPT_TEX, "manuscript path")
for (d in c(OUTPUT_DIR, DATA_RAW, DATA_PROC, PROV_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

qlib <- Sys.getenv("QAIDR_R_LIB", "")
if (nzchar(qlib) && dir.exists(qlib)) .libPaths(c(qlib, .libPaths()))

RUN_MODE <- Sys.getenv("QAIDR_RUN_MODE", "full")
stopifnot(RUN_MODE %in% c("full", "quick"))
SEED_BASE <- 20260715L

## Exact release pin. PRODUCTION_RERUN_DECISION.md records that substitutions
## are NOT acceptable for release comparison: R 4.6.1, QAIDR exactly 0.3.0,
## symbolicDA 0.6.2, coRanking 0.2.5. The previous floor (">= 0.2.0") is the
## mechanism that let 18 of 26 session records be produced under QAIDR 0.2.0
## while the manuscript claimed 0.3.0. It must fail closed instead.
suppressMessages(library(QAIDR))
QAIDR_REQUIRED_VERSION <- "0.3.0"
.qaidr_found <- utils::packageVersion("QAIDR")
if (.qaidr_found != package_version(QAIDR_REQUIRED_VERSION))
  stop("_common.r: QAIDR exactly ", QAIDR_REQUIRED_VERSION,
       " is required for release-comparable results; found ", .qaidr_found,
       ". Install the exact release version (do not substitute a newer one) ",
       "or set QAIDR_R_LIB to a library containing it.")

## Record the stack actually in force so every provenance sidecar can attest it.
RELEASE_STACK <- list(
  R          = paste(R.version$major, R.version$minor, sep = "."),
  QAIDR      = as.character(.qaidr_found),
  symbolicDA = tryCatch(as.character(utils::packageVersion("symbolicDA")),
                        error = function(e) NA_character_),
  coRanking  = tryCatch(as.character(utils::packageVersion("coRanking")),
                        error = function(e) NA_character_)
)

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
  cfile <- file.path(SCRIPT_DIR, "_common.r")
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
    ## Named, not positional. The previous form serialized to a bare JSON array,
    ## so a reader could not tell which version belonged to which package (and a
    ## package that failed to resolve silently shortened the array, shifting
    ## every later element). Cache validation now checks packages$QAIDR by name.
    packages = {
      .pkgs <- c("QAIDR", "dataSDA", "symbolicDA", "RSDA", "umap", "coRanking")
      stats::setNames(
        lapply(.pkgs, function(p)
          tryCatch(as.character(utils::packageVersion(p)),
                   error = function(e) NA_character_)),
        .pkgs)
    },
    release_stack_required = list(R = "4.6.1", QAIDR = QAIDR_REQUIRED_VERSION,
                                  symbolicDA = "0.6.2", coRanking = "0.2.5"),
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
