# Build IMMUTABLE, version-labelled package source manifests.
#
# Why this exists
# ---------------
# The single `package_source_manifest.json` was regenerated in place and stamped
# with whatever QAIDR happened to be installed at the time
# (b2_provenance_backfill.r:18). That is how outputs produced under QAIDR 0.2.0
# came to reference a manifest labelled 0.3.0, and it is why 50 sidecars carry a
# stale binding to it.
#
# The fix is one manifest PER VERSION, written once and never rewritten. Each
# artifact then binds to the manifest of the version that actually produced it.
# 0.2.0 outputs stay bound to the 0.2.0 manifest; they are NOT relabelled.
#
# Usage
#   Rscript make_versioned_manifests.r <version> <source_tree> [--force]
# e.g.
#   Rscript make_versioned_manifests.r 0.3.0 ../../            # current tree
#   Rscript make_versioned_manifests.r 0.2.0 /path/to/QAIDR-0.2.0-archive
#
# Refuses to overwrite an existing manifest unless --force is given, because the
# whole point is immutability.

SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L)
  stop("usage: make_versioned_manifests.r <version> <source_tree> [--force]")
ver <- args[1]
src <- normalizePath(args[2], winslash = "/", mustWork = TRUE)
force <- "--force" %in% args

## _common.r pins QAIDR to the release version, which we do NOT want to require
## here: this script must be runnable to describe an OLD version too. So resolve
## paths directly rather than sourcing it.
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", "..", "..", ".."),
                              winslash = "/", mustWork = TRUE)
PROV_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", "provenance"),
                          winslash = "/", mustWork = TRUE)

sha256_file <- function(f) digest::digest(f, file = TRUE, algo = "sha256")

## The manifest covers the package's own source: DESCRIPTION, NAMESPACE, and
## every R/ file. Case-insensitive extension match, because the historical tree
## and the current tree do not agree on ".R" vs ".r".
pkg_files <- c(file.path(src, "DESCRIPTION"),
               file.path(src, "NAMESPACE"),
               list.files(file.path(src, "R"), full.names = TRUE,
                          pattern = "\\.[rR]$"))
pkg_files <- pkg_files[file.exists(pkg_files)]
if (!length(pkg_files)) stop("no package source files found under ", src)

## Verify the tree really is the version claimed.
desc_ver <- read.dcf(file.path(src, "DESCRIPTION"), fields = "Version")[1, 1]
if (!identical(as.character(desc_ver), ver))
  stop("DESCRIPTION says version ", desc_ver, " but ", ver, " was requested")

rel <- function(f) sub(paste0("^", src, "/"), "", f, fixed = FALSE)

files <- setNames(lapply(pkg_files, sha256_file), vapply(pkg_files, rel, character(1)))
files <- files[order(names(files))]      # deterministic ordering

out <- file.path(PROV_DIR, sprintf("package_source_manifest_qaidr-%s.json", ver))
if (file.exists(out) && !force)
  stop("refusing to overwrite immutable manifest: ", basename(out),
       " (pass --force only if you are certain)")

man <- list(
  package = "QAIDR",
  version = ver,
  immutable = TRUE,
  note = paste("Immutable source manifest for QAIDR", ver,
               "-- describes the source that produced artifacts recorded as",
               "this version. Never regenerate in place; create a new",
               "version-labelled manifest instead."),
  source_tree = if (identical(ver, "0.2.0"))
                  "recovered from sadm_submission/qaidr-0.2.0-data-code-archive.zip"
                else "Rcode-QAIDR/QAIDR (working tree)",
  n_files = length(files),
  created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
  files = files
)

jsonlite::write_json(man, out, auto_unbox = TRUE, pretty = TRUE)

## self-hash sidecar so the manifest itself is tamper-evident
sc <- sub("\\.json$", "_sha256.json", out)
jsonlite::write_json(list(manifest = basename(out),
                          sha256 = sha256_file(out),
                          created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
                     sc, auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("wrote %s (%d files)\n", basename(out), length(files)))
cat(sprintf("      %s\n", basename(sc)))
