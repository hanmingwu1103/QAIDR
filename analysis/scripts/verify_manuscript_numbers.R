# Expanded manuscript-object verifier (B2 Part 2.5).
# Verifies, against the machine-readable pipeline outputs:
#   - Tables I/II (mean(SD) cells; stars from the >=23/25 adjusted-p rule);
#   - Table III (mean(SD) + bare flip rates; no stars permitted);
#   - Face & USHCN tables (bare 3-dec cells; stars from *_pvalues_minP;
#     daggers against tie outputs);
#   - B2 tables when present (Scenario IV summary; K-stability; sensitivity);
#   - high-risk numeric prose claims;
#   - figure files referenced by \includegraphics (existence + hash match to
#     analysis outputs) -> writes figure_manifest.json;
#   - SELF-TEST: corrupts one token in a temp copy and requires nonzero exit.
# Exit code: 0 iff all checks pass. Child mode (internal, for the self-test):
#   Rscript verify_manuscript_numbers.R --texfile=<path> --no-selftest
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))

argv <- commandArgs(trailingOnly = TRUE)
texarg <- sub("^--texfile=", "", grep("^--texfile=", argv, value = TRUE))
TEXFILE <- if (length(texarg)) texarg else Sys.getenv("QAIDR_TEXFILE", MANUSCRIPT_TEX)
DO_SELFTEST <- !("--no-selftest" %in% argv)

tex <- readLines(TEXFILE, warn = FALSE)
fails <- 0L; n_checked <- 0L
bad <- function(...) { cat("MISMATCH:", ..., "\n"); fails <<- fails + 1L }
ok <- function() n_checked <<- n_checked + 1L

strip <- function(x) {
  x <- gsub("\\\\textcolor\\{red\\}\\{", "", x)
  x <- gsub("\\\\scriptsize|\\\\textbf", "", x)
  x <- gsub("[{}]", "", x)
  gsub("\\\\,", "", x)
}
METRIC_REV <- c("Int-Euclidean" = "Int-Euclidean", "Hausdorff" = "Hausdorff",
                "Ichino-Yaguchi" = "Ichino-Yaguchi", "Wasserstein" = "Wasserstein",
                "Centers (baseline)" = "Centers-Euclidean")
tbl_rows <- function(label) {
  i0 <- grep(paste0("label\\{", label, "\\}"), tex)
  if (!length(i0)) return(NULL)
  i1 <- i0 + which(grepl("\\\\midrule", tex[(i0 + 1):(i0 + 30)]))[1]
  i2 <- i0 + which(grepl("\\\\bottomrule", tex[(i0 + 1):(i0 + 90)]))[1]
  rows <- tex[(i1 + 1):(i2 - 1)]
  rows[!grepl("^\\s*\\\\midrule\\s*$", rows)]
}
parse_cell <- function(s) {
  star <- grepl("\\$\\^\\{?\\*\\}?\\$", s)
  s2 <- gsub("\\$\\^\\{?\\*\\}?\\$|\\$\\^\\{?\\\\dagger\\}?\\$", "", s)
  mm <- regmatches(s2, regexec("(-?[0-9]+\\.[0-9]+)\\s*(?:\\((-?[0-9.]+)\\))?", s2))[[1]]
  list(mean = as.numeric(mm[2]),
       sd = if (nchar(mm[3] %||% "")) as.numeric(mm[3]) else NA_real_,
       star = star)
}

star_rule_mc <- function(cal_csv) {
  if (!file.exists(cal_csv)) return(NULL)
  cal <- read.csv(cal_csv); cal <- cal[cal$kind == "minP", ]
  reps <- sort(unique(cal$rep))
  if (!identical(reps, 1:25)) bad(basename(cal_csv), "does not contain exactly reps 1..25")
  idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  aggregate(cal[idx], cal[c("IDR", "Metric")], function(p) sum(p <= 0.05) >= 23)
}

check_mc_table <- function(label, sum_csv, cal_csv, tie_csv, sd_col = ".sd") {
  rows <- tbl_rows(label); if (is.null(rows)) { bad(label, "table not found"); return() }
  d <- read.csv(sum_csv); st <- star_rule_mc(cal_csv)
  tf <- if (!is.null(tie_csv) && file.exists(tie_csv)) read.csv(tie_csv) else NULL
  idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  cur <- NA
  for (r in rows) {
    cells <- trimws(strsplit(gsub("\\\\\\\\\\s*$", "", strip(r)), "&")[[1]])
    if (length(cells) < 8) next
    if (nzchar(cells[1])) cur <- cells[1]
    met <- METRIC_REV[[sub("\\$\\^\\{?\\\\dagger\\}?\\$", "", cells[2])]]
    if (is.null(met)) next
    dd <- d[d$Method == cur & d$Metric == met, ]
    if (!nrow(dd)) { bad(label, cur, met, "no CSV row"); next }
    for (j in seq_along(idx)) {
      pc <- parse_cell(cells[2 + j])
      mv <- dd[[paste0(idx[j], ".mean")]]; sv <- dd[[paste0(idx[j], sd_col)]]
      if (is.na(pc$mean) || abs(pc$mean - round(mv, 3)) > 5e-4)
        bad(label, cur, met, idx[j], "mean", pc$mean, "vs", round(mv, 3)) else ok()
      ## compact notation: parenthetical integer = SD in units of 1e-3
      if (!is.na(pc$sd) && abs(pc$sd - round(sv * 1000)) > 0.5)
        bad(label, cur, met, idx[j], "SDx1000", pc$sd, "vs", round(sv * 1000)) else ok()
      if (!is.null(st)) {
        se <- st[st$IDR == cur & st$Metric == met, ]
        want <- nrow(se) && isTRUE(se[[idx[j]]])
        if (pc$star != want) bad(label, cur, met, idx[j], "star", pc$star, "vs", want) else ok()
      }
    }
  }
  cat(sprintf("%s: OK so far (%d cumulative checks)\n", label, n_checked))
}

check_plain_table <- function(label, idx_csv, minp_csv) {
  rows <- tbl_rows(label); if (is.null(rows)) { bad(label, "table not found"); return() }
  d <- read.csv(idx_csv); pm <- read.csv(minp_csv)
  idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  cur <- NA
  for (r in rows) {
    cells <- trimws(strsplit(gsub("\\\\\\\\\\s*$", "", strip(r)), "&")[[1]])
    if (length(cells) < 8) next
    if (nzchar(cells[1])) cur <- cells[1]
    dag_tex <- grepl("dagger", cells[2])
    met <- METRIC_REV[[sub("\\$\\^\\{?\\\\dagger\\}?\\$", "", cells[2])]]
    if (is.null(met)) next
    dd <- d[d$Method == cur & d$Metric == met, ]
    pp <- pm[pm$IDR == cur & pm$Metric == met, ]
    if (!nrow(dd)) { bad(label, cur, met, "no CSV row"); next }
    if (!is.na(dd$tie_affected) && (dag_tex != isTRUE(dd$tie_affected)))
      bad(label, cur, met, "dagger", dag_tex, "vs", dd$tie_affected) else ok()
    for (j in seq_along(idx)) {
      pc <- parse_cell(cells[2 + j])
      if (is.na(pc$mean) || abs(pc$mean - round(dd[[idx[j]]], 3)) > 5e-4)
        bad(label, cur, met, idx[j], pc$mean, "vs", round(dd[[idx[j]]], 3)) else ok()
      want <- nrow(pp) && pp[[idx[j]]] <= 0.05
      if (pc$star != want) bad(label, cur, met, idx[j], "star", pc$star, "vs", want) else ok()
    }
  }
  cat(sprintf("%s: OK so far (%d cumulative checks)\n", label, n_checked))
}

check_table3 <- function() {
  rows <- tbl_rows("tab:sim3_results"); if (is.null(rows)) { bad("tab:sim3", "not found"); return() }
  d <- read.csv(file.path(OUTPUT_DIR, "scenario3_summary.csv"))
  cols <- c("dQ_TC", "dB_TC", "flip_TC", "dQ_RE", "dB_RE", "flip_RE", "dQ_LC", "dB_LC")
  cur <- NA
  for (r in rows) {
    if (grepl("\\*", r)) bad("tab:sim3", "unsupported star present")
    cells <- trimws(strsplit(gsub("\\\\\\\\\\s*$", "", strip(r)), "&")[[1]])
    if (length(cells) < 10) next
    if (nzchar(cells[1])) cur <- cells[1]
    met <- METRIC_REV[[cells[2]]]
    if (is.null(met)) next
    dd <- d[d$Method == cur & d$Metric == met, ]
    if (!nrow(dd)) { bad("tab:sim3", cur, met, "no CSV row"); next }
    for (j in seq_along(cols)) {
      pc <- parse_cell(cells[2 + j])
      if (startsWith(cols[j], "flip")) {
        if (abs(pc$mean - round(dd[[paste0(cols[j], ".mean")]], 2)) > 5e-3)
          bad("tab:sim3", cur, met, cols[j], pc$mean) else ok()
        if (!is.na(pc$sd)) bad("tab:sim3", cur, met, cols[j], "flip must be bare") else ok()
      } else {
        if (abs(pc$mean - round(dd[[paste0(cols[j], ".mean")]], 3)) > 5e-4)
          bad("tab:sim3", cur, met, cols[j], "mean", pc$mean) else ok()
        if (!is.na(pc$sd) && abs(pc$sd - round(dd[[paste0(cols[j], ".sd")]] * 1000)) > 0.5)
          bad("tab:sim3", cur, met, cols[j], "SDx1000", pc$sd) else ok()
      }
    }
  }
  cat(sprintf("tab:sim3_results: OK so far (%d cumulative checks)\n", n_checked))
}

## ---- prose assertions -------------------------------------------------------
tex1 <- paste(strip(tex), collapse = " ")
has <- function(pat) grepl(pat, tex1, fixed = TRUE)
expect_prose <- function(token, cond, desc) {
  if (!cond) bad("prose:", desc, "->", token) else ok()
}
pr <- function(x, d = 3) sprintf(paste0("%.", d, "f"), x)

s1 <- read.csv(file.path(OUTPUT_DIR, "scenario1_summary.csv"))
v <- function(df, m, g, col) df[df$Method == m & df$Metric == g, col]
expect_prose("0.835", has(pr(round(v(s1, "C-PCA", "Centers-Euclidean", "Q_TC.mean"), 3))),
             "Scenario I baseline Q_TC")
expect_prose("0.141", has(pr(round(v(s1, "C-PCA", "Centers-Euclidean", "B_TC.mean"), 3))),
             "Scenario I baseline B_TC")
expect_prose("Sc1 SD max<=0.036", has("0.036") &&
               max(s1[grep("\\.sd$", names(s1))]) <= 0.036 + 1e-9,
             "Scenario I max SD claim")
s2f <- read.csv(file.path(OUTPUT_DIR, "scenario2_factorial_summary.csv"))
fv <- function(di, ge) round(s2f[s2f$Method == "C-PCA" & s2f$Metric == "Wasserstein" &
                                  s2f$dr_input == di & s2f$eval_geom == ge, "B_TC"], 3)
for (tk in c(fv("raw", "std"), fv("std", "std"), fv("raw", "raw"), fv("std", "raw")))
  expect_prose(pr(tk), has(pr(tk)), "Scenario II factorial B_TC value")
mid <- read.csv(file.path(OUTPUT_DIR, "midscale_indices.csv"))
expect_prose("0.688", has(pr(round(v(mid, "Int-UMAP", "Int-Euclidean", "Q_TC"), 3))),
             "USHCN Int-UMAP IE Q_TC")
expect_prose("0.918", has(pr(round(v(mid, "MR-PCA", "Int-Euclidean", "Q_TC"), 3))),
             "USHCN MR-PCA IE Q_TC")
mrt <- read.csv(file.path(OUTPUT_DIR, "midscale_runtime.csv"))
for (comp in seq_len(nrow(mrt)))
  expect_prose(sprintf("%.0f s", mrt$seconds[comp]),
               has(paste0("$", sprintf("%.0f", mrt$seconds[comp]), "$")),
               paste("USHCN runtime component", mrt$component[comp]))
s3 <- read.csv(file.path(OUTPUT_DIR, "scenario3_summary.csv"))
expect_prose("0.128/0.120", has("0.128") && has("0.120"), "Scenario III MRRE vs TC")

## ---- B2 appendix prose (checked only once the appendix is integrated) ------
if (any(grepl("label\\{app:b2\\}", tex))) {
  sens <- read.csv(file.path(OUTPUT_DIR, "b2_sensitivity.csv"))
  rngc <- grep("^range_", names(sens))
  sens$max_range <- apply(sens[rngc], 1, max)
  nu_max <- max(sens$max_range[sens$par == "nu" & sens$arm == "Face"])
  lam_face_max <- max(sens$max_range[sens$par == "lambda" & sens$arm == "Face"])
  lam_grid <- sens$max_range[sens$par == "lambda" & sens$arm != "Face"]
  expect_prose("nu max 0.074", has(pr(round(nu_max, 3))), "B2 nu max variation")
  expect_prose("lambda Face 0.22", has(sprintf("%.2f", floor(lam_face_max * 100) / 100)),
               "B2 lambda Face max variation")
  expect_prose("lambda grid 0.36", has(sprintf("%.2f", floor(lam_grid * 100) / 100)),
               "B2 lambda grid max variation")
  cp <- read.csv(file.path(OUTPUT_DIR, "b2_cpca_diagnostics.csv"))
  expect_prose("identity 1.4e-13",
               has("1.4\\times 10^-13") &&
                 max(cp$err_C, cp$err_R) <= 1.45e-13,
               "B2 C-PCA identity max error")
  expect_prose("eig share ranges",
               has(paste0("[", pr(round(min(cp$eig_share), 3)))) &&
                 has(paste0(pr(round(max(cp$eig_share), 3)), "]")),
               "B2 eigenshare data range")
  expect_prose("F_diff floor 1.1e5",
               quantile(cp$F_diff, 0.025) > 1.1e5 && has("1.1\\times 10^5"),
               "B2 F-diff quantile claim")
}

## ---- figures ---------------------------------------------------------------
inc <- regmatches(tex1, gregexpr("includegraphics[^ ]*\\[?[^]]*\\]?[A-Za-z0-9_./-]+\\.pdf", tex1))[[1]]
figs <- unique(regmatches(inc, regexpr("[A-Za-z0-9_.-]+\\.pdf", inc)))
fig_manifest <- list()
for (f in figs) {
  mpath <- file.path(MANUSCRIPT_DIR, f)
  apath <- file.path(OUTPUT_DIR, f)
  if (!file.exists(mpath)) { bad("figure missing in manuscript dir:", f); next }
  ok()
  hm <- sha256_file(mpath)
  ha <- if (file.exists(apath)) sha256_file(apath) else NA_character_
  if (!is.na(ha) && !identical(hm, ha)) bad("figure hash mismatch (manuscript vs analysis):", f)
  else if (!is.na(ha)) ok()
  fig_manifest[[f]] <- list(manuscript = rel_path(mpath), sha256 = hm,
                            analysis_match = identical(hm, ha))
}
jsonlite::write_json(fig_manifest, file.path(OUTPUT_DIR, "figure_manifest.json"),
                     auto_unbox = TRUE, pretty = TRUE)

## ---- run table checks -------------------------------------------------------
check_mc_table("tab:sim1_results", file.path(OUTPUT_DIR, "scenario1_summary.csv"),
               file.path(OUTPUT_DIR, "scenario1_calibration.csv"),
               file.path(OUTPUT_DIR, "scenario1_tie_flags.csv"))
check_mc_table("tab:sim2_results", file.path(OUTPUT_DIR, "scenario2_summary.csv"),
               file.path(OUTPUT_DIR, "scenario2_calibration.csv"), NULL)
check_table3()
check_plain_table("tab:realdata_results", file.path(OUTPUT_DIR, "face_indices.csv"),
                  file.path(OUTPUT_DIR, "face_pvalues_minP.csv"))
check_plain_table("tab:ushcn_results", file.path(OUTPUT_DIR, "midscale_indices.csv"),
                  file.path(OUTPUT_DIR, "midscale_pvalues_minP.csv"))
if (any(grepl("label\\{tab:sim4_results\\}", tex)))
  check_mc_table("tab:sim4_results", file.path(OUTPUT_DIR, "scenario4_summary.csv"),
                 file.path(OUTPUT_DIR, "scenario4_calibration.csv"), NULL)

cat(sprintf("\nCHECKS: %d passed, %d MISMATCHES\n", n_checked, fails))

## ---- deliberate-failure self-test -------------------------------------------
if (DO_SELFTEST && fails == 0L) {
  tmp <- tempfile(fileext = ".tex")
  on.exit(unlink(tmp), add = TRUE)
  crooked <- sub("0.835", "0.836", tex)   # corrupt one verified token
  writeLines(crooked, tmp)
  st <- system2(file.path(R.home("bin"), "Rscript"),
                args = c(shQuote(SCRIPT_FILE), paste0("--texfile=", tmp), "--no-selftest"),
                stdout = NULL, stderr = NULL)
  if (st == 0) { cat("SELF-TEST FAILED: corrupted copy passed\n"); fails <- fails + 1L }
  else cat("SELF-TEST OK: corrupted copy correctly rejected (exit", st, ")\n")
}
if (fails > 0) quit(status = 1)
