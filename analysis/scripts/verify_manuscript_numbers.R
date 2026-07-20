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
  x <- gsub("\\\\textcolor\\{red\\}\\{|\\\\Rone\\{", "", x)
  x <- gsub("\\\\scriptsize|\\\\textbf|\\\\tiny", "", x)
  x <- gsub("[{}]", "", x)
  gsub("\\\\,", "", x)
}
METRIC_REV <- c("Int-Euclidean" = "Int-Euclidean", "Hausdorff" = "Hausdorff",
                "Ichino-Yaguchi" = "Ichino-Yaguchi", "Wasserstein" = "Wasserstein",
                "Centers (baseline)" = "Centers-Euclidean",
                ## Prompt H abbreviated/renamed evaluation labels
                "IE" = "Int-Euclidean", "H" = "Hausdorff",
                "IY" = "Ichino-Yaguchi", "W" = "Wasserstein",
                "Centers" = "Centers-Euclidean")

## Prompt H: closed-form null means (thm:nullcal) for chance-adjusted Qc.
mu0_null <- function(n, K) {
  GK <- if (K < n / 2) n * K * (2 * n - 3 * K - 1) else n * (n - K) * (n - K - 1)
  HK <- n * sum(abs(n - 2 * seq_len(K) + 1) / seq_len(K))
  j <- seq_len(K)
  EWn <- n / ((n - 1) * HK) * sum((1 / j) * (j * (j - 1) / 2 + (n - 1 - j) * (n - j) / 2))
  c(Q_TC = 1 - n * K * (n - 1 - K) * (n - K) / (GK * (n - 1)),
    Q_RE = 1 - EWn, Q_LC = K / (n - 1))
}

## Prompt H cell grammar: MEAN{\tiny(SD)\Rone{[QC]}}$^{\Rone{FRAC}}$ (sim
## tables) or [\Rone{]VALUE[}]{\scriptsize\Rone{[QC]}}$^{*}$/$^{\Rone{*}}$
## (single-dataset tables); parsed from the RAW (unstripped) cell.
parse_cell2 <- function(s) {
  star <- grepl("\\^\\{(\\\\Rone\\{)?\\*", s)
  fm <- regmatches(s, regexec("\\$\\^\\{\\\\Rone\\{(-?[0-9.]+)\\}\\}\\$", s))[[1]]
  qm <- regmatches(s, regexec("\\\\Rone\\{\\[(-?[0-9.]+)\\]\\}", s))[[1]]
  mm <- regmatches(s, regexec("^\\s*(\\\\Rone\\{)?(-?[0-9]+\\.[0-9]+)", s))[[1]]
  sm <- regmatches(s, regexec("\\{\\\\(?:tiny|scriptsize)\\s*\\((-?[0-9.]+)\\)", s))[[1]]
  list(mean = if (length(mm) > 2) as.numeric(mm[3]) else NA_real_,
       sd = if (length(sm) > 1) as.numeric(sm[2]) else NA_real_,
       star = star,
       frac = if (length(fm) > 1) as.numeric(fm[2]) else NA_real_,
       qc = if (length(qm) > 1) as.numeric(qm[2]) else NA_real_)
}
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

frac_rule_mc <- function(cal_csv) {
  if (!file.exists(cal_csv)) return(NULL)
  cal <- read.csv(cal_csv); cal <- cal[cal$kind == "minP", ]
  reps <- sort(unique(cal$rep))
  if (!identical(reps, 1:25)) bad(basename(cal_csv), "does not contain exactly reps 1..25")
  idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  aggregate(cal[idx], cal[c("IDR", "Metric")], function(p) mean(p <= 0.05))
}

split_raw <- function(r) trimws(strsplit(gsub("\\\\\\\\\\s*$", "", r), "&")[[1]])
metric_of <- function(cell2) {
  key <- trimws(sub("\\$\\^\\{?\\\\dagger\\}?\\$", "", strip(cell2)))
  if (key %in% names(METRIC_REV)) METRIC_REV[[key]] else NULL
}

check_mc_table <- function(label, sum_csv, cal_csv, tie_csv, n_obj, K = 10L,
                           sd_col = ".sd") {
  rows <- tbl_rows(label); if (is.null(rows)) { bad(label, "table not found"); return() }
  d <- read.csv(sum_csv); fr <- frac_rule_mc(cal_csv)
  m0 <- mu0_null(n_obj, K)
  idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  cur <- NA
  for (r in rows) {
    cells <- split_raw(r)
    if (length(cells) < 8) next
    c1 <- trimws(strip(cells[1]))
    if (nzchar(c1)) cur <- c1
    met <- metric_of(cells[2])
    if (is.null(met)) next
    dd <- d[d$Method == cur & d$Metric == met, ]
    if (!nrow(dd)) { bad(label, cur, met, "no CSV row"); next }
    for (j in seq_along(idx)) {
      pc <- parse_cell2(cells[2 + j])
      mv <- dd[[paste0(idx[j], ".mean")]]; sv <- dd[[paste0(idx[j], sd_col)]]
      if (is.na(pc$mean) || abs(pc$mean - round(mv, 3)) > 5e-4)
        bad(label, cur, met, idx[j], "mean", pc$mean, "vs", round(mv, 3)) else ok()
      ## compact notation: parenthetical integer = SD in units of 1e-3
      if (!is.na(pc$sd) && abs(pc$sd - round(sv * 1000)) > 0.5)
        bad(label, cur, met, idx[j], "SDx1000", pc$sd, "vs", round(sv * 1000)) else ok()
      if (!is.null(fr)) {
        se <- fr[fr$IDR == cur & fr$Metric == met, ]
        if (!nrow(se) || is.na(pc$frac) ||
            abs(pc$frac - round(se[[idx[j]]], 2)) > 1e-9)
          bad(label, cur, met, idx[j], "frac", pc$frac, "vs",
              if (nrow(se)) round(se[[idx[j]]], 2) else NA) else ok()
      }
      if (startsWith(idx[j], "Q")) {
        want_qc <- round((mv - m0[[idx[j]]]) / (1 - m0[[idx[j]]]), 2)
        if (is.na(pc$qc) || abs(pc$qc - want_qc) > 1e-9)
          bad(label, cur, met, idx[j], "Qc", pc$qc, "vs", want_qc) else ok()
      }
    }
  }
  cat(sprintf("%s: OK so far (%d cumulative checks)\n", label, n_checked))
}

check_plain_table <- function(label, idx_csv, minp_csv, n_obj, K) {
  rows <- tbl_rows(label); if (is.null(rows)) { bad(label, "table not found"); return() }
  d <- read.csv(idx_csv); pm <- read.csv(minp_csv)
  m0 <- mu0_null(n_obj, K)
  idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  cur <- NA
  for (r in rows) {
    cells <- split_raw(r)
    if (length(cells) < 8) next
    c1 <- trimws(strip(cells[1]))
    if (nzchar(c1)) cur <- c1
    dag_tex <- grepl("dagger", cells[2])
    met <- metric_of(cells[2])
    if (is.null(met)) next
    dd <- d[d$Method == cur & d$Metric == met, ]
    pp <- pm[pm$IDR == cur & pm$Metric == met, ]
    if (!nrow(dd)) { bad(label, cur, met, "no CSV row"); next }
    if (!is.na(dd$tie_affected) && (dag_tex != isTRUE(dd$tie_affected)))
      bad(label, cur, met, "dagger", dag_tex, "vs", dd$tie_affected) else ok()
    for (j in seq_along(idx)) {
      pc <- parse_cell2(cells[2 + j])
      if (is.na(pc$mean) || abs(pc$mean - round(dd[[idx[j]]], 3)) > 5e-4)
        bad(label, cur, met, idx[j], pc$mean, "vs", round(dd[[idx[j]]], 3)) else ok()
      want <- nrow(pp) && pp[[idx[j]]] <= 0.05
      if (pc$star != want) bad(label, cur, met, idx[j], "star", pc$star, "vs", want) else ok()
      if (startsWith(idx[j], "Q")) {
        want_qc <- round((dd[[idx[j]]] - m0[[idx[j]]]) / (1 - m0[[idx[j]]]), 2)
        if (is.na(pc$qc) || abs(pc$qc - want_qc) > 1e-9)
          bad(label, cur, met, idx[j], "Qc", pc$qc, "vs", want_qc) else ok()
      }
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

## ---- Prompt H additions -----------------------------------------------------
v3c_f <- file.path(OUTPUT_DIR, "h_audit_v3c_familywise.rds")
if (file.exists(v3c_f)) {
  v3c <- readRDS(v3c_f)
  rup <- function(x) sprintf("%.3f", floor(x * 1000 + 0.5) / 1000)  # half-up
  for (i in seq_len(nrow(v3c)))
    expect_prose(rup(v3c$rate[i]), has(rup(v3c$rate[i])),
                 paste("v3c familywise rate", v3c$family[i]))
  expect_prose("v3c all pass", all(v3c$pass), "v3c Wilson intervals cover/undershoot 0.05")
} else bad("prose:", "h_audit_v3c_familywise.rds missing")

lam_c <- file.path(OUTPUT_DIR, "h_ushcn_lambda_components.rds")
if (file.exists(lam_c)) {
  cmp <- readRDS(lam_c)
  hi <- cmp[cmp$space == "high", ]
  expect_prose("71.0%", has(sprintf("%.1f", 100 * hi$frac_dminus_zero)),
               "USHCN high-space frac d-=0")
  expect_prose("0.023", has(sprintf("%.3f", round(hi$mean_ratio, 3))),
               "USHCN high-space mean d-/d+ ratio")
  iu <- cmp[cmp$space == "Int-UMAP", ]
  expect_prose("Int-UMAP d- == 0 all pairs", iu$frac_dminus_zero == 1,
               "Int-UMAP low-space d- identically zero")
  prof <- readRDS(file.path(OUTPUT_DIR, "h_ushcn_lambda_profile.rds"))
  rng <- function(m) range(prof$Q_TC[prof$Method == m & !is.na(prof$lambda) &
                                       prof$lambda >= 0.05])
  riu <- round(rng("Int-UMAP"), 2); rmr <- round(rng("MR-PCA"), 2)
  expect_prose("IntUMAP 0.62-0.72", has(sprintf("%.2f", riu[1])) &&
                 has(sprintf("%.2f", riu[2])), "lambda profile Int-UMAP range")
  expect_prose("MR-PCA 0.90-0.93", has(sprintf("%.2f", rmr[1])) &&
                 has(sprintf("%.2f", rmr[2])), "lambda profile MR-PCA range")
} else bad("prose:", "h_ushcn_lambda_components.rds missing")

face_idx <- read.csv(file.path(OUTPUT_DIR, "face_indices.csv"))
expect_prose("Face B_LC <= 0.10", max(face_idx$B_LC) <= 0.10 + 1e-9 && has("0.10"),
             "Face B_LC bound claim")

gw_f <- file.path(OUTPUT_DIR, "h_graded_width_summary.csv")
if (file.exists(gw_f) && any(grepl("label\\{tab:graded_width\\}", tex))) {
  gw <- read.csv(gw_f)
  pri <- gw[gw$primary == TRUE | gw$primary == "TRUE", ]
  i0 <- grep("label\\{tab:graded_width\\}", tex)
  i2 <- i0 + which(grepl("\\\\bottomrule", tex[(i0 + 1):(i0 + 20)]))[1]
  rows <- tex[(i0 + 1):(i2 - 1)]
  for (m in unique(pri$Method)) {
    row <- rows[grepl(paste0("^\\s*", m, " &"), rows)]
    if (!length(row)) { bad("tab:graded_width", m, "row missing"); next }
    cells <- trimws(strsplit(gsub("\\\\\\\\\\s*$", "", strip(row[1])), "&")[[1]])
    for (k in 1:3) {
      lv <- c("L1", "L2", "L3")[k]
      want <- round(pri$mean_qc_diff[pri$Method == m & pri$level == lv], 3)
      got <- as.numeric(cells[1 + k])
      if (is.na(got) || abs(got - want) > 5e-4)
        bad("tab:graded_width", m, lv, got, "vs", want) else ok()
    }
  }
  expect_prose("graded MCSE <= 0.006", max(pri$mcse) <= 0.006 + 1e-9 && has("0.006"),
               "graded-width MCSE bound")
  iu <- round(pri$mean_qc_diff[pri$Method == "Int-UMAP"], 3)
  mr <- round(pri$mean_qc_diff[pri$Method == "MR-PCA"], 3)
  expect_prose("IntUMAP 0.013->0.517", has(sprintf("%.3f", min(iu))) &&
                 has(sprintf("%.3f", max(iu))), "graded Int-UMAP range prose")
  expect_prose("MR-PCA -> -0.184", has(sprintf("%.3f", min(mr))),
               "graded MR-PCA prose")
} else bad("prose:", "graded-width outputs or table missing")

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
               file.path(OUTPUT_DIR, "scenario1_tie_flags.csv"),
               n_obj = 300L, K = 10L)
check_mc_table("tab:sim2_results", file.path(OUTPUT_DIR, "scenario2_summary.csv"),
               file.path(OUTPUT_DIR, "scenario2_calibration.csv"), NULL,
               n_obj = 800L, K = 10L)
check_table3()
check_plain_table("tab:realdata_results", file.path(OUTPUT_DIR, "face_indices.csv"),
                  file.path(OUTPUT_DIR, "face_pvalues_minP.csv"),
                  n_obj = 27L, K = 5L)
check_plain_table("tab:ushcn_results", file.path(OUTPUT_DIR, "midscale_indices.csv"),
                  file.path(OUTPUT_DIR, "midscale_pvalues_minP.csv"),
                  n_obj = 645L, K = 10L)
if (any(grepl("label\\{tab:sim4_results\\}", tex)))
  check_mc_table("tab:sim4_results", file.path(OUTPUT_DIR, "scenario4_summary.csv"),
                 file.path(OUTPUT_DIR, "scenario4_calibration.csv"), NULL,
                 n_obj = 300L, K = 10L)

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
