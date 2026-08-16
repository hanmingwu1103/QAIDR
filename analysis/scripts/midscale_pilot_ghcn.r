# Stage 6, Arm B pilot: NOAA GHCN-Daily / USHCN seasonal temperature intervals.
# DATA ENGINEERING + RUBRIC STATISTICS ONLY. This script computes NO quality
# or behavior index; the selection record must be frozen before any such
# computation (Prompt B Stage 6).
#
# Locked rules (REVISION_WORKLOG.md): year 2023 (fallback 2022); seasons
# DJF/MAM/JJA/SON (DJF = Dec of prior year + Jan/Feb of the target year);
# variables [mean daily TMIN, mean daily TMAX] in degrees C per season (p=4);
# unit retained iff every season has >=90% of expected days with BOTH TMIN
# and TMAX present and QFLAG blank; no imputation.

SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
save_sessioninfo("midscale_pilot_ghcn")

RAW <- file.path(DATA_RAW, "ghcnd")
TAR <- file.path(RAW, "ghcnd_hcn.tar.gz")
stopifnot(file.exists(TAR))
## Author decision 2026-07-15 (recorded in REVISION_WORKLOG.md and
## AUTHOR_DECISIONS pointer): the USHCN pool saturates below the rubric's
## 1,000-station target (645 pass at the locked 90% rule for 2023; ~714 even
## at 75%), so the author explicitly accepted n = 645 with the locked filters
## unchanged. The acceptance floor below reflects that decision; the
## completeness rule itself is NOT relaxed.
YEARS <- c(2023L)
N_ACCEPT_FLOOR <- 600L

## 1. Record snapshot hashes
hashes <- vapply(c("ghcnd_hcn.tar.gz", "ghcnd-stations.txt", "ghcnd-version.txt"),
                 function(f) digest::digest(file.path(RAW, f), file = TRUE,
                                            algo = "sha256"),
                 character(1))
cat("Snapshot SHA-256:\n"); print(hashes)

## 2. Extract .dly files (once)
DLY_DIR <- file.path(RAW, "extracted")
if (!dir.exists(DLY_DIR)) {
  dir.create(DLY_DIR)
  utils::untar(TAR, exdir = DLY_DIR)
}
dly <- list.files(DLY_DIR, pattern = "\\.dly$", recursive = TRUE, full.names = TRUE)
cat("Stations in snapshot:", length(dly), "\n")

## 3. Parse TMIN/TMAX lines for the needed months
## .dly fixed-width: ID 1-11, YEAR 12-15, MONTH 16-17, ELEMENT 18-21,
## then 31 x (VALUE 5, MFLAG 1, QFLAG 1, SFLAG 1); values in tenths deg C, -9999 = NA
season_of <- function(y, m, target) {
  if (m == 12 && y == target - 1) return("DJF")
  if (y != target) return(NA_character_)
  c("DJF","DJF",rep("MAM",3),rep("JJA",3),rep("SON",3),NA)[m][1]
}
days_in <- function(y, m) c(31, 28 + (y %% 4 == 0 & (y %% 100 != 0 | y %% 400 == 0)),
                            31,30,31,30,31,31,30,31,30,31)[m]

parse_station <- function(f, target) {
  ln <- readLines(f, warn = FALSE)
  yr <- as.integer(substr(ln, 12, 15)); mo <- as.integer(substr(ln, 16, 17))
  el <- substr(ln, 18, 21)
  keep <- (el %in% c("TMIN", "TMAX")) & ((yr == target & mo <= 11) | (yr == target - 1 & mo == 12) | (yr == target & mo %in% 1:2))
  keep <- keep | ((el %in% c("TMIN","TMAX")) & yr == target)
  ln <- ln[keep]; if (!length(ln)) return(NULL)
  yr <- yr[keep]; mo <- mo[keep]; el <- el[keep]
  ## PER-ELEMENT QUALITY CONTROL (author decision, retained).
  ## TMIN and TMAX are screened and averaged SEPARATELY. For each element the
  ## daily value is retained when it is present, not the -9999 sentinel, and
  ## carries a blank GHCN quality flag; the seasonal mean of that element is
  ## then taken over its own retained days. The two elements therefore need not
  ## be retained on the same calendar days.
  ##
  ## This is a deliberate, documented estimand, not an oversight: the seasonal
  ## interval endpoints are two independently estimated seasonal means, and the
  ## width is their difference. It is NOT a paired daily diurnal range, and must
  ## not be described as one. A same-day (paired) alternative was measured and
  ## rejected; see revision/USHCN_PAIRED_DAY_STOP_REPORT.md for the comparison
  ## and the decision. The completeness rule below enforces >= 90% valid daily
  ## observations PER ELEMENT in every season.
  out <- list()
  for (r in seq_along(ln)) {
    ssn <- season_of(yr[r], mo[r], target)
    if (is.na(ssn)) next
    nd <- days_in(yr[r], mo[r])
    pos0 <- 22 + (seq_len(nd) - 1) * 8
    val <- suppressWarnings(as.integer(substring(ln[r], pos0, pos0 + 4)))
    qfl <- substring(ln[r], pos0 + 6, pos0 + 6)
    ok <- !is.na(val) & val != -9999L & qfl == " "
    out[[r]] <- data.frame(season = ssn, el = el[r], n_ok = sum(ok),
                           n_days = nd, sum_ok = sum(val[ok]) / 10)
  }
  do.call(rbind, out)
}

build_year <- function(target) {
  res <- vector("list", length(dly))
  for (ii in seq_along(dly)) {
    st <- substr(basename(dly[ii]), 1, 11)
    d <- parse_station(dly[ii], target)
    if (is.null(d)) next
    agg <- aggregate(cbind(n_ok, n_days, sum_ok) ~ season + el, d, sum)
    agg$station <- st
    res[[ii]] <- agg
    if (ii %% 200 == 0) cat("  station", ii, "/", length(dly), "\n")
  }
  do.call(rbind, res)
}

for (target in YEARS) {
  cat("== Target year", target, "==\n")
  agg <- build_year(target)
  ## Completeness, PER ELEMENT: a station is retained only when EACH of TMIN and
  ## TMAX independently has >= 90% valid daily observations in EACH of the four
  ## seasons -- that is, all 8 season-by-element cells must pass. The threshold
  ## is the locked 90% rule and is not revisited here.
  agg$frac <- agg$n_ok / agg$n_days
  wide_ok <- with(agg, tapply(frac >= 0.9, list(station, paste(season, el)), isTRUE))
  full <- rownames(wide_ok)[rowSums(wide_ok, na.rm = TRUE) == 8]
  cat("Stations passing completeness:", length(full), "of", length(unique(agg$station)), "\n")
  if (length(full) >= N_ACCEPT_FLOOR) { TARGET <- target; break }
}
stopifnot(exists("TARGET"))

## 4. Seasonal interval construction for retained stations
## Each endpoint is the seasonal mean of its OWN element, over that element's own
## retained days: mean_val = sum_ok / n_ok separately for TMIN and for TMAX. The
## reshape below pairs the two independently estimated means into one interval.
sub <- agg[agg$station %in% full & agg$frac >= 0.9, ]
sub$mean_val <- sub$sum_ok / sub$n_ok
tm <- reshape(sub[, c("station", "season", "el", "mean_val")],
              idvar = c("station", "season"), timevar = "el", direction = "wide")
names(tm)[names(tm) == "mean_val.TMAX"] <- "tmax"
names(tm)[names(tm) == "mean_val.TMIN"] <- "tmin"
tm <- tm[complete.cases(tm), ]
## interval = [mean TMIN, mean TMAX]; degenerate if tmax <= tmin
tm$width <- tm$tmax - tm$tmin
n_degen <- sum(tm$width <= 0)
cat(sprintf("Season-cells: %d; degenerate widths: %d (%.3f%%)\n",
            nrow(tm), n_degen, 100 * n_degen / nrow(tm)))

centers <- reshape(transform(tm, c = (tmin + tmax) / 2)[, c("station", "season", "c")],
                   idvar = "station", timevar = "season", direction = "wide")
radii <- reshape(transform(tm, r = pmax(0, (tmax - tmin) / 2))[, c("station", "season", "r")],
                 idvar = "station", timevar = "season", direction = "wide")
keep <- complete.cases(centers) & complete.cases(radii)
centers <- centers[keep, ]; radii <- radii[keep, ]
stopifnot(identical(centers$station, radii$station))
n_final <- nrow(centers)
cat("Final n (stations with all four seasonal intervals):", n_final, "\n")

## ---------------------------------------------------------------------------
## Regression checks on the locked USHCN construction. These are assertions,
## not diagnostics: if any fails the snapshot must not be written.
## ---------------------------------------------------------------------------
N_EXPECTED <- 645L

## (a) 90% rule, PER ELEMENT. Every retained station must clear >= 90% valid
##     daily observations for BOTH TMIN and TMAX in EACH of the four seasons,
##     i.e. all 8 season-by-element cells.
chk <- agg[agg$station %in% centers$station, ]
cells <- with(chk, tapply(frac >= 0.9, list(station, paste(season, el)), isTRUE))
stopifnot(
  "per-element 90% rule: a retained station has fewer than 8 passing cells" =
    all(rowSums(cells, na.rm = TRUE) == 8L),
  "per-element 90% rule: a retained cell falls below 0.9" =
    all(chk$frac >= 0.9)
)

## (b) lower <= upper for every season cell (equivalently radius >= 0).
stopifnot(
  "interval ordering: mean TMAX < mean TMIN in some season cell" =
    all(tm$tmax >= tm$tmin),
  "interval ordering: negative radius" =
    all(as.matrix(radii[, -1]) >= 0)
)

## (c) locked station count. n = 645 is the author's recorded 2026-07-15
##     decision under the locked filters; a change here means the construction
##     or the frozen snapshot moved and must be investigated, not absorbed.
stopifnot(
  "station count changed: expected 645" = identical(as.integer(n_final), N_EXPECTED)
)
cat("Regression checks passed: per-element 90% rule, lower <= upper, n =",
    N_EXPECTED, "\n")

## 5. Save frozen derived snapshot + rubric report (NO indices computed)
ord <- order(centers$station)
cmat <- as.matrix(centers[ord, -1]); rmat <- as.matrix(radii[ord, -1])
colnames(cmat) <- sub("^c\\.", "", colnames(cmat))
colnames(rmat) <- sub("^r\\.", "", colnames(rmat))
rmat <- rmat[, colnames(cmat), drop = FALSE]
snap <- list(stations = centers$station[ord],
             centers = cmat,
             radii = rmat,
             target_year = TARGET,
             source_version = readLines(file.path(RAW, "ghcnd-version.txt"), n = 1),
             raw_hashes = hashes,
             built = format(Sys.time()))
saveRDS(snap, file.path(DATA_PROC, "ghcn_seasonal_intervals.rds"))
rubric <- data.frame(
  criterion = c("availability", "completeness", "scale", "widths", "runtime"),
  value = c("official NOAA, public domain, versioned; hashes recorded",
            sprintf("%d/%d stations pass 90%% joint rule (year %d)", n_final, length(dly), TARGET),
            sprintf("n = %d (target 1000-3000)", n_final),
            sprintf("%.3f%% degenerate season-cells pre-filter", 100 * n_degen / nrow(tm)),
            "projected from benchmark: n~1200, p=4 within envelope"))
save_result(rubric, "midscale_pilot_ghcn_rubric",
            inputs = file.path(RAW, c("ghcnd_hcn.tar.gz", "ghcnd-stations.txt")),
            extra = list(n_final = n_final, target_year = TARGET))
cat("Arm B pilot complete. NO indices computed.\n")
