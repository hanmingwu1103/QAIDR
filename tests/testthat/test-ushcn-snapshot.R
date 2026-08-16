# Regression checks for the locked USHCN seasonal-interval construction.
#
# The USHCN endpoints are two INDEPENDENTLY estimated seasonal means: TMIN and
# TMAX are quality-screened and averaged separately, each over its own retained
# days, and the interval width is the difference between the separately
# estimated seasonal mean TMAX and mean TMIN. It is NOT a paired daily diurnal
# range. A same-day (paired) alternative was measured and rejected; see
# revision/USHCN_PAIRED_DAY_STOP_REPORT.md.
#
# The per-element 90% completeness rule itself is asserted inside
# analysis/scripts/midscale_pilot_ghcn.r, which has access to the daily counts.
# What is checkable from the frozen snapshot is the resulting invariant set:
# the station count, interval ordering, and structural consistency.

snapshot_path <- function() {
  # tests may run from tests/testthat, from the package root, or from a check dir
  cands <- c(
    file.path("..", "..", "analysis", "data-processed", "ghcn_seasonal_intervals.rds"),
    file.path("analysis", "data-processed", "ghcn_seasonal_intervals.rds"),
    file.path("..", "..", "..", "analysis", "data-processed", "ghcn_seasonal_intervals.rds")
  )
  hit <- cands[file.exists(cands)]
  if (length(hit)) normalizePath(hit[1]) else NA_character_
}

test_that("USHCN frozen snapshot retains exactly 645 stations", {
  p <- snapshot_path()
  skip_if(is.na(p), "frozen USHCN snapshot not available in this checkout")
  z <- readRDS(p)

  expect_true(all(c("stations", "centers", "radii") %in% names(z)))
  expect_identical(nrow(z$centers), 645L)
  expect_identical(nrow(z$radii), 645L)
  expect_identical(length(z$stations), 645L)
  expect_identical(nrow(z$centers), length(unique(z$stations)))
})

test_that("USHCN intervals satisfy lower <= upper in every season cell", {
  p <- snapshot_path()
  skip_if(is.na(p), "frozen USHCN snapshot not available in this checkout")
  z <- readRDS(p)

  # interval = [center - radius, center + radius]; lower <= upper iff radius >= 0
  expect_true(all(is.finite(z$radii)))
  expect_true(all(z$radii >= 0))

  lower <- z$centers - z$radii
  upper <- z$centers + z$radii
  expect_true(all(upper >= lower))
})

test_that("USHCN snapshot is structurally consistent (4 seasons, aligned)", {
  p <- snapshot_path()
  skip_if(is.na(p), "frozen USHCN snapshot not available in this checkout")
  z <- readRDS(p)

  expect_identical(ncol(z$centers), 4L)
  expect_identical(ncol(z$radii), 4L)
  expect_identical(dim(z$centers), dim(z$radii))
  expect_identical(colnames(z$centers), colnames(z$radii))
  expect_true(all(is.finite(z$centers)))
  expect_false(anyDuplicated(z$stations) > 0)
})
