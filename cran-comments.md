## Test environments

* Local: Windows 10 x64 (build 19045), R 4.6.0 (2026-04-24 ucrt),
  `R CMD check --as-cran` on the built source tarball, **including** the
  PDF manual and full vignette re-building.
* GitHub Actions `R-CMD-check` matrix (macOS release, Windows release,
  Ubuntu release, Ubuntu devel, Ubuntu oldrel-1). The most recent successful
  matrix run is 2026-07-20 on commit 6f22749c. That workflow builds with
  `--no-manual`, so the PDF manual is covered only by the local check above.

## R CMD check results

0 errors | 0 warnings | 2 notes

Both notes are expected:

* `checking CRAN incoming feasibility ... NOTE` — "New submission".
  This is the first CRAN submission of QAIDR.

* `checking HTML version of manual ... NOTE` — "Skipping checking HTML
  validation: no command 'tidy' found." HTML Tidy is not installed on the
  local machine. This reflects the local toolchain only, not the package;
  CRAN check machines provide `tidy`.

Package vignettes, re-building of vignette outputs, the PDF version of the
manual, and all `testthat` tests check OK.

## Downstream dependencies

There are currently no downstream CRAN dependencies.

## Correctness note

Version 0.3.0 replaces the familywise min-P calculation used in 0.2.0 with a
fully symmetric randomization construction. Adjusted p-values from 0.2.0 must
be regenerated; marginal p-values and descriptive indices are unchanged by
this correction. The exact algorithm, scope of weak FWER control, validation
tests, and migration note are documented prominently in `NEWS.md`. All
affected manuscript inference outputs were regenerated with 0.3.0.
