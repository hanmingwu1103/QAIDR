## Test environments

* Windows 10 x64 (build 19045), R 4.6.1
* GitHub Actions: R release on macOS, Ubuntu, and Windows; R devel and oldrel
  on Ubuntu (configured; the QAIDR 0.3.0 run will be recorded before CRAN
  submission)

## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new CRAN submission.
* The recorded QAIDR 0.3.0 local check used `--no-manual --ignore-vignettes`;
  the full GitHub Actions matrix and vignette/manual checks will be recorded
  before CRAN submission.

## Downstream dependencies

There are currently no downstream CRAN dependencies.

## Correctness note

Version 0.3.0 replaces the familywise min-P calculation used in 0.2.0 with a
fully symmetric randomization construction. Adjusted p-values from 0.2.0 must
be regenerated; marginal p-values and descriptive indices are unchanged by
this correction. The exact algorithm, scope of weak FWER control, validation
tests, and migration note are documented prominently in `NEWS.md`. All
affected manuscript inference outputs were regenerated with 0.3.0.
