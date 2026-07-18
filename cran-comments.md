## Test environments

* Windows 11 x64, R 4.6.1
* GitHub Actions: R release on macOS, Ubuntu, and Windows; R devel and oldrel
  on Ubuntu (all passed for release commit fcdcef1)

## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new CRAN submission.
* HTML validation was skipped because the local Windows check environment
  does not provide the external HTML Tidy executable. Vignettes and both the
  PDF and HTML manuals were otherwise built successfully.

## Downstream dependencies

There are currently no downstream CRAN dependencies.

## Correctness note

Version 0.2.0 corrects the MRRE and LCMC-behavior formulas in the earlier
GitHub-only 0.1.0 release. The correction is documented prominently in
`NEWS.md`; all manuscript analyses were regenerated with 0.2.0.
