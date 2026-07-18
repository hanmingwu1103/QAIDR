# Self-testing audit of the exceptional-value enumeration machinery
# (Prompt E 2026-07-17). Exercises the PRODUCTION code in b2_breaks_lib.R on
# the eight adjudicated cases: affine root; two distinct quadratic roots;
# exact double root (the Prompt D (1-5*nu)^2 tangency, both as raw
# coefficients and from real interval geometry); no real root; boundary
# root; identically-zero comparison; deduplication without raw replacement;
# union of original/embedded exceptional sets. Exits 0 iff all pass.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "b2_breaks_lib.R"))
suppressMessages(library(QAIDR))

coef_row <- function(space, id, q2, q1, q0) data.frame(
  space = space, anchor = id, j = 2L, jprime = 3L, q2 = q2, q1 = q1, q0 = q0,
  s2 = abs(q2), s1 = abs(q1), s0 = abs(q0))

affine  <- classify_polynomials(coef_row("high", 1, 0, 1, -0.25), c(0, 0.5))
two     <- classify_polynomials(coef_row("high", 2, 1, -0.5, 0.04), c(0, 0.5))
double  <- classify_polynomials(coef_row("high", 3, 25, -10, 1), c(0, 0.5))
none    <- classify_polynomials(coef_row("high", 4, 1, 0, 1), c(0, 0.5))
bound   <- classify_polynomials(coef_row("high", 5, 1, -0.2, 0), c(0, 0.5))
persist <- classify_polynomials(coef_row("high", 6, 0, 0, 0), c(0, 0.5))
dedup   <- classify_polynomials(rbind(
  coef_row("high", 7, 0, 1, -0.2),
  coef_row("high", 8, 0, 1, -(0.2 + 0.25 * TOL_ROOT))), c(0, 0.5))
hi <- classify_polynomials(rbind(coef_row("high",  9, 0, 1, -0.1),
                                 coef_row("high", 10, 0, 1, -0.3)), c(0, 0.5))
lo <- classify_polynomials(rbind(coef_row("low", 11, 0, 1, -0.1),
                                 coef_row("low", 12, 0, 1, -0.4)), c(0, 0.5))
both <- union_breaks(hi, lo)

stopifnot(
  all.equal(affine$internal, 0.25),
  max(abs(two$internal - c(0.1, 0.4))) < 1e-14,
  all.equal(double$tangency, 0.2),
  length(double$internal) == 1L,
  !length(none$internal),
  identical(bound$boundary, 0),
  all.equal(bound$internal, 0.2),
  nrow(persist$persistent) == 1L,
  length(dedup$internal) == 1L,
  dedup$internal_meta$n_source == 2L,
  dedup$internal_meta$raw_max > dedup$internal_meta$raw_min,
  max(abs(both$internal - c(0.1, 0.3, 0.4))) < 1e-14
)
cat("coefficient-level cases: PASS\n")

## Geometry-level tangency regression (Prompt D construction): the squared
## IY comparison difference at anchor 1 is (1 - 5 nu)^2; the enumeration
## must classify nu = 0.2 as one internal tangency, not miss it and not
## fabricate spurious ordinary roots from cancellation.
centers <- rbind(c(0, 0), c(6, sqrt(45) / 2), c(3 * sqrt(80) / 4, 0), c(50, 70))
radii   <- rbind(c(0, 0), c(0, sqrt(45) / 2), c(sqrt(80) / 4, 0), c(0, 0))
br <- nu_breaks(centers, radii, "high")
hit <- abs(br$tangency - 0.2) < 1e-12
stopifnot(any(hit))
near_spurious <- br$internal[abs(br$internal - 0.2) < 1e-6 & abs(br$internal - 0.2) > 1e-12]
stopifnot(length(near_spurious) == 0L)
D02 <- idist_ichino_yaguchi(centers, radii, nu = 0.2)
err <- tryCatch({ rank_matrix(D02, ties = "error"); FALSE }, error = function(e) TRUE)
stopifnot(err)
ord <- vapply(1:64, function(s) { set.seed(s)
  R <- rank_matrix(D02, ties = "random"); sign(R[1, 2] - R[1, 3]) }, numeric(1))
stopifnot(setequal(unique(ord), c(-1, 1)))
cat("geometry-level tangency regression: PASS\n")
cat("b2_breaks_selftest: ALL PASS\n")
