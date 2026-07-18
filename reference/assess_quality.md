# Assess quality of interval DR projections

Computes the adjudicated quality and behavior indices for each
combination of DR method and interval dissimilarity, optionally with
permutation inference. P-values use the finite-sample Monte Carlo
estimator `(c + 1) / (m + 1)` throughout. When `perm_test = TRUE`,
family-wise adjusted p-values are additionally computed by single-step
Westfall-Young min-P resampling using the shared joint permutation
draws; the multiplicity family is, per DR method and index family (T&C,
MRRE, LCMC) at the fixed K, the set of {quality, behavior} tests across
all requested metrics. Quality tests are one-sided (upper); behavior
tests are two-sided.

## Usage

``` r
assess_quality(
  x,
  projections,
  K = 5,
  metrics = c("Int-Euclidean", "Hausdorff", "Ichino-Yaguchi", "Wasserstein"),
  lambda = 0.5,
  nu = 0.5,
  baseline = TRUE,
  perm_test = FALSE,
  n_perm = 999,
  ties = c("error", "random"),
  seed = NULL
)
```

## Arguments

- x:

  An `interval_data` object (standardized or raw; the caller controls
  preprocessing).

- projections:

  An `idr_projections` object from
  [`run_idr()`](https://hanmingwu1103.github.io/QAIDR/reference/run_idr.md).

- K:

  Integer neighborhood size (default 5).

- metrics:

  Character vector of interval dissimilarities.

- lambda:

  Optimism index for the Interval Euclidean score (default 0.5).

- nu:

  Ichino-Yaguchi span weight (default 0.5).

- baseline:

  Logical; add the center-only Euclidean evaluation row (default
  `TRUE`).

- perm_test:

  Logical; perform permutation inference (default `FALSE`).

- n_perm:

  Integer number of permutations (default 999).

- ties:

  Tie policy passed to
  [`rank_matrix()`](https://hanmingwu1103.github.io/QAIDR/reference/rank_matrix.md)
  (`"error"` for primary analyses).

- seed:

  Optional integer seed for the permutation stream.

## Value

A `qaidr_assessment` object: `results` (indices), `pvalues`
(unadjusted), `pvalues_adj` (min-P adjusted), `null_draws` (joint
permutation draws, for reproducible multiplicity adjustment), `K`,
`params`.

## Details

A center-only Euclidean baseline (`baseline = TRUE`) evaluates every
projection with plain Euclidean distances between interval centers in
both spaces, quantifying what interval-aware evaluation adds.
