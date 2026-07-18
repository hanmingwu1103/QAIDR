# Compute quality/behavior index profiles over K

Ranks each space once per metric and evaluates the indices for every K
on the same co-ranking matrix (the co-ranking matrix does not depend on
K), reducing the cost from the previous per-K re-ranking.

## Usage

``` r
k_profiles(
  x,
  projections,
  K_max = NULL,
  metrics = c("Int-Euclidean", "Hausdorff", "Ichino-Yaguchi", "Wasserstein"),
  lambda = 0.5,
  nu = 0.5,
  baseline = TRUE,
  ties = c("error", "random")
)
```

## Arguments

- x:

  An `interval_data` object (standardized or raw; the caller controls
  preprocessing).

- projections:

  An `idr_projections` object from
  [`run_idr()`](https://hanmingwu1103.github.io/QAIDR/reference/run_idr.md).

- K_max:

  Maximum neighborhood size (default `n - 2`).

- metrics:

  Character vector of interval dissimilarities.

- lambda:

  Optimism index for the Interval Euclidean score (default 0.5).

- nu:

  Ichino-Yaguchi span weight (default 0.5).

- baseline:

  Logical; add the center-only Euclidean evaluation row (default
  `TRUE`).

- ties:

  Tie policy passed to
  [`rank_matrix()`](https://hanmingwu1103.github.io/QAIDR/reference/rank_matrix.md)
  (`"error"` for primary analyses).

## Value

Data frame with columns Method, Metric, K and the six indices.
