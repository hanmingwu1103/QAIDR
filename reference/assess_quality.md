# Assess quality of interval DR projections

Computes quality and behavior indices for each combination of DR method
and interval distance metric. Optionally performs permutation tests for
statistical significance.

## Usage

``` r
assess_quality(
  x,
  projections,
  K = 5,
  metrics = c("Int-Euclidean", "Hausdorff", "Ichino-Yaguchi", "Wasserstein"),
  perm_test = FALSE,
  n_perm = 1000
)
```

## Arguments

- x:

  An `interval_data` object (standardized).

- projections:

  An `idr_projections` object from
  [`run_idr()`](https://hanmingwu1103.github.io/QAIDR/reference/run_idr.md),
  or a named list with the same structure.

- K:

  Integer neighborhood size (default 5).

- metrics:

  Character vector of distance metrics to evaluate.

- perm_test:

  Logical; perform permutation tests (default `FALSE`).

- n_perm:

  Integer number of permutations (default 1000).

## Value

An object of class `qaidr_assessment` containing:

- results:

  Data frame with columns IDR, Metric, Q_TC, B_TC, Q_RE, B_RE, Q_LC,
  B_LC.

- pvalues:

  Data frame of p-values (if `perm_test = TRUE`).

- K:

  The neighborhood size used.

## Examples

``` r
if (FALSE) { # \dontrun{
data(cars_mm)
x <- standardize(cars_mm)
proj <- run_idr(x)
result <- assess_quality(x, proj, K = 5, perm_test = TRUE)
print(result)
} # }
```
