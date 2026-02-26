# Compute quality/behavior index profiles over K

Computes co-ranking indices for all neighborhood sizes from 1 to
`K_max`, for each combination of DR method and distance metric.

## Usage

``` r
k_profiles(
  x,
  projections,
  K_max = NULL,
  metrics = c("Int-Euclidean", "Hausdorff", "Ichino-Yaguchi", "Wasserstein")
)
```

## Arguments

- x:

  An `interval_data` object (standardized).

- projections:

  An `idr_projections` object.

- K_max:

  Maximum neighborhood size (default `nrow(x$centers) - 2`).

- metrics:

  Character vector of distance metrics.

## Value

A data frame with columns Method, Metric, K, Q_TC, B_TC, Q_RE, B_RE,
Q_LC, B_LC.

## Examples

``` r
if (FALSE) { # \dontrun{
data(cars_mm)
x <- standardize(cars_mm)
proj <- run_idr(x)
profiles <- k_profiles(x, proj, K_max = 10)
} # }
```
