# Compute interval distance matrix

Dispatcher that computes a distance matrix using a named interval
metric.

## Usage

``` r
idist(centers, radii = NULL, metric = "Wasserstein")
```

## Arguments

- centers:

  Numeric matrix of interval midpoints (n x p), or an `interval_data`
  object.

- radii:

  Numeric matrix of interval half-widths (n x p). Ignored if `centers`
  is an `interval_data` object.

- metric:

  Character string: one of `"Int-Euclidean"`, `"Hausdorff"`,
  `"Ichino-Yaguchi"`, or `"Wasserstein"`.

## Value

A symmetric n x n distance matrix.

## Examples

``` r
C <- matrix(rnorm(12), 4, 3)
R <- matrix(runif(12, 0.1, 0.5), 4, 3)
D <- idist(C, R, "Wasserstein")
```
