# Ichino-Yaguchi dissimilarity for intervals

Computes the Ichino-Yaguchi dissimilarity between interval-valued
observations based on intersection and union of intervals.

## Usage

``` r
idist_ichino_yaguchi(centers, radii, gamma = 0.5)
```

## Arguments

- centers:

  Numeric matrix of interval midpoints (n x p).

- radii:

  Numeric matrix of interval half-widths (n x p).

- gamma:

  Weighting parameter (default 0.5).

## Value

A symmetric n x n distance matrix.

## Examples

``` r
C <- matrix(rnorm(12), 4, 3)
R <- matrix(runif(12, 0.1, 0.5), 4, 3)
D <- idist_ichino_yaguchi(C, R)
```
