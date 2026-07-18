# Ichino-Yaguchi dissimilarity for intervals

Computes the Ichino-Yaguchi dissimilarity between interval-valued
observations based on the interval join (hull) and meet (intersection).

## Usage

``` r
idist_ichino_yaguchi(centers, radii, nu = 0.5, gamma = NULL)
```

## Arguments

- centers:

  Numeric matrix of interval midpoints (n x p).

- radii:

  Numeric matrix of interval half-widths (n x p).

- nu:

  Span weight in the canonical Ichino-Yaguchi range \\\[0, 0.5\]\\
  (default 0.5); named \\\nu\\ in the manuscript.

- gamma:

  Deprecated alias for `nu`.

## Value

A symmetric n x n distance matrix.

## Examples

``` r
C <- matrix(rnorm(12), 4, 3)
R <- matrix(runif(12, 0.1, 0.5), 4, 3)
D <- idist_ichino_yaguchi(C, R)
```
