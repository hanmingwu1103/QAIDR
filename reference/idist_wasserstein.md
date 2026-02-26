# L2-Wasserstein (Mallows) distance for intervals

Computes the L2-Wasserstein distance between interval-valued
observations, assuming uniform distributions within each interval.

## Usage

``` r
idist_wasserstein(centers, radii)
```

## Arguments

- centers:

  Numeric matrix of interval midpoints (n x p).

- radii:

  Numeric matrix of interval half-widths (n x p).

## Value

A symmetric n x n distance matrix.

## Examples

``` r
C <- matrix(rnorm(12), 4, 3)
R <- matrix(runif(12, 0.1, 0.5), 4, 3)
D <- idist_wasserstein(C, R)
```
