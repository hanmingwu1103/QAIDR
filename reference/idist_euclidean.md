# Interval Euclidean distance

Computes the interval Euclidean distance between interval-valued
observations, defined as a weighted combination of the maximum and
minimum distances between hyperrectangles.

## Usage

``` r
idist_euclidean(centers, radii, lambda = 0.5)
```

## Arguments

- centers:

  Numeric matrix of interval midpoints (n x p).

- radii:

  Numeric matrix of interval half-widths (n x p).

- lambda:

  Weight parameter in \\\[0, 1\]\\ (default 0.5).

## Value

A symmetric n x n distance matrix.

## Examples

``` r
C <- matrix(rnorm(12), 4, 3)
R <- matrix(runif(12, 0.1, 0.5), 4, 3)
D <- idist_euclidean(C, R)
```
