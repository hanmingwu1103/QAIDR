# Product-Hausdorff distance for intervals

Computes the coordinatewise product-Hausdorff dissimilarity between
interval-valued observations: the interval Hausdorff distance is taken
per coordinate and aggregated in \\\ell_2\\. This is the manuscript's
product-Hausdorff form; it is generally NOT the set-Hausdorff metric
between hyperrectangles under the Euclidean point metric. The public
metric key `"Hausdorff"` is unchanged.

## Usage

``` r
idist_hausdorff(centers, radii)
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
D <- idist_hausdorff(C, R)
```
