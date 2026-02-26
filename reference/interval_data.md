# Create an interval_data object

Constructs an `interval_data` object from center and radii matrices.

## Usage

``` r
interval_data(centers, radii, labels = NULL)
```

## Arguments

- centers:

  A numeric matrix of interval midpoints (n x p).

- radii:

  A numeric matrix of interval half-widths (n x p).

- labels:

  An optional factor or character vector of class labels (length n).

## Value

An object of class `interval_data`.

## Examples

``` r
C <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)
R <- matrix(c(0.1, 0.2, 0.3, 0.1, 0.2, 0.3), nrow = 3)
x <- interval_data(C, R)
```
