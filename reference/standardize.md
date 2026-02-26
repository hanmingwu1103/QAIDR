# Standardize interval data

Centers the midpoints and scales both midpoints and radii by the
midpoint standard deviation per variable. This is the standard symbolic
standardization procedure.

## Usage

``` r
standardize(x)
```

## Arguments

- x:

  An `interval_data` object.

## Value

A new `interval_data` object with standardized values.

## Examples

``` r
C <- matrix(rnorm(30), 10, 3)
R <- matrix(runif(30, 0.1, 1), 10, 3)
x <- interval_data(C, R)
xs <- standardize(x)
```
