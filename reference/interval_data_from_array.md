# Create interval_data from a 3D min-max array

Converts a 3D array (n x p x 2) where `[,,1]` is the lower bound and
`[,,2]` is the upper bound into an `interval_data` object.

## Usage

``` r
interval_data_from_array(arr, labels = NULL)
```

## Arguments

- arr:

  A 3D numeric array of dimensions (n, p, 2).

- labels:

  An optional factor or character vector of class labels.

## Value

An `interval_data` object.

## Examples

``` r
arr <- array(c(1, 2, 3, 4, 5, 6, 2, 4, 5, 6, 8, 9),
             dim = c(3, 2, 2))
x <- interval_data_from_array(arr)
```
