# Create interval_data from a min-max data frame

Reads a data frame with alternating min/max columns and extracts centers
and radii. This corresponds to the former `extract_CR_from_MM`.

## Usage

``` r
interval_data_from_mm(mm_df, int_min_at, y_at = NULL)
```

## Arguments

- mm_df:

  A data frame with min/max column pairs.

- int_min_at:

  Integer vector of column indices for interval minima.

- y_at:

  Integer or `NULL`. Column index for class labels.

## Value

An `interval_data` object.

## Examples

``` r
df <- data.frame(
  A.min = c(1, 3), A.max = c(2, 5),
  B.min = c(10, 20), B.max = c(12, 25),
  class = c("a", "b")
)
x <- interval_data_from_mm(df, int_min_at = c(1, 3), y_at = 5)
```
