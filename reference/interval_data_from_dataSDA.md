# Convert a dataSDA data set to interval_data

Converts interval columns used by the dataSDA package into the
center-radius representation used by QAIDR. Both `symbolic_interval`
columns (as in `cars.int`) and character columns formatted as
`"lower,upper"` (as in `face.iGAP`) are supported. Data are loaded from
dataSDA at run time and are not redistributed by QAIDR.

## Usage

``` r
interval_data_from_dataSDA(x, label_col = NULL, labels = NULL)
```

## Arguments

- x:

  A data frame or `symbolic_tbl` containing interval columns.

- label_col:

  Optional name or index of a column containing one class label per
  observation. If omitted, a non-interval column named `class` is used
  when present.

- labels:

  Optional factor or character vector of labels. Supply either `labels`
  or `label_col`, not both.

## Value

An `interval_data` object.

## Examples

``` r
if (requireNamespace("dataSDA", quietly = TRUE)) {
  data("face.iGAP", package = "dataSDA")
  person <- sub("[[:digit:]]+$", "", rownames(face.iGAP))
  face <- interval_data_from_dataSDA(face.iGAP, labels = person)
  print(face)
}
#> interval_data: 27 observations, 6 variables
#> Labels: FRA, HUS, INC, ISA, JPL, KHA, LOT, PHI, ROM
```
