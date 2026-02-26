# Face anthropometry interval data

Interval-valued anthropometric face measurements of 27 individuals
across 6 variables: AD, BC, AH, DH, EH, GH (distances between facial
landmarks).

## Usage

``` r
facedata_mm
```

## Format

An `interval_data` object with 27 observations and 6 variables:

- centers:

  27 x 6 matrix of interval midpoints.

- radii:

  27 x 6 matrix of interval half-widths.

- labels:

  Factor with 9 levels (3-letter person identifiers).

## Source

Billard, L. and Diday, E. (2006). *Symbolic Data Analysis: Conceptual
Statistics and Data Mining*. Wiley.
