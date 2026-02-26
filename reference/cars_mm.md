# Cars interval data

Interval-valued measurements of 27 car models across 4 variables: Price,
Engine Capacity (EngCap), Top Speed, and Acceleration. Each observation
is an interval representing the range of values within a car class.

## Usage

``` r
cars_mm
```

## Format

An `interval_data` object with 27 observations and 4 variables:

- centers:

  27 x 4 matrix of interval midpoints.

- radii:

  27 x 4 matrix of interval half-widths.

- labels:

  Factor with 4 levels: Berlina, Luxury, Sportive, Utilitarian.

## Source

Billard, L. and Diday, E. (2006). *Symbolic Data Analysis: Conceptual
Statistics and Data Mining*. Wiley.
