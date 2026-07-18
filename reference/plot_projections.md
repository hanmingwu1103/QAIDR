# Plot 2D interval projections

Creates a faceted plot of 2D interval projections from multiple DR
methods, displaying intervals as rectangles.

## Usage

``` r
plot_projections(projections, labels = NULL, obs_labels = NULL)
```

## Arguments

- projections:

  An `idr_projections` object or a named list.

- labels:

  A factor or character vector of class labels.

- obs_labels:

  Optional character vector of observation labels for text.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
C <- matrix(rnorm(80), 20, 4)
R <- matrix(runif(80, 0.05, 0.25), 20, 4)
x <- standardize(interval_data(C, R))
proj <- run_idr(x)
plot_projections(proj)
} # }
```
