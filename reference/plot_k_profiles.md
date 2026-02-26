# Plot quality/behavior index profiles over K

Creates a 2 x 3 grid of line plots (Quality and Behavior for T&C, MRRE,
LCMC) for a given distance metric.

## Usage

``` r
plot_k_profiles(profile_data, metric = NULL)
```

## Arguments

- profile_data:

  A data frame as returned by
  [`k_profiles()`](https://hanmingwu1103.github.io/QAIDR/reference/k_profiles.md).

- metric:

  Character string specifying which distance metric to plot. If `NULL`,
  returns a list of plots for all metrics.

## Value

A `ggplot` (or list of plots if `metric` is `NULL`).

## Examples

``` r
if (FALSE) { # \dontrun{
data(cars_mm)
x <- standardize(cars_mm)
proj <- run_idr(x)
profiles <- k_profiles(x, proj, K_max = 10)
plot_k_profiles(profiles, metric = "Wasserstein")
} # }
```
