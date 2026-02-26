# Generate interval Swiss roll data

Generates a synthetic interval-valued dataset embedded on a Swiss roll
manifold. Useful for benchmarking interval dimensionality reduction
methods.

## Usage

``` r
gen_interval_swissroll(
  n = 800,
  seed = 1,
  t_min = 1.5 * pi,
  t_max = 4.5 * pi,
  y_max = 70,
  noise_sd = 0.02,
  r_min = 0.2,
  r_max = 0.5
)
```

## Arguments

- n:

  Integer number of observations (default 800).

- seed:

  Integer random seed (default 1).

- t_min:

  Minimum angle parameter (default `1.5 * pi`).

- t_max:

  Maximum angle parameter (default `4.5 * pi`).

- y_max:

  Maximum value for the linear dimension (default 70).

- noise_sd:

  Standard deviation of Gaussian noise on centers (default 0.02).

- r_min:

  Minimum interval radius (default 0.2).

- r_max:

  Maximum interval radius (default 0.5).

## Value

An `interval_data` object with 3 variables (x, y, z).

## Examples

``` r
dat <- gen_interval_swissroll(n = 100, seed = 42)
print(dat)
#> interval_data: 100 observations, 3 variables
```
