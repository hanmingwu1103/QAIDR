# Run interval dimensionality reduction methods

Applies multiple interval DR methods (C-PCA, V-PCA, MR-PCA, SPCA, IMDS,
Int-UMAP) to an `interval_data` object and returns the projections.

## Usage

``` r
run_idr(
  x,
  methods = c("C-PCA", "V-PCA", "MR-PCA", "SPCA", "IMDS", "Int-UMAP"),
  labels = NULL,
  no_dim = 2,
  umap_config = umap_config_default(),
  verbose = TRUE
)
```

## Arguments

- x:

  An `interval_data` object (should be standardized).

- methods:

  Character vector of methods to run. Default includes all six:
  `"C-PCA"`, `"V-PCA"`, `"MR-PCA"`, `"SPCA"`, `"IMDS"`, `"Int-UMAP"`.

- labels:

  Optional factor of class labels. If `NULL`, uses `x$labels`.

- no_dim:

  Integer number of target dimensions (default 2).

- umap_config:

  List of UMAP hyperparameters. Default uses
  [`umap_config_default()`](https://hanmingwu1103.github.io/QAIDR/reference/umap_config_default.md).

- verbose:

  Logical; print progress messages (default `TRUE`).

## Value

An object of class `idr_projections`: a named list where each element
contains `C` (centers matrix), `R` (radii matrix), and `type` ("Point"
or "Interval").

## Examples

``` r
if (FALSE) { # \dontrun{
C <- matrix(rnorm(80), 20, 4)
R <- matrix(runif(80, 0.05, 0.25), 20, 4)
x <- standardize(interval_data(C, R))
proj <- run_idr(x)
} # }
```
