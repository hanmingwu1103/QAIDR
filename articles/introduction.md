# Introduction to QAIDR

## Overview

**QAIDR** (Quality Assessment for Interval-based Dimensionality
Reduction) provides tools for evaluating interval DR methods using
co-ranking matrices. It supports:

- **4 interval distance metrics**: Interval Euclidean, Hausdorff,
  Ichino-Yaguchi, L2-Wasserstein
- **3 quality/behavior index families**: T&C, MRRE, LCMC (6 indices
  total)
- **6 DR methods**: C-PCA, V-PCA, MR-PCA, SPCA, IMDS, Int-UMAP
- **Statistical tests**: permutation inference with family-wise
  adjustment
- **Baselines**: center-only assessment to isolate information from
  widths

## Quick Start

### 1. Creating Interval Data

Interval-valued data consists of observations where each variable is an
interval $`[a, b]`$, represented internally by its center $`(a+b)/2`$
and radius $`(b-a)/2`$.

``` r

library(QAIDR)

# From centers and radii matrices
C <- matrix(rnorm(40), 10, 4)
R <- matrix(runif(40, 0.1, 0.5), 10, 4)
x <- interval_data(C, R, labels = rep(c("A", "B"), each = 5))
print(x)
#> interval_data: 10 observations, 4 variables
#> Labels: A, B
```

You can also create interval data from common formats:

``` r

# From a min-max data frame
df <- data.frame(
  A.min = c(1, 3, 5), A.max = c(2, 5, 8),
  B.min = c(10, 20, 30), B.max = c(12, 25, 35),
  class = c("x", "y", "x")
)
x2 <- interval_data_from_mm(df, int_min_at = c(1, 3), y_at = 5)
print(x2)
#> interval_data: 3 observations, 2 variables
#> Labels: x, y

# From a 3D array (n x p x 2)
arr <- array(c(1, 2, 3, 4, 5, 6, 2, 4, 5, 6, 8, 9), dim = c(3, 2, 2))
x3 <- interval_data_from_array(arr)
```

### 2. Standardization

Standardize interval data so that each variable has zero-mean centers
and unit standard deviation:

``` r

xs <- standardize(x)
# Column means of centers are now 0
round(colMeans(xs$centers), 10)
#> [1] 0 0 0 0
# Column SDs of centers are now 1
round(apply(xs$centers, 2, sd), 10)
#> [1] 1 1 1 1
```

### 3. Computing Interval Distances

Four distance metrics are available for comparing interval-valued
observations:

``` r

D1 <- idist(xs, metric = "Wasserstein")
D2 <- idist(xs, metric = "Hausdorff")
D3 <- idist(xs, metric = "Int-Euclidean")
D4 <- idist(xs, metric = "Ichino-Yaguchi")

cat("Distance matrix dimensions:", dim(D1), "\n")
#> Distance matrix dimensions: 10 10
```

### 4. Co-ranking Indices

Given high-dimensional and low-dimensional distance matrices, compute
the six quality/behavior indices:

``` r

set.seed(42)
Y <- matrix(rnorm(20), 10, 2)
Dl <- as.matrix(dist(Y))
Dh <- idist(xs, metric = "Wasserstein")

indices <- coranking_indices(Dh, Dl, K = 3)
print(indices)
#>        Q_TC        B_TC        Q_RE        B_RE        Q_LC        B_LC 
#>  0.39333333 -0.06666667  0.42647059 -0.01294118  0.20000000  0.13333333 
#> attr(,"Q")
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    0    2    0    0    2    2    0    3    1
#>  [2,]    0    0    2    1    0    0    1    2    4
#>  [3,]    0    0    2    1    2    2    0    3    0
#>  [4,]    3    0    1    0    1    1    1    0    3
#>  [5,]    1    1    1    3    0    0    1    1    2
#>  [6,]    3    2    0    2    0    0    3    0    0
#>  [7,]    0    2    1    0    4    1    1    1    0
#>  [8,]    1    2    2    0    1    2    2    0    0
#>  [9,]    2    1    1    3    0    2    1    0    0
#> attr(,"class")
#> [1] "coranking" "matrix"    "array"    
#> attr(,"n_ties")
#> [1] 0
```

- **Quality indices** (Q_TC, Q_RE, Q_LC): range \[0, 1\], higher is
  better
- **Behavior indices** (B_TC, B_RE, B_LC): range \[-1, 1\], 0 means
  balanced; positive values indicate intrusion-dominated embeddings

### 5. Using `dataSDA` Datasets

QAIDR does not redistribute example datasets. The converter below reads
the interval encodings used by the CRAN package `dataSDA`.

``` r

if (requireNamespace("dataSDA", quietly = TRUE)) {
  data("cars.int", package = "dataSDA")
  cars <- interval_data_from_dataSDA(cars.int)
  print(cars)

  data("face.iGAP", package = "dataSDA")
  face_labels <- sub("[[:digit:]]+$", "", rownames(face.iGAP))
  face <- interval_data_from_dataSDA(face.iGAP, labels = face_labels)
  print(face)
} else {
  message("Install 'dataSDA' to run the dataset examples.")
}
#> interval_data: 27 observations, 4 variables
#> Labels: Utilitarian, Berlina, Sportive, Luxury
#> interval_data: 27 observations, 6 variables
#> Labels: FRA, HUS, INC, ISA, JPL, KHA, LOT, PHI, ROM
```

## Full Workflow

The typical QAIDR workflow is:

1.  Load or create `interval_data`
2.  Standardize with
    [`standardize()`](https://hanmingwu1103.github.io/QAIDR/reference/standardize.md)
3.  Run DR methods with
    [`run_idr()`](https://hanmingwu1103.github.io/QAIDR/reference/run_idr.md)
4.  Assess quality with
    [`assess_quality()`](https://hanmingwu1103.github.io/QAIDR/reference/assess_quality.md)
5.  Visualize with
    [`plot_projections()`](https://hanmingwu1103.github.io/QAIDR/reference/plot_projections.md)
    and
    [`plot_k_profiles()`](https://hanmingwu1103.github.io/QAIDR/reference/plot_k_profiles.md)

``` r

# Complete example (requires dataSDA, symbolicDA, RSDA, umap, and dplyr)
data("cars.int", package = "dataSDA")
cars <- interval_data_from_dataSDA(cars.int)
x <- standardize(cars)

proj <- run_idr(x)
result <- assess_quality(x, proj, K = 5, perm_test = TRUE, n_perm = 999,
                         baseline = TRUE)
print(result)

plot_projections(proj, labels = cars$labels)

profiles <- k_profiles(x, proj, K_max = 20)
plot_k_profiles(profiles, metric = "Wasserstein")
```

See the “Real Data Analysis” and “Simulation Study” vignettes for
compact tutorials. The complete seeded manuscript pipeline and frozen
outputs are documented in `README_REPRODUCIBILITY.md` in the GitHub
repository.
