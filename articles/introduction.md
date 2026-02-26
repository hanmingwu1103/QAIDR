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
- **Statistical tests**: Permutation-based significance testing

## Quick Start

### 1. Creating Interval Data

Interval-valued data consists of observations where each variable is an
interval $\lbrack a,b\rbrack$, represented internally by its center
$(a + b)/2$ and radius $(b - a)/2$.

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
#>         Q_TC         B_TC         Q_RE         B_RE         Q_LC         B_LC 
#>  0.393333333 -0.066666667  0.442941176 -0.003529412  0.200000000  0.666666667
```

- **Quality indices** (Q_TC, Q_RE, Q_LC): range \[0, 1\], higher is
  better
- **Behavior indices** (B_TC, B_RE, B_LC): range \[-1, 1\], 0 means
  balanced

### 5. Using Built-in Datasets

``` r
data(cars_mm)
print(cars_mm)
#> interval_data: 27 observations, 4 variables
#> Labels: Berlina, Luxury, Sportive, Utilitarian
summary(cars_mm)
#> interval_data: 27 observations, 4 variables
#> 
#> Center summary:
#>      Price            EngCap        TopSpeed      Acceleration   
#>  Min.   : 21342   Min.   :1173   Min.   :157.0   Min.   : 4.500  
#>  1st Qu.: 34420   1st Qu.:1681   1st Qu.:194.2   1st Qu.: 8.000  
#>  Median : 60900   Median :2388   Median :213.5   Median : 9.000  
#>  Mean   :102406   Mean   :2580   Mean   :216.0   Mean   : 9.667  
#>  3rd Qu.:169238   3rd Qu.:3378   3rd Qu.:235.0   3rd Qu.:11.500  
#>  Max.   :315992   Max.   :4530   Max.   :296.5   Max.   :14.500  
#> 
#> Radii summary:
#>      Price            EngCap          TopSpeed      Acceleration  
#>  Min.   :  2850   Min.   :  78.5   Min.   : 0.50   Min.   :0.000  
#>  1st Qu.:  6064   1st Qu.: 259.8   1st Qu.: 6.50   1st Qu.:0.500  
#>  Median : 11890   Median : 447.0   Median :12.00   Median :1.000  
#>  Mean   : 31869   Mean   : 591.7   Mean   :11.07   Mean   :1.333  
#>  3rd Qu.: 42689   3rd Qu.: 822.2   3rd Qu.:14.00   3rd Qu.:2.000  
#>  Max.   :160081   Max.   :1720.5   Max.   :24.50   Max.   :4.000  
#> 
#> Label distribution:
#> 
#>     Berlina      Luxury    Sportive Utilitarian 
#>           8           4           8           7

data(facedata_mm)
print(facedata_mm)
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
# Complete example (requires symbolicDA, RSDA)
data(cars_mm)
x <- standardize(cars_mm)

proj <- run_idr(x)
result <- assess_quality(x, proj, K = 5, perm_test = TRUE, n_perm = 1000)
print(result)

plot_projections(proj, labels = cars_mm$labels)

profiles <- k_profiles(x, proj, K_max = 20)
plot_k_profiles(profiles, metric = "Wasserstein")
```

See the “Real Data Analysis” and “Simulation Study” vignettes for
complete examples.
