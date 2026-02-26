
# QAIDR <img src="man/figures/logo.png" align="right" height="138" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/hanmingwu1103/QAIDR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hanmingwu1103/QAIDR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Quality Assessment for Interval-Based Dimensionality Reduction**

QAIDR provides tools for evaluating how well interval-based dimensionality
reduction (DR) methods preserve the structure of high-dimensional interval
data, using a co-ranking matrix framework.

## Features

- **4 interval distance metrics**: Interval Euclidean, Hausdorff,
  Ichino-Yaguchi, L2-Wasserstein
- **6 co-ranking indices**: Quality (Q) and Behavior (B) variants of
  Trustworthiness & Continuity, MRRE, and LCMC
- **6 DR method wrappers**: C-PCA, V-PCA, MR-PCA, SPCA, IMDS, Int-UMAP
- **Permutation tests** for statistical significance
- **Visualization**: 2D projection plots and K-neighbourhood profile plots

## Installation

Install the development version from GitHub:

```r
# install.packages("pak")
pak::pak("hanmingwu1103/QAIDR")
```

## Quick Start

```r
library(QAIDR)

# Load and standardize the built-in Cars dataset
data(cars_mm)
x <- standardize(cars_mm)

# Compute interval distances
D <- idist(x, metric = "Wasserstein")

# Run all 6 DR methods (requires symbolicDA, RSDA)
proj <- run_idr(x)

# Assess quality across all method-metric combinations
result <- assess_quality(x, proj, K = 5, perm_test = TRUE, n_perm = 1000)
print(result)

# Visualize
plot_projections(proj, labels = cars_mm$labels)

profiles <- k_profiles(x, proj)
plot_k_profiles(profiles, metric = "Wasserstein")
```

## Documentation

- [Package website](https://hanmingwu1103.github.io/QAIDR/)
- [`vignette("introduction")`](https://hanmingwu1103.github.io/QAIDR/articles/introduction.html) -- Basic workflow
- [`vignette("real-data-analysis")`](https://hanmingwu1103.github.io/QAIDR/articles/real-data-analysis.html) -- Reproduces the real data analysis
- [`vignette("simulation-study")`](https://hanmingwu1103.github.io/QAIDR/articles/simulation-study.html) -- Reproduces the simulation study

## License

MIT
