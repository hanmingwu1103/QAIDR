
# QAIDR

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
- **Center-only baselines** to quantify the information added by interval widths
- **Visualization**: 2D projection plots and K-neighbourhood profile plots

## Installation

Install the development version from GitHub:

```r
# install.packages("pak")
pak::pak("hanmingwu1103/QAIDR")

# Example interval datasets used below
install.packages("dataSDA")
```

## Quick Start

```r
library(QAIDR)

# Load the Cars data from dataSDA and convert its interval columns
data("cars.int", package = "dataSDA")
cars <- interval_data_from_dataSDA(cars.int)
x <- standardize(cars)

# Compute interval distances
D <- idist(x, metric = "Wasserstein")

# Run all 6 DR methods (requires symbolicDA, RSDA)
proj <- run_idr(x)

# Assess quality across all method-metric combinations
result <- assess_quality(x, proj, K = 5, perm_test = TRUE, n_perm = 999,
                         baseline = TRUE)
print(result)

# Visualize
plot_projections(proj, labels = cars$labels)

profiles <- k_profiles(x, proj)
plot_k_profiles(profiles, metric = "Wasserstein")
```

## Documentation

- [Package website](https://hanmingwu1103.github.io/QAIDR/)
- [`vignette("introduction")`](https://hanmingwu1103.github.io/QAIDR/articles/introduction.html) -- Basic workflow
- [`vignette("real-data-analysis")`](https://hanmingwu1103.github.io/QAIDR/articles/real-data-analysis.html) -- Cars and Face tutorials using `dataSDA`
- [`vignette("simulation-study")`](https://hanmingwu1103.github.io/QAIDR/articles/simulation-study.html) -- A compact simulation workflow
- [Manuscript reproducibility guide](README_REPRODUCIBILITY.md) -- seeded scripts and archived outputs for the complete study

## Citation

Paulo Canas Rodrigues and Han-Ming Wu (2026), *Interval-metric co-ranking:
quality assessment of dimensionality reduction for interval-valued data*.

For the package citation returned by R, run `citation("QAIDR")`.

## License

MIT
