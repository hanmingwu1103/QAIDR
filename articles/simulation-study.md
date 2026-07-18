# A Simulation Workflow with QAIDR

This vignette gives a compact template for designing and evaluating a
simulation with QAIDR. It is intentionally small and is not an exact
reproduction of the manuscript study. The complete study contains four
scenarios: width-structured groups, a nonlinear-manifold factorial
design, repeated endpoint perturbations, and an adversarial joint-stress
design with unequal clusters, center–width coupling, contamination,
forced ties, and rotated-frame misspecification.

The full seeded designs, 100-replication outputs, Monte Carlo
uncertainty, K-profile robustness checks, parameter-sensitivity
analysis, and calibration gates are available under `analysis/` in the
GitHub repository. See `README_REPRODUCIBILITY.md` before attempting the
full run.

## Generate interval-valued observations

The example below separates center structure from width structure. This
makes the center-only baseline scientifically interpretable: a
difference between an interval-metric row and the baseline measures
information that cannot be attributed to midpoints alone.

``` r

library(QAIDR)
set.seed(20260718)

n <- 120
p <- 5
group <- factor(rep(c("narrow", "wide"), each = n / 2))

centers <- matrix(rnorm(n * p), n, p)
radii <- matrix(runif(n * p, 0.05, 0.15), n, p)
radii[group == "wide", ] <- matrix(
  runif(sum(group == "wide") * p, 0.7, 1.0),
  sum(group == "wide"), p
)

colnames(centers) <- colnames(radii) <- paste0("V", seq_len(p))
x <- standardize(interval_data(centers, radii, labels = group))
```

## Define candidate embeddings

For an algorithm study, replace these illustrative projections with
outputs from
[`run_idr()`](https://hanmingwu1103.github.io/QAIDR/reference/run_idr.md)
or another DR implementation. The first projection retains both centers
and radii in the first two coordinates. The second discards widths,
while the third deliberately permutes the interval widths to create a
diagnostic failure.

``` r

set.seed(20260719)
bad_order <- sample.int(n)

projections <- structure(list(
  `Interval projection` = list(
    C = x$centers[, 1:2, drop = FALSE],
    R = x$radii[, 1:2, drop = FALSE],
    type = "Interval"
  ),
  `Center-only projection` = list(
    C = x$centers[, 1:2, drop = FALSE],
    R = matrix(0, n, 2),
    type = "Point"
  ),
  `Width-permuted projection` = list(
    C = x$centers[, 1:2, drop = FALSE],
    R = x$radii[bad_order, 1:2, drop = FALSE],
    type = "Interval"
  )
), class = "idr_projections")
```

## Assess quality, behavior, and the baseline

Use 999 permutations for inferential results. QAIDR applies the
finite-sample formula `(c + 1) / (m + 1)` and reports within-family
Westfall–Young min-P adjustments. For exploratory debugging, a smaller
number may be used only if the resulting p-values are not interpreted as
final evidence.

``` r

assessment <- assess_quality(
  x, projections, K = 10,
  perm_test = TRUE, n_perm = 999,
  seed = 20260720, baseline = TRUE
)
print(assessment)
assessment$pvalues_adj
```

The key comparison is not simply which row has the largest value. Check
whether the interval metrics distinguish the width-permuted projection
when the center-only baseline cannot, and report uncertainty across
independently generated datasets rather than treating one realization as
conclusive.

## Inspect neighborhood-scale sensitivity

``` r

profiles <- k_profiles(x, projections, K_max = 30, baseline = TRUE)
plot_k_profiles(profiles, metric = "Wasserstein")
```

Stable conclusions near the prespecified primary K are more persuasive
than a single optimized neighborhood size. Behavior values should also
be inspected: positive values indicate intrusion dominance, negative
values indicate extrusion dominance, and values near zero indicate
balance.

## Recommended reporting checklist

1.  Prespecify the data-generating mechanisms, primary K, metrics, and
    contrasts before viewing the index results.
2.  Include center-only baselines and negative controls that should show
    little or no interval advantage.
3.  Repeat independent datasets and report standard deviations or Monte
    Carlo standard errors.
4.  Audit exact and tolerance-detected ties; use seeded random
    resolution only when structural ties are unavoidable.
5.  Report K-profiles and sensitivity to the interval-dissimilarity
    parameters `lambda` and `nu` when conclusions depend on their
    defaults.
6.  Validate permutation calibration and multiplicity control under a
    null design before attaching significance markers.
