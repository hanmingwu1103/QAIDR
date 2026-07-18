# Package index

## Interval Data

Create and manipulate interval-valued data objects.

- [`interval_data()`](https://hanmingwu1103.github.io/QAIDR/reference/interval_data.md)
  : Create an interval_data object
- [`interval_data_from_mm()`](https://hanmingwu1103.github.io/QAIDR/reference/interval_data_from_mm.md)
  : Create interval_data from a min-max data frame
- [`interval_data_from_array()`](https://hanmingwu1103.github.io/QAIDR/reference/interval_data_from_array.md)
  : Create interval_data from a 3D min-max array
- [`interval_data_from_dataSDA()`](https://hanmingwu1103.github.io/QAIDR/reference/interval_data_from_dataSDA.md)
  : Convert a dataSDA data set to interval_data
- [`standardize()`](https://hanmingwu1103.github.io/QAIDR/reference/standardize.md)
  : Standardize interval data

## Distance Metrics

Interval distance functions for comparing interval-valued observations.

- [`idist()`](https://hanmingwu1103.github.io/QAIDR/reference/idist.md)
  : Compute interval distance matrix
- [`idist_euclidean()`](https://hanmingwu1103.github.io/QAIDR/reference/idist_euclidean.md)
  : Interval Euclidean distance
- [`idist_hausdorff()`](https://hanmingwu1103.github.io/QAIDR/reference/idist_hausdorff.md)
  : Product-Hausdorff distance for intervals
- [`idist_ichino_yaguchi()`](https://hanmingwu1103.github.io/QAIDR/reference/idist_ichino_yaguchi.md)
  : Ichino-Yaguchi dissimilarity for intervals
- [`idist_wasserstein()`](https://hanmingwu1103.github.io/QAIDR/reference/idist_wasserstein.md)
  : L2-Wasserstein (Mallows) distance for intervals

## Quality Assessment

Co-ranking indices, permutation tests, and K-profiles.

- [`coranking_indices()`](https://hanmingwu1103.github.io/QAIDR/reference/coranking_indices.md)
  : Compute co-ranking quality and behavior indices from distance
  matrices
- [`rank_matrix()`](https://hanmingwu1103.github.io/QAIDR/reference/rank_matrix.md)
  : Rank matrix with self-exclusion and tie detection
- [`coranking_matrix()`](https://hanmingwu1103.github.io/QAIDR/reference/coranking_matrix.md)
  : Co-ranking matrix from two rank matrices
- [`indices_from_coranking()`](https://hanmingwu1103.github.io/QAIDR/reference/indices_from_coranking.md)
  : Quality and behavior indices from a co-ranking matrix
- [`assess_quality()`](https://hanmingwu1103.github.io/QAIDR/reference/assess_quality.md)
  : Assess quality of interval DR projections
- [`perm_test()`](https://hanmingwu1103.github.io/QAIDR/reference/perm_test.md)
  : Permutation test for a single pair of distance matrices
- [`perm_null_indices()`](https://hanmingwu1103.github.io/QAIDR/reference/perm_null_indices.md)
  : Permutation null distribution of the indices
- [`perm_null_indices_reference()`](https://hanmingwu1103.github.io/QAIDR/reference/perm_null_indices_reference.md)
  : Reference (slow) permutation null via full co-ranking reconstruction
- [`k_profiles()`](https://hanmingwu1103.github.io/QAIDR/reference/k_profiles.md)
  : Compute quality/behavior index profiles over K

## Dimensionality Reduction

Wrappers for interval DR methods.

- [`run_idr()`](https://hanmingwu1103.github.io/QAIDR/reference/run_idr.md)
  : Run interval dimensionality reduction methods
- [`umap_config_default()`](https://hanmingwu1103.github.io/QAIDR/reference/umap_config_default.md)
  : Default UMAP configuration for interval UMAP
- [`umap_config_strong()`](https://hanmingwu1103.github.io/QAIDR/reference/umap_config_strong.md)
  : Strong UMAP configuration for interval UMAP
- [`gen_interval_swissroll()`](https://hanmingwu1103.github.io/QAIDR/reference/gen_interval_swissroll.md)
  : Generate interval Swiss roll data

## Visualization

Plotting functions for projections and index profiles.

- [`plot_projections()`](https://hanmingwu1103.github.io/QAIDR/reference/plot_projections.md)
  : Plot 2D interval projections
- [`plot_k_profiles()`](https://hanmingwu1103.github.io/QAIDR/reference/plot_k_profiles.md)
  : Plot quality/behavior index profiles over K
