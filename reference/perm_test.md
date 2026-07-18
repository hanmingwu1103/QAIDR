# Permutation test for a single pair of distance matrices

Shares the permutation engine of
[`assess_quality()`](https://hanmingwu1103.github.io/QAIDR/reference/assess_quality.md):
ranks are computed once and permuted (label-equivariance), the null uses
a single uniform bijection applied to both endpoints, and p-values use
`(c + 1) / (m + 1)` so the smallest attainable p-value is `1 / (m + 1)`.

## Usage

``` r
perm_test(
  D_high,
  D_low,
  K,
  n_perm = 999,
  ties = c("error", "random"),
  seed = NULL
)
```

## Arguments

- D_high:

  High-dimensional dissimilarity matrix.

- D_low:

  Low-dimensional dissimilarity matrix.

- K:

  Integer neighborhood size.

- n_perm:

  Number of permutations (default 999).

- ties:

  Tie policy (see
  [`rank_matrix()`](https://hanmingwu1103.github.io/QAIDR/reference/rank_matrix.md)).

- seed:

  Optional integer seed.

## Value

List with `vals` (observed indices), `pQ` (one-sided upper p-values for
quality indices), `pB` (two-sided p-values for behavior indices), and
`null_stats` (the m x 6 null draws).
