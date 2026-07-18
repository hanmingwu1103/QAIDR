# Permutation null distribution of the indices

Generates the joint permutation null for the six indices under the
uniform random-correspondence null: a single bijection \\\Pi\\ is drawn
uniformly and applied to both endpoints of the embedded ranks
(\\\gamma^\Pi\_{ij} = \gamma\_{\Pi(i)\Pi(j)}\\), while the
high-dimensional ranks stay fixed. Because ranks are label-equivariant,
no re-ranking is needed; each draw permutes the rows and columns of
`Rl`.

## Usage

``` r
perm_null_indices(Rh, Rl, K, m = 999, perms = NULL)
```

## Arguments

- Rh, Rl:

  Rank matrices from
  [`rank_matrix()`](https://hanmingwu1103.github.io/QAIDR/reference/rank_matrix.md).

- K:

  Integer neighborhood size.

- m:

  Integer number of permutations (default 999).

- perms:

  Optional m x n integer matrix whose rows are permutations of `1:n`.
  When supplied, these joint draws are used instead of fresh
  [`sample.int()`](https://rdrr.io/r/base/sample.html) draws. Required
  for Westfall-Young min-P families: every family member (e.g., every
  metric for one DR method) must be evaluated under the SAME permutation
  draws so that the joint null dependence is preserved.

## Value

An m x 6 matrix of null index values (columns as in
[`indices_from_coranking()`](https://hanmingwu1103.github.io/QAIDR/reference/indices_from_coranking.md)).
