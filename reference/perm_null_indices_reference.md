# Reference (slow) permutation null via full co-ranking reconstruction

Used to validate the fast pair-list engine in
[`perm_null_indices()`](https://hanmingwu1103.github.io/QAIDR/reference/perm_null_indices.md).

## Usage

``` r
perm_null_indices_reference(Rh, Rl, K, m = 999)
```

## Arguments

- Rh, Rl:

  Rank matrices from
  [`rank_matrix()`](https://hanmingwu1103.github.io/QAIDR/reference/rank_matrix.md).

- K:

  Integer neighborhood size.

- m:

  Integer number of permutations (default 999).
