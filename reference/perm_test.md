# Permutation test for a single DR projection

Performs a permutation test to assess statistical significance of
quality and behavior indices.

## Usage

``` r
perm_test(D_high, D_low, K, n_perm = 1000)
```

## Arguments

- D_high:

  High-dimensional distance matrix (n x n).

- D_low:

  Low-dimensional distance matrix (n x n).

- K:

  Integer neighborhood size.

- n_perm:

  Number of permutations (default 1000).

## Value

A list with elements:

- vals:

  Named vector of observed index values.

- pQ:

  P-values for quality indices (one-tailed).

- pB:

  P-values for behavior indices (two-tailed).

## Examples

``` r
set.seed(42)
Dh <- as.matrix(dist(matrix(rnorm(50), 10, 5)))
Dl <- as.matrix(dist(matrix(rnorm(20), 10, 2)))
pt <- perm_test(Dh, Dl, K = 3, n_perm = 99)
```
