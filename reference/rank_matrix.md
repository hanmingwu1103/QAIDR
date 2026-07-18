# Rank matrix with self-exclusion and tie detection

Computes, for each observation, the ranks of all other observations by
distance, excluding the self-distance before ranking. Ranks take values
in `1, ..., N-1`. Off-diagonal ties are detected within a scale-aware
tolerance; under the primary (default) policy `ties = "error"` any tie
aborts with an error, so that ordinary tie-free cells are computed on
strict rank systems as assumed by the co-ranking theory. The opt-in
`ties = "random"` resolves each tolerance-connected tie block by a
uniform random ordering of the whole block (seed the session RNG for
reproducibility); it is used only in explicitly labelled tie audits and
daggered structural-tie cells, where results are summarized over seeded
repeated resolutions, and is never used silently as if the ranks were
tie-free. Randomization covers every detected block, including near-ties
whose members are not exactly equal (base R's random tie method alone
would randomize only exact equals).

## Usage

``` r
rank_matrix(D, ties = c("error", "random"), tol = 1e-12)
```

## Arguments

- D:

  Numeric symmetric distance/dissimilarity matrix (n x n), non-negative,
  zero diagonal.

- ties:

  Character, one of `"error"` (primary) or `"random"` (documented
  fallback).

- tol:

  Numeric relative tolerance used for tie detection (default `1e-12`):
  two consecutive sorted distances are a tie when their gap is at most
  `tol` times the larger of the two. The default catches structural ties
  (duplicate objects, coincident interval geometry) while ignoring
  floating-point near-misses between genuinely distinct distances.

## Value

Integer matrix (n x n) of ranks with 0 on the diagonal and each row's
off-diagonal entries a permutation of `1:(n-1)`. Attribute `"n_ties"` is
the total number of adjacent tolerance links in the sorted off-diagonal
rows; a tolerance-connected block of size b contributes b - 1.

## Examples

``` r
D <- as.matrix(dist(matrix(rnorm(30), 10, 3)))
R <- rank_matrix(D)
```
