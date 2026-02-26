# Compute co-ranking quality and behavior indices

Given high-dimensional and low-dimensional distance matrices, computes
quality (Q) and behavior (B) indices based on the co-ranking matrix
framework. Three index families are computed: Trustworthiness &
Continuity (T&C), Mean Relative Rank Errors (MRRE), and Local Continuity
Meta-Criterion (LCMC).

## Usage

``` r
coranking_indices(Dh, Dl, K)
```

## Arguments

- Dh:

  Numeric matrix of high-dimensional pairwise distances (n x n).

- Dl:

  Numeric matrix of low-dimensional pairwise distances (n x n).

- K:

  Integer neighborhood size.

## Value

A named numeric vector with elements `Q_TC`, `B_TC`, `Q_RE`, `B_RE`,
`Q_LC`, `B_LC`.

## Examples

``` r
set.seed(42)
X <- matrix(rnorm(50), 10, 5)
Y <- matrix(rnorm(20), 10, 2)
Dh <- as.matrix(dist(X))
Dl <- as.matrix(dist(Y))
coranking_indices(Dh, Dl, K = 3)
#>        Q_TC        B_TC        Q_RE        B_RE        Q_LC        B_LC 
#> 0.653333333 0.066666667 0.670000000 0.008235294 0.500000000 0.266666667 
```
