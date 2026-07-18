# Compute co-ranking quality and behavior indices from distance matrices

Convenience wrapper: ranks both distance matrices (self-excluded,
tie-checked), forms the co-ranking matrix, and evaluates the six
indices.

## Usage

``` r
coranking_indices(Dh, Dl, K, ties = c("error", "random"), tol = 1e-12)
```

## Arguments

- Dh:

  Numeric matrix of high-dimensional pairwise dissimilarities.

- Dl:

  Numeric matrix of low-dimensional pairwise dissimilarities.

- K:

  Integer neighborhood size.

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

Named numeric vector `Q_TC, B_TC, Q_RE, B_RE, Q_LC, B_LC`, with the
co-ranking matrix attached as attribute `"Q"` and the tie counts as
attribute `"n_ties"`.

## Examples

``` r
set.seed(42)
X <- matrix(rnorm(50), 10, 5)
Y <- matrix(rnorm(20), 10, 2)
coranking_indices(as.matrix(dist(X)), as.matrix(dist(Y)), K = 3)
#>        Q_TC        B_TC        Q_RE        B_RE        Q_LC        B_LC 
#>  0.65333333  0.06666667  0.60470588  0.01176471  0.50000000 -0.03333333 
#> attr(,"Q")
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    1    1    2    0    1    1    1    3    0
#>  [2,]    2    3    2    2    0    0    1    0    0
#>  [3,]    2    2    0    1    2    2    0    0    1
#>  [4,]    0    1    1    2    1    2    2    1    0
#>  [5,]    1    1    0    1    1    2    2    0    2
#>  [6,]    0    0    2    1    2    1    3    0    1
#>  [7,]    2    0    0    3    0    2    0    1    2
#>  [8,]    1    1    3    0    1    0    1    2    1
#>  [9,]    1    1    0    0    2    0    0    3    3
#> attr(,"class")
#> [1] "coranking" "matrix"    "array"    
#> attr(,"n_ties")
#> [1] 0
```
