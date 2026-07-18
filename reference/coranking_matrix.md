# Co-ranking matrix from two rank matrices

Builds the (N-1) x (N-1) co-ranking matrix \\q\_{kl} = \|\\(i,j):
\rho\_{ij} = k, \gamma\_{ij} = l\\\|\\ from self-excluded rank matrices.
Row index k is the high-dimensional rank, column index l the
low-dimensional rank.

## Usage

``` r
coranking_matrix(Rh, Rl)
```

## Arguments

- Rh:

  Integer rank matrix for the high-dimensional space (from
  [`rank_matrix()`](https://hanmingwu1103.github.io/QAIDR/reference/rank_matrix.md)).

- Rl:

  Integer rank matrix for the low-dimensional space.

## Value

Integer (N-1) x (N-1) matrix of class `"coranking"` with
`sum(Q) == N * (N - 1)`.
