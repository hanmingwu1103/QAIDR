# Quality and behavior indices from a co-ranking matrix

Computes the six adjudicated indices from a co-ranking matrix at
neighborhood size K, following Lee & Verleysen (2009, Neurocomputing
72:1431-1443): Trustworthiness & Continuity (Eqs. 7-9, 26-27), Mean
Relative Rank Errors (Eqs. 10-12, 28-29), the Local Continuity
Meta-Criterion (Eq. 13), and the behavior index \\B\_{LC} \equiv B\_{NX}
= U_X - U_N\\ (Eqs. 22-23) computed on the bounded triangles of the K x
K upper-left block. Positive behavior values indicate
intrusion-dominated embeddings for all three families.

## Usage

``` r
indices_from_coranking(Q, K)
```

## Arguments

- Q:

  Co-ranking matrix from
  [`coranking_matrix()`](https://hanmingwu1103.github.io/QAIDR/reference/coranking_matrix.md).

- K:

  Integer neighborhood size, `1 <= K <= N - 2`.

## Value

Named numeric vector: `Q_TC, B_TC, Q_RE, B_RE, Q_LC, B_LC`.
