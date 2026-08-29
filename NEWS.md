# QAIDR 0.3.0 (2026-07-19)

Inference-correctness release. Adjusted p-values produced with QAIDR <= 0.2.0
are not valid familywise quantities and must be regenerated; unadjusted
(marginal) p-values and all descriptive indices are unchanged.

## Corrections

* `.p_minP()` (the familywise adjustment behind `pvalues_adj` in
  `assess_quality()`) is replaced by a fully symmetric single-step
  randomization min-P construction. The previous implementation compared an
  observed plus-one p-value on the `(m + 1)` denominator against null
  pseudo-p-values ranked on an `m` denominator; this asymmetry breaks the
  exchangeability of the observed row and can inflate the familywise error
  rate (an exact two-coordinate counterexample with `m = 19` gives FWER
  0.0883 at nominal 0.05). The corrected construction pools the observed row
  with all `m` joint draws into one `B = m + 1` row set, computes every
  row's marginal p-value with the common denominator `B` and inclusive
  comparisons, and adjusts by the ECDF of the per-row family minima. It has
  finite-sample weak FWER control at attainable levels under the complete
  random-correspondence null (row exchangeability), verified by an
  independent triple-loop oracle, complete-orbit enumeration at all
  attainable levels, and the counterexample above (new FWER 0.0108).
  Observed marginal p-values are byte-identical to `pvalues`
  (`(c + 1) / (m + 1)`); only `pvalues_adj` changes, and only upward or
  downward within the same family definition (families, tails, shared joint
  draws, and the API are unchanged).

## Tests

* New oracle suite `test-minp-symmetric.R` with independently generated
  fixtures (`minp_oracle_vectors.rds`): exact rational agreement on random
  tied arrays, edge cases (`m = 1`, duplicated rows, all-equal statistics,
  family sizes 1 and 4), complete-orbit weak-FWER enumeration at every
  attainable level, the `m = 19` counterexample, one-hypothesis reduction,
  observed-row identity with the marginal p-value, Q upper-tail versus |B|
  two-sided handling, metric-row relabeling equivariance, and a
  deliberate-failure self-test.

# QAIDR 0.2.0 (2026-07-18)

Correctness release implementing the adjudicated index specification from
Lee and Verleysen (2009), *Neurocomputing* 72:1431--1443. Results produced
with QAIDR 0.1.0 for the MRRE and LCMC-behavior families are not comparable
and should be regenerated.

## Corrections

* `rank_matrix(ties = "random")` now uniformly randomizes every
  tolerance-detected tie block as a whole (transitive adjacent near-tie
  chains included). Previously the detected blocks were passed to base R's
  `rank(..., ties.method = "random")`, which randomizes only exactly equal
  values, so detected near-ties that were not exactly equal were never
  randomized. Tie-free ranks are bit-for-bit unchanged and consume no RNG;
  published daggered cells (all caused by exact structural duplicates) are
  unaffected, as confirmed by the tie-impact audit.
* `.check_dist()` now accepts n = 3 observations, matching the null
  calibration theorem's domain (n >= 3, K <= n - 2).
* `idist_euclidean()` validates `lambda` directly (one finite numeric scalar
  in [0, 1]) instead of relying on the `idist()` dispatcher.
* `idist_hausdorff()` documentation now states explicitly that it computes
  the coordinatewise product-Hausdorff dissimilarity (L2 aggregation of
  per-coordinate interval Hausdorff distances), not the set-Hausdorff metric
  between hyperrectangles. The metric key `"Hausdorff"` is unchanged.
* `rank_matrix()` documentation now explains that the `n_ties` attribute
  counts adjacent tolerance links (a block of size b contributes b - 1) and
  describes the seeded structural-tie protocol used in the manuscript.

## Breaking corrections

* `coranking_indices()` MRRE now follows the printed Lee--Verleysen Eqs.
  (10)--(12): `W_n` sums over all pairs with low-space rank `l <= K` with
  weight `|k-l|/l`; `W_v` sums over all pairs with high-space rank `k <= K`
  with weight `|k-l|/k`. The 0.1.0 implementation omitted the upper-left
  block from both sums; all `Q_RE` and `B_RE` values change.
* `B_LC` is now Lee--Verleysen's behavior index `B_NX` (Eqs. 22--23):
  `B_LC = U_X - U_N` computed on the strict triangles of the K-by-K
  upper-left block. The 0.1.0 upper region was unbounded and therefore had a
  nonzero null mean. Positive values indicate intrusion-dominated embeddings.
* Ranks are computed after excluding the self-distance (`rank_matrix()`).
  The 0.1.0 rank-the-diagonal-and-subtract approach could silently drop pairs
  for duplicate observations and violate `sum(Q) = n(n-1)`.
* Index values are no longer clamped into their nominal ranges; out-of-range
  values raise a diagnostic warning instead.

## New behavior

* Primary analyses error on detected off-diagonal distance ties
  (`ties = "error"`). A documented, seeded fallback (`ties = "random"`)
  resolves tolerance-detected tie blocks uniformly at random.
* `perm_test()` and `assess_quality()` share one permutation engine and both
  use the finite-sample estimator `p = (c + 1)/(m + 1)`; the default is
  `n_perm = 999`.
* Single-step Westfall--Young min-P family-wise adjusted p-values
  (`pvalues_adj`) are provided within each index family (T&C, MRRE, LCMC) at
  fixed K, across metrics, per DR method, reusing the joint permutation draws.
* A center-only Euclidean baseline row (`baseline = TRUE`) is available in
  `assess_quality()` and `k_profiles()`.
* `idist()` exposes `lambda` (Interval Euclidean optimism index) and `nu`
  (Ichino--Yaguchi span weight); `idist_ichino_yaguchi()` renames `gamma` to
  `nu` (a deprecated alias is retained).
* `k_profiles()` computes the co-ranking matrix once per method-metric pair
  and evaluates all K on it.
* The Int-UMAP wrapper passes a configuration object to `RSDA::sym.umap()`
  and records the effective configuration.
* New exported building blocks: `rank_matrix()`, `coranking_matrix()`,
  `indices_from_coranking()`, and `perm_null_indices()`.
* New `interval_data_from_dataSDA()` converter loads the Cars and Face
  examples from the suggested CRAN package `dataSDA`; QAIDR no longer
  redistributes copies of those datasets.
* Package metadata now credits Paulo Canas Rodrigues and Han-Ming Wu as
  authors.
