# Prompt B2 preregistration (FROZEN 2026-07-16 after three-provider review)

Review record: Gemini adversarial review returned 4 required changes
(B_LC/MRRE in stability metrics; reduced contamination scale; removal of the
"well-performing linear methods" reporting restriction; refutable C-PCA
criteria). Codex audit returned detailed specification changes for D1–D4 and
identified two pre-existing defects repaired under this addendum (Table III
aggregation; stale prose claims). All changes are incorporated below; the
AMENDMENTS section lists them explicitly. This file is frozen at the hash
recorded in the B2 completion report before any new quality-index result is
computed.

All designs below are fixed before any quality- or behavior-index result of
the new analyses is computed or inspected. Seeds derive from the Phase B base
`SEED_BASE = 20260715` with the stated offsets. Existing Phase B outputs are
reused where valid; nothing already computed is re-selected.

## D1. K-profile evidence (Stage 5 item 5)

- Replication subset rule (all scenarios identically): replications r = 1..25
  (the same prespecified subset already used for permutation calibration).
- K grid rule, defined only from n and the primary K = 10: absolute grid
  K ∈ {5, 10, 20, 40} for every scenario (all satisfy 1 ≤ K ≤ n−2 for
  n ∈ {200, 300, 800}); the corresponding relative fractions K/n are reported
  alongside so scenarios of different n are comparable. The primary K stays 10
  regardless of profile appearance.
- Quantities: all 6 methods × (4 interval dissimilarities + center baseline)
  × all 6 corrected indices, per replication; summaries = mean and SD over the
  25 replications per (K, method, metric, index).
- Embeddings re-derived from the original per-replication seeds (identical
  RNG offsets to scenario{1,2,3}.R; Scenario III uses its unperturbed base
  data). The new adversarial scenario (D2) uses its own first 25 replications.
- Deliverable: one compact appendix table/figure per decision; stability
  statement compares method rankings at K = 5/20/40 vs K = 10 (Kendall tau of
  method ordering by Q_TC per metric, plus descriptive text). No re-selection
  of primary K.

## D2. Heterogeneous/adversarial scenario (Stage 5 item 7) — "Scenario IV"

One scenario jointly stressing all five required factors; n = 300, p = 5,
K = 10, 100 replications (seeds SEED_BASE + 7000 + r, r = 1..100), all 6
methods + center baseline, min-P inference on replications 1..25 with m = 999
(identical policy to Scenarios I–II).

Data-generating mechanism (per replication):
1. Unequal clusters: sizes (200, 80, 20); cluster centers drawn once per
   replication at mu_g ~ N(0, 2^2 I_5), g = 1..3; object centers
   c_i ~ N(mu_{g(i)}, I_5).
2. Center–width association: base radii r_is = 0.15·||c_i − mu_{g(i)}||_2 ·
   u_is with u_is ~ U(0.8, 1.2) — objects farther from their cluster core are
   wider (association strength fixed at 0.15).
3. Contamination: a fixed 5% of objects (rounded, sampled uniformly) have
   centers redrawn from U(−10, 10)^5 and radii multiplied by 5.
4. Partial degeneracy and ties: a fixed 10% of non-contaminated objects get
   all radii set to 0 (degenerate points); additionally the first 2% of
   objects are exact duplicates of the immediately preceding object
   (structural distance ties by construction; handled by the documented
   tie-audit protocol).
5. Metric misspecification: the five coordinates are mixed by a fixed random
   rotation Q_r (drawn once per replication) applied to centers only, so the
   coordinate-wise dissimilarities operate in a rotated frame relative to the
   generating axes.
Estimands and summaries: mean, SD, MC-SE, and 2.5%/97.5% quantiles of each
index over replications; method failure counts; tie-affected cell counts.
Interpretation rule fixed in advance: this scenario tests robustness of the
ASSESSMENT framework (does it still rank methods coherently and flag
distortions under stress), not a claim that any method wins; favorable and
unfavorable outcomes are reported symmetrically.

## D3. Empirical lambda/nu sensitivity (Stage 5 item 8)

- Fixed data (chosen before results): (a) Face dataset (n = 27, K = 5);
  (b) Scenario I replication r = 1 (seed SEED_BASE + 1001, standardized, with
  its seeded C-PCA embedding), K = 10.
- Face: EXACT breakpoint enumeration per Proposition prop:break — for lambda
  (Int-Euclidean): all anchor-pair crossings of the affine scores in [0,1];
  for nu (Ichino–Yaguchi): all roots of the quadratic comparison polynomials
  in [0, 0.5]. Indices evaluated at every interval midpoint plus endpoints
  {0, 1} (resp. {0, 0.5}); report the number of candidate breakpoints, the
  number of DISTINCT index values, the observed range per index, and the
  locations of material changes (|Δindex| > 0.01).
- Scenario I replicate (n = 300): exact enumeration is documented as
  infeasible (≈ n·C(n−1,2) ≈ 1.3×10^7 crossing candidates for lambda alone);
  per the preregistered fallback, a frozen uniform grid of 101 points on
  [0,1] for lambda and 51 points on [0,0.5] for nu is used, justified by the
  step-function structure (the grid lower-bounds the true variation; this is
  stated in the manuscript).
- Guidance rule fixed in advance: defaults lambda = 0.5, nu = 0.5 are judged
  adequate if all indices vary by < 0.02 across the enumerated/grid profile
  for the well-performing linear methods; otherwise the manuscript reports
  the sensitivity honestly and recommends profile reporting. No global
  robustness claim beyond the examined data.

## D4. C-PCA mechanism diagnostics (Stage 5 item 9)

- Data: the existing Scenario I subset r = 1..25 (re-derived embeddings, same
  seeds). Method: C-PCA only (plus V-PCA/SPCA descriptively if the same
  mechanism applies).
- Preregistered mechanism hypothesis H_mech: C-PCA fits its projection from
  interval CENTERS only (pure noise in Scenario I), and the projected
  interval half-widths are formed by propagating the original radii through
  the absolute loadings, radius_out(i,t) = sum_s |w_st| r_is; hence the
  width hierarchy survives in the embedding regardless of the (noise-fitted)
  loadings, which explains high interval-aware quality with uninformative
  centers.
- Diagnostics (per replication, aggregated with mean ± SD):
  (i) R^2 of regressing observed projected radii on |W|^T r (H_mech predicts
  R^2 ≈ 1 up to implementation constants);
  (ii) between/within-group F-ratio of projected mean radii (width-hierarchy
  preservation);
  (iii) center-signal check: variance explained by the top-2 PCs of centers
  (≈ 2/5 under pure noise) and the across-replication dispersion of |w_st|
  (loadings are arbitrary rotations);
  (iv) baseline contrast: same F-ratio computed on projected centers
  (H_mech predicts no group structure).
- Manuscript action rule: if (i) median R^2 > 0.95 and (ii) F-ratio large and
  stable, replace the speculative §5.1 sentence with the evidence-backed
  mechanism (red-wrapped); if not, weaken the manuscript claim to report the
  diagnostics honestly.

## AMENDMENTS (incorporated before freezing; supersede conflicting text above)

A-D1 (K profiles): the grid {5,10,20,40} is a primary-K multiplier rule
{0.5K0, K0, 2K0, 4K0}, K0 = 10 (NOT a relative-fraction design; K/n reported
descriptively). Uncertainty: mean, SD, MCSE, and 2.5%/97.5% quantiles over the
25 replications. Stability metric: replication-level Kendall tau (ties =
average-rank convention via cor(method="kendall")) between the method ordering
at each K and at K0, computed separately for EVERY index (Q_TC, B_TC, Q_RE,
B_RE, Q_LC, B_LC) and every evaluation (4 dissimilarities + baseline);
summarized as mean and 2.5%/97.5% quantiles over replications. Scenario II
retains raw-DR/standardized-evaluation. Scenario IV K-profile cells use the
multi-resolution tie audit (not single-resolution ties="random"). Baseline =
fifth evaluation column, not a seventh method.

A-D2 (Scenario IV, exact executable spec): per replication r (seed
SEED_BASE+7000+r; all draws from this stream in the stated order):
 (0) Reserve six non-overlapping donor/recipient index pairs: donors
     {1,3,5,7,9,11}, recipients {2,4,6,8,10,12} (recipient j is an exact copy
     of donor j-1 AFTER all transformations). These 12 indices are excluded
     from contamination and degeneracy.
 (1) Cluster means mu_g ~ N(0, 4·I5), sizes (200,80,20); centers
     c_i ~ N(mu_g(i), I5); base radii r_is = 0.15 · ||c_i - mu_g(i)||_2 ·
     u_is, u_is ~ U(0.8, 1.2) (purely multiplicative; no additive term).
 (2) Contamination: 15 objects sampled (without replacement) from indices
     13..300; their centers redrawn ~ U(-3, 3)^5 and radii multiplied by 2
     (reduced scale per Gemini to keep the stressor local, not a global
     distance compressor).
 (3) Degeneracy: 30 objects sampled from the remaining eligible pool
     (disjoint from contaminated and reserved): all radii set to 0.
 (4) Misspecification: seeded Haar rotation Q_r (QR decomposition of a 5x5
     N(0,1) matrix; R-sign convention Q <- Q %*% diag(sign(diag(R)))),
     applied to CENTERS ONLY. This is declared as a deliberately misspecified
     observed representation (radii remain axis-aligned); it is NOT a
     geometric rotation of the hyperrectangles, and the manuscript will state
     it as such.
 (5) Duplication: recipients overwritten as exact copies of their donors.
 Records kept per replication: contaminated set, degenerate set, donor/
 recipient pairs, method failures, high-space and low-space structural tie
 counts (separately), tie-affected cell fraction, multi-resolution spread.
 Frozen estimands: (E1) per method x evaluation: mean, SD, MCSE, 2.5%/97.5%
 quantiles of all six indices over 100 replications; (E2) interval-awareness
 contrast: mean over the four dissimilarities of Q_TC minus baseline Q_TC,
 per method, with replication SD; (E3) tie diagnostics as above; (E4) failure
 counts. Inference: min-P policy identical to Scenarios I-II (reps 1..25,
 m = 999, joint draws). A quality-blind runtime/tie-cost pilot (2
 replications; only wall time, memory floor, and tie counts inspected) runs
 before the full launch; if projected full cost exceeds 12 h, stop for author
 approval.

A-D3 (lambda/nu sensitivity): Face arm enumerates the UNION of breakpoints
from the high-space matrix and each method's low-space matrix (both rankings
change with the parameter), per method, all six methods; IY quadratics use
explicitly constructed coefficients with discriminant guards, roots clipped
to the domain, identically-zero comparisons recorded as persistent ties, and
root deduplication at tolerance 1e-10. Evaluation at one representative per
open interval (midpoint) plus domain endpoints; breakpoint values themselves
are tie cases and are characterized by their left/right interval values.
Distinct rankings counted via rank-matrix hash; distinct index vectors at
tolerance 1e-9. Scenario I arm: C-PCA ONLY (scope contradiction resolved);
the 101/51-point frozen grid is reported as DESCRIPTIVE SAMPLED SENSITIVITY
that lower-bounds true variation — it can demonstrate sensitivity but cannot
certify adequacy; the manuscript language must respect this asymmetry.
Guidance rule (amended per Gemini): defaults judged adequate only if the
MAXIMUM variation ACROSS ALL EVALUATED METHODS on the exact (Face) profiles
is < 0.02 per index; the simulation grid is reported descriptively either
way. No global robustness claim.

A-D4 (C-PCA mechanism): the implementation identity C_out = C_in W,
R_out = R_in |W| with W = eigen(cor(centers))$vectors[,1:2] (verified from
symbolicDA::PCA.SDA.r source) is CHECKED numerically via max_abs_error (plus
slope/intercept of the trivial fit), expected ~ machine precision, and is
labeled an implementation identity, NOT a mechanism test. Refutable mechanism
tests, frozen: (T1) independence of loadings from radii: per replication,
correlation between column-wise |w_st| profiles and the group-informative
radius profile; prereg expectation: distribution over 25 reps centered at 0
(95% replication interval covering 0) — a systematic association REFUTES the
center-only-fit claim; (T2) paired separation: pseudo-F (between/within group
mean squares) of projected radii vs projected centers; criterion: the 2.5%
quantile over 25 replications of (F_radii − F_centers) must exceed 0;
(T3) top-2 eigenvalue share of cor(centers) compared to a matched
isotropic-null reference (25 draws of n=300 x p=5 standard normal, seeds
SEED_BASE+7900+r): overlap of the two replication distributions reported;
H_mech does NOT require share ≈ 0.4 exactly. Sign/subspace non-identifiability
handled by using |W| and subspace-invariant statistics only. The mechanism
claim is limited to C-PCA; V-PCA/SPCA are mentioned only descriptively.
Action rule: if T2's criterion holds and T1 shows no systematic association,
replace the speculative manuscript sentence with the evidence-backed
mechanism; otherwise weaken the manuscript claim and report the diagnostics.

A-Repairs (non-design, executed under this addendum): Table III SDs
recomputed as between-dataset SDs of dataset means from
scenario3_per_perturbation.csv (current values pool all 10,000 rows and are
mislabeled); caption/prose relabeled accordingly. Stale prose corrected:
mid-scale runtimes (226/159/84 s; 9.4 min end-to-end incl. K profiles), Face
Q_RE lower bound 0.948, USHCN "linear methods" scoped to C-PCA/V-PCA/SPCA,
Face K-profile figure caption text aligned with what the plots actually show
(zero reference lines on behavior panels only). Stage 4 min-P gate re-run
under the corrected joint-draw engine with the full {Q,B} x (4 dissimilarities
+ baseline) family per index family, machine-readable record per the Codex
schema; the 2026-07-15 record is retained as deprecated audit history.

## Fixed non-design elements

- Familywise policy, tie policy, p-value formula: unchanged from Phase B.
- No primary-K change, no metric re-selection, no scenario re-runs beyond the
  preregistered subsets; Scenario II warning classification (Part 2.7) is a
  diagnostic re-run in quick mode and produces no manuscript numbers.
