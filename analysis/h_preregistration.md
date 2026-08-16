# Prompt H implementation contract and preregistration (R1)

STATUS: FROZEN before any new result is computed. The SHA-256 of this file
at freeze time is recorded in `H_PREREGISTRATION.freeze` and in the Prompt H
report. Author decisions H.1-H.10 (2026-07-19) govern.

## 1. Symmetric randomization min-P (replaces `.p_minP`)

### 1.1 Algorithm (exact)

Inputs per DR method and index family F in {TC, RE, LC}: the observed
statistic vector and the m = 999 stored joint random-correspondence draws
(shared across every family member), family membership fixed in advance as
every {Q, B} x evaluation-dissimilarity test currently passed jointly to
`.p_minP` (Q_F and B_F for each of the evaluated dissimilarities, including
the center baseline where it is evaluated). No member may be added or
removed after inspecting p-values.

Scale: larger-is-more-extreme. Quality indices enter as T = Q (one-sided
upper); behavior indices enter as T = |B| (two-sided via absolute value),
identically for observed and null rows.

With B = m + 1 rows (b = 0 the observed row, b = 1..m the draws):

    p[b,h]  = (1/B) * sum_{r=0..m} 1{ T[r,h] >= T[b,h] }     (inclusive)
    M[b]    = min_h p[b,h]
    padj[h] = (1/B) * sum_{b=0..m} 1{ M[b] <= p[0,h] }        (inclusive)

Properties required of the implementation: p[0,h] equals the conventional
plus-one Monte Carlo p-value (1 + #{s : T[s,h] >= T[0,h]})/(m+1); one common
denominator B for every row; inclusive comparisons throughout (conservative
under statistic ties); the old rank/m null-only array must not appear.

### 1.2 Claim (exact scope) and proof obligations

Claimed: finite-sample WEAK familywise error control at attainable levels
under the complete random-correspondence null, assuming (i) the B joint
statistic rows are exchangeable under that null, (ii) the family is fixed in
advance, (iii) the same joint permutations are used for all family members,
and (iv) any auxiliary tie-resolution randomness is independent of the
permutations and label-equivariant or jointly coupled. NOT claimed: strong
control, subset pivotality, or anything about partial nulls.

Claude-seat derivation (independent): let F_M(x) = (1/B) sum_b 1{M[b] <= x}.
Since F_M is nondecreasing and padj[h] = F_M(p[0,h]), we have
min_h padj[h] = F_M(min_h p[0,h]) = F_M(M[0]). Any rejection at level alpha
occurs iff F_M(M[0]) <= alpha, i.e., iff the inclusive lower rank R of M[0]
among {M[0], ..., M[m]} satisfies R <= alpha*B. Under the complete null the
rows are exchangeable, hence so are the M[b]; the inclusive rank of an
exchangeable coordinate is super-uniform on {1/B, ..., 1} (ties only
increase R), so P(R/B <= alpha) <= alpha, with equality at attainable
levels alpha in {1/B, ..., 1} in the absence of ties among the M[b]. QED.

Exchangeability with the actual sampler: the m permutations are drawn
i.i.d. uniformly WITH replacement, independently of data (may repeat; may
equal the identity); under H0 the observed correspondence is itself
uniform-permutation distributed and independent of the draws, so the B rows
are i.i.d.-conditionally-on-data and hence exchangeable. Duplicated rows and
identity draws are therefore admissible and handled by inclusive ties.

Gate: Codex must independently derive the same result and supply an oracle
implementation + test vectors; Gemini must adversarially search for
counterexamples (discreteness, statistic ties, absolute-value statistics,
repeated/identity permutations, sampled vs complete orbits, multi-metric
families). Two-provider acceptance of statement + proof + code before the
manuscript text is written. If the exchangeability argument fails against
the actual sampler: STOP (no silent Holm switch; Holm requires a new author
decision).

### 1.3 Oracle tests (all mandatory, see Prompt H Stage 3.1 list)

Exact independent-oracle agreement on random small arrays; one-hypothesis
reduction; common denominator + inclusive ties; Q vs |B| tails; shared
draws; relabeling equivariance; complete small-orbit weak-FWER enumeration
at ALL attainable alphas; the Prompt G m=19 two-coordinate counterexample
(old formula must fail, new must control); edge cases m=1, duplicated rows,
all-equal statistics, family size 1 and >2; deliberate-failure self-test.

### 1.4 Downstream regeneration (fixed list)

Regenerate with the corrected algorithm, production seeds, m=999, unchanged
family definitions: scenario1_calibration, scenario2_calibration,
scenario4_calibration (reps 1-25 each), face p-values, midscale (USHCN)
p-values, and every manuscript marker derived from `pvalues_adj`. Marginal
p-values and all descriptive indices must be byte-unchanged (assertion in
each regeneration script). Produce a complete old-vs-new adjusted-p/marker
delta table. Rerun the m=999, 2,000-null-replicate familywise calibration
audit (frozen design, fresh disjoint seeds: base 20260719 offsets) and
report per-family rates with Wilson 95% intervals; acceptance = statistical
compatibility with 0.05 (interval covers 0.05), plus exact orbit tests
passing. Simulation-table display rule (frozen now): significance markers in
simulation tables are replaced by rejection FRACTIONS across the 25
calibration replicates with exact Clopper-Pearson 95% intervals reported in
the calibration text; single-dataset tables (Face, USHCN) keep one star =
adjusted p <= 0.05 under the corrected procedure. Chance-adjusted quality
Qc = (Q - mu0)/(1 - mu0) (mu0 = closed-form null mean; negative allowed; not
a z-score) is added to the main tables in compact form alongside raw Q; the
format is fixed now, before any regenerated value is seen.

## 2. Empirical addition A: USHCN Interval-Euclidean lambda mechanism audit

- Data/stations/preprocessing/embeddings/seed: exactly the production USHCN
  analysis (frozen snapshot, standardize(), seed base 20260715 + 6000/6001,
  the stored six methods, K = 10). No refits; embeddings recomputed with the
  production stream only if required for byte-identical projections.
- Grid: lambda in seq(0, 1, by = 0.05) — 21 values, none omitted.
- Cells: EVERY prespecified USHCN Interval-Euclidean cell = all six methods
  (C-PCA, V-PCA, MR-PCA, SPCA, IMDS, Int-UMAP); the center-only baseline is
  lambda-invariant by construction and is recorded once as such.
- Outputs per (method, lambda): all six indices under the tie-audit
  protocol; plus fixed aggregate d^-/d^+ diagnostics computed over ALL
  station pairs (no anchor selection): mean and SD of d^-, of d^+, of the
  ratio d^-/d^+ (pairs with d^+ > 0), and the fraction of pairs with
  d^- = 0, in the high space and in each method's low space.
- Interpretation rule (fixed): descriptive mechanism/sensitivity audit; the
  manuscript's Int-UMAP sentence is qualified as implementation- and
  dissimilarity-specific REGARDLESS of outcome; the full profile ships in
  the outputs even if it weakens the narrative.
- Outputs: analysis/output/h_ushcn_lambda_profile.csv (+rds, meta),
  h_ushcn_lambda_components.csv (+rds, meta), sessionInfo, runtime.
- Script: analysis/scripts/h_ushcn_lambda.R, registered in make.R.

## 3. Empirical addition B: graded width-signal Scenario-I extension

- DGP: Scenario I's center geometry, n = 300, p = 5, three groups of 100,
  methods, K = 10, standardization, fitting procedure, d = 2 — unchanged.
- Width-signal levels (fixed before results, anchored to the current DGP
  scale): L1 weak: radii U(0.1,0.2) / U(0.3,0.4) / U(0.5,0.7);
  L2 moderate: U(0.1,0.2) / U(0.5,0.7) / U(1.5,1.9);
  L3 = the current Scenario-I reference: U(0.1,0.2) / U(1.0,1.2) / U(5,6).
  Rationale: geometric narrowing of the between-group radius separation
  toward overlap; L3 preserves continuity with the published scenario.
- Replications: 25 per level, seeds SEED_BASE_H = 20260719 + 8000 + r for
  r = 1..25 at L1, +8100+r at L2; L3 reuses the production Scenario-I
  replications 1..25 (seeds 20260715+1000+r) unchanged.
- Evaluations on the SAME embeddings: the four named dissimilarities
  (continuity with Scenario I), center-only Euclidean, and ONE
  endpoint/center-radius Euclidean evaluation (recorded as rank-equivalent
  to endpoint Euclidean up to the constant sqrt(2)). The algebraic
  identities are stated: endpoint = sqrt(2) x (c,r)-Euclidean
  (rank-identical); Wasserstein = (c, r/sqrt(3))-Euclidean (identical). No
  additional tuned radius weight (H.5).
- Primary estimand: per replication and method, the PAIRED difference
  Qc_interval - Qc_center of chance-adjusted T&C quality (Qc as in 1.4),
  where Qc_interval uses the Wasserstein evaluation (fixed choice: the only
  parameter-free interval dissimilarity with a proper-metric interpretation);
  summarized per level by mean and MC standard error over the 25
  replications. Secondary: the same contrast under the other
  dissimilarities and the (c,r) baseline; all six indices recorded.
- Reporting rule: all levels, all methods, all evaluations, failures, and
  null/negative findings are reported; convergence of interval-aware and
  center-only assessment at weak signal is an informative finding, not a
  failure. Presentation replaces/condenses current Scenario-I discussion
  (one compact figure or small table), no large new table.
- Outputs: analysis/output/h_graded_width_per_rep.csv, _summary.csv
  (+rds/meta), sessionInfo, runtime. Script:
  analysis/scripts/h_graded_width.R, registered in make.R.

## 4. Tie-policy verification (audit, not a study)

On the tie-affected production cells (Scenario I flagged cells, Scenario II
flagged cells, Scenario IV audited replications, USHCN daggered cells):
compare the existing 50 uniform strict refinements vs 500 refinements
(where runtime permits: USHCN + one scenario replication subset fixed as
S1 reps {1..10} flagged cells) vs a FIXED coupled-priority strict
refinement (priority = a seeded uniform permutation of object indices drawn
once per cell from seed 20260719+9000+cell_index, applied identically in
observed and permuted rank systems; label-equivariance verified by a
relabeling test). Midranks/fractional ranks are prohibited (H.6). Material
threshold (fixed): any displayed-value change at manuscript precision or
any conclusion change triggers prominent disclosure + three-provider
review; otherwise one manuscript sentence + archived table
(analysis/output/h_tie_policy_audit.csv). Script:
analysis/scripts/h_tie_policy_audit.R.

## 5. Explicit rejections (frozen)

No high-dimensional legacy-DR simulations; no third real dataset; no USHCN
multiple years or imputation models; no station map; no n = 5000 benchmark;
no cross-dissimilarity concordance study; no robust-scaling variants; no
within-data multi-seed UMAP study; no midrank/fractional co-ranking theory.

## 6. Page allocation (targets)

Additions: min-P proposition + algorithm (+0.4), taxonomy table (+0.5),
graded-width figure/table (+0.4), lambda-profile figure/sentences (+0.3),
factorial table (+0.3), algorithm box (+0.4), guidance/estimand sentences
(+0.5). Funded by: revision-history purge (-1.5), K-profile figure
consolidation Figs 1-4+6 -> 1 faceted figure (-2.0), Section 1/2.2
dedup + protocol dedup (-0.5), Face discussion condensation (-0.3).
Net target: <= 0 vs the 49-page submitted PDF; hard cap 48 pages.

## Amendment A1 (2026-07-19, adopted at the provider gate, BEFORE any
## gated result was computed)

Origin: three-provider gate outputs (Gemini: ACCEPT CONTRACT, 7 adversarial
attacks all held; Codex: AMEND-then-ACCEPT with the exact texts below, all
15 oracle self-checks TRUE). The frozen algorithm (1.1), family definition,
designs, levels, grids, seeds, estimands, display rules, and rejections are
UNCHANGED. A1 replaces one proof-wording paragraph, adds execution
guardrails, and corrects one rationale phrase. Original freeze hash
6612d73e...b204 remains on record; the post-A1 hash is appended to
`H_PREREGISTRATION.freeze`.

A1.1 (replaces the "Exchangeability with the actual sampler" paragraph in
1.2): Condition on the permutation-invariant orbit information, not on the
fully labelled observed correspondence. Under the complete
random-correspondence null, let the observed relative alignment Pi0 be
uniform on the permutation group. The code draws Delta_1..Delta_m
independently and uniformly with replacement, independently of Pi0, and its
effective rows correspond (up to composition convention) to Pi0,
Delta_1 Pi0, ..., Delta_m Pi0. Their joint probability at every ordered
(m+1)-tuple is |G|^-(m+1); hence the B effective alignments and joint
statistic rows are i.i.d. conditional on the orbit information and are
exchangeable. Identity and repeated draws are allowed and are handled by
inclusive ties. If ranks require auxiliary tie resolution, this conclusion
additionally requires the resolution variable to be independent and
label-equivariant in distribution, or a fixed priority/refinement to be
carried equivariantly with object identities through every permutation. A
fixed seed alone does not establish pathwise equivariance.

A1.2 (operational tie claim): validity in tie-affected inferential cells is
stated as unconditional over an independent uniform refinement (the current
`ties="random"` estimand), or via the contract-4 coupled object-priority
refinement carried through every permutation. Tied cells are not described
as conditionally exact for a fixed seed unless the relabelling/coupling
test passes. Tie-free cells need no qualification.

A1.3 (package identity guard): install source commit 6bb925f0 (or its
reviewed successor containing the min-P patch) into an isolated library and
set the library path to it before sourcing `_common.R` for any Prompt H
production run. Abort unless `find.package("QAIDR")` resolves inside that
library and `interval_data_from_dataSDA` exists in the namespace; the other
installed QAIDR 0.2.0 is not an admissible production artifact.

A1.4 (contract 2 execution): run the USHCN embedding reconstruction as a
standalone cold R process. After loading and standardizing the frozen
645-station snapshot, call `set.seed(20260715L + 6000L)` immediately before
one `run_methods_timed(xs, methods = METHODS)` call in canonical method
order; no warm-up or prior stochastic call is permitted. Every lambda-grid
cell that triggers the rank tie detector, especially lambda = 0, is
retained and evaluated through the 50-refinement tie-audit path with
`tie_seed = 20260715L + 6000L`, with tie indicators/spread saved. Compute
component summaries over exactly the unordered off-diagonal pairs i < j,
use `stats::sd`, and save both the total pair count and the eligible
d^+ > 0 count.

A1.5 (contract 3 execution): because production Scenario-I projections were
not stored, reconstruct L3 in a standalone cold R process, sequentially for
r = 1..25, before any L1/L2 fit, using the original data seed and
`set.seed(seed_r + 500000L)` immediately before the canonical full-method
call. Assert agreement of all pre-existing L3 method/evaluation/index cells
with `scenario1_per_rep` before computing the new endpoint/center-radius
row. For L1/L2, use the frozen data seed as seed_r and the unchanged
Scenario-I DR stream `set.seed(seed_r + 500000L)`; fit once and reuse that
embedding for every evaluation.

A1.6 (replaces the contract 3 rationale sentence): geometric narrowing of
the between-group radius separation relative to the common center geometry;
the L1 radius supports remain disjoint, so no support-overlap claim is
made. L3 preserves continuity with the published scenario.

A1.7 (contract 4 execution): define `cell_index` by lexicographic order of
the frozen cell key (scenario, replication, method, metric), before
inspecting any result. The seeded priority is attached to object
identities. Under every data relabelling or random-correspondence
permutation, permute the priority vector with those identities and use that
carried vector in both rank systems; never redraw or reassign priorities to
the new numeric positions. The relabelling test compares data and carried
priorities jointly.

## Amendment A2 (2026-07-19, author decision after the Stage 3 abort
## condition fired; see PROMPT_H_STAGE3_REPRODUCTION_FINDING.md)

The contract 1.4 promise "marginal p-values and descriptive indices must be
byte-unchanged" fired its abort condition: the stored production inference
outputs are NOT reproducible from the committed pipeline for Face and the
scenario calibrations (pre-existing defect introduced between the
2026-07-16 production runs and the 2026-07-18 CRAN-prep commit fcdcef1:
Face data-loader migration to dataSDA; epsilon-level numeric-environment
differences cascading through discrete tie/permutation streams). The
symmetric min-P change was exonerated by byte-identical cross-build probes.
USHCN reproduces its raw p-values exactly.

Author decision (2026-07-19): RE-BASELINE. The committed pipeline + QAIDR
0.3.0 + frozen data snapshots are the canonical R1 basis. Regenerated
outputs are accepted wholesale (raw + adjusted + affected descriptive
cells); every affected manuscript value/marker is updated in blue; the
re-baselining, the attribution evidence, and all star flips are disclosed
in the revision report and change ledger. Consequential edits:

- The regeneration assertion becomes: USHCN raw p-values byte-identical
  (verified); simulation/Face raw deltas RECORDED in full in the delta
  table instead of asserted zero.
- A1.5's L3 agreement assertion becomes: max |index delta| <= 1e-5 vs
  scenario1_per_rep (numerical-noise tolerance), with the observed maximum
  recorded; any delta above 1e-5 still aborts.
- The old-vs-new delta table now covers raw AND adjusted p-values.

## Amendment A3 (2026-07-19, before the tie-policy audit was run)

Scope discovery: the flagged-cell population is far larger than anticipated
(S1: 260 cells across 96 replications; S2: 498 cells across all 100
replications; S4: all 750 audited K=10 cells across 25 replications; USHCN:
4 daggered cells). Auditing every flagged cell requires ~221 embedding
refits (multi-day). Per the no-silent-caps rule the audit is scoped NOW,
before any audit result is computed, to fixed replication-index subsets
(outcome-independent): S1 replications 1..10 (with the 500-refinement arm,
per the original text), S2 replications 1..10, S4 replications 1..5 (all
their flagged K=10 cells), USHCN all 4 daggered cells (with the
500-refinement arm). Dropped counts are reported in the audit output.
Protocol note: the S4 production estimand shares one refinement stream
across methods within each (metric, resolution); the audit recomputes each
cell independently, so its 50-refinement arm is a same-protocol re-draw,
not a byte-reproduction of the stored S4 means; policies are compared
within the audit on identical streams. All other provisions of sec. 4 and
A1.7 are unchanged.

## 7. Abort conditions

Any oracle/orbit/FWER test failure; calibration audit incompatible with
0.05; exchangeability mismatch with the sampler; marginal p-values or
descriptive indices changed by regeneration; a frozen design violated or
run before this file's freeze hash was recorded; page cap exceeded without
author approval. On abort: stop, preserve evidence, report.
