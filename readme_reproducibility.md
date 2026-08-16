# Reproducing the QAIDR manuscript results

Every manuscript table, figure, numerical claim, and significance marker is
produced by a named, seeded script under `analysis/scripts/`. Machine-readable
CSV, RDS, and JSON outputs are saved under `analysis/output/` before LaTeX
rendering; tables are not hand-transcribed.

The repository includes the curated scripts, processed USHCN snapshot,
provenance records, and final outputs. The approximately 2.9 GB raw NOAA
download is deliberately excluded. Its official source, version, filtering
rules, and SHA-256 records are retained by the acquisition pipeline.

## Environment

The manuscript run used:

* R 4.6.1 (Windows UCRT)
* QAIDR 0.3.0
* dataSDA 0.2.6 (the `face.iGAP` object is loaded at run time)
* symbolicDA 0.6-2
* RSDA 3.2.5
* umap 0.2.10.0
* coRanking 0.2.5 (verification only)

Per-stage `sessionInfo()` files are stored in `analysis/output/`. JSON
sidecars record package versions, input hashes, output hashes, timestamps,
and elapsed times. `analysis/provenance/external_data_manifest.json` records
serialized SHA-256 digests for the `dataSDA` Face and Cars objects.

RSDA 3.2.5 has an important implementation caveat: `sym.umap()` ignores
arguments passed through `...` and overwrites `n_components` with the expanded
vertex dimension. QAIDR therefore passes a configuration object and records
the requested and effective settings.

## Running the pipeline

From the package repository root on PowerShell:

```powershell
$env:QAIDR_RUN_MODE = "quick"
Rscript analysis/make.R

$env:QAIDR_RUN_MODE = "full"
Rscript analysis/make.R
```

On a POSIX shell:

```sh
QAIDR_RUN_MODE=quick Rscript analysis/make.R
QAIDR_RUN_MODE=full Rscript analysis/make.R
```

Quick mode performs a validation pass and may reuse only full artifacts whose
contents and recorded hashes validate. Full mode recomputes every stage. A
sequential full rebuild requires roughly 70 or more hours on the reference
laptop; `analysis/make.R` prints the expected duration before each stage and
stops at the first failure.

## Manuscript object map

| Object | Main script | Primary outputs |
|:--|:--|:--|
| Scenario I | `scenario1.R` | `scenario1_summary.csv`, calibration and tie files |
| Scenario II | `scenario2.R` | `scenario2_summary.csv`, factorial, calibration and tie files |
| Scenario III | `scenario3.R`, `b2_fix_scenario3_aggregation.R` | `scenario3_per_perturbation.csv`, `scenario3_summary.csv` |
| Scenario IV | `b2_scenario4.R` | `scenario4_summary.csv`, calibration and baseline contrasts |
| Face table and figures | `face.R` | `face_indices.csv`, p-values, K-profiles and PDFs |
| USHCN analysis | `midscale_pilot_ghcn.R`, `midscale_ghcn.R` | processed snapshot, indices, p-values, runtime and K-profiles |
| Scalability | `benchmark.R` | `scalability_benchmark.csv` |
| K robustness | `b2_k_profiles.R`, `b2_k_profiles_merge.R` | profile summaries and Kendall stability |
| Lambda/nu sensitivity | `b2_sensitivity.R` | exact Face breakpoints, grids and endpoint audit |
| C-PCA mechanism | `b2_cpca_diagnostics.R` | `b2_cpca_diagnostics.csv` |
| Symmetric min-P regeneration | `regen_calibration.R`, `h_regen_s4_calibration.R` | regenerated calibration and adjusted-p outputs |
| Symmetric min-P calibration audit | `h_audit_v3c.R` | `h_audit_v3c_familywise.csv` |
| Old/new inference audit | `h_delta_table.R` | `h_delta_table.csv` |
| USHCN lambda mechanism | `h_ushcn_lambda.R` | lambda profiles, components, runtime and figure |
| Graded width signal | `h_graded_width.R` | per-replication results, summaries and runtime |
| Tie-policy audit | `h_tie_policy_audit.R` | `h_tie_policy_audit.csv` |
| Calibration gates | `stage4_gate2.R` | `stage4_gates.json`, per-replication calibration |
| LaTeX tables | `tables_to_tex.R` | `analysis/output/tex/*.tex` |
| Final verification | `verify_manuscript_numbers.R` | exit status and figure manifest |

The later `e10b_*` scripts and outputs are the final m = 999 multiplicity
audit used to confirm the published significance markers.

## Statistical protocol

* The T&C, MRRE, and LCMC families follow Lee and Verleysen (2009). MRRE uses
  the printed rank weights, and `B_LC` is the bounded `B_NX = U_X - U_N`.
  Positive behavior values mean intrusion dominance.
* Ranks exclude self-distances. Primary cells are tie-free; unavoidable
  structural ties use 50 seeded uniform resolutions, are flagged, and report
  their maximum resolution spread.
* Permutations apply one uniform label bijection to both endpoints of the
  embedded ranks. P-values use `(c + 1)/(m + 1)` with `m = 999`. QAIDR 0.3.0
  uses the symmetric single-step randomization min-P construction documented
  in `NEWS.md`: the observed row and all randomization rows share one
  denominator and inclusive comparisons. It provides finite-sample weak FWER
  control at attainable levels under the complete random-correspondence null
  for each prespecified index family across metrics at fixed K.
* The simulation study uses 100 independently seeded datasets per scenario
  and reports standard deviations and Monte Carlo standard errors. Calibration
  is evaluated on a prespecified subset of replications.
* Center-only Euclidean baselines measure what interval widths add beyond the
  midpoints.

## Provenance boundary

The saved numerical outputs and their sidecars are the canonical QAIDR 0.3.0
manuscript run. Version 0.3.0 changes familywise-adjusted p-values; all affected
inference outputs were regenerated, and `h_delta_table.csv` records the audit
against the previous implementation. The Face object in `dataSDA` was verified
to reproduce the former centers, radii, and labels exactly. Protected
rejected-submission files and raw NOAA downloads are not part of this
repository release.
