# Fixed blood transcriptomic signatures: version 1.1.2

This release accompanies *Fixed blood transcriptomic signatures show distinct longitudinal behavior and gene-contribution patterns across sepsis cohorts*. It contains the corrected v1.1.2 analysis outputs, patient-level derived data, source tables for Figures 1–6, portable figure-generation code, environment records, modified Knapp–Hartung sensitivity results, post-freeze confirmation results, and scripted clean-environment verification.

The public v1.0.1 archive is retained for provenance but is superseded for scientific use because its SIG001 implementation omitted the published 5/6 coefficient on the down-regulated module. Version 1.1.1 contained the corrected analyses but retained a redundant legacy Supplementary Figure S6 file in archived engineering materials. Version 1.1.2 removes that unused file and normalizes the submission packaging; no scientific result, numerical value, or interpretation changed between v1.1.1 and v1.1.2. Version 1.1.2 is the corrected submission release.

## Evidence boundary

The study evaluates published blood transcriptomic formulas under repeated sampling. It does not create a classifier, validate clinical monitoring or treatment response, transport original assay thresholds, deconvolve cell abundance, or establish biological mechanism.

## Public inputs

Original expression data are available from NCBI GEO under GSE236713, GSE57065, GSE95233, GSE54514, GSE110487, GSE8121 and GSE106878. The first six cohorts constituted the frozen primary analysis; GSE106878 was used only as a standalone post-freeze confirmation cohort. Source-study terms govern reuse. No direct identifiers or newly linked clinical data are included. Patient and sample identifiers in `data/derived_patient_level` are the pseudonymous identifiers used in the public deposits.

## Reproduce Figures 1–6

1. Use Python 3.12.
2. Install `environment/requirements_release.txt`.
3. From the release root, run `python scripts/run_reproduction.py`.

The script uses only release-relative paths and regenerates Figures 1–6. Scientific content is verified through byte-identical comparison of all 12 source-data CSVs; author-reviewed presentation-only layout differences are recorded explicitly rather than treated as numerical changes. The verifier also checks that TIFF files are RGB, 600 dpi and LZW-compressed, SVG files retain editable text, and PDFs exist. A successful run writes `reproducibility_evidence/reproduction_result_v1.1.2.json` with `PASS`.

## Reproduce the corrected analysis from GEO inputs

Portable Stage 3 and Stage 4 scripts are in `code/portable_analysis`. Download the six public inputs into the layout documented by `analysis_config.json`, set `SEPSIS_SIGNATURE_ANALYSIS_ROOT` to that analysis directory, install the recorded Python and R dependencies, and run the Stage 3 and Stage 4 entry points below. The release does not redistribute the large raw GEO files.

- `03_00_09_environment_lock/run_stage3.py`
- `03_00_09_environment_lock/verify_stage3.py`
- `04_environment/reproducibility_rerun_stage4.py`
- `04_environment/verify_stage4.py`

Version 1.1.2 replaces the superseded SIG001 implementation with the published formula `(up module) - (5/6) × (down module)` and regenerates every downstream score, paired change, cohort effect, meta-analysis, gene contribution, pathway result, table and figure that depends on it. The other seven formulas are unchanged. The complete formula registry and sensitivity evidence are retained in `tables` and `reproducibility_evidence`.

## Version 1.1.2 changes

- Corrected the SIG001 Sepsis MetaScore coefficient; v1.0.1 is scientifically superseded.
- Regenerated corrected patient-level scores, paired changes, gene contributions, cohort estimates, meta-analyses and source data.
- Added the complete signature implementation registry and independent formula checks.
- Added architecture-threshold sensitivity evidence and marked threshold-sensitive categorical labels.
- Replaced Figure 5 panel-b hatching with distinct solid colours and thin white segment boundaries; the FDR marker remains independently visible in panel c.
- Replaced residual project-governance wording in Figure 1 with reader-facing study-design wording; no data or numerical output changed.
- Refreshed the outward-facing public-dataset search date to 20 July 2026 while preserving the dated underlying search records.
- Expanded Table 2 with T24/T48 cohort counts and explicit T48 heterogeneity and gene-contribution metrics.
- Locked the figure environment to Python 3.12 and Matplotlib 3.11.0; the clean run reproduces all six submission figures and 12 source-data files.
- Added modified Knapp–Hartung sensitivity using q*=max(1,qHK). Four of the six primary standard-Hartung–Knapp combinations retained multiplicity-adjusted evidence; the primary analysis itself was not replaced.
- Formally adjudicated four updated-search candidate cohorts. GSE106878 is reported as a standalone 47-patient post-freeze T24 confirmation analysis; E-MEXP-3850, PRJEB111201 and phs003608 were excluded because a reproducible public expression input compatible with the frozen workflow was unavailable.
- Added exact Shapley/Owen attribution definitions, exact formula-source locations, and the public S25-S26 machine-readable tables.
- Corrected the public `data/source_tables` copies so that every Stage 3 and Stage 4 source table is byte-identical to the corrected reproduction project.
- Removed the redundant legacy Supplementary Figure S6 from archived engineering materials and normalized final release filenames without changing scientific content.
- Standardized reader-facing statistical notation to `ρ`, `α`, `β`, `γ`, `τ²`, `I²`, `ΔZ` and `φ`; machine-readable ASCII variable names and canonical MSigDB identifiers remain unchanged for reproducibility.

## Directory map

- `data/derived_patient_level`: corrected de-identified score, paired-change, gene-contribution and pathway-change matrices.
- `data/figure_source_data`: versioned CSV source data for Figures 1–6.
- `tables`: Supplementary Data, figure-source workbook, implementation registry and main Tables 1–2.
- `reference_outputs/main_figures_v1.1.2`: submission-reference figures used by the verifier.
- `code/portable_analysis`: release-relative corrected analysis scripts.
- `reproducibility_evidence`: numerical comparisons, clean-run results and sensitivity evidence.

## Version and citation

Release identifier: `srep-sepsis-temporal-signatures-v1.1.2`.

Public repository: `https://github.com/LiXinzhuo0425/srep-sepsis-temporal-signatures`.

The version-specific DOI is recorded in the GitHub release, Zenodo record and release manifest. Do not cite the superseded v1.0.1 version DOI as the corrected analysis.

## License

Original release code is available under the MIT License. Derived data, documentation and original figure/source-data materials are available under CC BY 4.0. Third-party source data and dependencies retain their original terms. See `LICENSE_NOTICE.md`, `LICENSE_CODE` and `LICENSE_DATA`.
