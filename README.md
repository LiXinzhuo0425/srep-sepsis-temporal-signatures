# Fixed blood transcriptomic signatures: code and data release v1.0.0

This release candidate accompanies *Fixed blood transcriptomic signatures show distinct temporal architectures across sepsis cohorts*. It contains the minimum derived dataset required to inspect the reported results, the source tables used for all main figures, portable figure-generation code, frozen analysis code and environment records, and independent reproducibility evidence.

## Evidence boundary

The release audits unchanged published formulas under repeated sampling. It does not create a classifier, validate clinical monitoring or treatment response, transport original assay thresholds, deconvolve cell abundance, or establish mechanism.

## Public inputs

The original expression data are available from NCBI GEO under GSE236713, GSE57065, GSE95233, GSE54514, GSE110487 and GSE8121. Source-study terms govern reuse. No direct identifiers or newly linked clinical data are included. Patient and sample identifiers in `data/derived_patient_level` are the pseudonymous identifiers already used in GEO.

## Reproduce the six main figures

1. Create a Python 3 environment.
2. Install `environment/requirements_release.txt`.
3. From this release root, run `python scripts/run_reproduction.py`.

The script uses only release-relative paths, regenerates Figures 1–6 plus Supplementary Figure S1, and compares all six generated PNGs with the versioned reference images. A successful run writes `reproduction_project/reproduction_result.json` with `PASS`.

## Reproduce the frozen analysis from GEO inputs

The original Stage 3 and Stage 4 scripts are preserved in `code/frozen_archive`. Environment-variable portable copies are in `code/portable_analysis`. Download the six public inputs into the layout documented by `analysis_config.json`, set `SEPSIS_SIGNATURE_ANALYSIS_ROOT` to that analysis directory, install the frozen Python/R dependencies, and run the Stage 3 and Stage 4 entry points listed below. The source-data release does not redistribute the large raw GEO files.

- `03_00_09_environment_lock/run_stage3.py`
- `03_00_09_environment_lock/verify_stage3.py`
- `04_environment/reproducibility_rerun_stage4.py`
- `04_environment/verify_stage4.py`

The completed clean-from-frozen-input reruns are retained in `reproducibility_evidence`. Stage 3 rebuilt raw-data preprocessing for GSE236713 and GSE110487 and reproduced scores, pairings, cohort effects and meta-analysis. Stage 4 forced pathway-cache regeneration and reproduced all 12 canonical outputs.

## Directory map

- `data/derived_patient_level`: de-identified derived score, paired-change, gene-contribution and pathway-change matrices.
- `data/source_tables`: reviewer-facing source tables used by the figure script.
- `tables`: Supplementary Data, source-data workbook and the two main tables.
- `reference_outputs`: versioned main-figure files used for comparison.
- `code/portable_analysis`: environment-variable portable analysis scripts.
- `code/frozen_archive`: unchanged historical audit copy of the frozen analysis code.
- `reproducibility_evidence`: run logs, numerical comparisons and verification summaries.

## Version and citation

Local release identifier: `srep-sepsis-temporal-signatures-v1.0.0`. The repository URL and persistent archive identifier remain intentionally unset until the authors publish the release. Do not cite a DOI until one has been issued by the selected repository.
