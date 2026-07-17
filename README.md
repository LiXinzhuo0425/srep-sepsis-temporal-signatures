# Fixed blood transcriptomic signatures: code and data release v1.0.1

This release accompanies *Fixed blood transcriptomic signatures show distinct temporal architectures across sepsis cohorts*. It contains the minimum derived dataset required to inspect the reported results, the source tables used for all main figures, portable figure-generation code, frozen analysis code and environment records, and scripted clean-environment reproducibility evidence.

## Evidence boundary

The release evaluates unchanged published formulas under repeated sampling. It does not create a classifier, validate clinical monitoring or treatment response, transport original assay thresholds, deconvolve cell abundance, or establish mechanism.

## Public inputs

The original expression data are available from NCBI GEO under GSE236713, GSE57065, GSE95233, GSE54514, GSE110487 and GSE8121. Source-study terms govern reuse. No direct identifiers or newly linked clinical data are included. Patient and sample identifiers in `data/derived_patient_level` are the pseudonymous identifiers already used in GEO.

## Reproduce the six main figures

1. Create a Python 3 environment.
2. Install `environment/requirements_release.txt`.
3. From this release root, run `python scripts/run_reproduction.py`.

The script uses only release-relative paths, regenerates Figures 1–6 plus Supplementary Figure S1, and compares all six generated PNGs with the versioned reference images. A successful run writes `reproduction_project/reproduction_result.json` with `PASS`.

## Reproduce the analysis from GEO inputs

Environment-variable portable copies of the Stage 3 and Stage 4 scripts are in `code/portable_analysis`. Download the six public inputs into the layout documented by `analysis_config.json`, set `SEPSIS_SIGNATURE_ANALYSIS_ROOT` to that analysis directory, install the frozen Python/R dependencies, and run the Stage 3 and Stage 4 entry points listed below. The source-data release does not redistribute the large raw GEO files. The exact historical path-bound source archive remains permanently available in the immutable v1.0.0 release; it was not duplicated in v1.0.1 so that this version contains no local usernames or absolute workstation paths.

- `03_00_09_environment_lock/run_stage3.py`
- `03_00_09_environment_lock/verify_stage3.py`
- `04_environment/reproducibility_rerun_stage4.py`
- `04_environment/verify_stage4.py`

The completed clean-from-frozen-input reruns are retained in `reproducibility_evidence`. Stage 3 rebuilt raw-data preprocessing for GSE236713 and GSE110487 and reproduced scores, pairings, cohort effects and meta-analysis. Stage 4 forced pathway-cache regeneration and reproduced all 12 canonical outputs. These are scripted verification records and are not represented as an independent human review.

## Version 1.0.1 changes

- Corrected the SIG034 source attribution to Zheng et al. (2021); the DOI remains `10.1016/j.immuni.2021.03.002`.
- Added the prespecified reporting rule for meta-analyses with fewer than three contributing cohorts. Rows with one cohort are reported as cohort effects, and two-cohort syntheses are descriptive; formal Hartung–Knapp intervals and P values are suppressed in the outward-facing tables.
- Replaced outward-facing T1–T4 labels with T24, T48, T72 and Day 5, while retaining the internal mapping needed to reproduce the code.
- Updated the Supplementary Data workbook and figure exports for the Stage 7 submission package. Figure 1 now uses a hatch as a redundant grayscale cue, and all TIFFs are RGB without transparency.
- No primary effect estimate, multiplicity decision, heterogeneity estimate, prediction interval, gene-architecture classification or pathway conclusion changed.

Raw historical outputs, including the small-cohort rows produced by the original frozen program, remain available in `code/portable_analysis` and the immutable v1.0.0 archive; the new reporting rule affects presentation, not the stored analysis history.

## Directory map

- `data/derived_patient_level`: de-identified derived score, paired-change, gene-contribution and pathway-change matrices.
- `data/source_tables`: reviewer-facing source tables used by the figure script.
- `tables`: Supplementary Data, source-data workbook and the two main tables.
- `reference_outputs`: versioned main-figure files used for comparison.
- `code/portable_analysis`: environment-variable portable analysis scripts.
- `reproducibility_evidence`: run logs, numerical comparisons and verification summaries.

## Version and citation

Release identifier: `srep-sepsis-temporal-signatures-v1.0.1`. The public repository is `https://github.com/LiXinzhuo0425/srep-sepsis-temporal-signatures`. The version DOI is `https://doi.org/10.5281/zenodo.21417810`; the concept DOI is `https://doi.org/10.5281/zenodo.21415496`.

## License

Original release code is available under the MIT License. Derived data, documentation, and original figure/source-data materials are available under CC BY 4.0. Third-party source data and dependencies retain their original terms. See `LICENSE_NOTICE.md`, `LICENSE_CODE`, and `LICENSE_DATA`.
