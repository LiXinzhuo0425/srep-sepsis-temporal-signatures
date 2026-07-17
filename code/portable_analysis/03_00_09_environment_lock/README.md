# Stage 3 frozen environment

This directory contains the pre-unblinding analysis implementation.

- Statistical and plotting backend: Python 3.9.6.
- RNA-seq preprocessing only: R 4.5.2 with project-local DESeq2 1.50.2.
- All figure drawing, previewing, exporting and visual QA is Python-only.
- Random seed: 20260716; patient bootstrap replicates: 2,000.
- `signatures.py` contains the eight fixed A2 functions and frozen direction coefficients.
- `run_stage3.py` performs sequential scoring, pairing, cohort estimates, meta-analysis and sensitivities.
- `make_figures.py` creates all publication figures from machine-readable source tables.
- `verify_stage3.py` and `reproducibility_rerun.py` are independent verification gates.

The validation cohorts cannot be scored until `STAGE3_ANALYSIS_UNLOCKED.md` exists. `run_stage3.py unblind` also enforces the frozen validation order.

## Commands

```text
python3 run_stage3.py pilot-smoke
python3 run_stage3.py prepare-rnaseq
python3 run_stage3.py unblind GSE95233
python3 run_stage3.py unblind GSE110487
python3 run_stage3.py unblind GSE236713
python3 run_stage3.py unblind GSE8121
python3 run_stage3.py analyze
python3 make_figures.py
python3 verify_stage3.py
python3 reproducibility_rerun.py
python3 finalize_stage3.py
```
