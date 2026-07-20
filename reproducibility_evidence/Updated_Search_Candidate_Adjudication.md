# Updated-search candidate-cohort adjudication

Decision date: 20 July 2026

The six-cohort primary analysis remains frozen. Candidate records identified by the updated search were assessed against the same outward eligibility requirements: a human sepsis case population, recoverable patient-level repeated measurements, a source-defined baseline and T24 or T48 window, openly reviewer-accessible expression values, and complete coverage of the 74 component genes required by all eight formulas. No candidate was added retrospectively to the primary meta-analysis.

## GSE106878 — INCLUDE_POST_FREEZE_CONFIRMATION

GSE106878 provides paired pretreatment and 24-hour whole-blood microarray profiles for 47 patients with septic shock from the CORTICUS trial. The public normalized matrix and platform annotation recover all 74 genes after the locked historical-symbol map. All eight scores were independently reconstructed with the v1.1.2 formulas. The overall T24 estimates, with 2,000 paired bootstrap replicates and seed 20260720, were:

- SIG001: -0.972464 delta Z (95% CI -1.138908 to -0.812034)
- SIG002: -1.483551 (-1.765863 to -1.201001)
- SIG003: -1.116821 (-1.367120 to -0.883496)
- SIG004: +0.917188 (+0.706763 to +1.113164)
- SIG022: +0.182822 (+0.028462 to +0.334997)
- SIG023: -0.207620 (-0.366096 to -0.034997)
- SIG033: +0.207774 (-0.032137 to +0.445220)
- SIG034: +0.186964 (-0.033315 to +0.432814)

The four main longitudinal directions were concordant in the placebo and hydrocortisone strata. Because this is a randomized treatment cohort identified after the primary freeze, these results are reported as a standalone post-freeze confirmation and are not pooled with the six frozen study families.

## E-MEXP-3850 — EXCLUDE

The record describes five children sampled repeatedly through 48 hours, but the current public deposit contains only IDF and SDRF metadata. Neither raw intensity files nor a processed expression matrix can be downloaded. Formula implementation and 74-gene coverage therefore cannot be verified reproducibly.

## PRJEB111201 — EXCLUDE

The study reports adult whole-blood RNA sequencing on days 1, 3 and 7, including ten complete day-1/day-3 pairs. Raw FASTQ files were public, but no processed subject-by-gene matrix or study-matched quantification workflow was available. Inclusion would require introducing and validating a new raw-read alignment and quantification branch after analysis freeze, rather than applying the frozen gene-level input workflow. The cohort was therefore documented but not analyzed.

## phs003608.v1.p1 — EXCLUDE

The Ghana cohort has eligible longitudinal sampling at 0, 6, 24, 48 and 72 hours. Subject-level bulk RNA-seq data are controlled-access in dbGaP. All seven associated Zenodo source-data workbooks were inspected; they contain figure source data, group-level differential-expression tables, classifier summaries and single-cell records, but not the subject-by-gene bulk matrix needed to calculate paired values for all eight signatures. The record is therefore excluded from the current public-data workflow.

## Change-control conclusion

The updated search changes neither the membership nor the estimand of the frozen primary analysis. GSE106878 adds a clearly labelled, standalone confirmation analysis. The other three records have clinically relevant longitudinal designs but fail the public, reproducible expression-input requirement for this submission. No wet-laboratory experiment or new exploratory endpoint was introduced.
