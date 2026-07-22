# v1.2.0 (2026-07-22)

Version 1.2.0 is a non-overwriting Stage 4 correction release. It supersedes v1.1.2 for gene-contribution architecture summaries while preserving the v1.1.2 record for provenance.

## Corrected analysis-set consistency

- A cross-material review found that the superseded Stage 4 architecture summaries for SIG002, SIG003, and SIG004 used five eligible T48 cohorts even though their primary score syntheses excluded GSE54514 and used four independent cohorts.
- Stage 4 was rerun using the exact signature-specific primary-independent cohort IDs, cohort count, and patient count from the Stage 3 primary row for every architecture component.
- The affected architecture set is `GSE110487;GSE57065;GSE8121;GSE95233` with k = 4 and n = 118.
- Dominance, cancellation, leading-gene identity, material-gene count, direction consistency, baseline labels, and all threshold-perturbation scenarios were recomputed from this single coherent set.

## Result impact

- SIG002 remains consistent multigene drift: dominance 0.510, cancellation 0.490, modal leading-gene agreement 1.00.
- SIG003 remains single-gene-dominant: dominance 0.679, cancellation 0.000, modal agreement 0.75; its label changes only in the +20% threshold scenario.
- SIG004 changes from single-gene-dominant to consistent multigene drift: dominance 0.589, cancellation 0.345, modal agreement 0.75.
- Three labels now change in at least two threshold scenarios: SIG002, SIG004, and SIG022.
- Stage 3 pooled effects, confidence intervals, prediction intervals, heterogeneity estimates, standard and modified Knapp–Hartung results, and the six primary standard-Hartung–Knapp findings are unchanged.

## Reproducibility guard

- Added an exact cohort-ID/count/patient-count regression check for all eight signatures.
- Added corrected machine-readable S11, S15, and S17 tables and matching Figure 4/Figure 6 source data.
- Regenerated submission figures and the reviewer-facing Table 2 from the corrected outputs.
