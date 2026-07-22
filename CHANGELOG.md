# Version 1.2.1 (2026-07-23)

This submission-matched maintenance release does not change the frozen primary estimands, cohort membership, pooled effects, confidence intervals, prediction intervals, heterogeneity estimates, or multiplicity conclusions.

## Reporting and reproducibility corrections

- Reused the primary-analysis results for the `BASELINE_SD` scaling control and restored the complete 556-record scaling table.
- Renamed healthy-control scaling as an available-cohort sensitivity and documented its T48 cohort-composition limitation.
- Removed superseded evidence-grade fields and replaced ambiguous stability codes with explicit threshold-sensitivity labels.
- Clarified GSE106878 as a post-freeze directional-replication analysis; treatment-arm summaries remain descriptive and use the common full-cohort baseline standard deviation.
- Revised Figure 6 to display separate baseline-anchored T24 and T48 landmark estimates.
- Regenerated Figures 1–6 in one locked Python/Matplotlib environment and synchronized the figure-source workbook.
- Redirected clean-run figure and source-data outputs to `reproducibility_evidence/generated_v1.2.1` so reruns do not overwrite retained historical project snapshots.
- Tightened recorded-render-environment detection to require the documented Python and package versions in addition to macOS arm64 and Arial; scientifically exact runs in other render environments are labelled `PASS_PLATFORM_RENDER_VARIATION` only after all 12 source files and every export check pass.
- Replaced the submission-reference layouts for Figures 1, 4, 5, and 6 with the author-reviewed PowerPoint/PDF versions, retaining the editable PowerPoint sources. One supplied Figure 4 label was reconciled to the frozen source (`RPGRIP1` to `ZDHHC19`); no numerical source data changed.
- Added the frozen 34-candidate Stage 2 eligibility, registry, and reproducibility-grade audit and documented the result-independent boundary that yielded the eight evaluated protein-coding A2 signatures.
- Defined the three S10 class-agreement analyses as all-cohort, strict-never-used, and prespecified non-pilot-only. The last corresponds to the frozen internal code label `BLINDED_VALIDATION_ONLY` and does not imply formal masking; MAD scaling remains a separate material-sign-reversal check.
- Renamed the S15 pathway field to `T48_lowest_FDR_primary_pathway` and made exact FDR ties deterministic by pathway name; no pathway estimate, interval, FDR, or selected row changed.
- Added the 74-gene HPA mapping, historical aliases, frozen assignment rules, unresolved reasons, HPA v25.1 source URLs, and source-file SHA-256 to Supplementary Data.
- Relabelled Supplementary Figure S2 as all-cohort, primary-independent, strict-never-used, and MAD-scaling estimates and synchronized all three figure formats and the Python generator.
- Re-laid out Supplementary Figure S1 direct labels with deterministic offsets and fine leader lines, then added locked-environment pixel and source-CSV verification to the clean reproduction workflow; no point coordinates or source values changed.
- Normalized manuscript, Supplementary Information, workbook, and archive terminology to version 1.2.1.
- Synchronized the field rename and deterministic tie-break with the canonical analysis tree; immutable historical releases, clean-room checks, and frozen archival snapshots were deliberately left unchanged.

The prior v1.2.0 release remains part of the immutable version history.
