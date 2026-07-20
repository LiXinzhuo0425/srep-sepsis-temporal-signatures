# Stage 4 figure contracts

Backend is frozen to Python/matplotlib. SVG with editable text is the primary export; PDF, 300-dpi PNG and 600-dpi TIFF are secondary exports. All figures use the same semantic palette and are generated from frozen CSV/Parquet sources.

## Figure 4A — total drift and leading gene contributions

- Conclusion: T48 drift differs across signatures, and the dominant mathematical contributors are signature-specific.
- Evidence: Stage 3 primary-independent random-effects delta-Z estimates plus Stage 4 primary-independent gene-contribution meta-estimates.
- Encoding: total effect with 95% CI on the left; the two largest absolute gene contributions on the right. Gene contributions are not stacked into the meta-analytic total because gene-specific random-effects weights differ.
- Review risk: avoid implying causal drivers; label contributions as exact mathematical attribution under the fixed formula.

## Figure 4B — dominance and cancellation architecture

- Conclusion: apparent score stability can arise from low component change or from internally cancelling gene changes; some drift is single-gene dominated.
- Evidence: frozen T48 median patient dominance ratio, cancellation index and absolute contribution sum.
- Encoding: x=dominance ratio, y=cancellation index, size=absolute contribution sum, colour=pre-coded architecture, direct signature labels.
- Review risk: architecture labels are rule-based descriptors, not a best-to-worst ranking.

## Figure 5 — gene, cell-source and pathway integration

- Conclusion: gene-level mathematical contributions have heterogeneous potential blood-cell sources, and only limited pathway coupling is consistently replicated.
- Evidence: HPA v25.1 blood-cell RNA annotations, primary-independent gene-contribution estimates, and independent-only T48 pathway coupling.
- Encoding: panel a shows the leading gene, contribution and broad lineage; panel b shows T48 Spearman ρ for nine prespecified Hallmark pathways, with an outline for FDR < 0.05.
- Review risk: bulk RNA cannot distinguish abundance shifts from within-cell transcription; HPA assignments are potential sources only.

## Figure 6 — context-of-use evidence matrix

- Conclusion: the study supports temporal-score description, offers limited biological anchoring, and does not support monitoring, prognosis or treatment decisions.
- Evidence: Stage 3 PASS plus Stage 4 module gates and claim boundaries.
- Encoding: scenario rows and evidence-domain columns with supported/limited/not-supported states; short boundary statements are printed alongside.
- Review risk: do not use “clinically useful”, “validated for monitoring” or threshold-oriented wording.

## Supplementary figures

- Eight T48 gene-contribution forest plots show every component gene, 95% CI and zero reference.
- The standalone pathway heatmap shows all 8×9 T48 primary couplings, not only FDR-positive cells.
- The standalone gene–cell-source map shows all signature genes with available HPA lineage evidence; unresolved genes remain explicit.
