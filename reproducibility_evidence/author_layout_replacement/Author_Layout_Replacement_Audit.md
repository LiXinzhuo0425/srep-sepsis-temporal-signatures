# Author layout replacement audit

## Scope

Figures 1, 4, 5, and 6 were replaced with the author-adjusted layouts. The supplied PowerPoint and PDF pairs were reviewed slide-by-slide. The PDF files were treated as the appearance masters, and the PowerPoint files were retained as editable author sources.

## Controlled content decision

- Figure 4 contained one non-layout discrepancy: the first SIG001 label was `RPGRIP1` in the supplied author file, whereas the frozen source ranks `CEACAM1` and `ZDHHC19` as the two largest absolute pooled T48 contributions. The replacement therefore corrects `RPGRIP1` to `ZDHHC19`.
- Figure 1 contains two non-numeric wording refinements: `exact gene decomposition` and `each study's locked baseline definition`.
- No numerical figure-source table, statistical estimate, confidence interval, patient count, signature set, or manuscript result was changed.

## Verification

- All four PowerPoint files passed slide overflow checks.
- The corrected Figure 4 PowerPoint passed the template-fidelity check with zero issues.
- Numerical-token multisets in the prior and replacement PDFs agree for Figures 1, 4, 5, and 6.
- Figure 4's `ZDHHC19` label is present, `RPGRIP1` is absent, and the label agrees with the frozen `Figure_4A_source_data.csv` ranking.
- Figure 6's 16 T24/T48 marker locations were measured from the rendered replacement and compared with `Figure_6_source_data.csv`; the maximum inferred-value difference was 0.0173 on the plotted scale, within the 0.03 rendering tolerance.
- Color and grayscale previews were visually inspected at full size. No clipping, overlap, missing glyph, or unreadable panel was identified.
- The four replaced submission TIFFs are RGB, 600 dpi, LZW-compressed, and 4,323 pixels wide (183 mm at 600 dpi); Figures 2 and 3 retain their code-generated dimensions.
- PDFs are one-page vector masters normalized to 183 mm width; fonts are embedded. SVG files retain live text elements.

## Deliverable status

Local package replacement passed the recorded checks. Public-release status is recorded in `release_manifest.json`; exact-version author approval and inspection of the journal-generated merged PDF remain separate external gates.
