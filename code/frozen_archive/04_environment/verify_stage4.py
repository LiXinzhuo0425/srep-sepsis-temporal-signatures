#!/usr/bin/env python3
"""Independent Stage 4 verification using frozen outputs and source records."""
from __future__ import annotations

import hashlib
import json
import math
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/Users/felix/Documents/New project/longitudinal_stage3")
SRC = ROOT / "04_source_data"
OUT = SRC / "04_22_independent_extension_verification.csv"
SUMMARY = SRC / "04_22_independent_extension_verification_summary.json"
TOL = 1e-8
rows: list[dict] = []


def record(check_id, category, detail, observed, expected, status, tolerance=""):
    rows.append(
        {
            "check_id": check_id,
            "category": category,
            "detail": detail,
            "observed": observed,
            "expected": expected,
            "tolerance": tolerance,
            "status": status,
        }
    )


def pass_if(value: bool) -> str:
    return "PASS" if bool(value) else "FAIL"


# 1. Exact reconstruction, independently summing the gene-level contribution table.
gene = pd.read_parquet(ROOT / "04_04_gene_level_longitudinal_dataset.parquet")
qc = pd.read_csv(SRC / "04_05_decomposition_reconstruction_qc.csv")
keys = ["dataset", "patient_id", "signature_id", "time_window"]
summed = gene.groupby(keys, as_index=False)["standardized_contribution"].sum().rename(
    columns={"standardized_contribution": "independent_reconstructed_delta_z"}
)
check = summed.merge(qc[keys + ["stage3_delta_z"]], on=keys, how="inner")
check["absolute_error"] = (check["independent_reconstructed_delta_z"] - check["stage3_delta_z"]).abs()
for sig, d in check.sort_values(keys).groupby("signature_id"):
    sampled = d.head(10)
    max_error = float(sampled["absolute_error"].max())
    record(
        f"RECON_{sig}",
        "Exact decomposition",
        "Independent sum for 10 deterministic patient-window records",
        f"n={len(sampled)}; max_error={max_error:.3e}",
        "10 records; max_error <= 1e-8",
        pass_if(len(sampled) == 10 and max_error <= TOL),
        "1e-8",
    )
global_max = float(check["absolute_error"].max())
record("RECON_ALL", "Exact decomposition", "All gene-level records reconstruct Stage 3 delta Z", f"n={len(check)}; max_error={global_max:.3e}", "n=3344; max_error <= 1e-8", pass_if(len(check) == 3344 and global_max <= TOL), "1e-8")


# 2. Formula methods and gene roster agreement.
coef = pd.read_csv(SRC / "04_03_signature_gene_coefficients.csv")
config = json.loads((ROOT / "04_environment/stage4_config.json").read_text())
for sig in config["signatures"]:
    expected_method = config["decomposition"][sig]
    coef_genes = set(coef.loc[coef.signature_id == sig, "gene"])
    data_genes = set(gene.loc[gene.signature_id == sig, "gene"])
    methods = set(gene.loc[gene.signature_id == sig, "attribution_method"])
    ok = coef_genes == data_genes and methods == {expected_method}
    record(f"FORMULA_{sig}", "Formula specification", "Frozen gene roster and attribution method", f"genes={len(data_genes)}; methods={sorted(methods)}", f"genes={len(coef_genes)}; method={expected_method}", pass_if(ok))


# 3. Protected source hashes.
lock = json.loads((ROOT / "04_environment/pathway_analysis_lock.json").read_text())
hash_files = {
    "msigdb_hallmark": ROOT / "04_pathway_data/gene_sets/h.all.v2026.1.Hs.symbols.gmt",
    "reactome_gmt_zip": ROOT / "04_pathway_data/gene_sets/ReactomePathways.gmt.zip",
    "gpl570": ROOT / "04_pathway_data/platform_annotations/GPL570.annot.gz",
    "gpl6947": ROOT / "04_pathway_data/platform_annotations/GPL6947.annot.gz",
    "gpl17077": ROOT / "04_pathway_data/platform_annotations/GPL17077.txt",
    "gencode_v25": ROOT / "04_pathway_data/platform_annotations/gencode.v25.annotation.gtf.gz",
}
for key, file in hash_files.items():
    digest = hashlib.sha256(file.read_bytes()).hexdigest()
    expected = lock["source_hashes"][key]
    record(f"HASH_{key.upper()}", "Source integrity", str(file.relative_to(ROOT)), digest, expected, pass_if(digest == expected))


# 4. Pathway gate, compatibility and prespecified window support.
coverage = pd.read_csv(SRC / "04_09_pathway_gene_coverage.csv")
primary = coverage[coverage.tier == "PRESET_PRIMARY"].copy()
support = primary.assign(ok=primary.compatible.eq("YES")).groupby("pathway")["ok"].sum()
record("PATHWAY_PRIMARY_SUPPORT", "Pathway feasibility", "Minimum compatible cohorts across nine primary Hallmarks", int(support.min()), ">=4", pass_if(len(support) == 9 and support.min() >= 4))
noncompatible = primary[primary.compatible.ne("YES")][["dataset", "pathway", "coverage_fraction"]]
observed_noncompatible = "; ".join(f"{r.dataset}:{r.pathway}:{r.coverage_fraction:.3f}" for r in noncompatible.itertuples())
record("PATHWAY_EXCLUSIONS", "Pathway feasibility", "Noncompatible primary dataset-pathway combinations", observed_noncompatible, "GSE54514:HALLMARK_COAGULATION:0.681", pass_if(len(noncompatible) == 1 and noncompatible.iloc[0].dataset == "GSE54514" and noncompatible.iloc[0].pathway == "HALLMARK_COAGULATION"))
windows = pd.read_csv(SRC / "04_09_pathway_time_window_support.csv")
window_map = windows.groupby("time_window")["dataset"].nunique().to_dict()
record("PATHWAY_WINDOWS", "Pathway feasibility", "Primary-window cohort support", json.dumps(window_map, sort_keys=True), "T1>=4 and T2>=4", pass_if(window_map.get("T1", 0) >= 4 and window_map.get("T2", 0) >= 4))


# 5. Independent fixed-effect reproduction of the central E2 coupling estimate.
cohort_coupling = pd.read_csv(SRC / "04_12_cohort_signature_pathway_coupling.csv")
f = cohort_coupling[
    (cohort_coupling.signature_id == "SIG034")
    & (cohort_coupling.time_window == "T2")
    & (cohort_coupling.pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
    & (cohort_coupling.primary_independent_rule == "INCLUDE_PRIMARY_INDEPENDENT")
]
w = 1.0 / np.square(f["fisher_z_se"].to_numpy(float))
pooled_z = float(np.sum(w * f["fisher_z"].to_numpy(float)) / np.sum(w))
pooled_rho = math.tanh(pooled_z)
meta = pd.read_csv(SRC / "04_12_meta_signature_pathway_coupling.csv")
target = meta[
    (meta.analysis_set == "INDEPENDENT_ONLY")
    & (meta.signature_id == "SIG034")
    & (meta.time_window == "T2")
    & (meta.pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
].iloc[0]
rho_error = abs(pooled_rho - float(target.pooled_spearman_rho))
record("META_SIG034_OXPHOS_T48", "Meta-analysis input", "Independent inverse-variance Fisher-z pooled rho", f"{pooled_rho:.12f}; cohorts={len(f)}", f"{target.pooled_spearman_rho:.12f}; cohorts=5", pass_if(len(f) == 5 and rho_error <= 1e-10 and float(target.tau2_fisher_z) == 0.0), "1e-10")


# Recalculate the BH family for the 72 independent-only T48 primary cells.
family = meta[
    (meta.analysis_set == "INDEPENDENT_ONLY")
    & (meta.time_window == "T2")
    & (meta.tier == "PRESET_PRIMARY")
    & (meta.pathway_method == "SINGSCORE")
].copy()
p = family["p_value"].to_numpy(float)
order = np.argsort(p)
adj = np.empty(len(p), float)
running = 1.0
for rank_index in range(len(p) - 1, -1, -1):
    idx = order[rank_index]
    candidate = p[idx] * len(p) / (rank_index + 1)
    running = min(running, candidate)
    adj[idx] = min(1.0, running)
family["independent_bh"] = adj
target_bh = family[(family.signature_id == "SIG034") & (family.pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION")].iloc[0]
fdr_error = abs(float(target_bh.independent_bh) - float(target_bh.fdr_within_analysis_window_tier))
record("FDR_SIG034_OXPHOS_T48", "Multiplicity", "Independent BH recalculation within T48 primary Hallmark family", f"family={len(family)}; FDR={target_bh.independent_bh:.12f}", f"family=72; FDR={target_bh.fdr_within_analysis_window_tier:.12f}", pass_if(len(family) == 72 and fdr_error <= 1e-10), "1e-10")


# 6. Cell-source annotation completeness and historical alias handling.
cell = pd.read_csv(SRC / "04_13_cell_source_annotation.csv")
alias_n = int(cell.mapping_status.astype(str).str.contains("HISTORICAL", case=False).sum())
record("CELL_ANNOTATION_COVERAGE", "Cell-source annotation", "Unique signature genes with an HPA mapping", cell.signature_gene.nunique(), 74, pass_if(cell.signature_gene.nunique() == 74 and cell.current_hpa_symbol.notna().all()))
record("CELL_ALIAS_RESOLUTION", "Cell-source annotation", "Historical symbols resolved", alias_n, 2, pass_if(alias_n == 2))
record("CELL_CLAIM_BOUNDARY", "Cell-source annotation", "Every row retains potential-source language", int(cell.interpretation.astype(str).str.contains("Potential|HPA reports", regex=True).sum()), 74, pass_if(cell.interpretation.astype(str).str.contains("Potential|HPA reports", regex=True).all()))


# 7. Clinical gate was applied without fitting associations.
clinical = pd.read_csv(SRC / "04_15_clinical_anchor_availability_audit.csv")
primary_clin = clinical[clinical.anchor_level.astype(str).str.contains("PRIMARY")]
max_primary_n = int(primary_clin.estimated_time_matched_paired_n.max())
status = pd.read_csv(SRC / "04_17_signature_clinical_anchor_analysis_status.csv").iloc[0]
record("CLINICAL_GATE", "Clinical anchoring", "Maximum time-matched Level 1 patient count", max_primary_n, ">=80 required; observed below gate", pass_if(max_primary_n < 80))
record("CLINICAL_NO_MODELS", "Clinical anchoring", "Formal models and P values generated", f"models={int(status.models_fitted)}; p_values={int(status.p_values_generated)}", "models=0; p_values=0", pass_if(int(status.models_fitted) == 0 and int(status.p_values_generated) == 0))


# 8. Integrated matrix and figure-export consistency.
integrated = pd.read_csv(SRC / "04_20_integrated_signature_interpretation_matrix.csv")
record("INTEGRATED_SIGNATURE_COUNT", "Integrated evidence", "One unique row per frozen signature", f"rows={len(integrated)}; unique={integrated.signature_id.nunique()}", "8 and 8", pass_if(len(integrated) == 8 and integrated.signature_id.nunique() == 8))
record("INTEGRATED_CLINICAL_BOUNDARY", "Integrated evidence", "Clinical anchor state remains gated", int(integrated.formal_clinical_anchor.eq("NOT_RUN_BY_PREDEFINED_GATE").sum()), 8, pass_if(integrated.formal_clinical_anchor.eq("NOT_RUN_BY_PREDEFINED_GATE").all()))
figqc = pd.read_csv(SRC / "04_figure_export_qc.csv")
svg = figqc[figqc.format.eq("svg")]
record("FIGURE_EXPORTS", "Figure reproducibility", "Nonempty PDF, PNG, TIFF and editable-text SVG exports", f"rows={len(figqc)}; svg_pass={int(svg.svg_text_status.eq('PASS').sum())}", "56 files; 14 SVG PASS", pass_if(len(figqc) == 56 and figqc.status.eq("PRESENT_NONEMPTY").all() and len(svg) == 14 and svg.svg_text_status.eq("PASS").all()))


# 9. Workbook package integrity and explicit absence of formula errors/formulas.
workbooks = sorted(p for p in ROOT.glob("04_*.xlsx") if p.name != "04_22_independent_extension_verification.xlsx")
valid_zip = 0
formula_cells = 0
for book in workbooks:
    if zipfile.is_zipfile(book):
        valid_zip += 1
        with zipfile.ZipFile(book) as zf:
            for name in zf.namelist():
                if name.startswith("xl/worksheets/") and name.endswith(".xml"):
                    formula_cells += len(re.findall(rb"<f(?:\s|>)", zf.read(name)))
record("WORKBOOK_INTEGRITY", "Deliverable integrity", "Artifact-tool workbooks are valid XLSX archives", f"valid={valid_zip}; total={len(workbooks)}", "15 valid workbooks before 04_22 creation", pass_if(len(workbooks) == 15 and valid_zip == 15))
record("WORKBOOK_FORMULAS", "Deliverable integrity", "Formula/error risk scan", formula_cells, "0 formula cells; all values are frozen outputs", pass_if(formula_cells == 0))


result = pd.DataFrame(rows)
result.to_csv(OUT, index=False)
summary = {
    "verification_date": "2026-07-16",
    "checks": int(len(result)),
    "passed": int((result.status == "PASS").sum()),
    "failed": int((result.status == "FAIL").sum()),
    "status": "PASS" if result.status.eq("PASS").all() else "FAIL",
    "max_reconstruction_error": global_max,
    "reproduced_sig034_oxphos_t48_rho": pooled_rho,
    "reproduced_sig034_oxphos_t48_fdr": float(target_bh.independent_bh),
}
SUMMARY.write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
