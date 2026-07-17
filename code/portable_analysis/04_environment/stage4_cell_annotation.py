#!/usr/bin/env python3
"""Annotate fixed-signature genes with HPA blood-cell RNA evidence and integrate contributions."""

from __future__ import annotations

import os

import hashlib
import json
import re
import zipfile
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ.get("SEPSIS_SIGNATURE_ANALYSIS_ROOT", Path.cwd())).resolve()
SOURCE = ROOT / "04_source_data"
HPA_ZIP = ROOT / "04_cell_annotation_sources/proteinatlas.tsv.zip"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_values(text) -> list[tuple[str, float]]:
    if pd.isna(text) or not str(text).strip():
        return []
    rows = []
    for item in str(text).split(";"):
        if ":" not in item:
            continue
        name, value = item.rsplit(":", 1)
        try:
            rows.append((name.strip(), float(value.strip())))
        except ValueError:
            continue
    return sorted(rows, key=lambda x: x[1], reverse=True)


def broad_lineage(cell: str | None) -> str:
    if not cell:
        return "BROAD_OR_UNRESOLVED"
    x = cell.lower()
    if any(k in x for k in ("neutrophil", "eosinophil", "basophil")):
        return "GRANULOCYTE"
    if "monocyte" in x or "macrophage" in x:
        return "MONOCYTE"
    if "dendritic" in x or re.search(r"\bdc\b", x):
        return "DENDRITIC"
    if "b-cell" in x or "b cell" in x or "plasma" in x:
        return "B_LYMPHOCYTE"
    if any(k in x for k in ("t-cell", "t cell", "t-reg", "cd4", "cd8")):
        return "T_LYMPHOCYTE"
    if "nk" in x or "natural killer" in x:
        return "NK_CELL"
    if "platelet" in x or "megakary" in x:
        return "PLATELET"
    if any(k in x for k in ("eryth", "reticulocyte")):
        return "ERYTHROID"
    return "OTHER_OR_UNRESOLVED"


def main() -> None:
    coefficient = pd.read_csv(SOURCE / "04_03_signature_gene_coefficients.csv")
    genes = sorted(set(coefficient["gene"]))
    usecols = [
        "Gene", "Gene synonym", "Ensembl", "RNA blood cell specificity", "RNA blood cell distribution",
        "RNA blood cell specificity score", "RNA blood cell specific nTPM",
        "RNA single cell type specificity", "RNA single cell type distribution", "RNA single cell type specific nCPM",
    ]
    with zipfile.ZipFile(HPA_ZIP) as zf:
        hpa = pd.read_csv(zf.open("proteinatlas.tsv"), sep="\t", usecols=usecols)
    by_gene = {row.Gene: row for row in hpa.itertuples(index=False)}
    alias = {}
    for row in hpa.itertuples(index=False):
        synonyms = getattr(row, "_1") if hasattr(row, "_1") else None
        # itertuples sanitizes the spaced column; a row-wise fallback below resolves aliases reliably.
        if pd.notna(synonyms):
            for s in str(synonyms).split(", "):
                alias.setdefault(s, row.Gene)
    # Resolve aliases explicitly from the original dataframe to avoid tuple-name ambiguity.
    for _, row in hpa[["Gene", "Gene synonym"]].dropna().iterrows():
        for s in str(row["Gene synonym"]).split(", "):
            alias.setdefault(s, row["Gene"])

    rows = []
    for gene in genes:
        current = gene if gene in by_gene else alias.get(gene)
        if current is None:
            rows.append({
                "signature_gene": gene, "current_hpa_symbol": "", "mapping_status": "NOT_FOUND",
                "hpa_blood_cell_specificity": "", "hpa_blood_cell_distribution": "",
                "primary_cell_type": "", "secondary_cell_types": "", "broad_lineage": "BROAD_OR_UNRESOLVED",
                "hpa_quantitative_evidence": "", "pbmc_sepsis_atlas_scope": "NOT_ASSESSABLE",
                "interpretation": "No HPA entry; no cell-source assignment forced.",
            })
            continue
        r = hpa[hpa["Gene"] == current].iloc[0]
        values = parse_values(r["RNA blood cell specific nTPM"])
        specificity = str(r["RNA blood cell specificity"]) if pd.notna(r["RNA blood cell specificity"]) else ""
        if values:
            primary_cell = values[0][0]
            secondary = [name for name, value in values[1:] if value >= 0.5 * values[0][1]]
            lineage = broad_lineage(primary_cell)
            quantitative = ";".join(f"{name}:{value:g}" for name, value in values)
            interpretation = "Potential cell source based on HPA blood-cell RNA enrichment; bulk change may also reflect abundance shifts."
        else:
            primary_cell = ""
            secondary = []
            lineage = "BROAD_OR_UNRESOLVED"
            quantitative = ""
            interpretation = "HPA reports low/broad immune-cell specificity; no dominant lineage assigned."
        pbmc_scope = (
            "DIRECTLY_WITHIN_PROFILED_PBMC_LINEAGE_SCOPE"
            if lineage in {"MONOCYTE", "DENDRITIC", "B_LYMPHOCYTE", "T_LYMPHOCYTE", "NK_CELL"}
            else "OUTSIDE_OR_LIMITED_IN_PBMC_ATLAS_SCOPE"
            if lineage in {"GRANULOCYTE", "PLATELET", "ERYTHROID"}
            else "UNRESOLVED"
        )
        rows.append({
            "signature_gene": gene, "current_hpa_symbol": current,
            "mapping_status": "EXACT" if current == gene else "HISTORICAL_SYMBOL_RESOLVED",
            "hpa_blood_cell_specificity": specificity,
            "hpa_blood_cell_distribution": str(r["RNA blood cell distribution"]) if pd.notna(r["RNA blood cell distribution"]) else "",
            "primary_cell_type": primary_cell, "secondary_cell_types": ";".join(secondary), "broad_lineage": lineage,
            "hpa_quantitative_evidence": quantitative, "pbmc_sepsis_atlas_scope": pbmc_scope,
            "interpretation": interpretation,
        })
    annotation = pd.DataFrame(rows)
    annotation["primary_reference"] = "Human Protein Atlas v25.1 proteinatlas.tsv, downloaded 2026-07-16"
    annotation["secondary_context_reference"] = "Reyes et al. Nature Medicine 2020 sepsis PBMC atlas (SCP548/SCP550); gene-level matrices not redistributed in author repository and not used to override HPA assignment"
    annotation["source_url"] = "https://www.proteinatlas.org/about/download"
    annotation["secondary_source_doi"] = "10.1038/s41591-020-0752-4"
    annotation.to_csv(SOURCE / "04_13_cell_source_annotation.csv", index=False)

    meta = pd.read_csv(SOURCE / "04_07_gene_contribution_meta_analysis.csv")
    meta = meta[meta["analysis_set"] == "PRIMARY_INDEPENDENT"].merge(annotation, left_on="gene", right_on="signature_gene", how="left")
    meta["absolute_pooled_contribution"] = meta["pooled_contribution"].abs()
    meta["signature_window_total_absolute"] = meta.groupby(["signature_id", "time_window"])["absolute_pooled_contribution"].transform("sum")
    meta["absolute_share"] = meta["absolute_pooled_contribution"] / meta["signature_window_total_absolute"]
    meta.to_csv(SOURCE / "04_14_gene_contribution_cell_source_detail.csv", index=False)

    lineage = meta.groupby(["signature_id", "time_window", "broad_lineage"], as_index=False).agg(
        net_pooled_contribution=("pooled_contribution", "sum"),
        absolute_pooled_contribution=("absolute_pooled_contribution", "sum"),
        gene_count=("gene", "nunique"),
    )
    lineage["absolute_share"] = lineage["absolute_pooled_contribution"] / lineage.groupby(["signature_id", "time_window"])["absolute_pooled_contribution"].transform("sum")
    lineage.to_csv(SOURCE / "04_14_gene_contribution_cell_source_map.csv", index=False)

    top = lineage.sort_values(["signature_id", "time_window", "absolute_share"], ascending=[True, True, False]).groupby(["signature_id", "time_window"]).head(1)
    top = top.rename(columns={"broad_lineage": "leading_cell_source", "absolute_share": "leading_cell_source_absolute_share"})
    top.to_csv(SOURCE / "04_14_leading_cell_source_by_signature.csv", index=False)

    summary = {
        "signature_gene_count": len(genes),
        "exact_or_alias_mapped": int((annotation["mapping_status"] != "NOT_FOUND").sum()),
        "historical_symbols_resolved": int((annotation["mapping_status"] == "HISTORICAL_SYMBOL_RESOLVED").sum()),
        "lineage_assigned": int((~annotation["broad_lineage"].isin(["BROAD_OR_UNRESOLVED", "OTHER_OR_UNRESOLVED"])).sum()),
        "hpa_version": "25.1",
        "hpa_file_sha256": sha256(HPA_ZIP),
        "claim_boundary": "Potential cell source only; no inference of within-cell transcription or deconvolved cell proportion.",
    }
    (SOURCE / "04_cell_annotation_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
