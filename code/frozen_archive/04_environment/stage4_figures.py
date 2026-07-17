#!/usr/bin/env python3
"""Generate publication-grade Stage 4 figures from frozen result tables."""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd


# Mandatory Nature-figure font and editable-SVG rules.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["legend.frameon"] = False
plt.rcParams["font.size"] = 8


ROOT = Path("/Users/felix/Documents/New project/longitudinal_stage3")
SOURCE = ROOT / "04_source_data"
OUT = ROOT / "04_figures"
FOREST = ROOT / "04_07_gene_contribution_forest_plots"
OUT.mkdir(exist_ok=True)
FOREST.mkdir(exist_ok=True)

SIG_ORDER = ["SIG001", "SIG002", "SIG003", "SIG004", "SIG022", "SIG023", "SIG033", "SIG034"]
NAMES = {
    "SIG001": "Sepsis MetaScore",
    "SIG002": "SeptiCyte LAB",
    "SIG003": "FAIM3:PLAC8",
    "SIG004": "sNIP",
    "SIG022": "Bacterial/Viral MetaScore",
    "SIG023": "Herberg 2-gene DRS",
    "SIG033": "Lin 7-gene mortality score",
    "SIG034": "Severe-or-Mild score",
}
SHORT_PATH = {
    "HALLMARK_INFLAMMATORY_RESPONSE": "Inflammatory",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNF/NF-kB",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "IFN-alpha",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN-gamma",
    "HALLMARK_COMPLEMENT": "Complement",
    "HALLMARK_COAGULATION": "Coagulation",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION": "OxPhos",
    "HALLMARK_APOPTOSIS": "Apoptosis",
    "HALLMARK_ALLOGRAFT_REJECTION": "Allograft rejection",
}
PATH_ORDER = list(SHORT_PATH)
LINEAGE_COLORS = {
    "GRANULOCYTE": "#B64342",
    "MONOCYTE": "#D8892B",
    "DENDRITIC": "#9A4D8E",
    "B_LYMPHOCYTE": "#3775BA",
    "T_LYMPHOCYTE": "#0F4D92",
    "NK_CELL": "#42949E",
    "PLATELET": "#8C6D31",
    "ERYTHROID": "#7A5195",
    "BROAD_OR_UNRESOLVED": "#8F8F8F",
    "OTHER_OR_UNRESOLVED": "#B8B8B8",
}
ARCH_COLORS = {
    "CONSISTENT_MULTIGENE_DRIFT": "#0F4D92",
    "SINGLE_GENE_DOMINANT_DRIFT": "#B64342",
    "COHORT_DEPENDENT_DRIFT": "#9A4D8E",
    "INTERNAL_CANCELLATION_STABILITY": "#42949E",
    "OVERALL_LOW_CHANGE_STABILITY": "#8BCF8B",
    "EVIDENCE_INSUFFICIENT": "#8F8F8F",
}


def panel_label(ax, label: str) -> None:
    ax.text(-0.08, 1.03, label, transform=ax.transAxes, fontsize=10, fontweight="bold", ha="left", va="bottom")


def export(fig, base: Path, dpi: int = 300) -> None:
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def data():
    stage3 = pd.read_csv(ROOT / "source_data/03_14_signature_evidence_grade.csv")
    gene = pd.read_csv(SOURCE / "04_07_gene_contribution_meta_analysis.csv")
    gene = gene[(gene["analysis_set"] == "PRIMARY_INDEPENDENT") & (gene["time_window"] == "T2")]
    arch = pd.read_csv(SOURCE / "04_08_signature_drift_architecture.csv")
    ann = pd.read_csv(SOURCE / "04_13_cell_source_annotation.csv")
    cell = pd.read_csv(SOURCE / "04_14_gene_contribution_cell_source_map.csv")
    cell = cell[cell["time_window"] == "T2"]
    coupling = pd.read_csv(SOURCE / "04_12_meta_signature_pathway_coupling.csv")
    coupling = coupling[(coupling["analysis_set"] == "INDEPENDENT_ONLY") & (coupling["time_window"] == "T2") & (coupling["tier"] == "PRESET_PRIMARY")]
    context = pd.read_csv(SOURCE / "04_18_context_of_use_evidence_matrix.csv")
    return stage3, gene, arch, ann, cell, coupling, context


def figure_4a(stage3: pd.DataFrame, gene: pd.DataFrame, ann: pd.DataFrame) -> None:
    t48 = stage3[stage3["time_window"] == "T2"].set_index("signature_id").loc[SIG_ORDER].reset_index()
    gene = gene.merge(ann[["signature_gene", "broad_lineage"]], left_on="gene", right_on="signature_gene", how="left")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 4.4), gridspec_kw={"width_ratios": [1.0, 1.25], "wspace": 0.48})
    y = np.arange(len(SIG_ORDER))[::-1]
    labels = [f"{s}  {NAMES[s]}" for s in SIG_ORDER]
    for yi, row in zip(y, t48.itertuples(index=False)):
        color = "#3775BA" if row.pooled_delta_z < -0.2 else "#B64342" if row.pooled_delta_z > 0.2 else "#767676"
        ax1.plot([row.ci95_lower, row.ci95_upper], [yi, yi], color=color, lw=1.6)
        ax1.scatter(row.pooled_delta_z, yi, color=color, s=24, zorder=3)
    ax1.axvline(0, color="#767676", ls="--", lw=0.8)
    ax1.set_yticks(y)
    ax1.set_yticklabels(labels, fontsize=7)
    ax1.set_xlabel("T48 pooled score change, delta Z")
    ax1.set_title("Total fixed-signature drift", loc="left", fontweight="bold")
    panel_label(ax1, "a")

    for yi, sig in zip(y, SIG_ORDER):
        g = gene[gene["signature_id"] == sig].copy().nlargest(2, "pooled_contribution", keep="all")
        # nlargest by absolute magnitude, not signed value.
        g = gene[gene["signature_id"] == sig].assign(abs_effect=lambda x: x["pooled_contribution"].abs()).nlargest(2, "abs_effect")
        offsets = np.linspace(-0.13, 0.13, len(g))
        for off, row in zip(offsets, g.itertuples(index=False)):
            lineage = getattr(row, "broad_lineage") if pd.notna(getattr(row, "broad_lineage")) else "BROAD_OR_UNRESOLVED"
            color = LINEAGE_COLORS.get(lineage, "#8F8F8F")
            ax2.plot([row.ci95_lower, row.ci95_upper], [yi + off, yi + off], color=color, lw=1.2)
            ax2.scatter(row.pooled_contribution, yi + off, color=color, s=20, zorder=3)
            ax2.annotate(row.gene, (row.pooled_contribution, yi + off), xytext=(3 if row.pooled_contribution >= 0 else -3, 0),
                         textcoords="offset points", ha="left" if row.pooled_contribution >= 0 else "right", va="center", fontsize=6.2)
    ax2.axvline(0, color="#767676", ls="--", lw=0.8)
    ax2.set_yticks(y)
    ax2.set_yticklabels([])
    ax2.set_xlabel("Exact gene contribution to delta Z")
    ax2.set_title("Two largest mathematical contributors", loc="left", fontweight="bold")
    panel_label(ax2, "b")
    fig.text(0.01, 0.01, "T48 primary-independent estimates. Gene-specific random-effects weights differ; contribution points are not stacked into the total meta-estimate.", fontsize=6.5, color="#555555")
    fig.subplots_adjust(bottom=0.16)
    export(fig, OUT / "Figure_4A_drift_and_gene_contributions")


def figure_4b(arch: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 4.7))
    for row in arch.itertuples(index=False):
        color = ARCH_COLORS.get(row.drift_architecture, "#8F8F8F")
        size = 55 + 70 * min(row.median_patient_absolute_contribution_sum, 2.5) / 2.5
        ax.scatter(row.median_patient_dominance_ratio, row.median_patient_cancellation_index, s=size, color=color, edgecolor="white", linewidth=0.8, zorder=3)
        label_offset = {
            "SIG003": (-5, 8),
            "SIG023": (5, 6),
        }.get(row.signature_id, (4, 3))
        label_ha = "right" if row.signature_id == "SIG003" else "left"
        ax.annotate(row.signature_id, (row.median_patient_dominance_ratio, row.median_patient_cancellation_index), xytext=label_offset, textcoords="offset points", ha=label_ha, fontsize=7, fontweight="bold")
    ax.axvline(0.60, color="#B8B8B8", ls="--", lw=0.8)
    ax.axhline(0.50, color="#B8B8B8", ls="--", lw=0.8)
    ax.set_xlim(-0.02, 0.82)
    ax.set_ylim(-0.03, 0.90)
    ax.set_xlabel("Median patient dominance ratio")
    ax.set_ylabel("Median patient cancellation index")
    ax.set_title("T48 drift architecture under frozen classification rules", loc="left", fontweight="bold")
    handles = [mpl.lines.Line2D([], [], marker="o", ls="", color=c, label=k.replace("_", " ").title(), markersize=6) for k, c in ARCH_COLORS.items() if k in set(arch["drift_architecture"])]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=6.5)
    fig.text(0.01, 0.01, "Point size encodes median absolute contribution sum. Threshold lines are descriptive and were frozen before classification.", fontsize=6.5, color="#555555")
    fig.subplots_adjust(right=0.73, bottom=0.15)
    export(fig, OUT / "Figure_4B_dominance_cancellation")


def coupling_matrix(coupling: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    rho = coupling.pivot(index="signature_id", columns="pathway", values="pooled_spearman_rho").reindex(index=SIG_ORDER, columns=PATH_ORDER)
    fdr = coupling.pivot(index="signature_id", columns="pathway", values="fdr_within_analysis_window_tier").reindex(index=SIG_ORDER, columns=PATH_ORDER)
    return rho.to_numpy(float), fdr.to_numpy(float)


def draw_heatmap(ax, coupling: pd.DataFrame, annotate=True) -> None:
    matrix, fdr = coupling_matrix(coupling)
    cmap = mpl.colormaps["RdBu_r"].copy()
    cmap.set_bad("#F2F2F2")
    im = ax.imshow(matrix, cmap=cmap, vmin=-0.7, vmax=0.7, aspect="auto")
    ax.set_xticks(range(len(PATH_ORDER)))
    ax.set_xticklabels([SHORT_PATH[p] for p in PATH_ORDER], rotation=42, ha="right", fontsize=6.5)
    ax.set_yticks(range(len(SIG_ORDER)))
    ax.set_yticklabels(SIG_ORDER, fontsize=7)
    ax.tick_params(length=0)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if not np.isfinite(matrix[i, j]):
                continue
            if annotate:
                color = "white" if abs(matrix[i, j]) > 0.42 else "black"
                ax.text(j, i, f"{matrix[i,j]:.2f}", ha="center", va="center", fontsize=5.6, color=color)
            if np.isfinite(fdr[i, j]) and fdr[i, j] < 0.05:
                ax.add_patch(mpl.patches.Rectangle((j - 0.48, i - 0.48), 0.96, 0.96, fill=False, ec="#FFD700", lw=1.8))
    return im


def pathway_heatmap(coupling: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    im = draw_heatmap(ax, coupling)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Independent-only pooled Spearman rho", fontsize=7)
    ax.set_title("T48 coupling between fixed-signature change and prespecified pathway change", loc="left", fontweight="bold")
    fig.text(0.01, 0.01, "Gold outline: BH-FDR < 0.05 within the independent-only T48 primary Hallmark family. All prespecified cells are shown.", fontsize=6.5, color="#555555")
    fig.subplots_adjust(bottom=0.26)
    export(fig, ROOT / "04_12_signature_pathway_heatmap")


def figure_5(gene: pd.DataFrame, ann: pd.DataFrame, cell: pd.DataFrame, coupling: pd.DataFrame) -> None:
    detail = gene.merge(ann[["signature_gene", "broad_lineage"]], left_on="gene", right_on="signature_gene", how="left")
    leading = detail.assign(abs_effect=lambda x: x["pooled_contribution"].abs()).sort_values(["signature_id", "abs_effect"], ascending=[True, False]).groupby("signature_id").head(1).set_index("signature_id").loc[SIG_ORDER].reset_index()
    fig = plt.figure(figsize=(7.2, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.35], width_ratios=[1.0, 1.15], hspace=0.50, wspace=0.45)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])
    y = np.arange(len(SIG_ORDER))[::-1]
    for yi, row in zip(y, leading.itertuples(index=False)):
        lineage = row.broad_lineage if pd.notna(row.broad_lineage) else "BROAD_OR_UNRESOLVED"
        color = LINEAGE_COLORS.get(lineage, "#8F8F8F")
        ax1.plot([row.ci95_lower, row.ci95_upper], [yi, yi], color=color, lw=1.3)
        ax1.scatter(row.pooled_contribution, yi, s=22, color=color)
        ax1.annotate(row.gene, (row.pooled_contribution, yi), xytext=(3 if row.pooled_contribution >= 0 else -3, 0), textcoords="offset points", ha="left" if row.pooled_contribution >= 0 else "right", va="center", fontsize=6.2)
    ax1.axvline(0, color="#767676", ls="--", lw=0.8)
    ax1.set_yticks(y)
    ax1.set_yticklabels(SIG_ORDER)
    ax1.set_xlabel("Leading gene contribution")
    ax1.set_title("Mathematical driver and potential source", loc="left", fontweight="bold")
    panel_label(ax1, "a")

    pivot = cell.pivot(index="signature_id", columns="broad_lineage", values="absolute_share").fillna(0).reindex(SIG_ORDER)
    lineages = [c for c in LINEAGE_COLORS if c in pivot.columns]
    left = np.zeros(len(pivot))
    for lineage in lineages:
        vals = pivot[lineage].to_numpy(float)
        ax2.barh(np.arange(len(SIG_ORDER))[::-1], vals, left=left, height=0.68, color=LINEAGE_COLORS[lineage], label=lineage.replace("_", " ").title())
        left += vals
    ax2.set_yticks(y)
    ax2.set_yticklabels([])
    ax2.set_xlim(0, 1)
    ax2.set_xlabel("Share of absolute pooled contribution")
    ax2.set_title("Potential blood-cell source composition", loc="left", fontweight="bold")
    ax2.legend(fontsize=5.8, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    panel_label(ax2, "b")

    im = draw_heatmap(ax3, coupling)
    cbar = fig.colorbar(im, ax=ax3, fraction=0.022, pad=0.02)
    cbar.set_label("Pooled rho", fontsize=7)
    ax3.set_title("Prespecified pathway coupling at T48", loc="left", fontweight="bold")
    panel_label(ax3, "c")
    fig.text(0.01, 0.005, "Cell sources are HPA-based potential sources; bulk RNA cannot separate abundance from within-cell expression. Gold outline denotes BH-FDR < 0.05.", fontsize=6.2, color="#555555")
    fig.subplots_adjust(bottom=0.18, top=0.96)
    export(fig, OUT / "Figure_5_gene_pathway_cell_integration")


def gene_cell_map(gene: pd.DataFrame, ann: pd.DataFrame) -> None:
    d = gene.merge(ann[["signature_gene", "broad_lineage"]], left_on="gene", right_on="signature_gene", how="left")
    d["abs_effect"] = d["pooled_contribution"].abs()
    d["share"] = d["abs_effect"] / d.groupby("signature_id")["abs_effect"].transform("sum")
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    rng = np.random.default_rng(20260716)
    for i, sig in enumerate(SIG_ORDER):
        g = d[d["signature_id"] == sig]
        jitter = rng.uniform(-0.22, 0.22, len(g))
        colors = [LINEAGE_COLORS.get(x, "#8F8F8F") for x in g["broad_lineage"].fillna("BROAD_OR_UNRESOLVED")]
        ax.scatter(g["pooled_contribution"], np.full(len(g), len(SIG_ORDER)-1-i) + jitter, s=15 + 280 * g["share"], c=colors, alpha=0.82, edgecolor="white", linewidth=0.4)
        top = g.nlargest(1, "abs_effect").iloc[0]
        place_right = top["pooled_contribution"] >= 0 or top["pooled_contribution"] < -0.50
        ax.annotate(top["gene"], (top["pooled_contribution"], len(SIG_ORDER)-1-i), xytext=(4 if place_right else -4, 0), textcoords="offset points", ha="left" if place_right else "right", va="center", fontsize=6.5, fontweight="bold")
    ax.axvline(0, color="#767676", ls="--", lw=0.8)
    ax.set_yticks(np.arange(len(SIG_ORDER))[::-1])
    ax.set_yticklabels(SIG_ORDER)
    ax.set_xlabel("T48 primary-independent pooled gene contribution")
    ax.set_title("Exact gene contribution and potential blood-cell source", loc="left", fontweight="bold")
    present = [x for x in LINEAGE_COLORS if x in set(d["broad_lineage"].dropna())]
    handles = [mpl.lines.Line2D([], [], marker="o", ls="", color=LINEAGE_COLORS[x], label=x.replace("_", " ").title(), markersize=5) for x in present]
    ax.legend(handles=handles, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=6)
    fig.text(0.01, 0.01, "Point size is the gene share of total absolute contribution within a signature. Labels mark the largest absolute contributor.", fontsize=6.5, color="#555555")
    fig.subplots_adjust(bottom=0.25)
    export(fig, ROOT / "04_14_gene_contribution_cell_source_map")


def figure_6(context: pd.DataFrame) -> None:
    support_map = {
        "INDIRECT_ONLY": ("Indirect", "#D9D9D9"),
        "SUPPORTED_FOR_SCORE_CHANGE_DESCRIPTION": ("Supported: score change", "#8BCF8B"),
        "SUPPORTED_WITH_ATTRITION_CAUTION": ("Supported with caution", "#AADCA9"),
        "LIMITED_BIOLOGICAL_ANCHOR_ONLY": ("Limited E2 anchor", "#F3C677"),
        "EXPLORATORY_CONTEXT_ONLY": ("Exploratory only", "#E9A6A1"),
        "NOT_SUPPORTED": ("Not supported", "#B64342"),
    }
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, len(context) - 0.5)
    ax.axis("off")
    ax.text(0.02, 1.025, "Use scenario", transform=ax.transAxes, fontweight="bold", fontsize=8, va="bottom")
    ax.text(0.34, 1.025, "Current evidence", transform=ax.transAxes, fontweight="bold", fontsize=8, va="bottom")
    ax.text(0.60, 1.025, "Boundary for interpretation", transform=ax.transAxes, fontweight="bold", fontsize=8, va="bottom")
    for idx, row in enumerate(context.iloc[::-1].itertuples(index=False)):
        y = idx
        ax.add_patch(mpl.patches.Rectangle((0.01, y-0.38), 0.97, 0.76, facecolor="#F8FAFC" if idx % 2 == 0 else "#EEF3F5", edgecolor="none"))
        ax.text(0.02, y, row.use_scenario, va="center", fontsize=7.2, fontweight="bold")
        label, color = support_map[row.current_support]
        ax.add_patch(mpl.patches.FancyBboxPatch((0.34, y-0.19), 0.21, 0.38, boxstyle="round,pad=0.02,rounding_size=0.03", facecolor=color, edgecolor="none"))
        ax.text(0.445, y, label, ha="center", va="center", fontsize=6.5, color="white" if row.current_support == "NOT_SUPPORTED" else "#272727", fontweight="bold")
        boundary = row.claim_boundary
        if len(boundary) > 72:
            words = boundary.split(); lines=[]; current=[]
            for word in words:
                if len(" ".join(current+[word])) > 42:
                    lines.append(" ".join(current)); current=[word]
                else: current.append(word)
            if current: lines.append(" ".join(current))
            boundary = "\n".join(lines[:3])
        ax.text(0.60, y, boundary, va="center", fontsize=6.4, color="#374151")
    fig.text(0.12, 0.965, "Context of use for repeated host-response RNA scores", ha="left", va="top", fontsize=10, fontweight="bold")
    fig.text(0.01, 0.01, "The current study supports temporal description, not assay thresholds, prognostic validation or treatment decisions.", fontsize=6.5, color="#555555")
    fig.subplots_adjust(bottom=0.08, top=0.84)
    export(fig, OUT / "Figure_6_context_of_use_evidence_boundary")


def contribution_forests(gene: pd.DataFrame) -> None:
    for sig in SIG_ORDER:
        g = gene[gene["signature_id"] == sig].sort_values("pooled_contribution")
        height = max(3.2, 1.3 + 0.22 * len(g))
        fig, ax = plt.subplots(figsize=(5.8, height))
        y = np.arange(len(g))
        colors = np.where(g["pooled_contribution"] >= 0, "#B64342", "#3775BA")
        for yi, row, color in zip(y, g.itertuples(index=False), colors):
            ax.plot([row.ci95_lower, row.ci95_upper], [yi, yi], color=color, lw=1.2)
            ax.scatter(row.pooled_contribution, yi, color=color, s=16)
        ax.axvline(0, color="#767676", ls="--", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(g["gene"], fontsize=6.5)
        ax.set_xlabel("Pooled exact contribution to T48 delta Z")
        ax.set_title(f"{sig} — {NAMES[sig]}\nPrimary-independent gene contributions", loc="left", fontweight="bold")
        fig.text(0.01, 0.012, f"All {len(g)} component genes are shown; random-effects 95% confidence intervals. No significance filter.", fontsize=6.2, color="#555555")
        bottom_margin = min(0.24, max(0.10, 0.75 / height))
        fig.subplots_adjust(left=0.25, bottom=bottom_margin, top=0.91)
        export(fig, FOREST / f"{sig}_T48_gene_contributions", dpi=300)


def qa_inventory() -> None:
    rows = []
    for path in sorted(list(OUT.glob("*.*")) + list(FOREST.glob("*.*")) + list(ROOT.glob("04_12_signature_pathway_heatmap.*")) + list(ROOT.glob("04_14_gene_contribution_cell_source_map.*"))):
        if path.suffix.lower() not in {".svg", ".pdf", ".png", ".tiff"}:
            continue
        row = {"file": str(path.relative_to(ROOT)), "format": path.suffix.lower().lstrip("."), "bytes": path.stat().st_size, "status": "PRESENT_NONEMPTY" if path.stat().st_size > 0 else "EMPTY"}
        if path.suffix.lower() == ".svg":
            text = path.read_text(encoding="utf-8", errors="replace")
            row["editable_text_nodes"] = text.count("<text")
            row["svg_text_status"] = "PASS" if row["editable_text_nodes"] > 0 else "FAIL"
        rows.append(row)
    pd.DataFrame(rows).to_csv(SOURCE / "04_figure_export_qc.csv", index=False)


def main() -> None:
    stage3, gene, arch, ann, cell, coupling, context = data()
    figure_4a(stage3, gene, ann)
    figure_4b(arch)
    pathway_heatmap(coupling)
    figure_5(gene, ann, cell, coupling)
    gene_cell_map(gene, ann)
    figure_6(context)
    contribution_forests(gene)
    qa_inventory()
    print(json.dumps({"main_figures": 4, "standalone_figures": 2, "gene_forests": 8}, indent=2))


if __name__ == "__main__":
    main()
