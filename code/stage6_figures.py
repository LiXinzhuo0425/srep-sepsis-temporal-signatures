#!/usr/bin/env python3
"""Generate Scientific Reports figures from frozen Stage 3 and corrected Stage 4 results.

Only presentation and reviewer-facing terminology change here; all plotted
scientific values are read from the frozen Stage 3/4 result tables.
"""

from __future__ import annotations

import json
import math
import os
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
from PIL import Image


PROJECT = Path(os.environ.get("SEPSIS_SIGNATURE_PROJECT_ROOT", Path.cwd())).resolve()
S3 = PROJECT / "longitudinal_stage3/source_data"
S4 = PROJECT / "longitudinal_stage3/04_source_data"
OUT = Path(os.environ.get(
    "STAGE6_FIGURE_OUT",
    PROJECT / "stage6_submission/06_02_revision/06_figures/main_v1_1",
)).resolve()
SOURCE_OUT = Path(os.environ.get(
    "STAGE6_SOURCE_OUT",
    PROJECT / "stage6_submission/06_02_revision/07_source_data_v1_1",
)).resolve()
SUPP_OUT = Path(os.environ.get(
    "STAGE6_SUPP_OUT",
    PROJECT / "stage5_manuscript/06_supplementary/figures",
)).resolve()
OUT.mkdir(parents=True, exist_ok=True)
SOURCE_OUT.mkdir(parents=True, exist_ok=True)
SUPP_OUT.mkdir(parents=True, exist_ok=True)

SIGS = ["SIG001", "SIG002", "SIG003", "SIG004", "SIG022", "SIG023", "SIG033", "SIG034"]
NAMES = {
    "SIG001": "Sepsis MetaScore",
    "SIG002": "SeptiCyte LAB",
    "SIG003": "FAIM3:PLAC8",
    "SIG004": "sNIP",
    "SIG022": "Bacterial/Viral MetaScore",
    "SIG023": "Herberg two-transcript DRS",
    "SIG033": "Lin seven-gene mortality score",
    "SIG034": "Severe-or-Mild score",
}
COHORTS = ["GSE236713", "GSE57065", "GSE95233", "GSE54514", "GSE110487", "GSE8121"]
COHORT_INFO = pd.DataFrame([
    ["GSE236713", "Adult sepsis", "PBL", "Agilent", 106, 93, 0, "Prespecified non-pilot"],
    ["GSE57065", "Adult septic shock", "Whole blood", "Affymetrix", 28, 28, 26, "Pilot"],
    ["GSE95233", "Adult septic shock", "Whole blood", "Affymetrix", 51, 20, 31, "Prespecified non-pilot"],
    ["GSE54514", "Adult sepsis", "Whole blood", "Illumina", 31, 31, 28, "Pilot"],
    ["GSE110487", "Adult septic shock", "Whole blood", "RNA-seq", 31, 0, 31, "Prespecified non-pilot"],
    ["GSE8121", "Pediatric septic shock", "Whole blood", "Affymetrix", 30, 0, 30, "Prespecified non-pilot"],
], columns=["cohort", "population", "sample_type", "platform", "longitudinal_n", "T24_n", "T48_n", "role"])

COL = {
    "blue": "#2F6B8A",
    "sky": "#74A9CF",
    "red": "#B24745",
    "orange": "#D98B3A",
    "purple": "#7A5AA6",
    "teal": "#3A8F8D",
    "green": "#4E8B57",
    "grey": "#777777",
    "light": "#EEF3F6",
    "ink": "#24323D",
}
COHORT_COLORS = {
    "GSE236713": "#2F6B8A", "GSE57065": "#70AD47", "GSE95233": "#B24745",
    "GSE54514": "#3A8F8D", "GSE110487": "#7A5AA6", "GSE8121": "#D19A38",
}
ARCH_COLORS = {
    "CONSISTENT_MULTIGENE_DRIFT": COL["blue"],
    "SINGLE_GENE_DOMINANT_DRIFT": COL["red"],
    "COHORT_DEPENDENT_DRIFT": COL["purple"],
    "INTERNAL_CANCELLATION_STABILITY": COL["teal"],
    "OVERALL_LOW_CHANGE_STABILITY": COL["green"],
}
ARCH_SHORT = {
    "CONSISTENT_MULTIGENE_DRIFT": "Multigene drift",
    "SINGLE_GENE_DOMINANT_DRIFT": "Single-gene dominant",
    "COHORT_DEPENDENT_DRIFT": "Cohort-dependent",
    "INTERNAL_CANCELLATION_STABILITY": "Internal cancellation",
    "OVERALL_LOW_CHANGE_STABILITY": "Low change",
}
ARCH_STYLE = {
    "CONSISTENT_MULTIGENE_DRIFT": {"marker": "^", "linestyle": ":"},
    "SINGLE_GENE_DOMINANT_DRIFT": {"marker": "s", "linestyle": "--"},
    "COHORT_DEPENDENT_DRIFT": {"marker": "o", "linestyle": "-"},
    "INTERNAL_CANCELLATION_STABILITY": {"marker": "D", "linestyle": "-."},
    "OVERALL_LOW_CHANGE_STABILITY": {"marker": "v", "linestyle": (0, (1, 1))},
}
STABILITY_COLORS = {
    "STABLE": "#2F6B8A",
    "BOUNDARY_SENSITIVE": "#8A8A8A",
    "UNSTABLE": "#D9822B",
}
STABILITY_LABELS = {
    "STABLE": "Stable across perturbations",
    "BOUNDARY_SENSITIVE": "Changed in one scenario",
    "UNSTABLE": "Changed in ≥2 scenarios",
}
SHORT_PATH = {
    "HALLMARK_INFLAMMATORY_RESPONSE": "Inflammatory",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNF/NF-kB",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "IFN-α",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN-γ",
    "HALLMARK_COMPLEMENT": "Complement",
    "HALLMARK_COAGULATION": "Coagulation",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION": "OxPhos",
    "HALLMARK_APOPTOSIS": "Apoptosis",
    "HALLMARK_ALLOGRAFT_REJECTION": "Allograft rejection",
}
PATH_ORDER = list(SHORT_PATH)
LINEAGE_COLORS = {
    "GRANULOCYTE": "#B64342", "MONOCYTE": "#D8892B", "DENDRITIC": "#9A4D8E",
    "B_LYMPHOCYTE": "#3775BA", "T_LYMPHOCYTE": "#0F4D92", "NK_CELL": "#42949E",
    "PLATELET": "#8C6D31", "ERYTHROID": "#7A5195",
    "BROAD_OR_UNRESOLVED": "#8F8F8F", "OTHER_OR_UNRESOLVED": "#B8B8B8",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "font.size": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})


def panel(ax, letter: str) -> None:
    ax.text(-0.10, 1.04, letter, transform=ax.transAxes, ha="left", va="bottom", fontsize=11, fontweight="bold")


def export(fig: plt.Figure, stem: str) -> None:
    base = OUT / stem
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".eps"), format="eps", bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    tiff_path = base.with_suffix(".tiff")
    fig.savefig(tiff_path, dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    with Image.open(tiff_path) as image:
        image.convert("RGB").save(tiff_path, dpi=(600, 600), compression="tiff_lzw")
    plt.close(fig)


def export_to(fig: plt.Figure, base: Path) -> None:
    """Export a figure outside the main-figure directory."""
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def public_labels(data: pd.DataFrame) -> pd.DataFrame:
    """Translate legacy frozen machine labels in reviewer-facing source data."""
    return data.replace({
        "PRESPECIFIED_NON_PILOT": "PRESPECIFIED_NON_PILOT",
        "PRESPECIFIED_NON_PILOT_ONLY": "PRESPECIFIED_NON_PILOT_ONLY",
        "PRESPECIFIED_VALIDATION": "PRESPECIFIED_NON_PILOT",
        "PRESPECIFIED_VALIDATION_ONLY": "PRESPECIFIED_NON_PILOT_ONLY",
        "Blinded validation": "Prespecified non-pilot",
        "Prespecified validation": "Prespecified non-pilot",
    })


def box(ax, x, y, w, h, text, face, edge="#C9D4DB", fs=8.5) -> None:
    patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.014,rounding_size=0.018",
                           facecolor=face, edgecolor=edge, linewidth=0.9)
    ax.add_patch(patch)
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs, color=COL["ink"])


def figure1() -> None:
    fig = plt.figure(figsize=(7.2, 6.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.05, 1.0], hspace=0.32)
    ax = fig.add_subplot(gs[0]); ax.set_axis_off(); panel(ax, "a")
    box(ax, 0.01, 0.62, 0.16, 0.22, "6 public\ncohorts", "#E8F1F5", fs=9.2)
    box(ax, 0.21, 0.62, 0.17, 0.22, "302 baseline\ncase patients", "#E8F1F5", fs=9.2)
    box(ax, 0.42, 0.62, 0.17, 0.22, "277 with ≥2\ntime windows", "#EAF3EA", fs=9.2)
    box(ax, 0.63, 0.62, 0.17, 0.22, "264 with T24\nor T48 pair", "#EAF3EA", fs=9.2)
    box(ax, 0.83, 0.62, 0.15, 0.22, "8 fixed\nsignatures", "#FFF3D5", fs=9.2)
    for x1, x2 in [(0.17, 0.21), (0.38, 0.42), (0.59, 0.63), (0.80, 0.84)]:
        ax.annotate("", xy=(x2, 0.73), xytext=(x1, 0.73), arrowprops=dict(arrowstyle="->", color="#64737D", lw=1.1))
    box(ax, 0.04, 0.20, 0.25, 0.20, "Fixed-score reconstruction\n48 cohort-signature pairs", "#EFF4F7")
    box(ax, 0.375, 0.20, 0.25, 0.20,
        "Eligible landmark pairs\nT24: n=172  |  T48: n=146\nPrimary k: 3–4  |  4–5",
        "#EFF4F7", fs=7.5)
    box(ax, 0.705, 0.20, 0.25, 0.20,
        "Meta-analysis, gene\ndecomposition, and annotation", "#EFF4F7", fs=7.9)
    for x1, x2 in [(0.29, 0.375), (0.625, 0.71)]:
        ax.annotate("", xy=(x2, 0.30), xytext=(x1, 0.30), arrowprops=dict(arrowstyle="->", color="#64737D", lw=1.1))
    ax.text(0.01, 0.97, "Longitudinal analysis design", fontsize=11, fontweight="bold", va="top")
    ax.text(0.01, 0.04, "T24 = 12–36 h; T48 = 36–60 h relative to each study-specific baseline definition.", fontsize=6.8, color="#5F6B76")

    ax2 = fig.add_subplot(gs[1]); panel(ax2, "b")
    y = np.arange(len(COHORT_INFO))[::-1]
    ax2.barh(y + 0.17, COHORT_INFO["T24_n"], height=0.31, color=COL["blue"], label="T24")
    ax2.barh(y - 0.17, COHORT_INFO["T48_n"], height=0.31, color=COL["red"], label="T48")
    ax2.set_yticks(y)
    ax2.set_yticklabels([f"{r.cohort}  {r.population}" for r in COHORT_INFO.itertuples(index=False)], fontsize=7)
    ax2.set_xlabel("Paired patients contributing to the primary window")
    ax2.set_xlim(0, 103)
    ax2.legend(ncol=2, loc="lower right")
    for yi, row in zip(y, COHORT_INFO.itertuples(index=False)):
        if row.T24_n: ax2.text(row.T24_n + 1.2, yi + 0.17, str(row.T24_n), va="center", fontsize=6.5)
        if row.T48_n: ax2.text(row.T48_n + 1.2, yi - 0.17, str(row.T48_n), va="center", fontsize=6.5)
    ax2.set_title("Cohort contributions to each primary window", loc="left", fontweight="bold")
    fig.suptitle("Study design and analysis flow", x=0.02, y=0.995, ha="left", fontsize=12, fontweight="bold")
    export(fig, "Figure_1_study_design")
    COHORT_INFO.to_csv(SOURCE_OUT / "Figure_1_source_data.csv", index=False)


def figure2() -> None:
    meta = pd.read_csv(S3 / "03_05_signature_level_meta_analysis.csv")
    data = meta[(meta.analysis_set == "PRIMARY_INDEPENDENT") & meta.time_window.isin(["T1", "T2"])].copy()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.8), sharey=True)
    for ax, window, title, color in zip(axes, ["T1", "T2"], ["T24 (12–36 h)", "T48 (36–60 h)"], [COL["blue"], COL["red"]]):
        sub = data[data.time_window == window].set_index("signature_id").reindex(SIGS)
        y = np.arange(len(SIGS))[::-1]
        ax.axvline(0, color="#777777", ls="--", lw=0.8)
        for yi, sig in zip(y, SIGS):
            row = sub.loc[sig]
            if np.isfinite(row.prediction_lower):
                ax.plot([row.prediction_lower, row.prediction_upper], [yi-0.15, yi-0.15], color=color, lw=0.9, alpha=0.43)
            ax.plot([row.ci95_lower, row.ci95_upper], [yi, yi], color=color, lw=2.0)
            ax.scatter(row.pooled_delta_z, yi, s=26, color=color, edgecolor="white", linewidth=0.5, zorder=3)
            ax.text(ax.get_xlim()[1] if False else row.pooled_delta_z, yi+0.20, f"{row.pooled_delta_z:.2f}", ha="center", fontsize=5.8, color=color)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{s}  {NAMES[s]}" for s in SIGS], fontsize=6.7)
        ax.set_xlabel("Pooled within-patient change (ΔZ)")
        ax.set_title(title, fontweight="bold")
        ax.text(0.02, 0.02, "Thick: 95% CI\nThin: 95% prediction interval", transform=ax.transAxes, fontsize=6.1, color="#5F6B76")
    panel(axes[0], "a"); panel(axes[1], "b")
    fig.suptitle("Pooled longitudinal change in eight fixed signatures", x=0.02, y=0.98, ha="left", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    export(fig, "Figure_2_pooled_longitudinal_change")
    data.to_csv(SOURCE_OUT / "Figure_2_source_data.csv", index=False)


def heat(ax, matrix, rows, cols, *, vmin, vmax, cmap, fmt, title, cbar_label, mask=None) -> None:
    cm = mpl.colormaps[cmap].copy(); cm.set_bad("#F2F2F2")
    im = ax.imshow(matrix, aspect="auto", vmin=vmin, vmax=vmax, cmap=cm)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, fontsize=7.6)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=7.6)
    ax.tick_params(length=0)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if not np.isfinite(matrix[i, j]): continue
            color = "white" if abs(matrix[i, j]) > 0.65 * max(abs(vmin), abs(vmax)) else "black"
            ax.text(j, i, format(matrix[i, j], fmt), ha="center", va="center", fontsize=6.7, color=color)
            if mask is not None and bool(mask[i, j]):
                ax.add_patch(mpl.patches.Rectangle((j-.47, i-.47), .94, .94, fill=False, ec="#F5C542", lw=1.5))
    ax.set_title(title, loc="left", fontweight="bold")
    cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.018)
    cb.set_label(cbar_label, fontsize=7.4); cb.ax.tick_params(labelsize=7.0)


def figure3() -> None:
    effects = pd.read_csv(S3 / "03_04_cohort_signature_primary_effects.csv")
    profile = pd.read_csv(S3 / "03_06_temporal_stability_profile.csv")
    fig = plt.figure(figsize=(7.2, 7.0))
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 0.85], hspace=0.46)
    for idx, (window, title) in enumerate([("T1", "T24 cohort-specific mean ΔZ"), ("T2", "T48 cohort-specific mean ΔZ")]):
        ax = fig.add_subplot(gs[idx]); panel(ax, "a" if idx == 0 else "b")
        sub = effects[effects.time_window == window].pivot(index="dataset", columns="signature_id", values="mean_delta_z").reindex(index=COHORTS, columns=SIGS)
        heat(ax, sub.to_numpy(float), COHORTS, SIGS, vmin=-1.6, vmax=1.6,
             cmap="RdBu_r", fmt=".2f", title=title, cbar_label="Mean ΔZ")
        # GSE54514 was retained for transparent cohort-level display but was
        # excluded from the primary independent synthesis for these signatures.
        row_i = COHORTS.index("GSE54514")
        for sig in ["SIG002", "SIG003", "SIG004"]:
            col_j = SIGS.index(sig)
            if np.isfinite(sub.loc["GSE54514", sig]):
                ax.add_patch(mpl.patches.Rectangle(
                    (col_j - 0.47, row_i - 0.47), 0.94, 0.94,
                    fill=False, ec="#222222", lw=1.25,
                ))
                ax.text(col_j + 0.34, row_i - 0.31, "†", ha="center",
                        va="center", fontsize=8.2, fontweight="bold")
    ax3 = fig.add_subplot(gs[2]); panel(ax3, "c")
    p = profile[profile.time_window == "T2"].set_index("signature_id").reindex(SIGS)
    y = np.arange(len(SIGS))[::-1]
    widths = (p.prediction_upper - p.prediction_lower).to_numpy(float)
    bars = ax3.barh(y, p.I2_percent, height=0.58, color=COL["sky"])
    ax3.set_yticks(y); ax3.set_yticklabels([f"{s}  {NAMES[s]}" for s in SIGS], fontsize=6.6)
    ax3.set_xlim(0, 105); ax3.set_xlabel("T48 I² (%)")
    ax3.axvline(50, color="#999999", lw=0.8, ls="--")
    for bar, i2, width in zip(bars, p.I2_percent, widths):
        ax3.text(min(i2 + 2, 97), bar.get_y()+bar.get_height()/2, f"{i2:.0f}%  |  PI width {width:.2f}", va="center", fontsize=6.1)
    ax3.set_title("T48 between-cohort heterogeneity and prediction-interval width", loc="left", fontweight="bold")
    ax3.text(0.99, 0.02, "Bar color is descriptive only.", transform=ax3.transAxes,
             ha="right", va="bottom", fontsize=5.8, color="#5F6B76")
    fig.suptitle("Cross-cohort transportability of longitudinal changes", x=0.02, y=0.995, ha="left", fontsize=12, fontweight="bold")
    export(fig, "Figure_3_cross_cohort_transportability")
    public_labels(effects).to_csv(SOURCE_OUT / "Figure_3A_B_source_data.csv", index=False)
    profile.to_csv(SOURCE_OUT / "Figure_3C_source_data.csv", index=False)


def figure4() -> None:
    gene = pd.read_csv(S4 / "04_07_gene_contribution_meta_analysis.csv")
    gene = gene[(gene.analysis_set == "PRIMARY_INDEPENDENT") & (gene.time_window == "T2")].copy()
    ann = pd.read_csv(S4 / "04_13_cell_source_annotation.csv")[["signature_gene", "broad_lineage"]]
    gene = gene.merge(ann, left_on="gene", right_on="signature_gene", how="left")
    stability = pd.read_csv(S3 / "10_02_classification_stability_matrix.csv")
    arch = pd.read_csv(S4 / "04_08_signature_drift_architecture.csv")
    arch = arch.merge(stability, on="signature_id", how="left").set_index("signature_id").reindex(SIGS).reset_index()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 5.0), gridspec_kw={"width_ratios": [1.30, 1], "wspace": 0.44})
    y = np.arange(len(SIGS))[::-1]
    for yi, sig in zip(y, SIGS):
        g = gene[gene.signature_id == sig].assign(abs_effect=lambda d: d.pooled_contribution.abs()).nlargest(2, "abs_effect")
        offsets = np.linspace(-0.13, 0.13, len(g))
        for off, row in zip(offsets, g.itertuples(index=False)):
            color = COL["red"] if row.pooled_contribution > 0 else COL["blue"]
            ax1.plot([row.ci95_lower, row.ci95_upper], [yi+off, yi+off], color=color, lw=1.3)
            ax1.scatter(row.pooled_contribution, yi+off, s=22, color=color, edgecolor="white", linewidth=0.4)
            ax1.annotate(row.gene, (row.pooled_contribution, yi+off), xytext=(3 if row.pooled_contribution >= 0 else -3, 0),
                         textcoords="offset points", ha="left" if row.pooled_contribution >= 0 else "right", va="center", fontsize=5.8)
    ax1.axvline(0, color="#777777", lw=0.8, ls="--")
    ax1.set_yticks(y); ax1.set_yticklabels([f"{s}  {NAMES[s]}" for s in SIGS], fontsize=6.4)
    ax1.set_xlabel("Exact gene contribution to T48 ΔZ")
    ax1.set_title("Two largest pooled mathematical contributions", loc="left", fontweight="bold")
    ax1.text(-0.14, 1.035, "a", transform=ax1.transAxes, ha="left", va="bottom", fontsize=11, fontweight="bold")

    for row in arch.itertuples(index=False):
        color = STABILITY_COLORS.get(row.sensitivity_status, COL["grey"])
        size = 55 + 80 * min(row.median_patient_absolute_contribution_sum, 2.5)/2.5
        ax2.scatter(row.median_patient_dominance_ratio, row.median_patient_cancellation_index, s=size, color=color, edgecolor="white", linewidth=0.7, zorder=3)
        dx, dy = {
            "SIG001": (4, 4), "SIG002": (5, 10), "SIG003": (-7, 10), "SIG004": (4, 4),
            "SIG022": (-7, 9), "SIG023": (5, 9), "SIG033": (4, 4), "SIG034": (4, 4),
        }.get(row.signature_id, (4, 3))
        label_ha = "right" if row.signature_id in {"SIG003", "SIG022"} else "left"
        ax2.annotate(row.signature_id, (row.median_patient_dominance_ratio, row.median_patient_cancellation_index), xytext=(dx, dy), textcoords="offset points", fontsize=6.6, fontweight="bold", ha=label_ha)
    ax2.axvline(0.60, color="#8A8A8A", lw=0.8, ls="--"); ax2.axhline(0.50, color="#8A8A8A", lw=0.8, ls="--")
    ax2.set_xlim(-0.03, 0.83); ax2.set_ylim(-0.04, 0.90)
    ax2.set_xlabel("Median cohort dominance"); ax2.set_ylabel("Median cohort cancellation")
    ax2.set_title("Continuous metrics and label sensitivity", loc="left", fontweight="bold")
    ax2.text(-0.14, 1.035, "b", transform=ax2.transAxes, ha="left", va="bottom", fontsize=11, fontweight="bold")
    handles = [
        mpl.lines.Line2D([], [], marker="o", ls="", color=color,
                         label=STABILITY_LABELS[status],
                         markersize=6)
        for status, color in STABILITY_COLORS.items()
    ]
    ax2.legend(handles=handles, fontsize=5.8, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=1)
    fig.text(0.01, 0.012, "Continuous metrics are primary. Each point is the median of cohort-specific patient medians from the same primary-independent set used for score synthesis.\nColors summarize label sensitivity under uniform ±10% and ±20% threshold perturbations; contributions are mathematical, not causal.", fontsize=5.8, color="#5F6B76", linespacing=1.18)
    fig.subplots_adjust(bottom=0.25, top=0.86)
    fig.suptitle("Gene-contribution patterns at T48", x=0.02, y=0.985, ha="left", fontsize=12, fontweight="bold")
    export(fig, "Figure_4_gene_contribution_patterns")
    gene.to_csv(SOURCE_OUT / "Figure_4A_source_data.csv", index=False)
    arch.to_csv(SOURCE_OUT / "Figure_4B_source_data.csv", index=False)


def coupling_matrix(coupling: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    rho = coupling.pivot(index="signature_id", columns="pathway", values="pooled_spearman_rho").reindex(index=SIGS, columns=PATH_ORDER)
    fdr = coupling.pivot(index="signature_id", columns="pathway", values="fdr_within_analysis_window_tier").reindex(index=SIGS, columns=PATH_ORDER)
    return rho.to_numpy(float), fdr.to_numpy(float)


def draw_pathway_heatmap(ax, coupling: pd.DataFrame):
    matrix, fdr = coupling_matrix(coupling)
    cmap = mpl.colormaps["RdBu_r"].copy(); cmap.set_bad("#F2F2F2")
    im = ax.imshow(matrix, cmap=cmap, vmin=-0.7, vmax=0.7, aspect="auto")
    ax.set_xticks(range(len(PATH_ORDER)))
    ax.set_xticklabels([SHORT_PATH[p] for p in PATH_ORDER], rotation=38, ha="right", fontsize=9.0)
    ax.set_yticks(range(len(SIGS))); ax.set_yticklabels(SIGS, fontsize=8.4)
    ax.tick_params(length=0)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if not np.isfinite(matrix[i, j]):
                continue
            color = "white" if abs(matrix[i, j]) > 0.42 else "black"
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=7.2, color=color)
            if np.isfinite(fdr[i, j]) and fdr[i, j] < 0.05:
                ax.add_patch(mpl.patches.Rectangle((j - 0.48, i - 0.48), 0.96, 0.96,
                                                   fill=False, ec="#FFD700", lw=1.8))
                ax.text(j + 0.34, i - 0.30, "*", ha="center", va="center",
                        fontsize=11.0, fontweight="bold", color="black", zorder=5)
    return im


def figure5() -> None:
    gene = pd.read_csv(S4 / "04_07_gene_contribution_meta_analysis.csv")
    gene = gene[(gene.analysis_set == "PRIMARY_INDEPENDENT") & (gene.time_window == "T2")].copy()
    ann = pd.read_csv(S4 / "04_13_cell_source_annotation.csv")
    cell = pd.read_csv(S4 / "04_14_gene_contribution_cell_source_map.csv")
    cell = cell[cell.time_window == "T2"].copy()
    coupling = pd.read_csv(S4 / "04_12_meta_signature_pathway_coupling.csv")
    coupling = coupling[(coupling.analysis_set == "INDEPENDENT_ONLY") &
                        (coupling.time_window == "T2") & (coupling.tier == "PRESET_PRIMARY")].copy()

    detail = gene.merge(ann[["signature_gene", "broad_lineage"]], left_on="gene", right_on="signature_gene", how="left")
    leading = (detail.assign(abs_effect=lambda d: d.pooled_contribution.abs())
               .sort_values(["signature_id", "abs_effect"], ascending=[True, False])
               .groupby("signature_id").head(1).set_index("signature_id").loc[SIGS].reset_index())
    fig = plt.figure(figsize=(7.2, 7.15))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.95, 1.60], width_ratios=[1.0, 1.15], hspace=0.54, wspace=0.44)
    ax1 = fig.add_subplot(gs[0, 0]); ax2 = fig.add_subplot(gs[0, 1]); ax3 = fig.add_subplot(gs[1, :])
    y = np.arange(len(SIGS))[::-1]
    for yi, row in zip(y, leading.itertuples(index=False)):
        lineage = row.broad_lineage if pd.notna(row.broad_lineage) else "BROAD_OR_UNRESOLVED"
        color = LINEAGE_COLORS.get(lineage, "#8F8F8F")
        ax1.plot([row.ci95_lower, row.ci95_upper], [yi, yi], color=color, lw=1.3)
        ax1.scatter(row.pooled_contribution, yi, s=22, color=color)
        ax1.annotate(row.gene, (row.pooled_contribution, yi), xytext=(3 if row.pooled_contribution >= 0 else -3, 0),
                     textcoords="offset points", ha="left" if row.pooled_contribution >= 0 else "right",
                     va="center", fontsize=6.7)
    ax1.axvline(0, color="#767676", ls="--", lw=0.8)
    ax1.set_yticks(y); ax1.set_yticklabels(SIGS, fontsize=7.2)
    ax1.set_xlabel("Leading gene contribution")
    ax1.set_title("Mathematical driver and potential source", loc="left", fontweight="bold"); panel(ax1, "a")

    pivot = cell.pivot(index="signature_id", columns="broad_lineage", values="absolute_share").fillna(0).reindex(SIGS)
    lineages = [c for c in LINEAGE_COLORS if c in pivot.columns]
    left = np.zeros(len(pivot))
    for lineage in lineages:
        vals = pivot[lineage].to_numpy(float)
        ax2.barh(y, vals, left=left, height=0.68, color=LINEAGE_COLORS[lineage],
                 edgecolor="white", linewidth=0.55,
                 label=lineage.replace("_", " ").title())
        left += vals
    ax2.set_yticks(y); ax2.set_yticklabels([]); ax2.set_xlim(0, 1)
    ax2.set_xlabel("Share of absolute pooled contribution")
    ax2.set_title("Mathematical contribution by potential source", loc="left", fontweight="bold")
    ax2.legend(fontsize=6.1, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.18)); panel(ax2, "b")

    im = draw_pathway_heatmap(ax3, coupling)
    cbar = fig.colorbar(im, ax=ax3, fraction=0.022, pad=0.02)
    cbar.set_label("Pooled Spearman ρ", fontsize=7.6); cbar.ax.tick_params(labelsize=7.0)
    ax3.set_title("Primary singscore pathway coupling at T48", loc="left", fontweight="bold"); panel(ax3, "c")
    fig.text(0.01, 0.017,
             "* Primary singscore BH-FDR < 0.05; the corresponding ssGSEA result had FDR=0.066.",
             fontsize=7.0, color="#444444")
    fig.subplots_adjust(bottom=0.18, top=0.97)
    export(fig, "Figure_5_gene_pathway_cell_integration")
    for name in [
        "04_07_gene_contribution_meta_analysis.csv", "04_12_meta_signature_pathway_coupling.csv",
        "04_12_meta_signature_pathway_coupling_ssgsea.csv", "04_13_cell_source_annotation.csv",
        "04_14_gene_contribution_cell_source_map.csv",
    ]:
        public_labels(pd.read_csv(S4 / name)).to_csv(SOURCE_OUT / f"Figure_5_{name}", index=False)


def figure6() -> None:
    meta = pd.read_csv(S3 / "03_05_signature_level_meta_analysis.csv")
    m = meta[(meta.analysis_set == "PRIMARY_INDEPENDENT") & meta.time_window.isin(["T1", "T2"])].pivot(index="signature_id", columns="time_window", values="pooled_delta_z").reindex(SIGS)
    arch = pd.read_csv(S4 / "04_08_signature_drift_architecture.csv").set_index("signature_id").reindex(SIGS)
    stability = pd.read_csv(S3 / "10_02_classification_stability_matrix.csv").set_index("signature_id").reindex(SIGS)
    fig = plt.figure(figsize=(7.2, 5.85))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.30, 1.08], wspace=0.38)
    ax = fig.add_subplot(gs[0]); panel(ax, "a")
    x = np.array([0, 1, 2])
    signature_colors = dict(zip(SIGS, mpl.colormaps["tab10"](np.linspace(0, 0.9, len(SIGS)))))
    marker_cycle = dict(zip(SIGS, ["o", "s", "^", "D", "v", "P", "X", "h"]))
    for sig in SIGS:
        vals = [0, m.loc[sig, "T1"], m.loc[sig, "T2"]]
        color = signature_colors[sig]
        status = stability.loc[sig, "sensitivity_status"]
        linestyle = {"STABLE": "-", "BOUNDARY_SENSITIVE": ":", "UNSTABLE": "--"}[status]
        ax.plot(x, vals, marker=marker_cycle[sig], linestyle=linestyle,
                lw=1.7, ms=4.3, color=color, alpha=0.95)
        label_dy = {"SIG001": 10, "SIG002": 2, "SIG003": -2, "SIG004": 2,
                    "SIG022": 2, "SIG023": 8, "SIG033": -8, "SIG034": -2}[sig]
        ax.annotate(sig, (2, vals[-1]), xytext=(5, label_dy), textcoords="offset points", va="center", fontsize=6.4, color=color, fontweight="bold")
    ax.axhline(0, color="#777777", lw=0.8, ls="--")
    ax.set_xticks(x); ax.set_xticklabels(["Baseline", "T24", "T48"])
    ax.set_ylabel("Pooled within-patient change (ΔZ)")
    ax.set_xlim(-0.08, 2.35)
    ax.set_ylim(-1.02, 0.78)
    ax.set_title("Pooled trajectories depend on sampling window", loc="left", fontweight="bold")
    handles = [
        mpl.lines.Line2D([], [], color="#444444", linestyle="-", lw=1.6,
                         label="Label stable across threshold perturbations"),
        mpl.lines.Line2D([], [], color="#444444", linestyle=":", lw=1.6,
                         label="Label changed in one scenario"),
        mpl.lines.Line2D([], [], color="#444444", linestyle="--", lw=1.6,
                         label="Label changed in ≥2 scenarios"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.47, -0.13),
              ncol=1, fontsize=5.9, columnspacing=1.0, handlelength=2.2)

    ax2 = fig.add_subplot(gs[1]); ax2.set_axis_off(); panel(ax2, "b")
    ax2.text(0.02, 0.98, "Elements of measurement context", transform=ax2.transAxes,
             va="top", fontweight="bold", fontsize=9.0)

    cards = [
        (0.02, 0.70, "Published formula", "Genes, coefficients\nand direction"),
        (0.55, 0.70, "Source baseline", "Study-specific\ninitial sample"),
        (0.02, 0.18, "Sampling window", "T24 and T48\nremain distinct"),
        (0.55, 0.18, "Cohort evidence", "Effect, CI and\nprediction interval"),
    ]
    for x0, y0, head, body in cards:
        ax2.add_patch(FancyBboxPatch((x0, y0), 0.42, 0.15,
                                     boxstyle="round,pad=0.012,rounding_size=0.018",
                                     facecolor="#F2F6F8", edgecolor="#C9D5DC", linewidth=0.8,
                                     transform=ax2.transAxes))
        ax2.text(x0 + 0.21, y0 + 0.104, head, transform=ax2.transAxes,
                 ha="center", va="center", fontweight="bold", fontsize=6.8)
        ax2.text(x0 + 0.21, y0 + 0.050, body, transform=ax2.transAxes,
                 ha="center", va="center", fontsize=5.9, linespacing=1.12)

    ax2.add_patch(FancyBboxPatch((0.25, 0.43), 0.50, 0.14,
                                 boxstyle="round,pad=0.014,rounding_size=0.022",
                                 facecolor="#E8F1F5", edgecolor=COL["blue"], linewidth=1.0,
                                 transform=ax2.transAxes))
    ax2.text(0.50, 0.505, "Longitudinal interpretation", transform=ax2.transAxes,
             ha="center", va="center", fontweight="bold", fontsize=7.3, color=COL["ink"])
    ax2.text(0.50, 0.458, "Score trajectory + gene-contribution pattern", transform=ax2.transAxes,
             ha="center", va="center", fontsize=5.7, color=COL["ink"])
    for start, end in [((0.23, 0.70), (0.36, 0.57)), ((0.76, 0.70), (0.64, 0.57)),
                       ((0.23, 0.33), (0.36, 0.43)), ((0.76, 0.33), (0.64, 0.43))]:
        ax2.annotate("", xy=end, xytext=start, xycoords=ax2.transAxes,
                     arrowprops=dict(arrowstyle="->", color="#6C7C86", lw=0.9))

    ax2.add_patch(FancyBboxPatch((0.02, 0.005), 0.95, 0.125,
                                 boxstyle="round,pad=0.012,rounding_size=0.015",
                                 facecolor="#FFF3D5", edgecolor="#D9C889", linewidth=0.8,
                                 transform=ax2.transAxes))
    ax2.text(0.045, 0.104, "Interpretive scope", transform=ax2.transAxes,
             va="center", fontweight="bold", fontsize=6.8)
    ax2.text(0.045, 0.050,
             "Repeated-measurement behavior only. Clinical monitoring, prognosis, treatment\nresponse, and causal mechanisms remain unestablished.",
             transform=ax2.transAxes, va="center", fontsize=5.8, linespacing=1.14)

    fig.suptitle("Measurement context of longitudinal transcriptomic scores",
                 x=0.02, y=0.985, ha="left", fontsize=12, fontweight="bold")
    fig.subplots_adjust(bottom=0.18, top=0.91)
    export(fig, "Figure_6_temporal_measurement_context")
    fig6_source = (m.reset_index()
                   .merge(arch[["drift_architecture"]].reset_index(), on="signature_id", how="left")
                   .merge(stability[["sensitivity_status"]].reset_index(), on="signature_id", how="left"))
    fig6_source["architecture_label"] = fig6_source["drift_architecture"].map(ARCH_SHORT)
    fig6_source.to_csv(SOURCE_OUT / "Figure_6_source_data.csv", index=False)


def supplementary_figure_1() -> None:
    meta = pd.read_csv(S3 / "03_05_signature_level_meta_analysis.csv")
    rows = []
    for window in ["T1", "T2"]:
        for sig in SIGS:
            pilot = meta[(meta.signature_id == sig) & (meta.time_window == window) & (meta.analysis_set == "PILOT_ONLY")]
            validation = meta[(meta.signature_id == sig) & (meta.time_window == window) &
                              meta.analysis_set.isin(["PRESPECIFIED_NON_PILOT_ONLY", "PRESPECIFIED_VALIDATION_ONLY"])]
            if not pilot.empty and not validation.empty:
                rows.append({"signature_id": sig, "time_window": window,
                             "pilot": pilot.iloc[0].pooled_delta_z,
                             "validation": validation.iloc[0].pooled_delta_z})
    data = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.7))
    limits = [float(np.nanmin(data[["pilot", "validation"]].to_numpy())),
              float(np.nanmax(data[["pilot", "validation"]].to_numpy()))]
    margin = max(0.2, (limits[1] - limits[0]) * 0.12); limits = [limits[0] - margin, limits[1] + margin]
    for ax, window in zip(axes, ["T1", "T2"]):
        sub = data[data.time_window == window]
        ax.plot(limits, limits, color="#A8A8A8", ls="--", lw=0.9)
        ax.axhline(0, color="#E0E0E0", lw=0.6); ax.axvline(0, color="#E0E0E0", lw=0.6)
        ax.scatter(sub.pilot, sub.validation, color=COL["red"], s=30)
        ordered = sub.sort_values("validation").copy(); label_y = ordered.validation.to_numpy(float).copy()
        for idx in range(1, len(label_y)): label_y[idx] = max(label_y[idx], label_y[idx - 1] + 0.095)
        overflow = label_y[-1] - (limits[1] - 0.04)
        if overflow > 0: label_y -= overflow
        for idx in range(len(label_y) - 2, -1, -1): label_y[idx] = min(label_y[idx], label_y[idx + 1] - 0.095)
        for (_, row), text_y in zip(ordered.iterrows(), label_y):
            ax.annotate(row.signature_id, xy=(row.pilot, row.validation), xytext=(row.pilot + 0.07, text_y),
                        fontsize=6.4, ha="left", va="center",
                        arrowprops=dict(arrowstyle="-", color="#7F8C99", lw=0.45, shrinkA=1.5, shrinkB=2.5))
        ax.set_xlim(limits); ax.set_ylim(limits); ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Pilot pooled ΔZ"); ax.set_ylabel("Prespecified non-pilot pooled ΔZ")
        ax.set_title("T24" if window == "T1" else "T48", fontweight="bold", fontsize=9)
    panel(axes[0], "a"); panel(axes[1], "b")
    fig.suptitle("Pilot and prespecified non-pilot concordance", x=0.02, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    export_to(fig, SUPP_OUT / "Supplementary_Figure_1_pilot_vs_validation")
    data.to_csv(SOURCE_OUT / "Supplementary_Figure_1_source_data.csv", index=False)


def copy_unchanged_supplementary_figures() -> None:
    source_dir = SUPP_OUT
    for number in range(2, 7):
        for path in source_dir.glob(f"Supplementary_Figure_{number}_*"):
            shutil.copy2(path, SUPP_OUT / path.name)


def main() -> None:
    functions = {"1": figure1, "2": figure2, "3": figure3,
                 "4": figure4, "5": figure5, "6": figure6}
    selected = {x.strip() for x in os.environ.get(
        "STAGE6_FIGURE_SELECTION", "1,2,3,4,5,6"
    ).split(",") if x.strip()}
    for number, function in functions.items():
        if number in selected:
            function()
    manifest = {
        "freeze_id": "SCIREP-ANALYSIS-v1.2.0-20260722",
        "scientific_source": "Frozen Stage 3 results and corrected v1.2.0 Stage 4 architecture summaries using signature-specific primary-independent cohort sets",
        "figures": [p.name for p in sorted(OUT.glob("Figure_*.*"))],
        "source_files": [p.name for p in sorted(SOURCE_OUT.glob("*"))],
        "backend": "Python/matplotlib only",
        "export": {"SVG": "editable text", "PDF": "fonttype 42", "TIFF": "600 dpi LZW", "PNG": "300 dpi preview"},
    }
    (OUT / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
