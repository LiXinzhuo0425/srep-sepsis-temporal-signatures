#!/usr/bin/env python3
"""Generate all Stage 3 publication figures with the frozen Python backend."""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd


ROOT = Path("/Users/felix/Documents/New project/longitudinal_stage3")
SOURCE = ROOT / "source_data"
MAIN = ROOT / "03_15_main_figures"
FORESTS = ROOT / "03_05_forest_plots"
ATTRITION = ROOT / "03_09_attrition_flowcharts"
TRAJECTORIES = ROOT / "03_12_longitudinal_trajectory_plots"
QA = ROOT / "qa/figures"

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 8
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["legend.frameon"] = False

COLORS = {
    "T1": "#3775BA",
    "T2": "#B64342",
    "pilot": "#A8A8A8",
    "validation": "#0F4D92",
    "all": "#767676",
    "primary": "#0F4D92",
    "strict": "#B64342",
    "mad": "#42949E",
    "neutral": "#D8D8D8",
    "ink": "#272727",
}
COHORT_COLORS = {
    "GSE236713": "#0F4D92", "GSE57065": "#8BCF8B", "GSE95233": "#B64342",
    "GSE54514": "#42949E", "GSE110487": "#9A4D8E", "GSE8121": "#CF9D3E",
}
SIGNATURES = ["SIG001", "SIG002", "SIG003", "SIG004", "SIG022", "SIG023", "SIG033", "SIG034"]


def panel_label(ax, label):
    ax.text(-0.08, 1.03, label, transform=ax.transAxes, ha="left", va="bottom", fontweight="bold", fontsize=10)


def save(fig, base: Path, *, height=None):
    base.parent.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def draw_box(ax, xy, width, height, text, face, edge="#D4DEE7", fontsize=8):
    box = FancyBboxPatch(xy, width, height, boxstyle="round,pad=0.02,rounding_size=0.02", facecolor=face, edgecolor=edge, linewidth=0.8)
    ax.add_patch(box)
    ax.text(xy[0] + width / 2, xy[1] + height / 2, text, ha="center", va="center", fontsize=fontsize)


def figure1_flow():
    cfg = json.loads((ROOT / "03_00_09_environment_lock/analysis_config.json").read_text())
    cohorts = pd.DataFrame(cfg["cohorts"])
    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], wspace=0.44, hspace=0.28)
    ax = fig.add_subplot(gs[0, :]); ax.set_axis_off(); panel_label(ax, "a")
    draw_box(ax, (0.02, 0.61), 0.22, 0.24, "6 independent\ncohorts", "#EAF1F6", fontsize=10)
    draw_box(ax, (0.30, 0.61), 0.22, 0.24, "302 baseline\ncase patients", "#EAF1F6", fontsize=10)
    draw_box(ax, (0.58, 0.61), 0.18, 0.24, "264 with\nT24/T48 pair", "#E8F3EC", fontsize=10)
    draw_box(ax, (0.82, 0.61), 0.16, 0.24, "8 fixed A2\nsignatures", "#FFF2CC", fontsize=10)
    for x1, x2 in ((0.24, 0.30), (0.52, 0.58), (0.76, 0.82)):
        ax.annotate("", xy=(x2, 0.73), xytext=(x1, 0.73), arrowprops=dict(arrowstyle="->", color="#5F6B76", lw=1.2))
    ax.text(0.02, 0.36, "Platforms", fontweight="bold")
    ax.text(0.15, 0.36, "Agilent  |  Affymetrix  |  Illumina  |  RNA-seq", color="#5F6B76")
    ax.text(0.02, 0.18, "Primary windows", fontweight="bold")
    ax.text(0.22, 0.18, "T0 to 24 h (T1)   and   T0 to 48 h (T2)", color="#5F6B76")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    axb = fig.add_subplot(gs[1, 0]); panel_label(axb, "b")
    plot = cohorts.set_index("dataset")[["patients_T0_T1", "patients_T0_T2"]]
    x = np.arange(len(plot)); width = 0.36
    axb.bar(x - width/2, plot["patients_T0_T1"], width, color=COLORS["T1"], label="T24")
    axb.bar(x + width/2, plot["patients_T0_T2"], width, color=COLORS["T2"], label="T48")
    axb.set_xticks(x); axb.set_xticklabels(plot.index, rotation=45, ha="right")
    axb.set_ylabel("Paired patients"); axb.legend(ncol=2, loc="upper right", bbox_to_anchor=(0.86, 1.0))

    axc = fig.add_subplot(gs[1, 1]); panel_label(axc, "c")
    groups = [4, 2, 1, 1]; labels = ["Sepsis/infection\ndiagnostic", "Bacterial/viral\ndiagnostic", "Prognosis", "Severity"]
    colors = ["#0F4D92", "#42949E", "#9A4D8E", "#CF9D3E"]
    axc.barh(np.arange(4)[::-1], groups, color=colors)
    axc.set_yticks(np.arange(4)[::-1]); axc.set_yticklabels(labels)
    axc.set_xlabel("Frozen signatures")
    axc.set_xlim(0, 4.6)
    for y, val in zip(np.arange(4)[::-1], groups): axc.text(val + 0.08, y, str(val), va="center")
    fig.suptitle("Cohort, patient and signature analysis flow", x=0.02, ha="left", fontsize=11, fontweight="bold")
    fig.subplots_adjust(top=0.91)
    save(fig, MAIN / "Figure_1_analysis_flow")


def figure2_primary_forest():
    meta = pd.read_csv(SOURCE / "03_05_signature_level_meta_analysis.csv")
    primary = meta[(meta["analysis_set"] == "PRIMARY_INDEPENDENT") & meta["time_window"].isin(["T1", "T2"])].copy()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.4), sharey=True)
    for ax, window, label in zip(axes, ["T1", "T2"], ["T0 to 24 h", "T0 to 48 h"]):
        data = primary[primary["time_window"] == window].set_index("signature_id").reindex(SIGNATURES)
        y = np.arange(len(SIGNATURES))[::-1]
        ax.axvline(0, color="#767676", ls="--", lw=0.9)
        for yi, sig in zip(y, SIGNATURES):
            row = data.loc[sig]
            if pd.isna(row["pooled_delta_z"]): continue
            ax.plot([row["ci95_lower"], row["ci95_upper"]], [yi, yi], color=COLORS[window], lw=1.6)
            ax.scatter(row["pooled_delta_z"], yi, color=COLORS[window], s=22, zorder=3)
            if math.isfinite(row["prediction_lower"]):
                ax.plot([row["prediction_lower"], row["prediction_upper"]], [yi-0.17, yi-0.17], color=COLORS[window], lw=0.8, alpha=0.45)
        ax.set_yticks(y); ax.set_yticklabels(SIGNATURES)
        ax.set_xlabel("Pooled within-patient change (deltaZ)")
        ax.set_title(label, fontsize=9, fontweight="bold")
        ax.text(0.02, 0.02, "Thick: 95% CI\nThin: 95% prediction interval", transform=ax.transAxes, fontsize=6.5, color="#5F6B76")
    panel_label(axes[0], "a"); panel_label(axes[1], "b")
    fig.suptitle("Primary within-patient score changes", x=0.02, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, MAIN / "Figure_2_primary_effects")


def figure3_stability_heatmap():
    profile = pd.read_csv(SOURCE / "03_06_temporal_stability_profile.csv")
    index = pd.MultiIndex.from_product([SIGNATURES, ["T1", "T2"]], names=["signature_id", "time_window"])
    p = profile.set_index(["signature_id", "time_window"]).reindex(index)
    drift = np.column_stack([p.xs("T1", level="time_window")["pooled_delta_z"], p.xs("T2", level="time_window")["pooled_delta_z"]])
    i2 = np.column_stack([p.xs("T1", level="time_window")["I2_percent"], p.xs("T2", level="time_window")["I2_percent"]])
    piw = np.column_stack([(p.xs("T1", level="time_window")["prediction_upper"] - p.xs("T1", level="time_window")["prediction_lower"]), (p.xs("T2", level="time_window")["prediction_upper"] - p.xs("T2", level="time_window")["prediction_lower"])])
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 4.5), gridspec_kw={"width_ratios": [1, 1, 1]})
    specs = [(drift, "Pooled deltaZ", "RdBu_r", -np.nanmax(abs(drift)), np.nanmax(abs(drift))), (i2, "I-squared (%)", "Blues", 0, 100), (piw, "Prediction interval width", "Oranges", 0, np.nanmax(piw))]
    for idx, (ax, (matrix, title, cmap, vmin, vmax)) in enumerate(zip(axes, specs)):
        im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["T24", "T48"])
        ax.set_yticks(np.arange(len(SIGNATURES))); ax.set_yticklabels(SIGNATURES if idx == 0 else [])
        ax.set_title(title, fontsize=8.5, fontweight="bold")
        for i in range(matrix.shape[0]):
            for j in range(2):
                if np.isfinite(matrix[i, j]): ax.text(j, i, f"{matrix[i,j]:.2f}" if idx != 1 else f"{matrix[i,j]:.0f}", ha="center", va="center", fontsize=6.5)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        panel_label(ax, chr(ord("a") + idx))
    fig.suptitle("Temporal stability and between-cohort heterogeneity", x=0.02, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, ROOT / "03_06_stability_heatmap")
    # Main-figure copy is generated from the same figure source and values.
    for suffix in (".svg", ".pdf", ".tiff", ".png"):
        target = MAIN / f"Figure_3_stability_heatmap{suffix}"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((ROOT / f"03_06_stability_heatmap{suffix}").read_bytes())


def figure4_pilot_validation():
    meta = pd.read_csv(SOURCE / "03_05_signature_level_meta_analysis.csv")
    rows = []
    for window in ["T1", "T2"]:
        for sig in SIGNATURES:
            pilot = meta[(meta.signature_id == sig) & (meta.time_window == window) & (meta.analysis_set == "PILOT_ONLY")]
            validation = meta[(meta.signature_id == sig) & (meta.time_window == window) & (meta.analysis_set == "BLINDED_VALIDATION_ONLY")]
            if not pilot.empty and not validation.empty:
                rows.append({"signature_id": sig, "time_window": window, "pilot": pilot.iloc[0].pooled_delta_z, "validation": validation.iloc[0].pooled_delta_z})
    data = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.5))
    limits = [float(np.nanmin(data[["pilot", "validation"]].to_numpy())), float(np.nanmax(data[["pilot", "validation"]].to_numpy()))]
    margin = max(0.2, (limits[1] - limits[0]) * 0.12); limits = [limits[0] - margin, limits[1] + margin]
    for ax, window in zip(axes, ["T1", "T2"]):
        sub = data[data.time_window == window]
        ax.plot(limits, limits, color="#A8A8A8", ls="--", lw=0.9)
        ax.axhline(0, color="#E0E0E0", lw=0.6); ax.axvline(0, color="#E0E0E0", lw=0.6)
        ax.scatter(sub.pilot, sub.validation, color=COLORS["validation"], s=28)
        ordered = sub.sort_values("validation").copy()
        label_y = ordered["validation"].to_numpy(float).copy()
        minimum_gap = 0.095
        for idx in range(1, len(label_y)):
            label_y[idx] = max(label_y[idx], label_y[idx - 1] + minimum_gap)
        overflow = label_y[-1] - (limits[1] - 0.04)
        if overflow > 0:
            label_y -= overflow
        for idx in range(len(label_y) - 2, -1, -1):
            label_y[idx] = min(label_y[idx], label_y[idx + 1] - minimum_gap)
        for (_, row), text_y in zip(ordered.iterrows(), label_y):
            text_x = row.pilot + 0.07
            ax.annotate(
                row.signature_id, xy=(row.pilot, row.validation), xytext=(text_x, text_y), textcoords="data",
                fontsize=6.1, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color="#7F8C99", lw=0.45, shrinkA=1.5, shrinkB=2.5),
            )
        ax.set_xlim(limits); ax.set_ylim(limits); ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Pilot pooled deltaZ"); ax.set_ylabel("Blinded-validation pooled deltaZ")
        ax.set_title("T24" if window == "T1" else "T48", fontweight="bold", fontsize=9)
    panel_label(axes[0], "a"); panel_label(axes[1], "b")
    fig.suptitle("Pilot and blinded-validation concordance", x=0.02, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, MAIN / "Figure_4_pilot_vs_validation")


def figure5_robustness():
    meta = pd.read_csv(SOURCE / "03_05_signature_level_meta_analysis.csv")
    scaling = pd.read_csv(SOURCE / "03_11_scaling_sensitivity_analysis.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.6), sharey=True)
    sets = [("ALL_COHORTS", "All", COLORS["all"]), ("PRIMARY_INDEPENDENT", "Overlap-excluded", COLORS["primary"]), ("STRICT_NEVER_USED", "Strict never-used", COLORS["strict"])]
    offsets = [-0.18, -0.06, 0.06, 0.18]
    for ax, window in zip(axes, ["T1", "T2"]):
        y = np.arange(len(SIGNATURES))[::-1]
        ax.axvline(0, color="#767676", ls="--", lw=0.8)
        for offset, (analysis_set, label, color) in zip(offsets[:3], sets):
            data = meta[(meta.time_window == window) & (meta.analysis_set == analysis_set)].set_index("signature_id").reindex(SIGNATURES)
            for yi, sig in zip(y, SIGNATURES):
                row = data.loc[sig]
                if pd.isna(row.pooled_delta_z): continue
                ax.plot([row.ci95_lower, row.ci95_upper], [yi+offset, yi+offset], color=color, lw=1.0)
                ax.scatter(row.pooled_delta_z, yi+offset, color=color, s=12, label=label if yi == y[0] else None)
        mad = scaling[(scaling.record_type == "META_PRIMARY_INDEPENDENT") & (scaling.time_window == window) & (scaling.scaling_method == "BASELINE_MAD")].set_index("signature_id").reindex(SIGNATURES)
        for yi, sig in zip(y, SIGNATURES):
            row = mad.loc[sig]
            if pd.isna(row.mean_change): continue
            ax.plot([row.ci95_lower, row.ci95_upper], [yi+offsets[3], yi+offsets[3]], color=COLORS["mad"], lw=1.0)
            ax.scatter(row.mean_change, yi+offsets[3], color=COLORS["mad"], marker="D", s=12, label="MAD scaling" if yi == y[0] else None)
        ax.set_yticks(y); ax.set_yticklabels(SIGNATURES)
        ax.set_xlabel("Pooled within-patient change")
        ax.set_title("T24" if window == "T1" else "T48", fontweight="bold", fontsize=9)
    axes[1].legend(loc="lower right", fontsize=6.5)
    panel_label(axes[0], "a"); panel_label(axes[1], "b")
    fig.suptitle("Robustness to overlap exclusions and scaling", x=0.02, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, MAIN / "Figure_5_robustness")


def figure6_trajectories():
    grades = pd.read_csv(SOURCE / "03_14_signature_evidence_grade.csv")
    changes = pd.read_parquet(ROOT / "03_03_longitudinal_analysis_datasets/all_windows.parquet")
    priority = ["HIGH_CONSISTENCY", "MODERATE_CONSISTENCY", "COHORT_DEPENDENT", "UNSTABLE", "EVIDENCE_INSUFFICIENT"]
    selected = []
    for grade in priority:
        candidates = grades[grades.evidence_grade == grade].groupby("signature_id")["primary_patient_n"].sum().sort_values(ascending=False)
        if len(candidates):
            sig = sorted(candidates[candidates == candidates.max()].index)[0]
            if sig not in selected: selected.append(sig)
        if len(selected) == 4: break
    for sig in SIGNATURES:
        if sig not in selected: selected.append(sig)
        if len(selected) == 4: break
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.8), sharex=True)
    for ax, sig in zip(axes.ravel(), selected):
        group = changes[(changes.signature_id == sig) & (changes.time_window.isin(["T1", "T2"]))]
        for dataset, cohort in group.groupby("dataset"):
            for _, row in cohort.iterrows():
                x = [0, 1 if row.time_window == "T1" else 2]
                ax.plot(x, [0, row.delta_z], color=COHORT_COLORS[dataset], alpha=0.12, lw=0.6)
            means = cohort.groupby("time_window").delta_z.mean()
            for window, xpos in (("T1", 1), ("T2", 2)):
                if window in means: ax.scatter(xpos, means[window], color=COHORT_COLORS[dataset], s=20, edgecolor="white", linewidth=0.4, zorder=4)
        ax.axhline(0, color="#767676", ls="--", lw=0.7)
        ax.set_title(sig, fontweight="bold", fontsize=9)
        ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["T0", "T24", "T48"])
        ax.set_ylabel("Within-patient deltaZ")
    for idx, ax in enumerate(axes.ravel()): panel_label(ax, chr(ord("a") + idx))
    handles = [plt.Line2D([0], [0], color=color, lw=2, label=dataset) for dataset, color in COHORT_COLORS.items()]
    fig.legend(handles=handles, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 0.95), fontsize=6.5)
    fig.suptitle("Patient-level trajectories from pre-coded evidence archetypes", x=0.02, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    save(fig, MAIN / "Figure_6_patient_trajectories")


def supporting_figures():
    effects = pd.read_csv(SOURCE / "03_04_cohort_signature_primary_effects.csv")
    meta = pd.read_csv(SOURCE / "03_05_signature_level_meta_analysis.csv")
    cohort_order = list(COHORT_COLORS)
    cohort_y = {dataset: y for dataset, y in zip(cohort_order, np.arange(len(cohort_order))[::-1])}
    for sig in SIGNATURES:
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.8), sharey=True)
        for ax, window in zip(axes, ["T1", "T2"]):
            sub = effects[(effects.signature_id == sig) & (effects.time_window == window)].copy().set_index("dataset")
            ax.axvline(0, color="#767676", ls="--", lw=0.8)
            for dataset in cohort_order:
                if dataset not in sub.index:
                    continue
                row = sub.loc[dataset]
                yi = cohort_y[dataset]
                ax.plot([row.ci95_lower, row.ci95_upper], [yi, yi], color="#7884B4", lw=1.2)
                ax.scatter(row.mean_delta_z, yi, color="#484878", s=18)
            pooled = meta[(meta.signature_id == sig) & (meta.time_window == window) & (meta.analysis_set == "PRIMARY_INDEPENDENT")]
            if not pooled.empty:
                row = pooled.iloc[0]; ax.axvspan(row.ci95_lower, row.ci95_upper, color="#F0C0CC", alpha=0.35); ax.axvline(row.pooled_delta_z, color="#B64342", lw=1.2)
            ax.set_yticks([cohort_y[d] for d in cohort_order]); ax.set_yticklabels(cohort_order)
            ax.set_ylim(-0.5, len(cohort_order) - 0.5)
            ax.set_title("T24" if window == "T1" else "T48", fontsize=9, fontweight="bold")
            ax.set_xlabel("deltaZ")
        fig.suptitle(f"{sig}: cohort-specific and pooled effects", x=0.02, ha="left", fontsize=10, fontweight="bold")
        fig.tight_layout(rect=(0, 0, 1, 0.93)); save(fig, FORESTS / f"{sig}_forest")

    attr = pd.read_csv(SOURCE / "03_09_attrition_and_missingness_analysis.csv")
    for dataset in sorted(attr.dataset.unique()):
        sub = attr[(attr.dataset == dataset) & (attr.signature_id == "SIG001")].set_index("time_window")
        fig, ax = plt.subplots(figsize=(4.2, 2.8))
        windows = [w for w in ["T1", "T2"] if w in sub.index]
        observed = [sub.loc[w, "observed_pair_n"] for w in windows]
        missing = [sub.loc[w, "missing_followup_n"] for w in windows]
        x = np.arange(len(windows))
        ax.bar(x, observed, color="#3775BA", label="Observed pair")
        ax.bar(x, missing, bottom=observed, color="#D8D8D8", label="Missing/unknown")
        ax.set_xticks(x); ax.set_xticklabels(["T24" if w == "T1" else "T48" for w in windows])
        ax.set_ylabel("Baseline-eligible patients"); ax.set_title(dataset, fontweight="bold")
        ax.legend(fontsize=7); fig.tight_layout(); save(fig, ATTRITION / f"{dataset}_attrition")

    changes = pd.read_parquet(ROOT / "03_03_longitudinal_analysis_datasets/all_windows.parquet")
    for sig in SIGNATURES:
        sub = changes[changes.signature_id == sig]
        summary = sub.groupby(["dataset", "time_window"]).delta_z.agg(["mean", "sem"]).reset_index()
        fig, ax = plt.subplots(figsize=(5.2, 3.4))
        for dataset, group in summary.groupby("dataset"):
            order = [w for w in ["T1", "T2", "T3", "T4"] if w in set(group.time_window)]
            g = group.set_index("time_window").loc[order]
            x = [TARGET for TARGET in [1,2,3,4] if ["T1","T2","T3","T4"][TARGET-1] in order]
            ax.errorbar(x, g["mean"], yerr=1.96*g["sem"], color=COHORT_COLORS[dataset], marker="o", lw=1, ms=3, label=dataset)
        ax.axhline(0, color="#767676", ls="--", lw=0.7); ax.set_xticks([1,2,3,4]); ax.set_xticklabels(["T24","T48","T72","Day 5"])
        ax.set_ylabel("Mean deltaZ (95% CI)"); ax.set_title(sig, fontweight="bold"); ax.legend(fontsize=5.5, ncol=2)
        fig.tight_layout(); save(fig, TRAJECTORIES / f"{sig}_extended_trajectory")


def main():
    for path in (MAIN, FORESTS, ATTRITION, TRAJECTORIES, QA): path.mkdir(parents=True, exist_ok=True)
    figure1_flow()
    figure2_primary_forest()
    figure3_stability_heatmap()
    figure4_pilot_validation()
    figure5_robustness()
    figure6_trajectories()
    supporting_figures()
    print("figures generated")


if __name__ == "__main__":
    main()
