#!/usr/bin/env python3
"""Frozen Stage 4 gene-contribution decomposition and meta-analysis.

This module reads only Stage 3 frozen expression, score, pairing and overlap
artifacts. It never refits a signature, changes a time window, selects genes,
or uses clinical outcomes. Nonlinear scores are decomposed with exact Shapley
or fixed-group Owen values, not with a local linear approximation.
"""

from __future__ import annotations

import os

import hashlib
import html
import itertools
import json
import math
import sys
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(os.environ.get("SEPSIS_SIGNATURE_ANALYSIS_ROOT", Path.cwd())).resolve()
ENV3 = ROOT / "03_00_09_environment_lock"
ENV4 = ROOT / "04_environment"
SOURCE = ROOT / "04_source_data"
sys.path.insert(0, str(ENV3))

from signatures import (  # noqa: E402
    LIN7_WEIGHTS,
    SCORE_DIRECTION,
    SMS_DOWN,
    SMS_UP,
    SOM_M1,
    SOM_M2,
    SOM_M3,
    SOM_M4,
    STAGE3_FUNCTIONS,
    STAGE3_REQUIRED_GENES,
    SWEENEY_BV_BACTERIAL,
    SWEENEY_BV_VIRAL,
)
from run_stage3 import reml_meta  # noqa: E402


CONFIG = json.loads((ENV4 / "stage4_config.json").read_text(encoding="utf-8"))
SEED = int(CONFIG["seed"])
BOOTSTRAPS = int(CONFIG["bootstraps"])
TOL = float(CONFIG["reconstruction_absolute_tolerance"])
PRIMARY_WINDOWS = tuple(CONFIG["primary_windows"])
ALL_WINDOWS = PRIMARY_WINDOWS + tuple(CONFIG["extended_windows"])

LINEAR_COEFFICIENTS = {
    "SIG002": {"PLAC8": 1.0, "LAMP1": 1.0, "PLA2G7": -1.0, "CEACAM4": -1.0},
    "SIG023": {"FAM89A": 1.0, "IFI44L": -1.0},
    "SIG033": {gene: float(weight) for gene, weight in LIN7_WEIGHTS.items()},
}

GENE_GROUPS = {
    "SIG001": {"UP_GEOMETRIC_MEAN": tuple(SMS_UP), "DOWN_GEOMETRIC_MEAN": tuple(SMS_DOWN)},
    "SIG002": {"LINEAR_PANEL": tuple(STAGE3_REQUIRED_GENES["SIG002"])},
    "SIG003": {"RATIO_PANEL": tuple(STAGE3_REQUIRED_GENES["SIG003"])},
    "SIG004": {"RATIO_PANEL": tuple(STAGE3_REQUIRED_GENES["SIG004"])},
    "SIG022": {"VIRAL_GEOMETRIC_MEAN": tuple(SWEENEY_BV_VIRAL), "BACTERIAL_GEOMETRIC_MEAN": tuple(SWEENEY_BV_BACTERIAL)},
    "SIG023": {"LINEAR_PANEL": tuple(STAGE3_REQUIRED_GENES["SIG023"])},
    "SIG033": {"LINEAR_PANEL": tuple(STAGE3_REQUIRED_GENES["SIG033"])},
    "SIG034": {"MODULE_1": tuple(SOM_M1), "MODULE_2": tuple(SOM_M2), "MODULE_3": tuple(SOM_M3), "MODULE_4": tuple(SOM_M4)},
}

STRUCTURE_ROWS = [
    {"signature_id": "SIG001", "formula_structure": "weighted difference of two geometric means", "exact_method": "group-wise exact Shapley", "formula": "GM(6 up genes) - (5/6) * GM(5 down genes)", "intercept_in_change": "NO", "linear_approximation": "PROHIBITED", "exactly_decomposable": "YES"},
    {"signature_id": "SIG002", "formula_structure": "linear signed sum", "exact_method": "analytic coefficient times expression change", "formula": "PLAC8 + LAMP1 - PLA2G7 - CEACAM4", "intercept_in_change": "NO", "linear_approximation": "NOT_APPLICABLE", "exactly_decomposable": "YES"},
    {"signature_id": "SIG003", "formula_structure": "two-gene ratio", "exact_method": "exact two-point Shapley", "formula": "FAIM3 / PLAC8; orientation multiplied by -1", "intercept_in_change": "NO", "linear_approximation": "PROHIBITED", "exactly_decomposable": "YES"},
    {"signature_id": "SIG004", "formula_structure": "ratio of expression difference", "exact_method": "exact two-point Shapley", "formula": "(NLRP1 - IDNK) / PLAC8", "intercept_in_change": "NO", "linear_approximation": "PROHIBITED", "exactly_decomposable": "YES"},
    {"signature_id": "SIG022", "formula_structure": "difference of two geometric means", "exact_method": "group-wise exact Shapley", "formula": "GM(3 viral genes) - GM(4 bacterial genes)", "intercept_in_change": "NO", "linear_approximation": "PROHIBITED", "exactly_decomposable": "YES"},
    {"signature_id": "SIG023", "formula_structure": "two-gene expression difference", "exact_method": "analytic coefficient times expression change", "formula": "FAM89A - IFI44L", "intercept_in_change": "NO", "linear_approximation": "NOT_APPLICABLE", "exactly_decomposable": "YES"},
    {"signature_id": "SIG033", "formula_structure": "linear weighted sum", "exact_method": "analytic coefficient times expression change", "formula": "sum of 7 published beta_g * X_g", "intercept_in_change": "NO", "linear_approximation": "NOT_APPLICABLE", "exactly_decomposable": "YES"},
    {"signature_id": "SIG034", "formula_structure": "ratio of sums of four module geometric means", "exact_method": "exact Owen value under the frozen four-module partition", "formula": "[GM(M1)+GM(M2)]/[GM(M3)+GM(M4)]", "intercept_in_change": "NO", "linear_approximation": "PROHIBITED", "exactly_decomposable": "YES"},
]


def stable_seed(*parts: str) -> int:
    digest = hashlib.sha256("|".join(parts).encode()).hexdigest()
    return (SEED + int(digest[:8], 16)) % (2**32 - 1)


def popcount(value: int) -> int:
    """Python-version-independent integer population count."""

    return bin(value).count("1")


@lru_cache(None)
def factorial_weight(n: int, subset_size: int) -> float:
    return math.factorial(subset_size) * math.factorial(n - subset_size - 1) / math.factorial(n)


def geometric_mean(values: list[float] | np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    if np.any(~np.isfinite(array)) or np.any(array <= 0):
        raise ValueError("geometric mean attribution requires finite positive expression")
    return float(np.exp(np.mean(np.log(array))))


def exact_shapley(genes: tuple[str, ...], x0: dict[str, float], x1: dict[str, float], function) -> dict[str, float]:
    """Exact two-point Shapley decomposition of f(x1)-f(x0)."""

    n = len(genes)
    values: dict[int, float] = {}
    for mask in range(1 << n):
        state = {gene: (x1[gene] if mask & (1 << index) else x0[gene]) for index, gene in enumerate(genes)}
        values[mask] = float(function(state))
    phi = {gene: 0.0 for gene in genes}
    for index, gene in enumerate(genes):
        bit = 1 << index
        for mask in range(1 << n):
            if mask & bit:
                continue
            phi[gene] += factorial_weight(n, popcount(mask)) * (values[mask | bit] - values[mask])
    return phi


def shapley_geometric_mean(genes: tuple[str, ...], x0: dict[str, float], x1: dict[str, float]) -> dict[str, float]:
    return exact_shapley(genes, x0, x1, lambda state: geometric_mean([state[g] for g in genes]))


def owen_sig034(x0: dict[str, float], x1: dict[str, float]) -> dict[str, float]:
    """Exact Owen values for the frozen four-module partition of SIG034."""

    groups = [tuple(SOM_M1), tuple(SOM_M2), tuple(SOM_M3), tuple(SOM_M4)]
    m = len(groups)
    baseline_modules = [geometric_mean([x0[g] for g in group]) for group in groups]
    follow_modules = [geometric_mean([x1[g] for g in group]) for group in groups]

    def score(modules: list[float]) -> float:
        denominator = modules[2] + modules[3]
        if denominator == 0:
            raise ValueError("SIG034 denominator is zero in an Owen coalition")
        return (modules[0] + modules[1]) / denominator

    out = {gene: 0.0 for group in groups for gene in group}
    for group_index, genes in enumerate(groups):
        n = len(genes)
        base_logs = np.log(np.asarray([x0[g] for g in genes], dtype=float))
        follow_logs = np.log(np.asarray([x1[g] for g in genes], dtype=float))
        delta_logs = follow_logs - base_logs
        gm_cache = np.empty(1 << n, dtype=float)
        base_sum = float(base_logs.sum())
        for inner_mask in range(1 << n):
            delta_sum = sum(delta_logs[idx] for idx in range(n) if inner_mask & (1 << idx))
            gm_cache[inner_mask] = math.exp((base_sum + float(delta_sum)) / n)

        other_groups = [idx for idx in range(m) if idx != group_index]
        for outer_mask in range(1 << (m - 1)):
            outer_size = popcount(outer_mask)
            outer_weight = factorial_weight(m, outer_size)
            modules = baseline_modules.copy()
            for pos, idx in enumerate(other_groups):
                if outer_mask & (1 << pos):
                    modules[idx] = follow_modules[idx]
            values = np.empty(1 << n, dtype=float)
            for inner_mask in range(1 << n):
                modules[group_index] = float(gm_cache[inner_mask])
                values[inner_mask] = score(modules)
            for gene_index, gene in enumerate(genes):
                bit = 1 << gene_index
                subtotal = 0.0
                for inner_mask in range(1 << n):
                    if inner_mask & bit:
                        continue
                    subtotal += factorial_weight(n, popcount(inner_mask)) * (values[inner_mask | bit] - values[inner_mask])
                out[gene] += outer_weight * subtotal
    return out


def decompose(signature: str, x0: dict[str, float], x1: dict[str, float]) -> dict[str, float]:
    genes = tuple(STAGE3_REQUIRED_GENES[signature])
    if signature in LINEAR_COEFFICIENTS:
        return {gene: LINEAR_COEFFICIENTS[signature][gene] * (x1[gene] - x0[gene]) for gene in genes}
    if signature == "SIG001":
        up = shapley_geometric_mean(tuple(SMS_UP), x0, x1)
        down = shapley_geometric_mean(tuple(SMS_DOWN), x0, x1)
        return {**up, **{gene: -(5.0 / 6.0) * value for gene, value in down.items()}}
    if signature == "SIG022":
        viral = shapley_geometric_mean(tuple(SWEENEY_BV_VIRAL), x0, x1)
        bacterial = shapley_geometric_mean(tuple(SWEENEY_BV_BACTERIAL), x0, x1)
        return {**viral, **{gene: -value for gene, value in bacterial.items()}}
    if signature == "SIG003":
        return exact_shapley(genes, x0, x1, lambda s: s["FAIM3"] / s["PLAC8"])
    if signature == "SIG004":
        return exact_shapley(genes, x0, x1, lambda s: (s["NLRP1"] - s["IDNK"]) / s["PLAC8"])
    if signature == "SIG034":
        return owen_sig034(x0, x1)
    raise KeyError(signature)


def gene_group(signature: str, gene: str) -> str:
    for group, genes in GENE_GROUPS[signature].items():
        if gene in genes:
            return group
    raise KeyError((signature, gene))


def coefficient(signature: str, gene: str) -> float:
    if signature in LINEAR_COEFFICIENTS:
        return LINEAR_COEFFICIENTS[signature][gene]
    return math.nan


def load_expressions() -> dict[str, pd.DataFrame]:
    out = {}
    for dataset in CONFIG["cohorts"]:
        frame = pd.read_parquet(ROOT / f"intermediate/required_expression_{dataset}.parquet")
        pivot = frame.pivot(index="gene", columns="sample_id", values="expression")
        out[dataset] = pivot
    return out


def overlap_map() -> dict[tuple[str, str], str]:
    score = pd.read_parquet(ROOT / "03_02_full_signature_score_matrix.parquet")
    unique = score[["signature_id", "dataset", "development_overlap_status"]].drop_duplicates()
    return {(r.signature_id, r.dataset): r.development_overlap_status for r in unique.itertuples()}


def build_gene_level() -> tuple[pd.DataFrame, pd.DataFrame]:
    changes = pd.read_parquet(ROOT / "03_03_longitudinal_analysis_datasets/all_windows.parquet")
    changes = changes[changes["time_window"].isin(ALL_WINDOWS)].copy()
    expressions = load_expressions()
    overlaps = overlap_map()
    rows: list[dict] = []
    qc_rows: list[dict] = []

    for change in changes.itertuples(index=False):
        signature = change.signature_id
        genes = tuple(STAGE3_REQUIRED_GENES[signature])
        matrix = expressions[change.dataset]
        missing = [g for g in genes if g not in matrix.index]
        if missing:
            raise RuntimeError(f"{change.dataset} {signature}: missing genes {missing}")
        x0 = {g: float(matrix.at[g, change.baseline_sample_id]) for g in genes}
        x1 = {g: float(matrix.at[g, change.followup_sample_id]) for g in genes}
        raw = decompose(signature, x0, x1)
        direction = int(SCORE_DIRECTION[signature])
        raw_sum = float(sum(raw.values()))
        oriented_sum = raw_sum * direction
        delta_sum = oriented_sum / float(change.baseline_sd)
        raw_error = abs(oriented_sum - float(change.raw_change))
        z_error = abs(delta_sum - float(change.delta_z))
        qc_rows.append({
            "dataset": change.dataset,
            "patient_id": change.patient_id,
            "signature_id": signature,
            "time_window": change.time_window,
            "gene_n": len(genes),
            "oriented_reconstructed_raw_change": oriented_sum,
            "stage3_raw_change": float(change.raw_change),
            "reconstructed_delta_z": delta_sum,
            "stage3_delta_z": float(change.delta_z),
            "raw_absolute_error": raw_error,
            "delta_z_absolute_error": z_error,
            "tolerance": TOL,
            "status": "PASS" if raw_error <= TOL and z_error <= TOL else "FAIL",
        })
        for gene in genes:
            raw_contribution = float(raw[gene])
            oriented_contribution = raw_contribution * direction
            rows.append({
                "dataset": change.dataset,
                "patient_id": change.patient_id,
                "time_window": change.time_window,
                "signature_id": signature,
                "gene": gene,
                "gene_group": gene_group(signature, gene),
                "baseline_sample_id": change.baseline_sample_id,
                "followup_sample_id": change.followup_sample_id,
                "baseline_expression": x0[gene],
                "followup_expression": x1[gene],
                "expression_change": x1[gene] - x0[gene],
                "original_coefficient": coefficient(signature, gene),
                "direction_coefficient": direction,
                "raw_contribution": raw_contribution,
                "oriented_raw_contribution": oriented_contribution,
                "baseline_score_sd": float(change.baseline_sd),
                "standardized_contribution": oriented_contribution / float(change.baseline_sd),
                "attribution_method": CONFIG["decomposition"][signature],
                "development_overlap_status": overlaps[(signature, change.dataset)],
                "analysis_role": change.role,
            })

    gene_level = pd.DataFrame(rows)
    qc = pd.DataFrame(qc_rows)
    if not (qc["status"] == "PASS").all():
        failed = qc[qc["status"] != "PASS"].head(10)
        raise RuntimeError("decomposition reconstruction failed:\n" + failed.to_string(index=False))
    return gene_level, qc


def bootstrap_summary(values: np.ndarray, *parts: str) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    if len(values) == 1:
        return math.nan, float(values[0]), float(values[0])
    rng = np.random.default_rng(stable_seed(*parts))
    means = values[rng.integers(0, len(values), size=(BOOTSTRAPS, len(values)))].mean(axis=1)
    return float(means.std(ddof=1)), float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def cohort_summaries(gene_level: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    patient_metric_rows = []
    for keys, group in gene_level.groupby(["dataset", "patient_id", "signature_id", "time_window"], sort=False):
        values = group["standardized_contribution"].to_numpy(float)
        abs_values = np.abs(values)
        abs_sum = float(abs_values.sum())
        index = int(np.argmax(abs_values)) if len(values) else 0
        patient_metric_rows.append({
            "dataset": keys[0], "patient_id": keys[1], "signature_id": keys[2], "time_window": keys[3],
            "total_delta_z": float(values.sum()),
            "absolute_contribution_sum": abs_sum,
            "dominance_ratio": float(abs_values[index] / abs_sum) if abs_sum > 0 else 0.0,
            "cancellation_index": float(1.0 - abs(values.sum()) / abs_sum) if abs_sum > 0 else 0.0,
            "dominant_gene": str(group.iloc[index]["gene"]),
        })
    patient_metrics = pd.DataFrame(patient_metric_rows)

    rows = []
    for keys, group in gene_level.groupby(["dataset", "signature_id", "time_window", "gene", "gene_group"], sort=False):
        values = group["standardized_contribution"].to_numpy(float)
        se, lower, upper = bootstrap_summary(values, *map(str, keys))
        rows.append({
            "dataset": keys[0], "signature_id": keys[1], "time_window": keys[2], "gene": keys[3], "gene_group": keys[4],
            "paired_n": len(values), "mean_standardized_contribution": float(values.mean()),
            "median_standardized_contribution": float(np.median(values)), "bootstrap_se": se,
            "ci95_lower": lower, "ci95_upper": upper, "positive_fraction": float(np.mean(values > 0)),
            "negative_fraction": float(np.mean(values < 0)), "mean_absolute_contribution": float(np.mean(np.abs(values))),
            "development_overlap_status": group["development_overlap_status"].iloc[0],
            "analysis_role": group["analysis_role"].iloc[0],
        })
    cohort_gene = pd.DataFrame(rows)
    cohort_gene["analysis_role"] = cohort_gene["analysis_role"].replace(
        {"BLINDED_VALIDATION": "PRESPECIFIED_NON_PILOT"}
    )

    metric_summary = patient_metrics.groupby(["dataset", "signature_id", "time_window"], as_index=False).agg(
        paired_n=("patient_id", "size"),
        mean_total_delta_z=("total_delta_z", "mean"),
        median_absolute_contribution_sum=("absolute_contribution_sum", "median"),
        median_dominance_ratio=("dominance_ratio", "median"),
        median_cancellation_index=("cancellation_index", "median"),
    )
    dominant = patient_metrics.groupby(["dataset", "signature_id", "time_window"])["dominant_gene"].agg(lambda x: Counter(x).most_common(1)[0])
    dominant = dominant.reset_index(name="dominant_tuple")
    dominant[["dominant_gene", "dominant_gene_patient_n"]] = pd.DataFrame(dominant.pop("dominant_tuple").tolist(), index=dominant.index)
    cohort_metrics = metric_summary.merge(dominant, on=["dataset", "signature_id", "time_window"], how="left")
    cohort_metrics["dominant_gene_patient_fraction"] = cohort_metrics["dominant_gene_patient_n"] / cohort_metrics["paired_n"]
    cohort_context = gene_level[
        [
            "dataset",
            "signature_id",
            "time_window",
            "development_overlap_status",
            "analysis_role",
        ]
    ].drop_duplicates()
    if cohort_context.duplicated(["dataset", "signature_id", "time_window"]).any():
        raise RuntimeError("cohort context is not unique at dataset-signature-window grain")
    cohort_metrics = cohort_metrics.merge(
        cohort_context,
        on=["dataset", "signature_id", "time_window"],
        how="left",
        validate="one_to_one",
    )
    cohort_metrics["analysis_role"] = cohort_metrics["analysis_role"].replace(
        {"BLINDED_VALIDATION": "PRESPECIFIED_NON_PILOT"}
    )
    return cohort_gene, cohort_metrics


def analysis_subset(frame: pd.DataFrame, analysis_set: str) -> pd.DataFrame:
    if analysis_set == "ALL_COHORTS":
        return frame
    if analysis_set == "PRESPECIFIED_NON_PILOT_ONLY":
        return frame[frame["analysis_role"] == "PRESPECIFIED_NON_PILOT"]
    excluded_primary = {"DEVELOPMENT_OVERLAP", "POSSIBLE_SAME_MARS_PROGRAM", "POSSIBLE_SAME_PROGRAM"}
    excluded_strict = excluded_primary | {"PRIOR_EXTERNAL_VALIDATION", "POSSIBLE_PRIOR_EXTERNAL_BENCHMARK"}
    excluded = excluded_primary if analysis_set == "PRIMARY_INDEPENDENT" else excluded_strict
    return frame[~frame["development_overlap_status"].isin(excluded)]


def bh(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    out = np.full(len(p), np.nan)
    finite = np.isfinite(p)
    idx = np.where(finite)[0]
    if not len(idx):
        return out
    vals = p[idx]
    order = np.argsort(vals)
    ranked = vals[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    out[idx[order]] = adjusted
    return out


def gene_meta(cohort_gene: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta_rows = []
    loo_rows = []
    for keys, base in cohort_gene.groupby(["signature_id", "time_window", "gene", "gene_group"], sort=False):
        for analysis_set in ("ALL_COHORTS", "PRESPECIFIED_NON_PILOT_ONLY", "PRIMARY_INDEPENDENT", "STRICT_NEVER_USED"):
            subset = analysis_subset(base, analysis_set)
            subset = subset[np.isfinite(subset["bootstrap_se"]) & (subset["bootstrap_se"] > 0)]
            if subset.empty:
                continue
            result = reml_meta(subset["mean_standardized_contribution"].to_numpy(), subset["bootstrap_se"].to_numpy())
            pooled_sign = np.sign(result["mu"])
            consistency = float(np.mean(np.sign(subset["mean_standardized_contribution"]) == pooled_sign)) if pooled_sign else float(np.mean(subset["mean_standardized_contribution"] == 0))
            meta_rows.append({
                "signature_id": keys[0], "time_window": keys[1], "gene": keys[2], "gene_group": keys[3],
                "analysis_set": analysis_set, "cohort_n": len(subset), "patient_n": int(subset["paired_n"].sum()),
                "pooled_contribution": result["mu"], "pooled_se": result["se"], "ci95_lower": result["lower"],
                "ci95_upper": result["upper"], "tau2": result["tau2"], "I2_percent": result["I2"],
                "prediction_lower": result["prediction_lower"], "prediction_upper": result["prediction_upper"],
                "direction_consistency": consistency, "p_value": result["p_value"], "cohorts": ";".join(subset["dataset"]),
            })
            if analysis_set == "PRIMARY_INDEPENDENT" and len(subset) >= 3:
                for omitted in subset["dataset"]:
                    reduced = subset[subset["dataset"] != omitted]
                    rr = reml_meta(reduced["mean_standardized_contribution"].to_numpy(), reduced["bootstrap_se"].to_numpy())
                    loo_rows.append({
                        "signature_id": keys[0], "time_window": keys[1], "gene": keys[2], "omitted_dataset": omitted,
                        "cohort_n": len(reduced), "pooled_contribution": rr["mu"], "ci95_lower": rr["lower"],
                        "ci95_upper": rr["upper"], "tau2": rr["tau2"], "I2_percent": rr["I2"],
                    })
    meta = pd.DataFrame(meta_rows)
    for (signature, window), index in meta[(meta["analysis_set"] == "PRIMARY_INDEPENDENT") & meta["time_window"].isin(PRIMARY_WINDOWS)].groupby(["signature_id", "time_window"]).groups.items():
        meta.loc[index, "fdr_within_signature_window"] = bh(meta.loc[index, "p_value"].to_numpy())
    return meta, pd.DataFrame(loo_rows)


def drift_architecture(cohort_metrics: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    stage3 = pd.read_csv(ROOT / "source_data/03_05_signature_level_meta_analysis.csv")
    stage3 = stage3[stage3["analysis_set"] == "PRIMARY_INDEPENDENT"]
    th = CONFIG["drift_architecture_thresholds"]
    rows = []
    for signature in CONFIG["signatures"]:
        selected_window = None
        selected_stage3 = None
        for window in th["classification_time_priority"]:
            stage3_candidates = stage3[
                (stage3["signature_id"] == signature)
                & (stage3["time_window"] == window)
            ]
            if (
                len(stage3_candidates) == 1
                and int(stage3_candidates.iloc[0]["cohort_n"]) >= int(th["minimum_cohorts"])
            ):
                selected_window = window
                selected_stage3 = stage3_candidates.iloc[0]
                break
        if selected_window is None:
            rows.append({"signature_id": signature, "classification_window": "NONE", "drift_architecture": "EVIDENCE_INSUFFICIENT", "classification_reason": "fewer than three cohorts"})
            continue
        s3 = selected_stage3
        primary_cohorts = tuple(str(s3["cohorts"]).split(";"))
        if len(primary_cohorts) != int(s3["cohort_n"]):
            raise RuntimeError(
                f"{signature} {selected_window}: Stage 3 cohort list length does not match cohort_n"
            )
        cm = cohort_metrics[
            (cohort_metrics["signature_id"] == signature)
            & (cohort_metrics["time_window"] == selected_window)
            & (cohort_metrics["dataset"].isin(primary_cohorts))
        ].copy()
        observed_cohorts = set(cm["dataset"])
        expected_cohorts = set(primary_cohorts)
        if observed_cohorts != expected_cohorts or len(cm) != len(primary_cohorts):
            raise RuntimeError(
                f"{signature} {selected_window}: architecture cohorts {sorted(observed_cohorts)} "
                f"do not match Stage 3 PRIMARY_INDEPENDENT cohorts {sorted(expected_cohorts)}"
            )
        if int(cm["paired_n"].sum()) != int(s3["patient_n"]):
            raise RuntimeError(
                f"{signature} {selected_window}: architecture patient count "
                f"{int(cm['paired_n'].sum())} does not match Stage 3 patient_n {int(s3['patient_n'])}"
            )
        gm = meta[(meta["signature_id"] == signature) & (meta["time_window"] == selected_window) & (meta["analysis_set"] == "PRIMARY_INDEPENDENT")].copy()
        gene_abs = gm.assign(abs_pooled=gm["pooled_contribution"].abs()).sort_values("abs_pooled", ascending=False)
        total_gene_abs = float(gene_abs["abs_pooled"].sum())
        leading_gene = str(gene_abs.iloc[0]["gene"])
        leading_share = float(gene_abs.iloc[0]["abs_pooled"] / total_gene_abs) if total_gene_abs > 0 else 0.0
        gene_shares_ge = int(np.sum(gene_abs["abs_pooled"] / total_gene_abs >= float(th["multigene_share"]))) if total_gene_abs > 0 else 0
        top_by_cohort = cm["dominant_gene"].value_counts(normalize=True)
        modal_dominant_gene = str(top_by_cohort.index[0]) if len(top_by_cohort) else ""
        top_identity_fraction = float(top_by_cohort.iloc[0]) if len(top_by_cohort) else math.nan
        median_dominance = float(cm["median_dominance_ratio"].median())
        median_cancellation = float(cm["median_cancellation_index"].median())
        median_abs_sum = float(cm["median_absolute_contribution_sum"].median())
        pooled_delta = float(s3["pooled_delta_z"])
        direction_consistency = float(s3["direction_consistency"])
        if float(s3["I2_percent"]) >= float(th["cohort_dependence_I2_percent"]) or top_identity_fraction < float(th["dominant_gene_cohort_identity_fraction"]):
            label = "COHORT_DEPENDENT_DRIFT"
            reason = "high total-score heterogeneity or low modal dominant-gene agreement across cohorts"
        elif median_dominance >= float(th["single_gene_dominance_ratio"]) and top_identity_fraction >= float(th["dominant_gene_cohort_identity_fraction"]):
            label = "SINGLE_GENE_DOMINANT_DRIFT"
            reason = "dominance ratio and modal dominant-gene agreement thresholds met"
        elif abs(pooled_delta) < float(th["stable_absolute_pooled_delta_z"]) and median_cancellation >= float(th["cancellation_index"]) and median_abs_sum >= float(th["nontrivial_absolute_contribution_sum"]):
            label = "INTERNAL_CANCELLATION_STABILITY"
            reason = "small total drift despite nontrivial opposing gene contributions"
        elif abs(pooled_delta) >= float(th["stable_absolute_pooled_delta_z"]) and gene_shares_ge >= int(th["multigene_minimum_count"]) and direction_consistency >= float(th["total_direction_consistency"]):
            label = "CONSISTENT_MULTIGENE_DRIFT"
            reason = "multiple genes contribute materially and total direction is consistent"
        elif abs(pooled_delta) < float(th["stable_absolute_pooled_delta_z"]) and median_abs_sum < float(th["nontrivial_absolute_contribution_sum"]) and median_cancellation < float(th["cancellation_index"]):
            label = "OVERALL_LOW_CHANGE_STABILITY"
            reason = "small total and component-level change"
        else:
            label = "EVIDENCE_INSUFFICIENT"
            reason = "no frozen architecture rule was fully met"
        rows.append({
            "signature_id": signature, "classification_window": selected_window, "drift_architecture": label,
            "classification_reason": reason, "stage3_pooled_delta_z": pooled_delta, "stage3_I2_percent": float(s3["I2_percent"]),
            "stage3_prediction_lower": float(s3["prediction_lower"]), "stage3_prediction_upper": float(s3["prediction_upper"]),
            "median_patient_dominance_ratio": median_dominance, "median_patient_cancellation_index": median_cancellation,
            "median_patient_absolute_contribution_sum": median_abs_sum, "leading_gene": leading_gene,
            "leading_gene_pooled_absolute_share": leading_share,
            "modal_cohort_dominant_gene": modal_dominant_gene,
            "modal_dominant_gene_agreement_fraction": top_identity_fraction,
            "architecture_cohort_n": int(len(cm)),
            "architecture_patient_n": int(cm["paired_n"].sum()),
            "architecture_cohorts": ";".join(primary_cohorts),
            "cohort_summary_aggregation": "median_of_cohort_specific_patient_medians",
            "material_gene_count": gene_shares_ge, "stage3_direction_consistency": direction_consistency,
        })
    return pd.DataFrame(rows)


def write_specification() -> None:
    specs = pd.DataFrame(STRUCTURE_ROWS)
    specs["gene_n"] = specs["signature_id"].map(lambda s: len(STAGE3_REQUIRED_GENES[s]))
    specs["direction_coefficient"] = specs["signature_id"].map(SCORE_DIRECTION)
    specs["reconstruction_tolerance"] = TOL
    specs.to_csv(SOURCE / "04_03_signature_decomposition_specification.csv", index=False)
    gene_rows = []
    for signature in CONFIG["signatures"]:
        for gene in STAGE3_REQUIRED_GENES[signature]:
            gene_rows.append({
                "signature_id": signature, "gene": gene, "gene_group": gene_group(signature, gene),
                "original_coefficient": coefficient(signature, gene), "direction_coefficient": SCORE_DIRECTION[signature],
                "attribution_method": CONFIG["decomposition"][signature],
            })
    pd.DataFrame(gene_rows).to_csv(SOURCE / "04_03_signature_gene_coefficients.csv", index=False)


def write_qc_html(qc: pd.DataFrame) -> None:
    summary = qc.groupby(["signature_id", "time_window"], as_index=False).agg(
        records=("patient_id", "size"), passed=("status", lambda x: int(np.sum(x == "PASS"))),
        max_raw_error=("raw_absolute_error", "max"), max_delta_z_error=("delta_z_absolute_error", "max"),
    )
    table = summary.to_html(index=False, float_format=lambda x: f"{x:.3e}", border=0)
    document = f"""<!doctype html><html><head><meta charset='utf-8'><title>Stage 4 decomposition reconstruction QC</title>
<style>body{{font-family:Arial,sans-serif;margin:32px;color:#172b3a}}table{{border-collapse:collapse;width:100%}}th{{background:#1f4e78;color:white;padding:8px}}td{{padding:6px;border-bottom:1px solid #d9e2f3}}.pass{{color:#137333;font-weight:700}}</style></head>
<body><h1>Stage 4 decomposition reconstruction QC</h1><p class='pass'>PASS: {int((qc.status == 'PASS').sum())}/{len(qc)} patient-signature-window records reconstructed within absolute tolerance {TOL:g}.</p>
<p>Both oriented raw-score change and Stage 3 deltaZ were reconstructed. Nonlinear scores used exact Shapley or fixed-partition Owen values; no linear approximation was used.</p>{table}</body></html>"""
    (ROOT / "04_05_decomposition_reconstruction_qc.html").write_text(document, encoding="utf-8")


def main() -> None:
    SOURCE.mkdir(parents=True, exist_ok=True)
    write_specification()
    gene_level, qc = build_gene_level()
    gene_level.to_parquet(ROOT / "04_04_gene_level_longitudinal_dataset.parquet", index=False)
    qc.to_csv(SOURCE / "04_05_decomposition_reconstruction_qc.csv", index=False)
    qc_summary = qc.groupby(["dataset", "signature_id", "time_window"], as_index=False).agg(
        record_n=("patient_id", "size"), pass_n=("status", lambda x: int(np.sum(x == "PASS"))),
        max_raw_absolute_error=("raw_absolute_error", "max"), max_delta_z_absolute_error=("delta_z_absolute_error", "max"),
    )
    qc_summary["status"] = np.where(qc_summary["record_n"] == qc_summary["pass_n"], "PASS", "FAIL")
    qc_summary.to_csv(SOURCE / "04_04_gene_level_data_qc.csv", index=False)
    write_qc_html(qc)

    cohort_gene, cohort_metrics = cohort_summaries(gene_level)
    cohort_gene.to_csv(SOURCE / "04_06_cohort_level_gene_contributions.csv", index=False)
    cohort_metrics.to_csv(SOURCE / "04_06_cohort_level_signature_architecture_metrics.csv", index=False)
    meta, loo = gene_meta(cohort_gene)
    meta.to_csv(SOURCE / "04_07_gene_contribution_meta_analysis.csv", index=False)
    loo.to_csv(SOURCE / "04_07_gene_contribution_leave_one_out.csv", index=False)
    architecture = drift_architecture(cohort_metrics, meta)
    architecture.to_csv(SOURCE / "04_08_signature_drift_architecture.csv", index=False)

    summary = {
        "gene_level_rows": len(gene_level),
        "patient_signature_window_records": len(qc),
        "reconstruction_passed": int((qc.status == "PASS").sum()),
        "maximum_raw_error": float(qc.raw_absolute_error.max()),
        "maximum_delta_z_error": float(qc.delta_z_absolute_error.max()),
        "cohort_gene_rows": len(cohort_gene),
        "gene_meta_rows": len(meta),
        "signatures_classified": int(architecture.signature_id.nunique()),
    }
    (SOURCE / "04_decomposition_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
