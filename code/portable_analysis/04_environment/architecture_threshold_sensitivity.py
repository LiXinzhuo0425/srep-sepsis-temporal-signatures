from __future__ import annotations

import json
import math
from copy import deepcopy
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("SEPSIS_SIGNATURE_ANALYSIS_ROOT", Path.cwd())).resolve()
SOURCE = ROOT / "04_source_data"
ENV4 = ROOT / "04_environment"
ANALYSIS = SOURCE
ANALYSIS.mkdir(parents=True, exist_ok=True)

CONFIG = json.loads((ENV4 / "stage4_config.json").read_text(encoding="utf-8"))
BASE_THRESHOLDS = CONFIG["drift_architecture_thresholds"]
SIGNATURES = CONFIG["signatures"]

COHORT_METRICS = pd.read_csv(SOURCE / "04_06_cohort_level_signature_architecture_metrics.csv")
GENE_META = pd.read_csv(SOURCE / "04_07_gene_contribution_meta_analysis.csv")
STAGE3 = pd.read_csv(ROOT / "source_data/03_05_signature_level_meta_analysis.csv")
STORED_BASELINE = pd.read_csv(SOURCE / "04_08_signature_drift_architecture.csv")

GENE_META = GENE_META[GENE_META["analysis_set"] == "PRIMARY_INDEPENDENT"].copy()
STAGE3 = STAGE3[STAGE3["analysis_set"] == "PRIMARY_INDEPENDENT"].copy()

SCENARIOS = {
    "BASELINE": 1.00,
    "A_MINUS_10_PERCENT": 0.90,
    "B_PLUS_10_PERCENT": 1.10,
    "C_MINUS_20_PERCENT": 0.80,
    "D_PLUS_20_PERCENT": 1.20,
}

NUMERIC_THRESHOLD_KEYS = [
    "minimum_cohorts",
    "cohort_dependence_I2_percent",
    "dominant_gene_cohort_identity_fraction",
    "single_gene_dominance_ratio",
    "stable_absolute_pooled_delta_z",
    "cancellation_index",
    "nontrivial_absolute_contribution_sum",
    "multigene_share",
    "multigene_minimum_count",
    "total_direction_consistency",
]


def thresholds_for(multiplier: float) -> dict:
    """Uniformly scale every numeric classification cutoff.

    Count gates remain real-valued after scaling. Because observed counts are
    integers, this is equivalent to the next attainable integer when a gate is
    increased (e.g. 3 x 1.10 = 3.3 requires at least 4 cohorts). This avoids
    silently cancelling the requested perturbation by integer rounding.
    """

    thresholds = deepcopy(BASE_THRESHOLDS)
    for key in NUMERIC_THRESHOLD_KEYS:
        thresholds[key] = float(BASE_THRESHOLDS[key]) * multiplier
    return thresholds


def metric_row(signature: str, thresholds: dict) -> dict:
    selected_window = None
    selected_stage3 = None
    contributing_cohorts = 0
    for window in thresholds["classification_time_priority"]:
        candidates = STAGE3[
            (STAGE3["signature_id"] == signature)
            & (STAGE3["time_window"] == window)
        ]
        if len(candidates) == 1 and float(candidates.iloc[0]["cohort_n"]) >= float(
            thresholds["minimum_cohorts"]
        ):
            selected_window = window
            selected_stage3 = candidates.iloc[0]
            contributing_cohorts = int(selected_stage3["cohort_n"])
            break

    if selected_window is None:
        return {
            "signature_id": signature,
            "classification_window": "NONE",
            "contributing_cohorts": 0,
            "drift_architecture": "EVIDENCE_INSUFFICIENT",
            "classification_reason": "fewer cohorts than the scenario-specific minimum",
            "rule_triggered": "minimum_cohorts",
        }

    s3 = selected_stage3
    primary_cohorts = tuple(str(s3["cohorts"]).split(";"))
    if len(primary_cohorts) != contributing_cohorts:
        raise RuntimeError(
            f"{signature} {selected_window}: Stage 3 cohort list length does not match cohort_n"
        )
    cm = COHORT_METRICS[
        (COHORT_METRICS["signature_id"] == signature)
        & (COHORT_METRICS["time_window"] == selected_window)
        & (COHORT_METRICS["dataset"].isin(primary_cohorts))
    ].copy()
    if set(cm["dataset"]) != set(primary_cohorts) or len(cm) != contributing_cohorts:
        raise RuntimeError(
            f"{signature} {selected_window}: architecture cohort set does not match "
            "the Stage 3 PRIMARY_INDEPENDENT cohort set"
        )
    if int(cm["paired_n"].sum()) != int(s3["patient_n"]):
        raise RuntimeError(
            f"{signature} {selected_window}: architecture patient count does not match Stage 3"
        )
    gm = GENE_META[
        (GENE_META["signature_id"] == signature)
        & (GENE_META["time_window"] == selected_window)
    ].copy()
    gene_abs = gm.assign(abs_pooled=gm["pooled_contribution"].abs()).sort_values(
        "abs_pooled", ascending=False
    )
    total_gene_abs = float(gene_abs["abs_pooled"].sum())
    leading_gene = str(gene_abs.iloc[0]["gene"])
    leading_share = (
        float(gene_abs.iloc[0]["abs_pooled"] / total_gene_abs)
        if total_gene_abs > 0
        else 0.0
    )
    material_gene_count = (
        int(
            np.sum(
                gene_abs["abs_pooled"] / total_gene_abs
                >= float(thresholds["multigene_share"])
            )
        )
        if total_gene_abs > 0
        else 0
    )
    top_by_cohort = cm["dominant_gene"].value_counts(normalize=True)
    modal_dominant_gene = str(top_by_cohort.index[0]) if len(top_by_cohort) else ""
    top_identity_fraction = (
        float(top_by_cohort.iloc[0]) if len(top_by_cohort) else math.nan
    )
    median_dominance = float(cm["median_dominance_ratio"].median())
    median_cancellation = float(cm["median_cancellation_index"].median())
    median_abs_sum = float(cm["median_absolute_contribution_sum"].median())
    pooled_delta = float(s3["pooled_delta_z"])
    direction_consistency = float(s3["direction_consistency"])
    i2_percent = float(s3["I2_percent"])

    cohort_dependent = (
        i2_percent >= float(thresholds["cohort_dependence_I2_percent"])
        or top_identity_fraction
        < float(thresholds["dominant_gene_cohort_identity_fraction"])
    )
    single_gene = (
        median_dominance >= float(thresholds["single_gene_dominance_ratio"])
        and top_identity_fraction
        >= float(thresholds["dominant_gene_cohort_identity_fraction"])
    )
    internal_cancellation = (
        abs(pooled_delta) < float(thresholds["stable_absolute_pooled_delta_z"])
        and median_cancellation >= float(thresholds["cancellation_index"])
        and median_abs_sum
        >= float(thresholds["nontrivial_absolute_contribution_sum"])
    )
    consistent_multigene = (
        abs(pooled_delta) >= float(thresholds["stable_absolute_pooled_delta_z"])
        and material_gene_count >= float(thresholds["multigene_minimum_count"])
        and direction_consistency
        >= float(thresholds["total_direction_consistency"])
    )
    low_change = (
        abs(pooled_delta) < float(thresholds["stable_absolute_pooled_delta_z"])
        and median_abs_sum
        < float(thresholds["nontrivial_absolute_contribution_sum"])
        and median_cancellation < float(thresholds["cancellation_index"])
    )

    if cohort_dependent:
        label = "COHORT_DEPENDENT_DRIFT"
        reason = "high total-score heterogeneity or low modal dominant-gene agreement across cohorts"
        triggered = "1_cohort_dependent"
    elif single_gene:
        label = "SINGLE_GENE_DOMINANT_DRIFT"
        reason = "dominance ratio and modal dominant-gene agreement thresholds met"
        triggered = "2_single_gene"
    elif internal_cancellation:
        label = "INTERNAL_CANCELLATION_STABILITY"
        reason = "small total drift despite nontrivial opposing gene contributions"
        triggered = "3_internal_cancellation"
    elif consistent_multigene:
        label = "CONSISTENT_MULTIGENE_DRIFT"
        reason = "multiple genes contribute materially and total direction is consistent"
        triggered = "4_consistent_multigene"
    elif low_change:
        label = "OVERALL_LOW_CHANGE_STABILITY"
        reason = "small total and component-level change"
        triggered = "5_low_change"
    else:
        label = "EVIDENCE_INSUFFICIENT"
        reason = "no scenario-specific architecture rule was fully met"
        triggered = "6_fallback_insufficient"

    return {
        "signature_id": signature,
        "classification_window": selected_window,
        "contributing_cohorts": contributing_cohorts,
        "contributing_patients": int(cm["paired_n"].sum()),
        "contributing_cohort_ids": ";".join(primary_cohorts),
        "architecture_cohort_n": contributing_cohorts,
        "architecture_patient_n": int(cm["paired_n"].sum()),
        "drift_architecture": label,
        "classification_reason": reason,
        "rule_triggered": triggered,
        "stage3_pooled_delta_z": pooled_delta,
        "stage3_I2_percent": i2_percent,
        "stage3_prediction_lower": float(s3["prediction_lower"]),
        "stage3_prediction_upper": float(s3["prediction_upper"]),
        "median_patient_dominance_ratio": median_dominance,
        "median_patient_cancellation_index": median_cancellation,
        "median_patient_absolute_contribution_sum": median_abs_sum,
        "leading_gene": leading_gene,
        "leading_gene_pooled_absolute_share": leading_share,
        "modal_cohort_dominant_gene": modal_dominant_gene,
        "modal_dominant_gene_agreement_fraction": top_identity_fraction,
        "cohort_summary_aggregation": "median_of_cohort_specific_patient_medians",
        "material_gene_count": material_gene_count,
        "stage3_direction_consistency": direction_consistency,
        "flag_cohort_dependent": cohort_dependent,
        "flag_single_gene": single_gene,
        "flag_internal_cancellation": internal_cancellation,
        "flag_consistent_multigene": consistent_multigene,
        "flag_low_change": low_change,
    }


scenario_rows = []
threshold_rows = []
for scenario, multiplier in SCENARIOS.items():
    th = thresholds_for(multiplier)
    threshold_rows.append(
        {
            "scenario": scenario,
            "multiplier": multiplier,
            **{key: th[key] for key in NUMERIC_THRESHOLD_KEYS},
        }
    )
    for signature in SIGNATURES:
        scenario_rows.append(
            {
                "scenario": scenario,
                "multiplier": multiplier,
                **metric_row(signature, th),
            }
        )

scenario_results = pd.DataFrame(scenario_rows)
threshold_table = pd.DataFrame(threshold_rows)

baseline = scenario_results[scenario_results["scenario"] == "BASELINE"].copy()
baseline_compare = STORED_BASELINE.merge(
    baseline,
    on="signature_id",
    suffixes=("_stored", "_reconstructed"),
    how="outer",
)
baseline_compare["class_match"] = (
    baseline_compare["drift_architecture_stored"]
    == baseline_compare["drift_architecture_reconstructed"]
)
baseline_compare["window_match"] = (
    baseline_compare["classification_window_stored"]
    == baseline_compare["classification_window_reconstructed"]
)
numeric_fields = [
    "stage3_pooled_delta_z",
    "stage3_I2_percent",
    "stage3_prediction_lower",
    "stage3_prediction_upper",
    "median_patient_dominance_ratio",
    "median_patient_cancellation_index",
    "median_patient_absolute_contribution_sum",
    "leading_gene_pooled_absolute_share",
    "modal_dominant_gene_agreement_fraction",
    "architecture_cohort_n",
    "architecture_patient_n",
    "material_gene_count",
    "stage3_direction_consistency",
]
for field in numeric_fields:
    baseline_compare[f"{field}_abs_diff"] = (
        baseline_compare[f"{field}_stored"]
        - baseline_compare[f"{field}_reconstructed"]
    ).abs()
baseline_compare["numeric_match"] = baseline_compare[
    [f"{field}_abs_diff" for field in numeric_fields]
].max(axis=1) <= 1e-12
baseline_compare["overall_match"] = (
    baseline_compare["class_match"]
    & baseline_compare["window_match"]
    & baseline_compare["numeric_match"]
)

class_matrix = scenario_results.pivot(
    index="signature_id", columns="scenario", values="drift_architecture"
).reindex(SIGNATURES)
for col in SCENARIOS:
    if col not in class_matrix.columns:
        class_matrix[col] = np.nan
class_matrix = class_matrix[list(SCENARIOS)]
class_matrix["n_changed_A_to_D"] = class_matrix[
    [
        "A_MINUS_10_PERCENT",
        "B_PLUS_10_PERCENT",
        "C_MINUS_20_PERCENT",
        "D_PLUS_20_PERCENT",
    ]
].ne(class_matrix["BASELINE"], axis=0).sum(axis=1)
class_matrix["sensitivity_status"] = np.select(
    [
        class_matrix["n_changed_A_to_D"] == 0,
        class_matrix["n_changed_A_to_D"] == 1,
    ],
    ["STABLE", "BOUNDARY_SENSITIVE"],
    default="UNSTABLE",
)
class_matrix = class_matrix.reset_index()

continuous_metrics = baseline[
    [
        "signature_id",
        "classification_window",
        "contributing_cohorts",
        "contributing_patients",
        "contributing_cohort_ids",
        "stage3_pooled_delta_z",
        "stage3_I2_percent",
        "stage3_prediction_lower",
        "stage3_prediction_upper",
        "median_patient_dominance_ratio",
        "median_patient_cancellation_index",
        "median_patient_absolute_contribution_sum",
        "leading_gene",
        "leading_gene_pooled_absolute_share",
        "modal_cohort_dominant_gene",
        "modal_dominant_gene_agreement_fraction",
        "cohort_summary_aggregation",
        "material_gene_count",
        "stage3_direction_consistency",
    ]
].copy()
continuous_metrics.insert(0, "scenario", "E_CONTINUOUS_ONLY")
continuous_metrics["discrete_classification"] = "NOT_APPLIED"

rule_registry = pd.DataFrame(
    [
        {
            "priority": 0,
            "label": "EVIDENCE_INSUFFICIENT",
            "rule": "No T48 window with at least 3 contributing independent study families; then try T24. If neither qualifies, stop.",
            "baseline_thresholds": "minimum_cohorts=3; classification_time_priority=T48 then T24",
            "role": "Window eligibility gate",
        },
        {
            "priority": 1,
            "label": "COHORT_DEPENDENT_DRIFT",
            "rule": "I2 >= 75% OR modal dominant-gene agreement fraction < 0.60.",
            "baseline_thresholds": "I2=75%; modal agreement fraction=0.60",
            "role": "First classification rule; overrides all later labels",
        },
        {
            "priority": 2,
            "label": "SINGLE_GENE_DOMINANT_DRIFT",
            "rule": "Median of cohort-specific patient-level dominance medians >= 0.60 AND modal dominant-gene agreement fraction >= 0.60.",
            "baseline_thresholds": "dominance=0.60; modal agreement fraction=0.60",
            "role": "Applied only after cohort-dependence rule fails",
        },
        {
            "priority": 3,
            "label": "INTERNAL_CANCELLATION_STABILITY",
            "rule": "abs(pooled delta Z) < 0.20 AND median cancellation >= 0.50 AND median absolute contribution sum >= 0.40.",
            "baseline_thresholds": "abs delta Z=0.20; cancellation=0.50; absolute sum=0.40",
            "role": "Applied only after priorities 1-2 fail",
        },
        {
            "priority": 4,
            "label": "CONSISTENT_MULTIGENE_DRIFT",
            "rule": "abs(pooled delta Z) >= 0.20 AND at least 2 genes each have >=15% of pooled absolute contribution AND direction consistency >=0.75.",
            "baseline_thresholds": "abs delta Z=0.20; gene share=0.15; gene count=2; direction=0.75",
            "role": "Applied only after priorities 1-3 fail",
        },
        {
            "priority": 5,
            "label": "OVERALL_LOW_CHANGE_STABILITY",
            "rule": "abs(pooled delta Z) < 0.20 AND median absolute contribution sum < 0.40 AND median cancellation < 0.50.",
            "baseline_thresholds": "abs delta Z=0.20; absolute sum=0.40; cancellation=0.50",
            "role": "Applied only after priorities 1-4 fail",
        },
        {
            "priority": 6,
            "label": "EVIDENCE_INSUFFICIENT",
            "rule": "No preceding architecture rule is fully met.",
            "baseline_thresholds": "No additional cutoff",
            "role": "Fallback after all eligible-window rules",
        },
    ]
)
rule_registry["classification_scope"] = (
    "T48 primary; the exact Stage 3 PRIMARY_INDEPENDENT cohort set for cohort summaries, "
    "pooled total effects and pooled gene effects"
)
rule_registry["interpretive_status"] = (
    "Descriptive label; continuous metrics are primary"
)

provenance = pd.DataFrame(
    [
        {
            "artifact": "04_01_stage4_scope_freeze_v1.0.md",
            "version_or_commit": "Stage 4 scope freeze v1.0",
            "frozen_or_committed_at": "2026-07-16 11:40:18 +0800",
            "sha256_or_git_blob": "7a16a3cad78747b5ddefa2cdfe8d241b5419deaf2eb4d011a49236262ecefcd9",
            "evidentiary_role": "Human-readable pre-result classification rules and priority",
        },
        {
            "artifact": "04_environment/stage4_config.json",
            "version_or_commit": "stage4_extension_v1.0",
            "frozen_or_committed_at": "2026-07-16 11:40:18 +0800",
            "sha256_or_git_blob": "5dd1534de29dabb2b548f0d8bf41e490a7894eb183e98387e3f36d084244d5f9",
            "evidentiary_role": "Machine-readable thresholds used by code",
        },
        {
            "artifact": "Public code/data release v1.0.1",
            "version_or_commit": "67b43348fa3fb2bec6ac7ec4313eae0fad505d66",
            "frozen_or_committed_at": "2026-07-18 01:29:12 +0800",
            "sha256_or_git_blob": "stage4_config git blob 9602e94c92173c3f807f6c9143f18afcf90ab055",
            "evidentiary_role": "Public release containing the byte-identical threshold configuration",
        },
        {
            "artifact": "03_00_07_statistical_analysis_plan_frozen_v1.0.docx",
            "version_or_commit": "Stage 3 SAP v1.0",
            "frozen_or_committed_at": "2026-07-16 09:10:44 +0800",
            "sha256_or_git_blob": "8ce694da9e383ffe5234d648b60bfaec676f440620616aca59bfefd01ccbef8e",
            "evidentiary_role": "Parent SAP for total-score temporal classification; gene-contribution architecture rules were added in the separately frozen Stage 4 scope",
        },
    ]
)

summary = {
    "analysis": "architecture classification threshold sensitivity",
    "baseline_reconstruction_pass": bool(baseline_compare["overall_match"].all()),
    "baseline_reconstruction_pass_n": int(baseline_compare["overall_match"].sum()),
    "baseline_reconstruction_total_n": int(len(baseline_compare)),
    "stable_signature_n": int((class_matrix["sensitivity_status"] == "STABLE").sum()),
    "boundary_sensitive_signature_n": int(
        (class_matrix["sensitivity_status"] == "BOUNDARY_SENSITIVE").sum()
    ),
    "unstable_signature_n": int(
        (class_matrix["sensitivity_status"] == "UNSTABLE").sum()
    ),
    "signatures_with_any_change": class_matrix.loc[
        class_matrix["n_changed_A_to_D"] > 0, "signature_id"
    ].tolist(),
    "signatures_with_two_or_more_changes": class_matrix.loc[
        class_matrix["n_changed_A_to_D"] >= 2, "signature_id"
    ].tolist(),
    "title_weakening_trigger": bool(
        (class_matrix["n_changed_A_to_D"] >= 2).sum() >= 2
    ),
    "threshold_perturbation_note": (
        "Every numeric cutoff was multiplied uniformly. Count gates were compared "
        "as real-valued cutoffs, so increases move to the next attainable integer."
    ),
    "classification_priority": [
        "cohort-dependent",
        "single-gene-dominant",
        "internal cancellation",
        "consistent multigene",
        "overall low change",
        "insufficient evidence",
    ],
}

threshold_table.to_csv(ANALYSIS / "10_02_threshold_scenarios.csv", index=False)
scenario_results.to_csv(ANALYSIS / "10_02_scenario_classifications.csv", index=False)
baseline_compare.to_csv(ANALYSIS / "10_02_baseline_reconstruction_check.csv", index=False)
class_matrix.to_csv(ANALYSIS / "10_02_classification_stability_matrix.csv", index=False)
continuous_metrics.to_csv(ANALYSIS / "10_02_continuous_metrics_scenario_E.csv", index=False)
rule_registry.to_csv(ANALYSIS / "10_02_rule_registry.csv", index=False)
provenance.to_csv(ANALYSIS / "10_02_rule_provenance.csv", index=False)
(ANALYSIS / "10_02_summary.json").write_text(
    json.dumps(summary, indent=2), encoding="utf-8"
)

print(json.dumps(summary, indent=2))
print("\nClassification stability matrix")
print(class_matrix.to_string(index=False))
