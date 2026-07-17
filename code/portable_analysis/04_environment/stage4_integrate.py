#!/usr/bin/env python3
"""Build Stage 4 claim-boundary, context-of-use and integrated interpretation source tables."""

from __future__ import annotations

import os

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("SEPSIS_SIGNATURE_ANALYSIS_ROOT", Path.cwd())).resolve()
SOURCE = ROOT / "04_source_data"
SIGNATURES = ["SIG001", "SIG002", "SIG003", "SIG004", "SIG022", "SIG023", "SIG033", "SIG034"]


def evidence_levels() -> None:
    rows = [
        ["E1", "Exact contribution decomposition under the frozen formula", "Explain the mathematical source of score change", "Causal driver; within-cell regulation; clinical utility", "CORE_REQUIRED", "PASS"],
        ["E2", "Multicohort association with prespecified pathway change", "State that score change is consistent with or coupled to a biological process", "State that a pathway drives score change or treatment response", "CONDITIONAL", "PASS_LIMITED"],
        ["E3", "External blood/single-cell reference annotation", "Describe a potential cell source of component genes", "Infer cell abundance or within-cell transcription from bulk RNA", "REQUIRED_LIGHTWEIGHT", "PASS"],
        ["E4", "Association with synchronous clinical-severity change", "Limited clinical anchoring if the formal gate is met", "Prognostic validation, monitoring utility or treatment guidance", "CONDITIONAL", "STOPPED_BY_GATE"],
        ["E5", "Causal mechanism or clinical decision utility", "None in the current study", "Any causal or decision-support claim", "PROHIBITED", "NOT_ASSESSED"],
    ]
    pd.DataFrame(rows, columns=["evidence_level", "evidence_type", "permitted_claim", "prohibited_claim", "stage4_role", "stage4_status"]).to_csv(
        SOURCE / "04_02_evidence_level_and_claim_boundary.csv", index=False
    )


def context_of_use() -> None:
    rows = [
        ["Baseline adjunctive diagnosis", "INDIRECT_ONLY", "The current study preserves original diagnostic context but does not re-establish thresholds.", "Assay-specific threshold recovery, calibration and prospective consecutive recruitment", "Do not describe as newly validated baseline diagnostic performance."],
        ["Repeat testing at 24 h", "SUPPORTED_FOR_SCORE_CHANGE_DESCRIPTION", "Within-patient score stability or drift can be quantified at T24.", "Interpretation bands and prospective clinical anchors", "Sampling time must be stated with every score."],
        ["Repeat testing at 48 h", "SUPPORTED_WITH_ATTRITION_CAUTION", "Within-patient drift and cohort heterogeneity can be quantified at T48.", "Prospective handling of informative dropout and interpretation bands", "T48 evidence is more vulnerable to dropout and recovery-state selection."],
        ["Treatment-response monitoring", "LIMITED_BIOLOGICAL_ANCHOR_ONLY", "Some score changes may track prespecified pathway changes.", "Replicated synchronous clinical anchor, assay repeatability and prespecified response criterion", "Do not use pathway coupling as proof of treatment responsiveness."],
        ["Prognostic assessment", "EXPLORATORY_CONTEXT_ONLY", "Fixed outcome labels can describe cohort composition.", "Prospective prediction model, calibration and independent validation", "No mortality association is interpreted as prognostic validation."],
        ["Clinical decision support", "NOT_SUPPORTED", "No treatment or triage recommendation.", "Impact study, explicit action threshold, benefit-harm analysis and external validation", "Current evidence cannot guide treatment."],
    ]
    pd.DataFrame(rows, columns=["use_scenario", "current_support", "permitted_conclusion", "missing_evidence", "claim_boundary"]).to_csv(
        SOURCE / "04_18_context_of_use_evidence_matrix.csv", index=False
    )


def integrated_matrix() -> None:
    names = {r["signature_id"]: r["name"] for r in json.loads((ROOT / "03_00_09_environment_lock/analysis_config.json").read_text())["signatures"]}
    stage3 = pd.read_csv(ROOT / "source_data/03_14_signature_evidence_grade.csv")
    arch = pd.read_csv(SOURCE / "04_08_signature_drift_architecture.csv")
    cell = pd.read_csv(SOURCE / "04_14_leading_cell_source_by_signature.csv")
    coupling = pd.read_csv(SOURCE / "04_12_meta_signature_pathway_coupling.csv")
    coupling = coupling[(coupling["analysis_set"] == "INDEPENDENT_ONLY") & (coupling["tier"] == "PRESET_PRIMARY")]
    coupling_ssgsea = pd.read_csv(SOURCE / "04_12_meta_signature_pathway_coupling_ssgsea.csv")
    coupling_ssgsea = coupling_ssgsea[(coupling_ssgsea["analysis_set"] == "INDEPENDENT_ONLY") & (coupling_ssgsea["tier"] == "PRESET_PRIMARY")]
    interpretation = {
        "SIG001": "COHORT_DEPENDENT_DRIFT",
        "SIG002": "CONSISTENT_TEMPORAL_DRIFT",
        "SIG003": "CONSISTENT_TEMPORAL_DRIFT_SINGLE_GENE_DOMINANT",
        "SIG004": "CONSISTENT_TEMPORAL_DRIFT_SINGLE_GENE_DOMINANT",
        "SIG022": "COHORT_DEPENDENT_DRIFT",
        "SIG023": "RELATIVE_SCORE_STABILITY_WITH_SINGLE_GENE_DOMINANCE",
        "SIG033": "INTERNAL_CANCELLATION_AND_COHORT_DEPENDENT_STABILITY",
        "SIG034": "TIME_STABLE_WITH_INTERPRETABLE_CANCELLATION_AND_LIMITED_E2_ANCHOR",
    }
    rows = []
    for sig in SIGNATURES:
        s1 = stage3[(stage3["signature_id"] == sig) & (stage3["time_window"] == "T1")].iloc[0]
        s2 = stage3[(stage3["signature_id"] == sig) & (stage3["time_window"] == "T2")].iloc[0]
        a = arch[arch["signature_id"] == sig].iloc[0]
        c = cell[(cell["signature_id"] == sig) & (cell["time_window"] == "T2")].iloc[0]
        q = coupling[(coupling["signature_id"] == sig) & (coupling["time_window"] == "T2")].sort_values("fdr_within_analysis_window_tier")
        best = q.iloc[0]
        q2 = coupling_ssgsea[(coupling_ssgsea["signature_id"] == sig) & (coupling_ssgsea["time_window"] == "T2") & (coupling_ssgsea["pathway"] == best["pathway"])]
        sens = q2.iloc[0] if len(q2) else None
        rows.append({
            "signature_id": sig,
            "signature_name": names[sig],
            "T24_pooled_delta_z": s1["pooled_delta_z"],
            "T24_ci": f"{s1['ci95_lower']:.3f} to {s1['ci95_upper']:.3f}",
            "T24_stability_class": s1["primary_stability_class"],
            "T48_pooled_delta_z": s2["pooled_delta_z"],
            "T48_ci": f"{s2['ci95_lower']:.3f} to {s2['ci95_upper']:.3f}",
            "T48_prediction_interval": f"{s2['prediction_lower']:.3f} to {s2['prediction_upper']:.3f}",
            "T48_I2_percent": s2["I2_percent"],
            "T48_evidence_grade": s2["evidence_grade"],
            "drift_architecture": a["drift_architecture"],
            "leading_gene": a["leading_gene"],
            "leading_gene_pooled_absolute_share": a["leading_gene_pooled_absolute_share"],
            "median_dominance_ratio": a["median_patient_dominance_ratio"],
            "median_cancellation_index": a["median_patient_cancellation_index"],
            "leading_potential_cell_source": c["leading_cell_source"],
            "cell_source_absolute_share": c["leading_cell_source_absolute_share"],
            "T48_strongest_primary_pathway": best["pathway"],
            "T48_pathway_rho": best["pooled_spearman_rho"],
            "T48_pathway_ci": f"{best['ci_low']:.3f} to {best['ci_high']:.3f}",
            "T48_pathway_prediction_interval": f"{best['prediction_low']:.3f} to {best['prediction_high']:.3f}",
            "T48_pathway_FDR": best["fdr_within_analysis_window_tier"],
            "T48_pathway_passes_FDR": "YES" if best["fdr_within_analysis_window_tier"] < 0.05 else "NO",
            "ssGSEA_same_pathway_rho": sens["pooled_spearman_rho"] if sens is not None else np.nan,
            "ssGSEA_same_pathway_FDR": sens["fdr_within_analysis_window_tier"] if sens is not None else np.nan,
            "formal_clinical_anchor": "NOT_RUN_BY_PREDEFINED_GATE",
            "development_overlap_robustness": "ASSESSED_IN_STAGE3_AND_GENE_CONTRIBUTION_META",
            "final_interpretation_level": interpretation[sig],
            "claim_boundary": "Mathematical contribution and association only; no causal, monitoring, prognostic or treatment claim.",
        })
    pd.DataFrame(rows).to_csv(SOURCE / "04_20_integrated_signature_interpretation_matrix.csv", index=False)


def multiplicity_robustness() -> None:
    rows = [
        ["Gene contribution decomposition", "All component genes reported; within-signature/window FDR supplied for inferential summaries", "Exact reconstruction; primary-independent meta; leave-one-dataset-out; development-overlap exclusion", "PASS", "No gene selected for entry based on significance"],
        ["Pathway longitudinal changes", "BH-FDR separately by T24/T48 and pathway tier", "Prediction intervals; leave-one-dataset-out; singscore versus ssGSEA", "PASS", "No prespecified primary pathway reached FDR<0.05 for its overall mean change"],
        ["Signature–pathway coupling", "BH-FDR by analysis set, time window and tier", "Independent-only cohorts; prediction intervals; ssGSEA sensitivity", "PASS_LIMITED", "One primary T48 coupling passed FDR: SIG034–oxidative phosphorylation"],
        ["Cell-source annotation", "Descriptive; no hypothesis test", "HPA exact/alias mapping; PBMC-atlas scope check; unresolved genes retained", "PASS", "Potential source only; no deconvolution"],
        ["Clinical anchoring", "No testing because the availability gate failed", "Static severity and fixed outcomes explicitly rejected as substitutes", "STOPPED_BY_GATE", "Zero clinical association P values generated"],
        ["Time-window families", "T24 and T48 kept separate", "No replacement by T72/day5", "PASS", "Extended windows do not alter primary interpretation"],
        ["Scale sensitivity", "Stage 3 baseline-SD primary; MAD and raw/rank sensitivity retained", "No result-driven scale selection", "PASS", "Stage 4 contributions inherit the frozen denominator"],
        ["Influential datasets", "No P-value gate", "Leave-one-dataset-out and prediction intervals", "PASS", "High heterogeneity is retained in interpretation"],
        ["Single-gene dominance", "Pre-coded dominance and cancellation thresholds", "Dominant-gene identity across cohorts; all genes retained", "PASS", "Dominance is descriptive, not a basis for dropping signatures"],
    ]
    pd.DataFrame(rows, columns=["analysis_family", "multiplicity_rule", "robustness_checks", "status", "interpretation_note"]).to_csv(
        SOURCE / "04_21_multiplicity_and_robustness_audit.csv", index=False
    )


def main() -> None:
    evidence_levels()
    context_of_use()
    integrated_matrix()
    multiplicity_robustness()
    print(json.dumps({"evidence_rows": 5, "context_rows": 6, "integrated_signatures": 8, "robustness_rows": 9}, indent=2))


if __name__ == "__main__":
    main()
