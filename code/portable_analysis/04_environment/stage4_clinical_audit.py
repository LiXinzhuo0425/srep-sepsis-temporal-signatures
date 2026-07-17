#!/usr/bin/env python3
"""Audit serial clinical anchors without fitting any association model."""

from __future__ import annotations

import os

import json
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ.get("SEPSIS_SIGNATURE_ANALYSIS_ROOT", Path.cwd())).resolve()
SOURCE = ROOT / "04_source_data"


def main() -> None:
    rows = [
        ["GSE236713", "Delta SOFA", "PRIMARY_LEVEL_1", "Mentioned in series summary / study methods", "NOT_DEPOSITED_SAMPLE_LEVEL", 0, "No serial sample-level SOFA values in GEO matrix, characteristics, RAW archive or supplementary listing", "ABANDON_FORMAL"],
        ["GSE57065", "Delta SOFA", "PRIMARY_LEVEL_1", "Not present", "UNAVAILABLE", 0, "No SOFA field", "ABANDON_FORMAL"],
        ["GSE95233", "Delta SOFA", "PRIMARY_LEVEL_1", "Not present", "UNAVAILABLE", 0, "No SOFA field", "ABANDON_FORMAL"],
        ["GSE54514", "Delta SOFA", "PRIMARY_LEVEL_1", "Not present", "UNAVAILABLE", 0, "No SOFA field", "ABANDON_FORMAL"],
        ["GSE110487", "Delta SOFA", "PRIMARY_LEVEL_1", "Not present in deposited workbook", "UNAVAILABLE", 0, "Expression-only workbook", "ABANDON_FORMAL"],
        ["GSE8121", "Delta SOFA", "PRIMARY_LEVEL_1", "Not present", "UNAVAILABLE", 0, "No SOFA field", "ABANDON_FORMAL"],
        ["GSE236713", "Shock or vasopressor state change", "PRIMARY_LEVEL_1", "Not present", "UNAVAILABLE", 0, "No serial vasoactive or shock-state field", "ABANDON_FORMAL"],
        ["GSE57065", "Shock or vasopressor state change", "PRIMARY_LEVEL_1", "Cohort is septic shock", "BASELINE_CASE_DEFINITION_ONLY", 0, "No serial resolution/persistence field", "ABANDON_FORMAL"],
        ["GSE95233", "Shock or vasopressor state change", "PRIMARY_LEVEL_1", "Cohort is septic shock", "BASELINE_CASE_DEFINITION_ONLY", 0, "No serial resolution/persistence field", "ABANDON_FORMAL"],
        ["GSE54514", "Shock or vasopressor state change", "PRIMARY_LEVEL_1", "Not present", "UNAVAILABLE", 0, "No serial field", "ABANDON_FORMAL"],
        ["GSE110487", "Shock or vasopressor state change", "PRIMARY_LEVEL_1", "Cohort is septic shock", "BASELINE_CASE_DEFINITION_ONLY", 0, "No serial field", "ABANDON_FORMAL"],
        ["GSE8121", "Shock or vasopressor state change", "PRIMARY_LEVEL_1", "Cohort is septic shock", "BASELINE_CASE_DEFINITION_ONLY", 0, "No serial field", "ABANDON_FORMAL"],
        ["ALL", "Delta lactate", "SECONDARY_LEVEL_2", "Not deposited in a time-matched sample-level form", "UNAVAILABLE", 0, "No eligible cohort", "ABANDON_FORMAL"],
        ["ALL", "Delta CRP or PCT", "SECONDARY_LEVEL_2", "Not deposited in a time-matched sample-level form", "UNAVAILABLE", 0, "No eligible cohort", "ABANDON_FORMAL"],
        ["ALL", "Mechanical ventilation state change", "SECONDARY_LEVEL_2", "Not deposited in a time-matched sample-level form", "UNAVAILABLE", 0, "No eligible cohort", "ABANDON_FORMAL"],
        ["GSE54514", "APACHE II change", "NONQUALIFYING_SEVERITY_FIELD", "Sample characteristic", "STATIC_BASELINE_VALUE_REPEATED", 0, "Within-patient APACHE II is unchanged across days; not a longitudinal anchor", "DESCRIPTIVE_ONLY"],
        ["GSE57065", "SAPS II change", "NONQUALIFYING_SEVERITY_FIELD", "Sample characteristic", "STATIC_HIGH_LOW_CATEGORY_REPEATED", 0, "Within-patient category is unchanged at 0, 24 and 48 h", "DESCRIPTIVE_ONLY"],
        ["GSE54514", "Neutrophil proportion change", "NONPRESET_LAB_FIELD", "Sample characteristic", "SERIAL_AVAILABLE_ONE_COHORT", 31, "One cohort only; not a prespecified severity anchor; no new cell-proportion correction analysis", "DESCRIPTIVE_SUPPLEMENT_ONLY"],
        ["GSE236713", "Mortality", "EXPLORATORY_OUTCOME", "Sample characteristic", "AVAILABLE_FIXED_OUTCOME", 93, "Not a synchronous clinical change; cannot validate prognosis", "EXPLORATORY_DESCRIPTION_ONLY"],
        ["GSE54514", "Mortality", "EXPLORATORY_OUTCOME", "Disease status", "AVAILABLE_FIXED_OUTCOME", 31, "Not a synchronous clinical change; cannot validate prognosis", "EXPLORATORY_DESCRIPTION_ONLY"],
        ["GSE95233", "Mortality", "EXPLORATORY_OUTCOME", "Sample characteristic", "AVAILABLE_FIXED_OUTCOME", 31, "Not a synchronous clinical change; cannot validate prognosis", "EXPLORATORY_DESCRIPTION_ONLY"],
    ]
    columns = ["dataset", "clinical_variable", "anchor_level", "public_source", "availability_status", "estimated_time_matched_paired_n", "audit_reason", "decision"]
    audit = pd.DataFrame(rows, columns=columns)
    audit.to_csv(SOURCE / "04_15_clinical_anchor_availability_audit.csv", index=False)

    decision = {
        "formal_gate": "FAIL",
        "formal_module_decision": "STOP_FORMAL_CLINICAL_ANCHORING",
        "eligible_level_1_anchor": None,
        "independent_cohorts_with_same_level_1_anchor": 0,
        "cumulative_time_matched_level_1_patients": 0,
        "required_cohorts": 2,
        "recommended_minimum_patients": 80,
        "reason": "No deposited serial SOFA, vasopressor/shock-state, lactate, CRP or PCT variable is available in at least two independent cohorts. Static APACHE II/SAPS II fields and mortality outcomes do not satisfy the frozen synchronous-change definition.",
        "analysis_04_17": "NOT_RUN_BY_PREDEFINED_GATE",
        "allowed_residual_use": "Availability report and claim-boundary matrix only; fixed mortality may be described as exploratory outcome availability but not analysed as prognostic validation.",
    }
    (SOURCE / "04_15_clinical_anchor_gate_decision.json").write_text(json.dumps(decision, indent=2) + "\n", encoding="utf-8")

    stopped = pd.DataFrame([{
        "status": "NOT_RUN_BY_PREDEFINED_GATE",
        "reason": decision["reason"],
        "main_anchor": "NONE",
        "secondary_anchors": "NONE",
        "models_fitted": 0,
        "p_values_generated": 0,
        "claim_boundary": "No clinical-change association or prognostic claim is made.",
    }])
    stopped.to_csv(SOURCE / "04_17_signature_clinical_anchor_analysis_status.csv", index=False)
    print(json.dumps(decision, indent=2))


if __name__ == "__main__":
    main()
