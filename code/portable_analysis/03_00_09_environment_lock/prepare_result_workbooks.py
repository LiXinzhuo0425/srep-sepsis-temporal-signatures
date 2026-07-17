#!/usr/bin/env python3
"""Convert Stage 3 CSV result sources to JSON for artifact-tool workbook authoring."""

from __future__ import annotations

import os

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("SEPSIS_SIGNATURE_ANALYSIS_ROOT", Path.cwd())).resolve()
SOURCE = ROOT / "source_data"

FILES = {
    "score_qc": "03_02_score_generation_qc.csv",
    "pairing_flow": "03_03_pairing_flow_report.csv",
    "primary_effects": "03_04_cohort_signature_primary_effects.csv",
    "model_diagnostics": "03_04_model_diagnostics.csv",
    "meta": "03_05_signature_level_meta_analysis.csv",
    "loo": "03_05_leave_one_dataset_out.csv",
    "stability": "03_06_temporal_stability_profile.csv",
    "within_use": "03_07_within_intended_use_comparison.csv",
    "development": "03_08_development_cohort_exclusion_analysis.csv",
    "attrition": "03_09_attrition_and_missingness_analysis.csv",
    "heterogeneity": "03_10_platform_and_population_heterogeneity.csv",
    "scaling": "03_11_scaling_sensitivity_analysis.csv",
    "extended": "03_12_extended_timepoint_analysis.csv",
    "auroc": "03_13_time_dependent_discrimination_analysis.csv",
    "evidence_grade": "03_14_signature_evidence_grade.csv",
    "verification": "03_16_independent_result_verification.csv",
}


def clean(frame: pd.DataFrame):
    frame = frame.replace([np.inf, -np.inf], np.nan)
    rows = []
    for record in frame.to_dict(orient="records"):
        clean_record = {}
        for key, value in record.items():
            if pd.isna(value):
                clean_record[key] = None
            elif isinstance(value, np.generic):
                clean_record[key] = value.item()
            else:
                clean_record[key] = value
        rows.append(clean_record)
    return {"columns": list(frame.columns), "rows": rows}


def main():
    payload = {}
    for name, file_name in FILES.items():
        path = SOURCE / file_name
        if not path.exists():
            raise FileNotFoundError(path)
        payload[name] = clean(pd.read_csv(path))

    cfg = json.loads((ROOT / "03_00_09_environment_lock/analysis_config.json").read_text())
    payload["cohort_characteristics"] = {"columns": list(pd.DataFrame(cfg["cohorts"]).columns), "rows": clean(pd.DataFrame(cfg["cohorts"]))["rows"]}
    payload["signature_characteristics"] = {"columns": list(pd.DataFrame(cfg["signatures"]).columns), "rows": clean(pd.DataFrame(cfg["signatures"]))["rows"]}
    out = ROOT / "work/result_workbook_sources.json"
    out.write_text(json.dumps(payload, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"tables": len(payload), "output": str(out)}, indent=2))


if __name__ == "__main__":
    main()
