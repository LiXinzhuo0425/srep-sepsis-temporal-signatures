#!/usr/bin/env python3
"""Verify that Stage 4 architecture metrics use the Stage 3 primary cohort set.

This is a release-blocking regression check for analysis v1.2.0.  It checks
cohort membership, cohort and patient counts, attainable modal-gene agreement
fractions, and the expected corrected labels for the affected signatures.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("SEPSIS_SIGNATURE_ANALYSIS_ROOT", Path.cwd())).resolve()
S3 = ROOT / "source_data" / "03_05_signature_level_meta_analysis.csv"
ARCH = ROOT / "04_source_data" / "04_08_signature_drift_architecture.csv"
SCENARIOS = ROOT / "04_source_data" / "10_02_scenario_classifications.csv"
STABILITY = ROOT / "04_source_data" / "10_02_classification_stability_matrix.csv"
OUT = ROOT / "04_source_data" / "10_16_architecture_primary_set_verification.json"


def cohort_set(value: object) -> set[str]:
    return {item for item in str(value).split(";") if item}


def main() -> None:
    stage3 = pd.read_csv(S3)
    stage3 = stage3[
        (stage3["analysis_set"] == "PRIMARY_INDEPENDENT")
        & stage3["time_window"].isin(["T1", "T2"])
    ].copy()
    architecture = pd.read_csv(ARCH)
    scenarios = pd.read_csv(SCENARIOS)
    baseline = scenarios[scenarios["scenario"] == "BASELINE"].copy()
    stability = pd.read_csv(STABILITY)

    failures: list[str] = []
    checked: list[dict[str, object]] = []
    for row in architecture.itertuples(index=False):
        meta = stage3[
            (stage3["signature_id"] == row.signature_id)
            & (stage3["time_window"] == row.classification_window)
        ]
        if len(meta) != 1:
            failures.append(f"{row.signature_id}: expected one Stage 3 primary row, found {len(meta)}")
            continue
        meta_row = meta.iloc[0]
        expected_cohorts = cohort_set(meta_row["cohorts"])
        observed_cohorts = cohort_set(row.architecture_cohorts)
        if observed_cohorts != expected_cohorts:
            failures.append(
                f"{row.signature_id}: architecture cohorts {sorted(observed_cohorts)} "
                f"!= Stage 3 primary cohorts {sorted(expected_cohorts)}"
            )
        if int(row.architecture_cohort_n) != int(meta_row["cohort_n"]):
            failures.append(f"{row.signature_id}: cohort_n mismatch")
        if int(row.architecture_patient_n) != int(meta_row["patient_n"]):
            failures.append(f"{row.signature_id}: patient_n mismatch")
        agreement = float(row.modal_dominant_gene_agreement_fraction)
        attainable = agreement * int(row.architecture_cohort_n)
        if not np.isclose(attainable, round(attainable), atol=1e-10):
            failures.append(
                f"{row.signature_id}: modal agreement {agreement} is not attainable "
                f"with k={int(row.architecture_cohort_n)}"
            )

        scenario_rows = scenarios[scenarios["signature_id"] == row.signature_id]
        for srow in scenario_rows.itertuples(index=False):
            if cohort_set(srow.contributing_cohort_ids) != expected_cohorts:
                failures.append(f"{row.signature_id}/{srow.scenario}: scenario cohort set mismatch")
            if int(srow.contributing_cohorts) != int(meta_row["cohort_n"]):
                failures.append(f"{row.signature_id}/{srow.scenario}: scenario cohort_n mismatch")
            if int(srow.contributing_patients) != int(meta_row["patient_n"]):
                failures.append(f"{row.signature_id}/{srow.scenario}: scenario patient_n mismatch")

        checked.append(
            {
                "signature_id": row.signature_id,
                "window": row.classification_window,
                "cohort_n": int(row.architecture_cohort_n),
                "patient_n": int(row.architecture_patient_n),
                "cohorts": sorted(observed_cohorts),
                "modal_agreement": agreement,
                "label": row.drift_architecture,
            }
        )

    affected = architecture.set_index("signature_id")
    expected = {
        "SIG002": (4, 118, "CONSISTENT_MULTIGENE_DRIFT"),
        "SIG003": (4, 118, "SINGLE_GENE_DOMINANT_DRIFT"),
        "SIG004": (4, 118, "CONSISTENT_MULTIGENE_DRIFT"),
    }
    for signature_id, (cohort_n, patient_n, label) in expected.items():
        row = affected.loc[signature_id]
        if (int(row.architecture_cohort_n), int(row.architecture_patient_n), row.drift_architecture) != (
            cohort_n,
            patient_n,
            label,
        ):
            failures.append(f"{signature_id}: corrected count or label does not match the v1.2.0 contract")
        if "GSE54514" in cohort_set(row.architecture_cohorts):
            failures.append(f"{signature_id}: excluded GSE54514 remained in architecture set")

    baseline_labels = baseline.set_index("signature_id")["drift_architecture"]
    for signature_id, row in affected.iterrows():
        if baseline_labels.get(signature_id) != row.drift_architecture:
            failures.append(f"{signature_id}: baseline sensitivity label differs from architecture table")

    status = stability.set_index("signature_id")["sensitivity_status"].to_dict()
    expected_status = {
        "SIG002": "UNSTABLE",
        "SIG003": "BOUNDARY_SENSITIVE",
        "SIG004": "UNSTABLE",
        "SIG022": "UNSTABLE",
    }
    for signature_id, expected_value in expected_status.items():
        if status.get(signature_id) != expected_value:
            failures.append(
                f"{signature_id}: sensitivity status {status.get(signature_id)!r} != {expected_value!r}"
            )

    report = {
        "analysis_release": "v1.2.0",
        "check": "Stage 4 architecture metrics use the exact Stage 3 PRIMARY_INDEPENDENT cohort set",
        "status": "PASS" if not failures else "FAIL",
        "checked_signatures": checked,
        "failures": failures,
    }
    OUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
