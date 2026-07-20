#!/usr/bin/env python3
"""Modified Knapp-Hartung sensitivity for the eight primary comparisons.

The prespecified standard Hartung-Knapp analysis remains primary. This script
reuses the exact cohort-level inputs and REML tau-squared values recorded in
Supplementary Data S18, constraining the Knapp-Hartung variance factor to at
least one. Holm adjustment is then applied across SIG001-SIG004 x T24/T48.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PRIMARY_SIGNATURES = ("SIG001", "SIG002", "SIG003", "SIG004")
PRIMARY_WINDOWS = ("T24", "T48")


def holm(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    running = 0.0
    m = len(p)
    for rank, index in enumerate(order):
        running = max(running, (m - rank) * p[index])
        adjusted[index] = min(1.0, running)
    return adjusted


def compute(group: pd.DataFrame) -> dict[str, float | int | str]:
    group = group.copy()
    y = group["cohort_mean_delta_z"].to_numpy(dtype=float)
    vi = group["sampling_variance"].to_numpy(dtype=float)
    tau2_values = group["tau2"].dropna().to_numpy(dtype=float)
    if len(tau2_values) == 0:
        raise ValueError("Missing tau2 in S18 input")
    tau2 = float(tau2_values[0])
    if not np.allclose(tau2_values, tau2, rtol=0, atol=1e-12):
        raise ValueError("Non-unique tau2 within a primary comparison")

    k = len(group)
    if k < 3:
        raise ValueError("Modified HK sensitivity requires at least three cohorts")
    w = 1.0 / (vi + tau2)
    mu = float(np.sum(w * y) / np.sum(w))
    q_hk = float(np.sum(w * (y - mu) ** 2) / (k - 1))
    standard_hk_se = float(np.sqrt(max(q_hk, 1e-12) / np.sum(w)))
    conventional_re_se = float(np.sqrt(1.0 / np.sum(w)))
    modified_hk_se = float(np.sqrt(max(q_hk, 1.0) / np.sum(w)))
    crit = float(stats.t.ppf(0.975, df=k - 1))
    t_value = mu / modified_hk_se
    p_value = float(2 * stats.t.sf(abs(t_value), df=k - 1))
    lower = mu - crit * modified_hk_se
    upper = mu + crit * modified_hk_se

    return {
        "signature_id": str(group["signature_id"].iloc[0]),
        "time_window": str(group["time_window"].iloc[0]),
        "analysis_set": "Primary independent set",
        "cohort_n": k,
        "patient_n": int(group["paired_n"].sum()),
        "cohorts": ";".join(group["dataset"].astype(str)),
        "pooled_delta_z": mu,
        "tau2": tau2,
        "q_hk": q_hk,
        "conventional_re_se": conventional_re_se,
        "standard_hk_se": standard_hk_se,
        "modified_hk_se": modified_hk_se,
        "modified_hk_ci95_lower": lower,
        "modified_hk_ci95_upper": upper,
        "modified_hk_p_value": p_value,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_workbook", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()

    data = pd.read_excel(args.input_workbook, sheet_name="S18_MetaInputs")
    data = data[
        data["signature_id"].isin(PRIMARY_SIGNATURES)
        & data["time_window"].isin(PRIMARY_WINDOWS)
        & (data["analysis_set"] == "PRIMARY_INDEPENDENT")
    ].copy()
    rows = []
    for signature in PRIMARY_SIGNATURES:
        for window in PRIMARY_WINDOWS:
            subset = data[(data["signature_id"] == signature) & (data["time_window"] == window)]
            if subset.empty:
                raise ValueError(f"Missing primary input: {signature} {window}")
            rows.append(compute(subset))

    result = pd.DataFrame(rows)
    result["modified_hk_holm_adjusted_p"] = holm(result["modified_hk_p_value"].to_numpy())
    result["modified_hk_adjusted_evidence"] = np.where(
        result["modified_hk_holm_adjusted_p"] < 0.05,
        "Retained",
        "Not retained",
    )
    result["interpretation"] = np.where(
        result["modified_hk_holm_adjusted_p"] < 0.05,
        "Multiplicity-adjusted evidence retained under modified Knapp-Hartung sensitivity",
        "Multiplicity-adjusted evidence not retained under modified Knapp-Hartung sensitivity",
    )
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_csv, index=False, float_format="%.15g")
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
