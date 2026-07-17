#!/usr/bin/env python3
"""Independent deterministic checks for Stage 3 pairing, scoring and synthesis."""

from __future__ import annotations

import os

import json
import math
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize


ROOT = Path(os.environ.get("SEPSIS_SIGNATURE_ANALYSIS_ROOT", Path.cwd())).resolve()
ENV = ROOT / "03_00_09_environment_lock"
sys.path.insert(0, str(ENV))
from run_stage3 import COHORTS, STAGE3_FUNCTIONS, manual_score  # noqa: E402


def independent_reml(y, se):
    y = np.asarray(y, float)
    vi = np.maximum(np.asarray(se, float) ** 2, 1e-12)

    def nll(log_tau_plus):
        tau = max(0.0, math.exp(log_tau_plus) - 1e-12)
        weights = 1.0 / (vi + tau)
        mean = np.dot(weights, y) / weights.sum()
        return 0.5 * (np.log(vi + tau).sum() + math.log(weights.sum()) + np.dot(weights, (y - mean) ** 2))

    solution = optimize.minimize_scalar(nll, bounds=(math.log(1e-12), math.log(20.0)), method="bounded")
    tau = max(0.0, math.exp(solution.x) - 1e-12)
    if nll(math.log(1e-12)) <= nll(solution.x) + 1e-9:
        tau = 0.0
    weights = 1.0 / (vi + tau)
    mean = float(np.dot(weights, y) / weights.sum())
    return mean, tau


def main():
    rng = random.Random(20260716)
    checks = []
    cfg = json.loads((ENV / "analysis_config.json").read_text())
    population = pd.DataFrame(cfg["population"])
    changes = pd.read_parquet(ROOT / "03_03_longitudinal_analysis_datasets/all_windows.parquet")
    scores = pd.read_parquet(ROOT / "03_02_full_signature_score_matrix.parquet")
    effects = pd.read_csv(ROOT / "source_data/03_04_cohort_signature_primary_effects.csv")
    meta = pd.read_csv(ROOT / "source_data/03_05_signature_level_meta_analysis.csv")

    primary_pairs = changes[changes.time_window.isin(["T1", "T2"])][["dataset", "patient_id"]].drop_duplicates()
    sample_pairs = primary_pairs.sample(n=min(20, len(primary_pairs)), random_state=20260716)
    for _, item in sample_pairs.iterrows():
        frozen = population[(population.dataset == item.dataset) & (population.patient_id.astype(str) == str(item.patient_id))].iloc[0]
        observed = changes[(changes.dataset == item.dataset) & (changes.patient_id.astype(str) == str(item.patient_id)) & (changes.signature_id == "SIG001")]
        status = True
        for _, row in observed.iterrows():
            expected = frozen[f"{row.time_window}_sample_id"]
            status &= row.baseline_sample_id == frozen.T0_sample_id and row.followup_sample_id == expected
        checks.append({"check_type": "PAIRING", "dataset": item.dataset, "item": str(item.patient_id), "expected": "frozen T0/follow-up IDs", "observed": "matched" if status else "mismatch", "absolute_difference": 0 if status else 1, "status": "PASS" if status else "FAIL"})

    expression = pd.concat([pd.read_parquet(ROOT / f"intermediate/required_expression_{dataset}.parquet") for dataset in COHORTS], ignore_index=True)
    expression_index = expression.set_index(["dataset", "sample_id", "gene"])["expression"]
    score_index = scores.set_index(["dataset", "sample_id", "signature_id"])["raw_score"]
    for signature in STAGE3_FUNCTIONS:
        candidates = scores[scores.signature_id == signature][["dataset", "sample_id"]].drop_duplicates()
        chosen = candidates.sample(n=5, random_state=20260716 + int(signature[3:]))
        for _, item in chosen.iterrows():
            genes = expression[(expression.dataset == item.dataset) & (expression.sample_id == item.sample_id)]
            values = dict(zip(genes.gene, genes.expression))
            independent = float(manual_score(signature, values))
            stored = float(score_index.loc[(item.dataset, item.sample_id, signature)])
            diff = abs(independent - stored)
            passed = diff <= 1e-10 * max(1.0, abs(independent), abs(stored))
            checks.append({"check_type": "SCORE_MANUAL", "dataset": item.dataset, "item": f"{signature}:{item.sample_id}", "expected": independent, "observed": stored, "absolute_difference": diff, "status": "PASS" if passed else "FAIL"})

    for _, row in meta[(meta.analysis_set == "PRIMARY_INDEPENDENT") & (meta.time_window.isin(["T1", "T2"]))].iterrows():
        datasets = row.cohorts.split(";")
        subset = effects[(effects.signature_id == row.signature_id) & (effects.time_window == row.time_window) & (effects.dataset.isin(datasets))]
        mean, tau = independent_reml(subset.mean_delta_z, subset.bootstrap_se)
        diff = abs(mean - row.pooled_delta_z)
        passed = diff < 1e-8 and abs(tau - row.tau2) < 1e-6
        checks.append({"check_type": "META_INPUT_AND_REML", "dataset": "MULTICOhort", "item": f"{row.signature_id}:{row.time_window}", "expected": mean, "observed": row.pooled_delta_z, "absolute_difference": diff, "status": "PASS" if passed else "FAIL"})

    for signature in ("SIG002", "SIG003", "SIG004"):
        rows = meta[(meta.signature_id == signature) & (meta.analysis_set == "PRIMARY_INDEPENDENT")]
        passed = all("GSE54514" not in str(value).split(";") for value in rows.cohorts)
        checks.append({"check_type": "DEVELOPMENT_OVERLAP_EXCLUSION", "dataset": "GSE54514", "item": signature, "expected": "excluded", "observed": "excluded" if passed else "included", "absolute_difference": 0 if passed else 1, "status": "PASS" if passed else "FAIL"})

    for path in sorted((ROOT / "03_15_main_figures").glob("Figure_*.svg")):
        passed = path.stat().st_size > 1000 and "<svg" in path.read_text(encoding="utf-8", errors="ignore")[:1000]
        checks.append({"check_type": "FIGURE_ARTIFACT", "dataset": "ALL", "item": path.name, "expected": "nonempty SVG", "observed": path.stat().st_size, "absolute_difference": 0 if passed else 1, "status": "PASS" if passed else "FAIL"})

    result = pd.DataFrame(checks)
    out = ROOT / "source_data/03_16_independent_result_verification.csv"
    result.to_csv(out, index=False)
    summary = {"checks": len(result), "passed": int((result.status == "PASS").sum()), "failed": int((result.status != "PASS").sum()), "status": "PASS" if (result.status == "PASS").all() else "FAIL"}
    (ROOT / "source_data/03_16_verification_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if summary["failed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
