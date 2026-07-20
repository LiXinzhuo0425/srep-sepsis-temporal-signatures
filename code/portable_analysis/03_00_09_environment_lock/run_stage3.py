#!/usr/bin/env python3
"""Frozen Stage 3 scoring, unblinding and longitudinal analysis pipeline."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import random
import shutil
import statistics
import subprocess
import sys
import warnings
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats
from sklearn.metrics import roc_auc_score
import statsmodels.formula.api as smf


ROOT = Path(os.environ.get("SEPSIS_SIGNATURE_ANALYSIS_ROOT", Path.cwd())).resolve()
ENV = ROOT / "03_00_09_environment_lock"
sys.path.insert(0, str(ENV))
from signatures import (  # noqa: E402
    SCORE_DIRECTION,
    STAGE3_FUNCTIONS,
    STAGE3_REQUIRED_GENES,
    TARGET_POSITIVE_STATE,
)

INPUT = ROOT / "inputs"
INTERMEDIATE = ROOT / "intermediate"
SOURCE_DATA = ROOT / "source_data"
UNBLIND_REPORTS = ROOT / "03_01_unblinding_reports"
SEED = 20260716
BOOTSTRAPS = 2000
PRIMARY_WINDOWS = ("T1", "T2")
EXTENDED_WINDOWS = ("T3", "T4")
ALL_WINDOWS = PRIMARY_WINDOWS + EXTENDED_WINDOWS
TARGET_HOURS = {"T0": 0.0, "T1": 24.0, "T2": 48.0, "T3": 72.0, "T4": 120.0}
COHORTS = ("GSE236713", "GSE57065", "GSE95233", "GSE54514", "GSE110487", "GSE8121")
PILOT = {"GSE57065", "GSE54514"}
VALIDATION_ORDER = ("GSE95233", "GSE110487", "GSE236713", "GSE8121")
DIAGNOSTIC_SIGNATURES = ("SIG001", "SIG002", "SIG003", "SIG004")
CONTROL_LABELS = {"HEALTHY", "CONTROL"}

DATA_FILES = {
    "GSE236713": INPUT / "expression/GSE236713_series_matrix.txt.gz",
    "GSE57065": INPUT / "expression/GSE57065_series_matrix.txt.gz",
    "GSE95233": INPUT / "expression/GSE95233_series_matrix.txt.gz",
    "GSE54514": INPUT / "expression/GSE54514_series_matrix.txt.gz",
    "GSE110487": INPUT / "expression/GSE110487_rawcounts.xlsx",
    "GSE8121": INPUT / "expression/GSE8121_series_matrix.txt.gz",
}

GSE236713_RAW_TAR = INPUT / "expression/GSE236713_RAW.tar"
GSE236713_RAW_REQUIRED = INTERMEDIATE / "GSE236713_raw75_log2_required.parquet"

COHORT_GROUPS = {
    "GSE236713": {"platform_class": "MICROARRAY", "platform": "AGILENT", "age_group": "ADULT", "sample_group": "PBL", "clinical_group": "BROADER_SEPSIS", "case_grade": "S2", "role": "PRESPECIFIED_NON_PILOT"},
    "GSE57065": {"platform_class": "MICROARRAY", "platform": "AFFYMETRIX", "age_group": "ADULT", "sample_group": "WHOLE_BLOOD", "clinical_group": "SEPTIC_SHOCK", "case_grade": "S2", "role": "PILOT"},
    "GSE95233": {"platform_class": "MICROARRAY", "platform": "AFFYMETRIX", "age_group": "ADULT", "sample_group": "WHOLE_BLOOD", "clinical_group": "SEPTIC_SHOCK", "case_grade": "S2", "role": "PRESPECIFIED_NON_PILOT"},
    "GSE54514": {"platform_class": "MICROARRAY", "platform": "ILLUMINA", "age_group": "ADULT", "sample_group": "WHOLE_BLOOD", "clinical_group": "BROADER_SEPSIS", "case_grade": "S2", "role": "PILOT"},
    "GSE110487": {"platform_class": "RNASEQ", "platform": "RNASEQ", "age_group": "ADULT", "sample_group": "WHOLE_BLOOD", "clinical_group": "SEPTIC_SHOCK", "case_grade": "S2", "role": "PRESPECIFIED_NON_PILOT"},
    "GSE8121": {"platform_class": "MICROARRAY", "platform": "AFFYMETRIX", "age_group": "PEDIATRIC", "sample_group": "WHOLE_BLOOD", "clinical_group": "SEPTIC_SHOCK", "case_grade": "S2", "role": "PRESPECIFIED_NON_PILOT"},
}


def ensure_dirs() -> None:
    for path in (INTERMEDIATE, SOURCE_DATA, UNBLIND_REPORTS):
        path.mkdir(parents=True, exist_ok=True)


def config() -> dict:
    locked = ENV / "analysis_config.json"
    fallback = ROOT / "work/freeze_sources.json"
    return json.loads((locked if locked.exists() else fallback).read_text(encoding="utf-8"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], columns: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        columns = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def manifest_df() -> pd.DataFrame:
    return pd.read_csv(INPUT / "metadata/longitudinal_sample_manifest.csv", dtype=str, keep_default_na=False)


def feature_map_df() -> pd.DataFrame:
    return pd.read_csv(INPUT / "metadata/platform_feature_map.csv", dtype=str, keep_default_na=False)


def overlap_lookup() -> dict[tuple[str, str], str]:
    return {(row["signature_id"], row["dataset"]): row["overlap_status"] for row in config()["overlap"]}


def load_microarray(dataset: str) -> pd.DataFrame:
    fmap = feature_map_df()
    fmap = fmap[(fmap["dataset"] == dataset) & (fmap["included_in_lock"] == "YES")]
    feature_to_gene = dict(zip(fmap["feature_id"].astype(str), fmap["canonical_symbol"].astype(str)))
    sums: dict[str, np.ndarray] = {}
    counts: dict[str, int] = defaultdict(int)
    sample_ids: list[str] = []
    in_table = False
    with gzip.open(DATA_FILES[dataset], "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                sample_ids = [value.strip('"') for value in next(handle).rstrip("\r\n").split("\t")[1:]]
                continue
            if line == "!series_matrix_table_end":
                break
            if not in_table or not line:
                continue
            fields = line.split("\t")
            feature = fields[0].strip('"')
            gene = feature_to_gene.get(feature)
            if not gene:
                continue
            values = np.asarray([float(value.strip('"')) for value in fields[1:]], dtype=float)
            if len(values) != len(sample_ids):
                raise RuntimeError(f"{dataset} {feature}: expression length mismatch")
            sums[gene] = sums.get(gene, np.zeros(len(values), dtype=float)) + values
            counts[gene] += 1
    expression = {gene: values / counts[gene] for gene, values in sums.items()}
    df = pd.DataFrame(expression, index=sample_ids).T
    df.index.name = "gene"
    return df


def prepare_gse236713_raw() -> Path:
    """Create the locked-probe, result-blind raw-data rescue matrix."""
    ensure_dirs()
    if GSE236713_RAW_REQUIRED.exists():
        return GSE236713_RAW_REQUIRED
    if not GSE236713_RAW_TAR.exists():
        raise RuntimeError("GSE236713 raw tar is required for scale-compatible scoring")
    cmd = [
        sys.executable,
        str(ENV / "preprocess_gse236713_raw.py"),
        str(GSE236713_RAW_TAR),
        str(INPUT / "metadata/platform_feature_map.csv"),
        str(GSE236713_RAW_REQUIRED),
        str(INTERMEDIATE / "GSE236713_raw_preprocessing_qc.json"),
    ]
    completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
    (INTERMEDIATE / "GSE236713_raw_preprocessing_log.txt").write_text(
        completed.stdout + completed.stderr, encoding="utf-8"
    )
    return GSE236713_RAW_REQUIRED


def load_gse236713_raw() -> pd.DataFrame:
    raw = pd.read_parquet(prepare_gse236713_raw())
    fmap = feature_map_df()
    fmap = fmap[(fmap["dataset"] == "GSE236713") & (fmap["included_in_lock"] == "YES")]
    mapped = raw.merge(
        fmap[["feature_id", "canonical_symbol"]], on="feature_id", how="inner", validate="many_to_one"
    )
    # Frozen rule: arithmetic mean across all retained probes for a gene.
    expression = mapped.groupby(["canonical_symbol", "sample_id"], sort=False)["expression"].mean().unstack("sample_id")
    expression.index.name = "gene"
    return expression


def prepare_gse110487_vst() -> Path:
    ensure_dirs()
    output = INTERMEDIATE / "GSE110487_required_gene_vst.tsv"
    if output.exists():
        return output
    raw = pd.read_excel(DATA_FILES["GSE110487"], sheet_name="rawcounts")
    raw["Geneid"] = raw["Geneid"].astype(str).str.replace(r"\..*$", "", regex=True)
    fmap = feature_map_df()
    fmap = fmap[(fmap["dataset"] == "GSE110487") & (fmap["included_in_lock"] == "YES")]
    mapping = dict(zip(fmap["feature_id"].str.replace(r"\..*$", "", regex=True), fmap["canonical_symbol"]))
    labels = raw["Geneid"].map(mapping).fillna(raw["Geneid"])
    counts = raw.drop(columns=["Geneid"]).apply(pd.to_numeric, errors="raise")
    counts.insert(0, "grouped_gene", labels)
    counts = counts.groupby("grouped_gene", sort=False).sum()
    if not np.allclose(counts.to_numpy(), np.round(counts.to_numpy())):
        raise RuntimeError("GSE110487 raw counts are not integers")
    counts = counts.astype(np.int64)
    grouped_path = INTERMEDIATE / "GSE110487_grouped_counts.tsv"
    counts.to_csv(grouped_path, sep="\t", index=True)
    required = sorted({gene for genes in STAGE3_REQUIRED_GENES.values() for gene in genes})
    required_path = INTERMEDIATE / "GSE110487_required_genes.txt"
    required_path.write_text("\n".join(required) + "\n", encoding="utf-8")
    cmd = [
        "Rscript", str(ENV / "preprocess_gse110487.R"), str(grouped_path), str(required_path), str(output)
    ]
    completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
    (INTERMEDIATE / "GSE110487_vst_log.txt").write_text(completed.stdout + completed.stderr, encoding="utf-8")
    return output


def load_gse110487() -> pd.DataFrame:
    vst = pd.read_csv(prepare_gse110487_vst(), sep="\t", index_col=0)
    manifest = manifest_df()
    rows = manifest[manifest["dataset"] == "GSE110487"]
    title_to_gsm = dict(zip(rows["title"], rows["sample_id"]))
    missing = [column for column in vst.columns if column not in title_to_gsm]
    if missing:
        raise RuntimeError(f"GSE110487 VST columns absent from manifest title map: {missing[:5]}")
    vst = vst.rename(columns=title_to_gsm)
    return vst


def load_expression(dataset: str) -> pd.DataFrame:
    if dataset == "GSE110487":
        return load_gse110487()
    if dataset == "GSE236713":
        return load_gse236713_raw()
    return load_microarray(dataset)


def geometric_mean(values: list[float]) -> float:
    if any(value <= 0 for value in values):
        raise ValueError("manual geometric mean requires positive values")
    return math.exp(sum(math.log(value) for value in values) / len(values))


def manual_score(signature: str, expression: dict[str, float]) -> float:
    if signature == "SIG001":
        return geometric_mean([expression[g] for g in ("CEACAM1", "ZDHHC19", "NMRK1", "GNA15", "BATF", "C3AR1")]) - (5.0 / 6.0) * geometric_mean([expression[g] for g in ("FAM214A", "TGFBI", "MTCH1", "RPGRIP1", "HLA-DPB1")])
    if signature == "SIG002":
        return expression["PLAC8"] + expression["LAMP1"] - expression["PLA2G7"] - expression["CEACAM4"]
    if signature == "SIG003":
        return expression["FAIM3"] / expression["PLAC8"]
    if signature == "SIG004":
        return (expression["NLRP1"] - expression["IDNK"]) / expression["PLAC8"]
    if signature == "SIG022":
        return geometric_mean([expression[g] for g in ("IFI27", "JUP", "LAX1")]) - geometric_mean([expression[g] for g in ("HK3", "TNIP1", "GPAA1", "CTSB")])
    if signature == "SIG023":
        return expression["FAM89A"] - expression["IFI44L"]
    if signature == "SIG033":
        weights = {"ADRB2": -0.4102, "CTSG": 0.1825, "CX3CR1": -0.1810, "CXCR6": 0.8549, "IL4R": -0.4270, "LTB": -0.5605, "TMSB10": -0.6836}
        return sum(weights[g] * expression[g] for g in weights)
    if signature == "SIG034":
        m1 = ("NQO2", "SLPI", "ORM1", "KLHL2", "ANXA3", "TXN", "AQP9", "BCL6", "DOK3", "PFKFB4", "TYK2")
        m2 = ("BCL2L11", "BCAT1", "BTBD7", "CEP55", "HMMR", "PRC1", "KIF15", "CAMP", "CEACAM8", "DEFA4", "LCN2", "CTSG", "AZU1")
        m3 = ("MAFB", "OASL", "UBE2L6", "VAMP5", "CCL2", "NAPA", "ATG3", "VRK2", "TMEM123", "CASP7")
        m4 = ("DOK2", "HLA-DPB1", "BUB3", "SMYD2", "SIDT1", "EXOC2", "TRIB2", "KLRB1")
        return (geometric_mean([expression[g] for g in m1]) + geometric_mean([expression[g] for g in m2])) / (geometric_mean([expression[g] for g in m3]) + geometric_mean([expression[g] for g in m4]))
    raise KeyError(signature)


def score_dataset(dataset: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    expression = load_expression(dataset)
    required_union = sorted({gene for genes in STAGE3_REQUIRED_GENES.values() for gene in genes})
    missing_union = sorted(set(required_union) - set(expression.index))
    if missing_union:
        raise RuntimeError(f"{dataset}: missing required genes {missing_union}")
    manifest = manifest_df()
    cohort_manifest = manifest[manifest["dataset"] == dataset].copy()
    sample_set = set(expression.columns)
    absent = sorted(set(cohort_manifest["sample_id"]) - sample_set)
    if absent:
        raise RuntimeError(f"{dataset}: manifest samples absent from expression matrix {absent[:5]}")
    overlaps = overlap_lookup()
    scores: list[dict] = []
    manual: list[dict] = []
    for _, row in cohort_manifest.iterrows():
        sample_id = row["sample_id"]
        values = expression[sample_id].to_dict()
        for signature, function in STAGE3_FUNCTIONS.items():
            raw_score = float(function(values))
            if not math.isfinite(raw_score):
                raise RuntimeError(f"{dataset} {sample_id} {signature}: nonfinite score")
            direction = SCORE_DIRECTION[signature]
            scores.append({
                "dataset": dataset,
                "patient_id": row["patient_id"],
                "sample_id": sample_id,
                "raw_time_label": row["raw_time_label"],
                "time_window": row["unified_time_window"],
                "relative_hours": row["relative_hours"],
                "disease_status": row["disease_status"],
                "outcome": row["outcome"],
                "main_longitudinal_case": row["main_longitudinal_case"],
                "signature_id": signature,
                "raw_score": raw_score,
                "direction_coefficient": direction,
                "oriented_score": direction * raw_score,
                "target_positive_state": TARGET_POSITIVE_STATE[signature],
                "gene_coverage_status": "COMPLETE",
                "calculation_status": "SUCCESS",
                "development_overlap_status": overlaps[(signature, dataset)],
                "analysis_role": COHORT_GROUPS[dataset]["role"],
            })
    score_df = pd.DataFrame(scores)
    qc_rows = []
    for signature, group in score_df.groupby("signature_id", sort=True):
        first_ids = sorted(group["sample_id"].unique())[:5]
        check_pass = 0
        for sample_id in first_ids:
            values = expression[sample_id].to_dict()
            code_value = float(STAGE3_FUNCTIONS[signature](values))
            manual_value = float(manual_score(signature, values))
            difference = abs(code_value - manual_value)
            status = "PASS" if difference <= 1e-10 * max(1.0, abs(code_value), abs(manual_value)) else "FAIL"
            check_pass += status == "PASS"
            manual.append({"dataset": dataset, "sample_id": sample_id, "signature_id": signature, "code_value": code_value, "independent_value": manual_value, "absolute_difference": difference, "status": status})
        qc_rows.append({
            "dataset": dataset,
            "signature_id": signature,
            "matrix_samples": expression.shape[1],
            "manifest_samples": cohort_manifest.shape[0],
            "scored_samples": len(group),
            "score_missing_n": int(group["oriented_score"].isna().sum()),
            "score_unique_n": int(group["oriented_score"].nunique()),
            "constant_score": "YES" if group["oriented_score"].nunique() < 2 else "NO",
            "score_min": float(group["oriented_score"].min()),
            "score_max": float(group["oriented_score"].max()),
            "manual_checks_passed": int(check_pass),
            "manual_checks_total": len(first_ids),
            "status": "PASS" if check_pass == len(first_ids) and group["oriented_score"].nunique() >= 2 else "FAIL",
        })
    expression_long = expression.loc[required_union].rename_axis("gene").reset_index().melt(id_vars="gene", var_name="sample_id", value_name="expression")
    expression_long.insert(0, "dataset", dataset)
    return score_df, pd.DataFrame(qc_rows), pd.DataFrame(manual)


def unblind(dataset: str, *, out_root: Path = ROOT) -> None:
    ensure_dirs()
    if dataset not in COHORTS:
        raise ValueError(dataset)
    if dataset not in PILOT:
        if not (ROOT / "STAGE3_ANALYSIS_UNLOCKED.md").exists():
            raise RuntimeError("Stage 3 analysis is not unlocked")
        index = VALIDATION_ORDER.index(dataset)
        for prior in VALIDATION_ORDER[:index]:
            if not (UNBLIND_REPORTS / f"{VALIDATION_ORDER.index(prior)+1:02d}_{prior}_unblinding_report.md").exists():
                raise RuntimeError(f"Sequential unblinding violation: {prior} must be completed first")
    score_df, qc_df, manual_df = score_dataset(dataset)
    if (qc_df["status"] != "PASS").any() or (manual_df["status"] != "PASS").any():
        raise RuntimeError(f"{dataset}: technical score QC failed")
    score_df.to_parquet(INTERMEDIATE / f"scores_{dataset}.parquet", index=False)
    qc_df.to_csv(INTERMEDIATE / f"score_qc_{dataset}.csv", index=False)
    manual_df.to_csv(INTERMEDIATE / f"manual_checks_{dataset}.csv", index=False)
    expression = load_expression(dataset)
    required_union = sorted({gene for genes in STAGE3_REQUIRED_GENES.values() for gene in genes})
    exp_long = expression.loc[required_union].rename_axis("gene").reset_index().melt(id_vars="gene", var_name="sample_id", value_name="expression")
    exp_long.insert(0, "dataset", dataset)
    exp_long.to_parquet(INTERMEDIATE / f"required_expression_{dataset}.parquet", index=False)
    cohort_scores = score_df[score_df["main_longitudinal_case"] == "YES"]
    patients = cohort_scores["patient_id"].nunique()
    samples = cohort_scores["sample_id"].nunique()
    pairs = config()["cohorts"]
    cohort = next(row for row in pairs if row["dataset"] == dataset)
    order = 0 if dataset in PILOT else VALIDATION_ORDER.index(dataset) + 1
    report = f"""# {dataset} unblinding report

**Execution order:** {order if order else 'pilot, previously unblinded'}  
**Execution date:** 2026-07-16  
**Input SHA256:** `{sha256(DATA_FILES[dataset])}`  
{f'**Raw rescue input SHA256:** `{sha256(GSE236713_RAW_TAR)}`  ' if dataset == 'GSE236713' else ''}
**Code SHA256:** `{sha256(ENV / 'run_stage3.py')}`

## Technical QC before proceeding

- Matrix samples scored: {score_df['sample_id'].nunique()}
- Longitudinal case samples: {samples}
- Independent case patients: {patients}
- Frozen expected T0-T24 pairs: {cohort['patients_T0_T1']}
- Frozen expected T0-T48 pairs: {cohort['patients_T0_T2']}
- Signature combinations passing score generation: {(qc_df['status'] == 'PASS').sum()}/8
- Constant-score combinations: {(qc_df['constant_score'] == 'YES').sum()}
- Manual arithmetic checks: {(manual_df['status'] == 'PASS').sum()}/{len(manual_df)}
- Missing-gene substitution: none
- Formula, direction or time-window changes after scoring: none
- Performance quantities generated at this step: none (no deltaZ, AUROC, ranking or trajectory plot)

**Decision:** TECHNICAL QC PASS; sequential analysis may proceed.
"""
    name = f"{order:02d}_{dataset}_unblinding_report.md" if order else f"pilot_{dataset}_technical_rescore_report.md"
    (UNBLIND_REPORTS / name).write_text(report, encoding="utf-8")
    state_path = INTERMEDIATE / "unblinding_state.json"
    state = json.loads(state_path.read_text()) if state_path.exists() else {"completed": []}
    if dataset not in state["completed"]:
        state["completed"].append(dataset)
    state["latest"] = dataset
    state_path.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"dataset": dataset, "qc": "PASS", "scores": len(score_df), "patients": patients}, indent=2))


def percentile(values: np.ndarray, q: float) -> float:
    return float(np.quantile(values, q, method="linear"))


def seed_for(*parts: str) -> int:
    return SEED + sum((index + 1) * ord(char) for index, char in enumerate("|".join(parts))) % 1_000_000


def bootstrap_mean(values: np.ndarray, *parts: str) -> tuple[float, float, float]:
    rng = np.random.default_rng(seed_for(*parts))
    n = len(values)
    draws = values[rng.integers(0, n, size=(BOOTSTRAPS, n))].mean(axis=1)
    return float(draws.std(ddof=1)), percentile(draws, 0.025), percentile(draws, 0.975)


def near_zero_scale(values: np.ndarray) -> bool:
    threshold = max(1e-8, 100 * np.finfo(float).eps * max(1.0, float(np.median(np.abs(values)))))
    return float(np.std(values, ddof=1)) <= threshold


def build_analysis_data(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    cfg = config()
    population = pd.DataFrame(cfg["population"])
    included = population[population["freeze_status"] == "INCLUDED_FROZEN"].copy()
    score_index = scores.set_index(["dataset", "sample_id", "signature_id"])
    baseline_scales: list[dict] = []
    change_rows: list[dict] = []
    manifest = manifest_df().set_index(["dataset", "sample_id"])
    for dataset in COHORTS:
        cohort_pop = population[population["dataset"] == dataset]
        for signature in STAGE3_FUNCTIONS:
            baseline_values = []
            for sample_id in cohort_pop["T0_sample_id"]:
                if not sample_id:
                    continue
                try:
                    baseline_values.append(float(score_index.loc[(dataset, sample_id, signature), "oriented_score"]))
                except KeyError:
                    pass
            base = np.asarray(baseline_values, dtype=float)
            if len(base) < 2 or near_zero_scale(base):
                sd = math.nan
                status = "BASELINE_SCALE_UNUSABLE"
            else:
                sd = float(np.std(base, ddof=1))
                status = "PASS"
            median = float(np.median(base)) if len(base) else math.nan
            mad = float(1.4826 * np.median(np.abs(base - median))) if len(base) else math.nan
            baseline_scales.append({"dataset": dataset, "signature_id": signature, "baseline_n": len(base), "baseline_sd": sd, "baseline_mad_scale": mad, "status": status})
            if status != "PASS":
                continue
            for _, patient in included[included["dataset"] == dataset].iterrows():
                t0 = patient["T0_sample_id"]
                if not t0:
                    continue
                try:
                    baseline_score = float(score_index.loc[(dataset, t0, signature), "oriented_score"])
                except KeyError:
                    continue
                for window in ALL_WINDOWS:
                    follow_id = patient[f"{window}_sample_id"]
                    if not follow_id:
                        continue
                    follow_score = float(score_index.loc[(dataset, follow_id, signature), "oriented_score"])
                    raw_change = follow_score - baseline_score
                    meta_row = manifest.loc[(dataset, follow_id)]
                    change_rows.append({
                        "dataset": dataset,
                        "patient_id": patient["patient_id"],
                        "signature_id": signature,
                        "time_window": window,
                        "baseline_sample_id": t0,
                        "followup_sample_id": follow_id,
                        "baseline_score": baseline_score,
                        "followup_score": follow_score,
                        "raw_change": raw_change,
                        "baseline_sd": sd,
                        "baseline_mad_scale": mad,
                        "delta_z": raw_change / sd,
                        "delta_z_mad": raw_change / mad if mad > 1e-8 else math.nan,
                        "time_interval_hours": meta_row["relative_hours"],
                        "outcome": meta_row["outcome"],
                        "analysis_eligible": "YES",
                        "exclusion_reason": "",
                        **COHORT_GROUPS[dataset],
                    })
    changes = pd.DataFrame(change_rows)
    scales = pd.DataFrame(baseline_scales)
    primary_patients = changes[changes["time_window"].isin(PRIMARY_WINDOWS)][["dataset", "patient_id"]].drop_duplicates()
    if len(primary_patients) != 264:
        raise RuntimeError(f"Frozen primary paired-patient count mismatch: observed {len(primary_patients)} expected 264")
    return changes, scales


def paired_effects(changes: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict] = []
    diagnostics: list[dict] = []
    for (dataset, signature, window), group in changes.groupby(["dataset", "signature_id", "time_window"], sort=True):
        values = group["delta_z"].to_numpy(float)
        se, lower, upper = bootstrap_mean(values, dataset, signature, window, "primary")
        if len(values) > 1:
            paired_corr = float(np.corrcoef(group["baseline_score"], group["followup_score"])[0, 1])
        else:
            paired_corr = math.nan
        loo = [abs((values.sum() - values[i]) / (len(values) - 1) - values.mean()) for i in range(len(values))] if len(values) > 2 else [math.nan]
        winsor = np.clip(values, np.quantile(values, 0.025), np.quantile(values, 0.975)) if len(values) >= 4 else values
        model_status = "NOT_RUN"
        model_effect = math.nan
        model_se = math.nan
        try:
            long = pd.concat([
                pd.DataFrame({"patient_id": group["patient_id"], "q_scaled": group["baseline_score"] / group["baseline_sd"], "followup": 0}),
                pd.DataFrame({"patient_id": group["patient_id"], "q_scaled": group["followup_score"] / group["baseline_sd"], "followup": 1}),
            ], ignore_index=True)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                fit = smf.mixedlm("q_scaled ~ followup", long, groups=long["patient_id"]).fit(reml=True, method="lbfgs", maxiter=500, disp=False)
            model_effect = float(fit.params["followup"])
            model_se = float(fit.bse["followup"])
            model_status = "CONVERGED" if bool(getattr(fit, "converged", False)) else "NONCONVERGED_FALLBACK_BOOTSTRAP"
            diagnostics.append({"dataset": dataset, "signature_id": signature, "time_window": window, "model": "MixedLM", "status": model_status, "warning_count": len(caught), "warning_text": " | ".join(str(item.message) for item in caught)[:1000]})
        except Exception as exc:
            model_status = "FAILED_FALLBACK_BOOTSTRAP"
            diagnostics.append({"dataset": dataset, "signature_id": signature, "time_window": window, "model": "MixedLM", "status": model_status, "warning_count": 0, "warning_text": str(exc)[:1000]})
        rows.append({
            "dataset": dataset,
            "signature_id": signature,
            "time_window": window,
            "paired_n": len(group),
            "mean_delta_z": float(values.mean()),
            "median_delta_z": float(np.median(values)),
            "q1_delta_z": float(np.quantile(values, 0.25)),
            "q3_delta_z": float(np.quantile(values, 0.75)),
            "bootstrap_se": se,
            "ci95_lower": lower,
            "ci95_upper": upper,
            "mean_raw_change": float(group["raw_change"].mean()),
            "proportion_increasing": float(np.mean(values > 0)),
            "paired_correlation": paired_corr,
            "winsorized_mean_delta_z": float(winsor.mean()),
            "max_leave_one_out_mean_shift": float(np.nanmax(loo)) if not np.all(np.isnan(loo)) else math.nan,
            "mixed_model_effect": model_effect,
            "mixed_model_se": model_se,
            "model_status": model_status,
            "analysis_role": COHORT_GROUPS[dataset]["role"],
        })
    return pd.DataFrame(rows), pd.DataFrame(diagnostics)


def reml_meta(y: np.ndarray, se: np.ndarray) -> dict[str, float]:
    y = np.asarray(y, dtype=float)
    vi = np.maximum(np.asarray(se, dtype=float) ** 2, 1e-12)
    k = len(y)
    if k < 2:
        return {"k": k, "mu": float(y[0]) if k else math.nan, "se": math.nan, "lower": math.nan, "upper": math.nan, "tau2": math.nan, "I2": math.nan, "prediction_lower": math.nan, "prediction_upper": math.nan, "p_value": math.nan}

    def objective(tau2: float) -> float:
        w = 1.0 / (vi + tau2)
        mu = float(np.sum(w * y) / np.sum(w))
        return 0.5 * (float(np.sum(np.log(vi + tau2))) + math.log(float(np.sum(w))) + float(np.sum(w * (y - mu) ** 2)))

    upper_bound = max(10.0, float(np.var(y, ddof=1) * 20 + np.max(vi) * 20))
    opt = optimize.minimize_scalar(objective, bounds=(0.0, upper_bound), method="bounded", options={"xatol": 1e-12})
    tau2 = max(0.0, float(opt.x))
    if objective(0.0) <= objective(tau2) + 1e-10:
        tau2 = 0.0
    w = 1.0 / (vi + tau2)
    mu = float(np.sum(w * y) / np.sum(w))
    q_hk = float(np.sum(w * (y - mu) ** 2) / (k - 1))
    hk_se = math.sqrt(max(q_hk, 1e-12) / float(np.sum(w)))
    crit = float(stats.t.ppf(0.975, df=k - 1))
    lower, upper = mu - crit * hk_se, mu + crit * hk_se
    t_value = mu / hk_se if hk_se > 0 else math.nan
    p_value = float(2 * stats.t.sf(abs(t_value), df=k - 1)) if math.isfinite(t_value) else math.nan
    wf = 1.0 / vi
    fixed = float(np.sum(wf * y) / np.sum(wf))
    q = float(np.sum(wf * (y - fixed) ** 2))
    i2 = max(0.0, (q - (k - 1)) / q * 100.0) if q > 0 else 0.0
    if k >= 3:
        pred_crit = float(stats.t.ppf(0.975, df=k - 2))
        pred_se = math.sqrt(tau2 + hk_se**2)
        pl, pu = mu - pred_crit * pred_se, mu + pred_crit * pred_se
    else:
        pl, pu = math.nan, math.nan
    return {"k": k, "mu": mu, "se": hk_se, "lower": lower, "upper": upper, "tau2": tau2, "I2": i2, "prediction_lower": pl, "prediction_upper": pu, "p_value": p_value}


def analysis_set_mask(effects: pd.DataFrame, signature: str, analysis_set: str) -> pd.Series:
    mask = effects["signature_id"] == signature
    if analysis_set == "PILOT_ONLY":
        mask &= effects["analysis_role"] == "PILOT"
    if analysis_set == "PRESPECIFIED_NON_PILOT_ONLY":
        mask &= effects["analysis_role"] == "PRESPECIFIED_NON_PILOT"
    if analysis_set in {"PRIMARY_INDEPENDENT", "STRICT_NEVER_USED"}:
        overlaps = overlap_lookup()
        excluded_primary = {"DEVELOPMENT_OVERLAP", "POSSIBLE_SAME_MARS_PROGRAM", "POSSIBLE_SAME_PROGRAM"}
        excluded_strict = excluded_primary | {"PRIOR_EXTERNAL_VALIDATION", "POSSIBLE_PRIOR_EXTERNAL_BENCHMARK"}
        excluded = excluded_strict if analysis_set == "STRICT_NEVER_USED" else excluded_primary
        mask &= ~effects["dataset"].map(lambda dataset: overlaps[(signature, dataset)] in excluded)
    return mask


def meta_analyses(effects: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta_rows: list[dict] = []
    loo_rows: list[dict] = []
    for signature in STAGE3_FUNCTIONS:
        for window in ALL_WINDOWS:
            base = effects[(effects["signature_id"] == signature) & (effects["time_window"] == window)]
            for analysis_set in ("ALL_COHORTS", "PILOT_ONLY", "PRESPECIFIED_NON_PILOT_ONLY", "PRIMARY_INDEPENDENT", "STRICT_NEVER_USED"):
                subset = base[analysis_set_mask(base, signature, analysis_set)] if analysis_set != "ALL_COHORTS" else base
                subset = subset[np.isfinite(subset["bootstrap_se"]) & (subset["bootstrap_se"] > 0)]
                if subset.empty:
                    continue
                meta = reml_meta(subset["mean_delta_z"].to_numpy(), subset["bootstrap_se"].to_numpy())
                direction = np.sign(subset["mean_delta_z"].to_numpy())
                pooled_sign = np.sign(meta["mu"])
                consistency = float(np.mean(direction == pooled_sign)) if pooled_sign != 0 else float(np.mean(direction == 0))
                meta_rows.append({
                    "signature_id": signature,
                    "time_window": window,
                    "analysis_set": analysis_set,
                    "cohort_n": len(subset),
                    "patient_n": int(subset["paired_n"].sum()),
                    "pooled_delta_z": meta["mu"],
                    "pooled_se": meta["se"],
                    "ci95_lower": meta["lower"],
                    "ci95_upper": meta["upper"],
                    "tau2": meta["tau2"],
                    "I2_percent": meta["I2"],
                    "prediction_lower": meta["prediction_lower"],
                    "prediction_upper": meta["prediction_upper"],
                    "direction_consistency": consistency,
                    "p_value": meta["p_value"],
                    "cohorts": ";".join(subset["dataset"]),
                })
                if len(subset) >= 3:
                    for omitted in subset["dataset"]:
                        reduced = subset[subset["dataset"] != omitted]
                        lometa = reml_meta(reduced["mean_delta_z"].to_numpy(), reduced["bootstrap_se"].to_numpy())
                        loo_rows.append({"signature_id": signature, "time_window": window, "analysis_set": analysis_set, "omitted_dataset": omitted, "cohort_n": len(reduced), "pooled_delta_z": lometa["mu"], "ci95_lower": lometa["lower"], "ci95_upper": lometa["upper"], "tau2": lometa["tau2"], "I2_percent": lometa["I2"]})
    meta_df = pd.DataFrame(meta_rows)
    primary_mask = meta_df["analysis_set"] == "PRIMARY_INDEPENDENT"
    family = primary_mask & meta_df["signature_id"].isin(DIAGNOSTIC_SIGNATURES) & meta_df["time_window"].isin(PRIMARY_WINDOWS)
    if family.any():
        meta_df.loc[family, "holm_adjusted_p"] = holm(meta_df.loc[family, "p_value"].to_numpy())
    return meta_df, pd.DataFrame(loo_rows)


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


def stability_class(row: pd.Series | dict) -> str:
    k = int(row["cohort_n"])
    n = int(row["patient_n"])
    if k < 3 or n < 60:
        return "EVIDENCE_INSUFFICIENT"
    pred_width = float(row["prediction_upper"] - row["prediction_lower"]) if math.isfinite(float(row["prediction_upper"])) else math.inf
    if float(row["I2_percent"]) >= 75 or pred_width >= 2.0:
        return "HIGHLY_HETEROGENEOUS"
    mu = float(row["pooled_delta_z"])
    if abs(mu) >= 0.25 and float(row["ci95_lower"]) * float(row["ci95_upper"]) > 0 and float(row["direction_consistency"]) >= 0.75:
        return "CONSISTENT_DRIFT"
    if abs(mu) < 0.25 and float(row["ci95_lower"]) >= -0.50 and float(row["ci95_upper"]) <= 0.50 and pred_width <= 1.0:
        return "RELATIVE_STABILITY"
    return "COHORT_DEPENDENT_DRIFT"


def stability_profiles(meta: pd.DataFrame, effects: pd.DataFrame) -> pd.DataFrame:
    rows = []
    primary = meta[meta["analysis_set"] == "PRIMARY_INDEPENDENT"].copy()
    for _, row in primary.iterrows():
        item = row.to_dict()
        item["stability_class"] = stability_class(row)
        subset = effects[(effects["signature_id"] == row["signature_id"]) & (effects["time_window"] == row["time_window"])]
        item["mean_patient_absolute_delta_z"] = float(np.average(subset["mean_delta_z"].abs(), weights=subset["paired_n"])) if len(subset) else math.nan
        rows.append(item)
    return pd.DataFrame(rows)


def within_use_comparison(changes: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for window in PRIMARY_WINDOWS:
        subset = changes[changes["time_window"] == window]
        for sig_a, sig_b in combinations(DIAGNOSTIC_SIGNATURES, 2):
            cohort_rows = []
            for dataset in COHORTS:
                a = subset[(subset["dataset"] == dataset) & (subset["signature_id"] == sig_a)][["patient_id", "delta_z"]].rename(columns={"delta_z": "a"})
                b = subset[(subset["dataset"] == dataset) & (subset["signature_id"] == sig_b)][["patient_id", "delta_z"]].rename(columns={"delta_z": "b"})
                common = a.merge(b, on="patient_id")
                if common.empty:
                    continue
                diff = np.abs(common["a"].to_numpy()) - np.abs(common["b"].to_numpy())
                se, lower, upper = bootstrap_mean(diff, dataset, sig_a, sig_b, window, "within_use")
                overlap = overlap_lookup()
                primary_excluded = any(overlap[(sig, dataset)] in {"DEVELOPMENT_OVERLAP", "POSSIBLE_SAME_MARS_PROGRAM", "POSSIBLE_SAME_PROGRAM"} for sig in (sig_a, sig_b))
                cohort_rows.append({"dataset": dataset, "signature_a": sig_a, "signature_b": sig_b, "time_window": window, "paired_n": len(common), "mean_absolute_drift_difference_a_minus_b": float(diff.mean()), "bootstrap_se": se, "ci95_lower": lower, "ci95_upper": upper, "primary_independent_eligible": "NO" if primary_excluded else "YES"})
            rows.extend(cohort_rows)
            eligible = pd.DataFrame(cohort_rows)
            if not eligible.empty:
                eligible = eligible[eligible["primary_independent_eligible"] == "YES"]
            if len(eligible) >= 2:
                meta = reml_meta(eligible["mean_absolute_drift_difference_a_minus_b"].to_numpy(), eligible["bootstrap_se"].to_numpy())
                rows.append({"dataset": "POOLED_PRIMARY_INDEPENDENT", "signature_a": sig_a, "signature_b": sig_b, "time_window": window, "paired_n": int(eligible["paired_n"].sum()), "mean_absolute_drift_difference_a_minus_b": meta["mu"], "bootstrap_se": meta["se"], "ci95_lower": meta["lower"], "ci95_upper": meta["upper"], "primary_independent_eligible": "YES", "p_value": meta["p_value"], "tau2": meta["tau2"], "I2_percent": meta["I2"]})
    result = pd.DataFrame(rows)
    for window in PRIMARY_WINDOWS:
        mask = (result["dataset"] == "POOLED_PRIMARY_INDEPENDENT") & (result["time_window"] == window)
        if mask.any():
            result.loc[mask, "holm_adjusted_p"] = holm(result.loc[mask, "p_value"].to_numpy())
    return result


def attrition_analysis(scores: pd.DataFrame, changes: pd.DataFrame) -> pd.DataFrame:
    population = pd.DataFrame(config()["population"])
    score_index = scores.set_index(["dataset", "sample_id", "signature_id"])
    rows = []
    for dataset in COHORTS:
        cohort = population[population["dataset"] == dataset]
        baseline = cohort[cohort["T0_sample_id"] != ""]
        for window in PRIMARY_WINDOWS:
            observed = baseline[baseline[f"{window}_sample_id"] != ""]
            missing = baseline[baseline[f"{window}_sample_id"] == ""]
            for signature in STAGE3_FUNCTIONS:
                obs_values = [float(score_index.loc[(dataset, sid, signature), "oriented_score"]) for sid in observed["T0_sample_id"]]
                miss_values = [float(score_index.loc[(dataset, sid, signature), "oriented_score"]) for sid in missing["T0_sample_id"]]
                difference = float(np.mean(obs_values) - np.mean(miss_values)) if obs_values and miss_values else math.nan
                rows.append({
                    "dataset": dataset,
                    "time_window": window,
                    "signature_id": signature,
                    "baseline_eligible_n": len(baseline),
                    "observed_pair_n": len(observed),
                    "missing_followup_n": len(missing),
                    "missing_percent": len(missing) / len(baseline) if len(baseline) else math.nan,
                    "death_before_window_n": 0,
                    "early_discharge_before_window_n": 0,
                    "technical_failure_n": 0,
                    "unknown_reason_n": len(missing),
                    "baseline_score_observed_mean": float(np.mean(obs_values)) if obs_values else math.nan,
                    "baseline_score_missing_mean": float(np.mean(miss_values)) if miss_values else math.nan,
                    "baseline_score_difference": difference,
                    "ipw_status": "NOT_ESTIMABLE_MISSING_AGE_SEX_BASELINE_SEVERITY",
                    "post_death_imputation": "PROHIBITED_NOT_PERFORMED",
                    "reason_note": "Public metadata do not establish event timing before the landmark; missing reasons remain UNKNOWN unless explicit.",
                })
    return pd.DataFrame(rows)


def scaling_analysis(scores: pd.DataFrame, changes: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    cohort_rows = []
    score_controls = scores[scores["disease_status"].isin(CONTROL_LABELS)]
    for (dataset, signature, window), group in changes.groupby(["dataset", "signature_id", "time_window"], sort=True):
        methods = {"BASELINE_SD": group["delta_z"].to_numpy(float)}
        if np.isfinite(group["delta_z_mad"]).all():
            methods["BASELINE_MAD"] = group["delta_z_mad"].to_numpy(float)
        controls = score_controls[(score_controls["dataset"] == dataset) & (score_controls["signature_id"] == signature)]["oriented_score"].to_numpy(float)
        if len(controls) >= 10 and np.std(controls, ddof=1) > 1e-8:
            methods["HEALTHY_CONTROL_SD"] = group["raw_change"].to_numpy(float) / np.std(controls, ddof=1)
        baseline_ref = group["baseline_score"].to_numpy(float)
        sorted_base = np.sort(baseline_ref)
        baseline_pct = np.searchsorted(sorted_base, group["baseline_score"].to_numpy(float), side="right") / len(sorted_base)
        follow_pct = np.searchsorted(sorted_base, group["followup_score"].to_numpy(float), side="right") / len(sorted_base)
        methods["BASELINE_PERCENTILE"] = follow_pct - baseline_pct
        for method, values in methods.items():
            se, lower, upper = bootstrap_mean(np.asarray(values), dataset, signature, window, method)
            cohort_rows.append({"record_type": "COHORT", "dataset": dataset, "signature_id": signature, "time_window": window, "scaling_method": method, "paired_n": len(values), "mean_change": float(np.mean(values)), "bootstrap_se": se, "ci95_lower": lower, "ci95_upper": upper})
        cohort_rows.append({"record_type": "COHORT", "dataset": dataset, "signature_id": signature, "time_window": window, "scaling_method": "RAW_SCORE", "paired_n": len(group), "mean_change": float(group["raw_change"].mean()), "bootstrap_se": math.nan, "ci95_lower": math.nan, "ci95_upper": math.nan})
    cohort_df = pd.DataFrame(cohort_rows)
    meta_rows = []
    for (signature, window, method), group in cohort_df[(cohort_df["scaling_method"] != "RAW_SCORE") & np.isfinite(cohort_df["bootstrap_se"])].groupby(["signature_id", "time_window", "scaling_method"]):
        mask = analysis_set_mask(group.rename(columns={"signature_id": "signature_id"}), signature, "PRIMARY_INDEPENDENT")
        eligible = group[mask]
        if len(eligible) < 2:
            continue
        meta = reml_meta(eligible["mean_change"].to_numpy(), eligible["bootstrap_se"].to_numpy())
        meta_rows.append({"record_type": "META_PRIMARY_INDEPENDENT", "dataset": "POOLED", "signature_id": signature, "time_window": window, "scaling_method": method, "paired_n": int(eligible["paired_n"].sum()), "cohort_n": len(eligible), "mean_change": meta["mu"], "bootstrap_se": meta["se"], "ci95_lower": meta["lower"], "ci95_upper": meta["upper"], "tau2": meta["tau2"], "I2_percent": meta["I2"], "prediction_lower": meta["prediction_lower"], "prediction_upper": meta["prediction_upper"]})
    return cohort_df, pd.DataFrame(meta_rows)


def heterogeneity_analysis(effects: pd.DataFrame) -> pd.DataFrame:
    rows = []
    subgroup_specs = {
        "platform_class": ["MICROARRAY", "RNASEQ"],
        "platform": ["AFFYMETRIX", "AGILENT", "ILLUMINA", "RNASEQ"],
        "age_group": ["ADULT", "PEDIATRIC"],
        "sample_group": ["WHOLE_BLOOD", "PBL"],
        "clinical_group": ["SEPTIC_SHOCK", "BROADER_SEPSIS"],
        "case_grade": ["S1", "S2"],
    }
    enriched = effects.copy()
    for column in subgroup_specs:
        enriched[column] = enriched["dataset"].map(lambda d: COHORT_GROUPS[d][column])
    for signature in STAGE3_FUNCTIONS:
        for window in PRIMARY_WINDOWS:
            base = enriched[(enriched["signature_id"] == signature) & (enriched["time_window"] == window)]
            for factor, levels in subgroup_specs.items():
                level_results = []
                for level in levels:
                    group = base[base[factor] == level]
                    if group.empty:
                        rows.append({"signature_id": signature, "time_window": window, "factor": factor, "level": level, "cohort_n": 0, "patient_n": 0, "analysis_status": "NOT_AVAILABLE"})
                        continue
                    meta = reml_meta(group["mean_delta_z"].to_numpy(), group["bootstrap_se"].to_numpy())
                    item = {"signature_id": signature, "time_window": window, "factor": factor, "level": level, "cohort_n": len(group), "patient_n": int(group["paired_n"].sum()), "pooled_delta_z": meta["mu"], "ci95_lower": meta["lower"], "ci95_upper": meta["upper"], "tau2": meta["tau2"], "I2_percent": meta["I2"], "analysis_status": "FORMAL_SUBGROUP_META" if len(group) >= 2 else "DESCRIPTIVE_SINGLE_COHORT"}
                    rows.append(item)
                    level_results.append((level, meta, len(group)))
                if len(level_results) == 2 and all(item[2] >= 2 for item in level_results):
                    (l1, m1, _), (l2, m2, _) = level_results
                    se_diff = math.sqrt(m1["se"] ** 2 + m2["se"] ** 2)
                    p = float(2 * stats.norm.sf(abs((m1["mu"] - m2["mu"]) / se_diff))) if se_diff > 0 else math.nan
                    rows.append({"signature_id": signature, "time_window": window, "factor": factor, "level": f"INTERACTION_{l1}_VS_{l2}", "cohort_n": sum(item[2] for item in level_results), "patient_n": math.nan, "pooled_delta_z": m1["mu"] - m2["mu"], "ci95_lower": (m1["mu"] - m2["mu"]) - 1.96 * se_diff, "ci95_upper": (m1["mu"] - m2["mu"]) + 1.96 * se_diff, "interaction_p": p, "analysis_status": "FORMAL_INTERACTION"})
    return pd.DataFrame(rows)


def auroc_analysis(scores: pd.DataFrame, changes: pd.DataFrame) -> pd.DataFrame:
    rows = []
    score_index = scores.set_index(["dataset", "sample_id", "signature_id"])
    controls_all = scores[scores["disease_status"].isin(CONTROL_LABELS)]
    for dataset in COHORTS:
        controls_ds = controls_all[controls_all["dataset"] == dataset]
        if controls_ds["sample_id"].nunique() < 10:
            continue
        for signature in DIAGNOSTIC_SIGNATURES:
            control_scores = controls_ds[controls_ds["signature_id"] == signature].drop_duplicates("sample_id")["oriented_score"].to_numpy(float)
            for window in PRIMARY_WINDOWS:
                group = changes[(changes["dataset"] == dataset) & (changes["signature_id"] == signature) & (changes["time_window"] == window)]
                if len(group) < 10:
                    continue
                baseline = group["baseline_score"].to_numpy(float)
                followup = group["followup_score"].to_numpy(float)
                labels = np.r_[np.ones(len(group)), np.zeros(len(control_scores))]
                auc0 = roc_auc_score(labels, np.r_[baseline, control_scores])
                auct = roc_auc_score(labels, np.r_[followup, control_scores])
                rng = np.random.default_rng(seed_for(dataset, signature, window, "auroc"))
                deltas = []
                for _ in range(BOOTSTRAPS):
                    pi = rng.integers(0, len(group), len(group))
                    ci = rng.integers(0, len(control_scores), len(control_scores))
                    c = control_scores[ci]
                    y = np.r_[np.ones(len(pi)), np.zeros(len(ci))]
                    deltas.append(roc_auc_score(y, np.r_[followup[pi], c]) - roc_auc_score(y, np.r_[baseline[pi], c]))
                rows.append({"dataset": dataset, "signature_id": signature, "time_window": window, "paired_case_n": len(group), "fixed_control_n": len(control_scores), "baseline_auroc": auc0, "followup_auroc": auct, "delta_auroc": auct - auc0, "delta_ci95_lower": float(np.quantile(deltas, 0.025)), "delta_ci95_upper": float(np.quantile(deltas, 0.975)), "analysis_status": "EXPLORATORY_FIXED_CONTROLS", "interpretation_lock": "Sampling-window dependence during biological evolution; not pure technical degradation"})
    return pd.DataFrame(rows)


def evidence_grades(meta: pd.DataFrame, profiles: pd.DataFrame, scaling_meta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for signature in STAGE3_FUNCTIONS:
        for window in PRIMARY_WINDOWS:
            m = meta[(meta["signature_id"] == signature) & (meta["time_window"] == window)]
            primary = m[m["analysis_set"] == "PRIMARY_INDEPENDENT"]
            if primary.empty:
                continue
            p = primary.iloc[0]
            primary_class = stability_class(p)
            sensitivities = []
            effects = []
            for name in ("ALL_COHORTS", "STRICT_NEVER_USED", "PRESPECIFIED_NON_PILOT_ONLY"):
                row = m[m["analysis_set"] == name]
                if not row.empty:
                    sensitivities.append(stability_class(row.iloc[0]))
                    effects.append(float(row.iloc[0]["pooled_delta_z"]))
            mad = scaling_meta[(scaling_meta["signature_id"] == signature) & (scaling_meta["time_window"] == window) & (scaling_meta["scaling_method"] == "BASELINE_MAD")]
            if not mad.empty:
                effects.append(float(mad.iloc[0]["mean_change"]))
            sign_reversal = len({int(np.sign(value)) for value in effects if abs(value) >= 0.1}) > 1
            agreement = sum(item == primary_class for item in sensitivities)
            if int(p["cohort_n"]) < 3 or int(p["patient_n"]) < 60:
                grade = "EVIDENCE_INSUFFICIENT"
            elif sign_reversal and agreement < 2:
                grade = "UNSTABLE"
            elif primary_class in {"HIGHLY_HETEROGENEOUS", "COHORT_DEPENDENT_DRIFT"}:
                grade = "COHORT_DEPENDENT"
            elif int(p["cohort_n"]) >= 4 and float(p["I2_percent"]) < 50 and agreement == len(sensitivities):
                grade = "HIGH_CONSISTENCY"
            elif agreement >= 2 and float(p["I2_percent"]) < 75:
                grade = "MODERATE_CONSISTENCY"
            else:
                grade = "COHORT_DEPENDENT"
            rows.append({"signature_id": signature, "time_window": window, "primary_stability_class": primary_class, "evidence_grade": grade, "primary_cohort_n": int(p["cohort_n"]), "primary_patient_n": int(p["patient_n"]), "pooled_delta_z": float(p["pooled_delta_z"]), "ci95_lower": float(p["ci95_lower"]), "ci95_upper": float(p["ci95_upper"]), "prediction_lower": float(p["prediction_lower"]), "prediction_upper": float(p["prediction_upper"]), "I2_percent": float(p["I2_percent"]), "sensitivity_class_agreement_n": agreement, "sensitivity_class_evaluable_n": len(sensitivities), "material_sign_reversal": "YES" if sign_reversal else "NO", "evidence_note": "Grade follows frozen multi-criterion rule; no single P-value gate."})
    return pd.DataFrame(rows)


def model_diagnostics_html(diagnostics: pd.DataFrame, path: Path) -> None:
    html = "<html><head><meta charset='utf-8'><style>body{font-family:Arial;margin:32px;color:#17212b}table{border-collapse:collapse;font-size:12px}th{background:#18324a;color:white}th,td{padding:6px 8px;border:1px solid #d4dee7;vertical-align:top}.pass{color:#176b3a}</style></head><body>"
    html += "<h1>Stage 3 model diagnostics</h1><p>Mixed models are companion estimates; the frozen paired bootstrap remains available when convergence fails.</p>"
    html += diagnostics.to_html(index=False, escape=True)
    html += "</body></html>"
    path.write_text(html, encoding="utf-8")


def analyze(out_root: Path = ROOT) -> None:
    ensure_dirs()
    if not (ROOT / "STAGE3_ANALYSIS_UNLOCKED.md").exists():
        raise RuntimeError("Stage 3 analysis is not unlocked")
    missing = [dataset for dataset in COHORTS if not (INTERMEDIATE / f"scores_{dataset}.parquet").exists()]
    if missing:
        raise RuntimeError(f"Datasets not scored: {missing}")
    scores = pd.concat([pd.read_parquet(INTERMEDIATE / f"scores_{dataset}.parquet") for dataset in COHORTS], ignore_index=True)
    qc = pd.concat([pd.read_csv(INTERMEDIATE / f"score_qc_{dataset}.csv") for dataset in COHORTS], ignore_index=True)
    if (qc["status"] != "PASS").any():
        raise RuntimeError("Score QC contains failures")
    scores.to_parquet(out_root / "03_02_full_signature_score_matrix.parquet", index=False)
    qc.to_csv(SOURCE_DATA / "03_02_score_generation_qc.csv", index=False)
    changes, scales = build_analysis_data(scores)
    analysis_dir = out_root / "03_03_longitudinal_analysis_datasets"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    changes.to_parquet(analysis_dir / "all_windows.parquet", index=False)
    for window, label in (("T1", "T0_to_T24"), ("T2", "T0_to_T48"), ("T3", "T0_to_T72"), ("T4", "T0_to_day5")):
        changes[changes["time_window"] == window].to_parquet(analysis_dir / f"{label}.parquet", index=False)
    scales.to_csv(SOURCE_DATA / "baseline_scales.csv", index=False)
    flow = changes.groupby(["dataset", "time_window"])["patient_id"].nunique().reset_index(name="paired_patients")
    flow.to_csv(SOURCE_DATA / "03_03_pairing_flow_report.csv", index=False)

    effects, diagnostics = paired_effects(changes)
    effects.to_csv(SOURCE_DATA / "03_04_cohort_signature_primary_effects.csv", index=False)
    model_diagnostics_html(diagnostics, out_root / "03_04_model_diagnostics.html")
    diagnostics.to_csv(SOURCE_DATA / "03_04_model_diagnostics.csv", index=False)

    meta, loo = meta_analyses(effects)
    meta.to_csv(SOURCE_DATA / "03_05_signature_level_meta_analysis.csv", index=False)
    loo.to_csv(SOURCE_DATA / "03_05_leave_one_dataset_out.csv", index=False)
    profiles = stability_profiles(meta, effects)
    profiles.to_csv(SOURCE_DATA / "03_06_temporal_stability_profile.csv", index=False)
    within = within_use_comparison(changes)
    within.to_csv(SOURCE_DATA / "03_07_within_intended_use_comparison.csv", index=False)

    development = meta[meta["analysis_set"].isin(["ALL_COHORTS", "PILOT_ONLY", "PRIMARY_INDEPENDENT", "STRICT_NEVER_USED", "PRESPECIFIED_NON_PILOT_ONLY"])].copy()
    development.to_csv(SOURCE_DATA / "03_08_development_cohort_exclusion_analysis.csv", index=False)
    attrition = attrition_analysis(scores, changes)
    attrition.to_csv(SOURCE_DATA / "03_09_attrition_and_missingness_analysis.csv", index=False)
    heterogeneity = heterogeneity_analysis(effects)
    heterogeneity.to_csv(SOURCE_DATA / "03_10_platform_and_population_heterogeneity.csv", index=False)
    scaling_cohort, scaling_meta = scaling_analysis(scores, changes)
    pd.concat([scaling_cohort, scaling_meta], ignore_index=True, sort=False).to_csv(SOURCE_DATA / "03_11_scaling_sensitivity_analysis.csv", index=False)
    extended = effects[effects["time_window"].isin(EXTENDED_WINDOWS)].copy()
    extended.to_csv(SOURCE_DATA / "03_12_extended_timepoint_analysis.csv", index=False)
    auroc = auroc_analysis(scores, changes)
    auroc.to_csv(SOURCE_DATA / "03_13_time_dependent_discrimination_analysis.csv", index=False)
    grades = evidence_grades(meta, profiles, scaling_meta)
    grades.to_csv(SOURCE_DATA / "03_14_signature_evidence_grade.csv", index=False)

    summary = {
        "generated_at": "2026-07-16",
        "cohorts": len(COHORTS),
        "signatures": len(STAGE3_FUNCTIONS),
        "dataset_signature_combinations": len(COHORTS) * len(STAGE3_FUNCTIONS),
        "score_rows": len(scores),
        "primary_unique_paired_patients": int(changes[changes["time_window"].isin(PRIMARY_WINDOWS)][["dataset", "patient_id"]].drop_duplicates().shape[0]),
        "T24_change_rows": int((changes["time_window"] == "T1").sum()),
        "T48_change_rows": int((changes["time_window"] == "T2").sum()),
        "cohort_effect_rows": len(effects),
        "primary_meta_rows": int(((meta["analysis_set"] == "PRIMARY_INDEPENDENT") & meta["time_window"].isin(PRIMARY_WINDOWS)).sum()),
        "model_converged_rows": int((diagnostics["status"] == "CONVERGED").sum()),
        "model_fallback_rows": int((diagnostics["status"] != "CONVERGED").sum()),
        "fixed_threshold_analysis": False,
        "calibration_or_brier": False,
        "retraining": False,
    }
    (SOURCE_DATA / "stage3_analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


def pilot_smoke() -> None:
    ensure_dirs()
    for dataset in ("GSE57065", "GSE54514"):
        unblind(dataset)
    print("pilot smoke complete")


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    p_unblind = sub.add_parser("unblind")
    p_unblind.add_argument("dataset", choices=COHORTS)
    sub.add_parser("pilot-smoke")
    sub.add_parser("prepare-rnaseq")
    sub.add_parser("analyze")
    args = parser.parse_args()
    if args.command == "unblind":
        unblind(args.dataset)
    elif args.command == "pilot-smoke":
        pilot_smoke()
    elif args.command == "prepare-rnaseq":
        print(prepare_gse110487_vst())
    elif args.command == "analyze":
        analyze()


if __name__ == "__main__":
    main()
