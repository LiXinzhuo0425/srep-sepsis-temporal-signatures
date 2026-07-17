#!/usr/bin/env python3
"""Frozen Stage 4 pathway feasibility, longitudinal change and coupling analyses."""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(os.environ.get("SEPSIS_SIGNATURE_ANALYSIS_ROOT", Path.cwd())).resolve()
ENV = ROOT / "04_environment"
PATHDATA = ROOT / "04_pathway_data"
ANN = PATHDATA / "platform_annotations"
GENESETS = PATHDATA / "gene_sets"
SOURCE = ROOT / "04_source_data"
INTERMEDIATE = ROOT / "intermediate"
INPUT = ROOT / "inputs"
STAGE3_ENV = ROOT / "03_00_09_environment_lock"
sys.path.insert(0, str(STAGE3_ENV))
from run_stage3 import reml_meta  # noqa: E402


SEED = 20260716
BOOTSTRAPS = 2000
FORCE_REBUILD = os.environ.get("STAGE4_FORCE_REBUILD", "0") == "1"
COHORTS = ("GSE236713", "GSE57065", "GSE95233", "GSE54514", "GSE110487", "GSE8121")
PRIMARY_WINDOWS = ("T1", "T2")
WINDOW_LABEL = {"T1": "T24", "T2": "T48"}
SERIES = {d: INPUT / "expression" / f"{d}_series_matrix.txt.gz" for d in COHORTS if d != "GSE110487"}
PLATFORM = {
    "GSE236713": "GPL17077",
    "GSE57065": "GPL570",
    "GSE95233": "GPL570",
    "GSE54514": "GPL6947",
    "GSE110487": "RNASEQ_GENCODE25",
    "GSE8121": "GPL570",
}

PRIMARY_HALLMARK = {
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_COAGULATION",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_APOPTOSIS",
    "HALLMARK_ALLOGRAFT_REJECTION",
}
REACTOME = {
    "Neutrophil degranulation",
    "Cross-presentation of particulate exogenous antigens (phagosomes)",
    "TCR signaling",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def bh_fdr(values: pd.Series) -> pd.Series:
    x = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    out = np.full(x.shape, np.nan)
    good = np.isfinite(x)
    if not good.any():
        return pd.Series(out, index=values.index)
    p = x[good]
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    restored = np.empty_like(q)
    restored[order] = q
    out[good] = restored
    return pd.Series(out, index=values.index)


def parse_gmt(path: Path, source: str) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            rows.append({"pathway": parts[0], "description": parts[1], "genes": set(parts[2:]), "source": source})
    return rows


def frozen_gene_sets() -> list[dict]:
    hallmark = parse_gmt(GENESETS / "h.all.v2026.1.Hs.symbols.gmt", "MSigDB_HALLMARK_2026.1.Hs")
    for row in hallmark:
        row["tier"] = "PRESET_PRIMARY" if row["pathway"] in PRIMARY_HALLMARK else "HALLMARK_SECONDARY"
    reactome_all = parse_gmt(GENESETS / "reactome_current/ReactomePathways.gmt", "REACTOME_CURRENT_2026-07-16")
    reactome = [row for row in reactome_all if row["pathway"] in REACTOME]
    for row in reactome:
        row["tier"] = "REACTOME_SUPPLEMENT"
    if len(hallmark) != 50 or len([r for r in hallmark if r["tier"] == "PRESET_PRIMARY"]) != 9 or len(reactome) != 3:
        raise RuntimeError("Frozen pathway collection is incomplete")
    return hallmark + reactome


def save_pathway_lock(gene_sets: list[dict]) -> None:
    lock = {
        "version": "stage4_pathway_v1.0",
        "freeze_date": "2026-07-16",
        "generated_before_pathway_scores": True,
        "primary_method": "centered_singscore_mean_percentile_rank_minus_0.5",
        "sensitivity_method": "fixed_alpha_0.25_ssgsea_rank_walk",
        "coverage_threshold_fraction": 0.70,
        "coverage_threshold_count": 10,
        "primary_pathway_minimum_compatible_cohorts": 4,
        "primary_windows": list(PRIMARY_WINDOWS),
        "within_cohort_standardization": "baseline_patient_SD",
        "baseline_sd_floor": 1e-8,
        "cohort_summary": "paired_patient_change_mean_and_bootstrap_CI",
        "cross_cohort_summary": "REML_random_effects_with_prediction_interval",
        "coupling": "within_cohort_Spearman_then_Fisher_z_REML",
        "multiplicity": "BH_FDR_separate_by_window_and_tier",
        "hallmark_count": 50,
        "preset_primary_hallmark_count": 9,
        "reactome_supplement_count": 3,
        "source_hashes": {
            "msigdb_hallmark": sha256(GENESETS / "h.all.v2026.1.Hs.symbols.gmt"),
            "reactome_gmt_zip": sha256(GENESETS / "ReactomePathways.gmt.zip"),
            "gpl570": sha256(ANN / "GPL570.annot.gz"),
            "gpl6947": sha256(ANN / "GPL6947.annot.gz"),
            "gpl17077": sha256(ANN / "GPL17077.txt"),
            "gencode_v25": sha256(ANN / "gencode.v25.annotation.gtf.gz"),
        },
        "gene_sets": [{"pathway": r["pathway"], "source": r["source"], "tier": r["tier"], "gene_count": len(r["genes"])} for r in gene_sets],
    }
    path = ENV / "pathway_analysis_lock.json"
    if path.exists():
        old = json.loads(path.read_text())
        old.pop("runtime_generated_at", None)
        if old != lock:
            raise RuntimeError("Existing pathway lock differs from current frozen specification")
    else:
        path.write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8")


def parse_geo_annotation(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f):
            if line.startswith("ID\t"):
                skip = i
                break
        else:
            raise RuntimeError(f"No annotation table in {path}")
    df = pd.read_csv(path, sep="\t", compression="gzip", skiprows=skip, dtype=str, keep_default_na=False)
    df = df[["ID", "Gene symbol"]].rename(columns={"ID": "feature_id", "Gene symbol": "gene"})
    df["gene"] = df["gene"].str.strip()
    df = df[(df["gene"] != "") & ~df["gene"].str.contains("///", regex=False)]
    df = df[~df["gene"].str.contains(r"\s|;|,", regex=True)]
    return df.drop_duplicates("feature_id")


def parse_gpl17077() -> pd.DataFrame:
    path = ANN / "GPL17077.txt"
    with path.open(encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f):
            if line.startswith("ID\tSPOT_ID\t"):
                skip = i
                break
        else:
            raise RuntimeError("No GPL17077 table")
    df = pd.read_csv(path, sep="\t", skiprows=skip, dtype=str, keep_default_na=False, nrows=62976)
    df = df[["ID", "GENE_SYMBOL", "CONTROL_TYPE"]].rename(columns={"ID": "feature_id", "GENE_SYMBOL": "gene"})
    df["gene"] = df["gene"].str.strip()
    df = df[(df["CONTROL_TYPE"].str.upper() == "FALSE") & (df["gene"] != "")]
    df = df[~df["gene"].str.contains(r"\s|;|,|///", regex=True)]
    return df[["feature_id", "gene"]].drop_duplicates("feature_id")


def platform_maps() -> dict[str, pd.DataFrame]:
    return {
        "GPL570": parse_geo_annotation(ANN / "GPL570.annot.gz"),
        "GPL6947": parse_geo_annotation(ANN / "GPL6947.annot.gz"),
        "GPL17077": parse_gpl17077(),
    }


def parse_gencode() -> dict[str, str]:
    mapping: dict[str, str] = {}
    gene_re = re.compile(r'gene_id "([^"]+)";')
    name_re = re.compile(r'gene_name "([^"]+)";')
    with gzip.open(ANN / "gencode.v25.annotation.gtf.gz", "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            gene = gene_re.search(fields[8])
            name = name_re.search(fields[8])
            if gene and name:
                mapping[gene.group(1).split(".")[0]] = name.group(1)
    return mapping


def load_series_matrix(dataset: str, amap: pd.DataFrame) -> pd.DataFrame:
    fmap = dict(zip(amap["feature_id"], amap["gene"]))
    sums: dict[str, np.ndarray] = {}
    counts: dict[str, int] = defaultdict(int)
    sample_ids: list[str] = []
    in_table = False
    with gzip.open(SERIES[dataset], "rt", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\r\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                sample_ids = [v.strip('"') for v in next(f).rstrip("\r\n").split("\t")[1:]]
                continue
            if line == "!series_matrix_table_end":
                break
            if not in_table or not line:
                continue
            fields = line.split("\t")
            gene = fmap.get(fields[0].strip('"'))
            if not gene:
                continue
            try:
                values = np.asarray([float(v.strip('"')) for v in fields[1:]], dtype=np.float32)
            except ValueError:
                continue
            if len(values) != len(sample_ids) or not np.all(np.isfinite(values)):
                continue
            if gene in sums:
                sums[gene] += values
            else:
                sums[gene] = values.copy()
            counts[gene] += 1
    expression = pd.DataFrame({g: v / counts[g] for g, v in sums.items()}, index=sample_ids).T
    expression.index.name = "gene"
    return expression


def prepare_gse110487(gencode: dict[str, str]) -> pd.DataFrame:
    output = PATHDATA / "full_expression_GSE110487.parquet"
    if output.exists() and not FORCE_REBUILD:
        return pd.read_parquet(output)
    raw = pd.read_excel(INPUT / "expression/GSE110487_rawcounts.xlsx", sheet_name="rawcounts")
    ids = raw.pop("Geneid").astype(str).str.replace(r"\..*$", "", regex=True)
    symbols = ids.map(gencode)
    counts = raw.apply(pd.to_numeric, errors="raise")
    counts.insert(0, "gene", symbols)
    counts = counts[counts["gene"].notna() & (counts["gene"] != "")]
    counts = counts.groupby("gene", sort=False).sum()
    if not np.allclose(counts.to_numpy(), np.round(counts.to_numpy())):
        raise RuntimeError("GSE110487 counts are not integer-valued")
    grouped = PATHDATA / "GSE110487_full_grouped_counts.tsv"
    counts.astype(np.int64).to_csv(grouped, sep="\t")
    vst_tsv = PATHDATA / "GSE110487_full_vst.tsv"
    cmd = ["Rscript", str(ENV / "preprocess_gse110487_full.R"), str(grouped), str(vst_tsv)]
    completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
    (PATHDATA / "GSE110487_full_vst.log").write_text(completed.stdout + completed.stderr, encoding="utf-8")
    vst = pd.read_csv(vst_tsv, sep="\t", index_col=0)
    manifest = pd.read_csv(INPUT / "metadata/longitudinal_sample_manifest.csv", dtype=str, keep_default_na=False)
    rows = manifest[manifest["dataset"] == "GSE110487"]
    title_to_gsm = dict(zip(rows["title"], rows["sample_id"]))
    vst = vst.rename(columns=title_to_gsm)
    vst.to_parquet(output)
    return vst


def load_full_expression(dataset: str, maps: dict[str, pd.DataFrame], gencode: dict[str, str]) -> pd.DataFrame:
    output = PATHDATA / f"full_expression_{dataset}.parquet"
    if output.exists() and not FORCE_REBUILD:
        return pd.read_parquet(output)
    if dataset == "GSE110487":
        return prepare_gse110487(gencode)
    expression = load_series_matrix(dataset, maps[PLATFORM[dataset]])
    expression.to_parquet(output)
    return expression


def centered_singscore(expr: pd.DataFrame, genes: set[str]) -> pd.Series:
    present = sorted(set(expr.index) & genes)
    ranks = expr.rank(axis=0, method="average", pct=True)
    return ranks.loc[present].mean(axis=0) - 0.5


def ssgsea_rank_walk(expr: pd.DataFrame, genes: set[str], alpha: float = 0.25) -> pd.Series:
    present = set(expr.index) & genes
    n = len(expr)
    nh = len(present)
    if nh == 0 or nh == n:
        return pd.Series(np.nan, index=expr.columns)
    out = {}
    for sample in expr.columns:
        order = expr[sample].sort_values(ascending=False).index.to_numpy()
        hit = np.fromiter((g in present for g in order), dtype=bool, count=n)
        descending_rank = np.arange(n, 0, -1, dtype=float)
        weights = np.where(hit, descending_rank ** alpha, 0.0)
        hit_step = weights / weights.sum()
        miss_step = np.where(~hit, 1.0 / (n - nh), 0.0)
        running = np.cumsum(hit_step - miss_step)
        out[sample] = float(running.sum() / n)
    return pd.Series(out)


def score_pathways(expr: pd.DataFrame, gene_sets: list[dict], dataset: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    coverage_rows = []
    score_rows = []
    gene_universe = set(expr.index)
    for gs in gene_sets:
        present = gene_universe & gs["genes"]
        frac = len(present) / len(gs["genes"])
        compatible = frac >= 0.70 and len(present) >= 10
        coverage_rows.append({
            "dataset": dataset, "platform": PLATFORM[dataset], "pathway": gs["pathway"], "tier": gs["tier"],
            "source": gs["source"], "reference_gene_count": len(gs["genes"]), "mapped_gene_count": len(present),
            "coverage_fraction": frac, "coverage_threshold_fraction": 0.70, "coverage_threshold_count": 10,
            "compatible": "YES" if compatible else "NO",
        })
        if not compatible:
            continue
        primary = centered_singscore(expr, gs["genes"])
        sensitivity = ssgsea_rank_walk(expr, gs["genes"])
        for sample in expr.columns:
            score_rows.append({
                "dataset": dataset, "sample_id": sample, "pathway": gs["pathway"], "tier": gs["tier"],
                "singscore": float(primary[sample]), "ssgsea_rank_walk": float(sensitivity[sample]),
                "mapped_gene_count": len(present), "coverage_fraction": frac,
            })
    return pd.DataFrame(coverage_rows), pd.DataFrame(score_rows)


def paired_pathway_changes(scores: pd.DataFrame) -> pd.DataFrame:
    pairs = pd.read_parquet(ROOT / "03_03_longitudinal_analysis_datasets/all_windows.parquet")
    pairs = pairs[pairs["time_window"].isin(PRIMARY_WINDOWS)]
    unique_pairs = pairs[["dataset", "patient_id", "time_window", "baseline_sample_id", "followup_sample_id"]].drop_duplicates()
    wide = scores.set_index(["dataset", "sample_id", "pathway", "tier"])[["singscore", "ssgsea_rank_walk"]]
    rows = []
    for r in unique_pairs.itertuples(index=False):
        subset = scores[(scores["dataset"] == r.dataset) & (scores["sample_id"].isin([r.baseline_sample_id, r.followup_sample_id]))]
        for (pathway, tier), g in subset.groupby(["pathway", "tier"]):
            values = g.set_index("sample_id")
            if r.baseline_sample_id not in values.index or r.followup_sample_id not in values.index:
                continue
            rows.append({
                "dataset": r.dataset, "patient_id": r.patient_id, "time_window": r.time_window,
                "pathway": pathway, "tier": tier,
                "baseline_singscore": values.loc[r.baseline_sample_id, "singscore"],
                "followup_singscore": values.loc[r.followup_sample_id, "singscore"],
                "raw_change": values.loc[r.followup_sample_id, "singscore"] - values.loc[r.baseline_sample_id, "singscore"],
                "baseline_ssgsea": values.loc[r.baseline_sample_id, "ssgsea_rank_walk"],
                "followup_ssgsea": values.loc[r.followup_sample_id, "ssgsea_rank_walk"],
                "ssgsea_change": values.loc[r.followup_sample_id, "ssgsea_rank_walk"] - values.loc[r.baseline_sample_id, "ssgsea_rank_walk"],
            })
    out = pd.DataFrame(rows)
    baseline_scales = out.groupby(["dataset", "pathway"])["baseline_singscore"].std(ddof=1).rename("baseline_sd")
    out = out.merge(baseline_scales, on=["dataset", "pathway"], how="left")
    out["delta_z"] = out["raw_change"] / out["baseline_sd"]
    out.loc[out["baseline_sd"] <= 1e-8, "delta_z"] = np.nan
    ssgsea_scales = out.groupby(["dataset", "pathway"])["baseline_ssgsea"].std(ddof=1).rename("baseline_ssgsea_sd")
    out = out.merge(ssgsea_scales, on=["dataset", "pathway"], how="left")
    out["ssgsea_delta_z"] = out["ssgsea_change"] / out["baseline_ssgsea_sd"]
    out.loc[out["baseline_ssgsea_sd"] <= 1e-8, "ssgsea_delta_z"] = np.nan
    return out


def bootstrap_ci(values: np.ndarray, seed_key: str) -> tuple[float, float]:
    digest = int(hashlib.sha256(seed_key.encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(SEED ^ digest)
    n = len(values)
    draws = rng.choice(values, size=(BOOTSTRAPS, n), replace=True).mean(axis=1)
    return tuple(np.quantile(draws, [0.025, 0.975]))


def cohort_pathway_effects(changes: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, g in changes.groupby(["dataset", "time_window", "pathway", "tier"]):
        values = g["delta_z"].dropna().to_numpy(dtype=float)
        sensitivity_values = g["ssgsea_delta_z"].dropna().to_numpy(dtype=float)
        if len(values) < 4:
            continue
        lo, hi = bootstrap_ci(values, "|".join(keys))
        rows.append({
            "dataset": keys[0], "time_window": keys[1], "time_label": WINDOW_LABEL[keys[1]],
            "pathway": keys[2], "tier": keys[3], "paired_n": len(values),
            "mean_delta_z": float(values.mean()), "median_delta_z": float(np.median(values)),
            "se_delta_z": float(values.std(ddof=1) / math.sqrt(len(values))),
            "bootstrap_ci_low": lo, "bootstrap_ci_high": hi,
            "proportion_increasing": float((values > 0).mean()),
            "mean_ssgsea_delta_z": float(np.mean(sensitivity_values)) if len(sensitivity_values) else math.nan,
            "se_ssgsea_delta_z": float(np.std(sensitivity_values, ddof=1) / math.sqrt(len(sensitivity_values))) if len(sensitivity_values) > 1 else math.nan,
        })
    result = pd.DataFrame(rows)
    result["p_value"] = 2 * stats.norm.sf(np.abs(result["mean_delta_z"] / result["se_delta_z"]))
    result["fdr_within_dataset_window_tier"] = result.groupby(["dataset", "time_window", "tier"], group_keys=False)["p_value"].apply(bh_fdr)
    return result


def meta_pathway_effects(effects: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, g in effects.groupby(["time_window", "pathway", "tier"]):
        g = g[np.isfinite(g["se_delta_z"]) & (g["se_delta_z"] > 0)]
        if len(g) < 2:
            continue
        m = reml_meta(g["mean_delta_z"].to_numpy(float), g["se_delta_z"].to_numpy(float))
        rows.append({
            "time_window": keys[0], "time_label": WINDOW_LABEL[keys[0]], "pathway": keys[1], "tier": keys[2],
            "cohort_count": len(g), "total_paired_n": int(g["paired_n"].sum()),
            "pooled_delta_z": m["mu"], "ci_low": m["lower"], "ci_high": m["upper"],
            "prediction_low": m["prediction_lower"], "prediction_high": m["prediction_upper"],
            "tau2": m["tau2"], "i2_percent": m["I2"],
            "direction_consistency": float(max((g["mean_delta_z"] > 0).mean(), (g["mean_delta_z"] < 0).mean())),
            "p_value": m["p_value"],
        })
    result = pd.DataFrame(rows)
    result["fdr_within_window_tier"] = result.groupby(["time_window", "tier"], group_keys=False)["p_value"].apply(bh_fdr)
    return result


def signature_pathway_coupling(changes: pd.DataFrame, pathway_change_column: str = "delta_z", method: str = "SINGSCORE") -> tuple[pd.DataFrame, pd.DataFrame]:
    stage3 = pd.read_parquet(ROOT / "03_03_longitudinal_analysis_datasets/all_windows.parquet")
    stage3 = stage3[stage3["time_window"].isin(PRIMARY_WINDOWS)]
    merged = stage3.merge(
        changes[["dataset", "patient_id", "time_window", "pathway", "tier", pathway_change_column]],
        on=["dataset", "patient_id", "time_window"], how="inner", validate="many_to_many",
    )
    if pathway_change_column == "delta_z":
        merged = merged.rename(columns={"delta_z_x": "signature_delta_z", "delta_z_y": "pathway_delta_z"})
    else:
        merged = merged.rename(columns={"delta_z": "signature_delta_z", pathway_change_column: "pathway_delta_z"})
    overlap = json.loads((STAGE3_ENV / "analysis_config.json").read_text())["overlap"]
    overlap_lookup = {(r["signature_id"], r["dataset"]): r for r in overlap}
    cohort_rows = []
    for keys, g in merged.groupby(["dataset", "signature_id", "time_window", "pathway", "tier"]):
        g = g[["signature_delta_z", "pathway_delta_z"]].dropna()
        if len(g) < 8 or g.nunique().min() < 2:
            continue
        rho, p = stats.spearmanr(g["signature_delta_z"], g["pathway_delta_z"])
        rho = float(np.clip(rho, -0.999999, 0.999999))
        ov = overlap_lookup[(keys[1], keys[0])]
        cohort_rows.append({
            "dataset": keys[0], "signature_id": keys[1], "time_window": keys[2], "time_label": WINDOW_LABEL[keys[2]],
            "pathway": keys[3], "tier": keys[4], "pathway_method": method, "paired_n": len(g), "spearman_rho": rho, "p_value": float(p),
            "fisher_z": float(np.arctanh(rho)), "fisher_z_se": float(1 / math.sqrt(len(g) - 3)),
            "overlap_status": ov["overlap_status"], "primary_independent_rule": ov["primary_independent_rule"],
        })
    cohort = pd.DataFrame(cohort_rows)
    cohort["fdr_within_dataset_signature_window_tier"] = cohort.groupby(["dataset", "signature_id", "time_window", "tier"], group_keys=False)["p_value"].apply(bh_fdr)
    meta_rows = []
    for independent_only in (False, True):
        source = cohort.copy()
        if independent_only:
            source = source[source["overlap_status"] == "NO_KNOWN_OVERLAP"]
        for keys, g in source.groupby(["signature_id", "time_window", "pathway", "tier"]):
            if len(g) < 2:
                continue
            m = reml_meta(g["fisher_z"].to_numpy(float), g["fisher_z_se"].to_numpy(float))
            meta_rows.append({
                "analysis_set": "INDEPENDENT_ONLY" if independent_only else "ALL_COHORTS",
                "signature_id": keys[0], "time_window": keys[1], "time_label": WINDOW_LABEL[keys[1]],
                "pathway": keys[2], "tier": keys[3], "pathway_method": method, "cohort_count": len(g), "total_paired_n": int(g["paired_n"].sum()),
                "pooled_spearman_rho": float(np.tanh(m["mu"])),
                "ci_low": float(np.tanh(m["lower"])), "ci_high": float(np.tanh(m["upper"])),
                "prediction_low": float(np.tanh(m["prediction_lower"])), "prediction_high": float(np.tanh(m["prediction_upper"])),
                "tau2_fisher_z": m["tau2"], "i2_percent": m["I2"],
                "direction_consistency": float(max((g["spearman_rho"] > 0).mean(), (g["spearman_rho"] < 0).mean())),
                "p_value": m["p_value"],
            })
    meta = pd.DataFrame(meta_rows)
    meta["fdr_within_analysis_window_tier"] = meta.groupby(["analysis_set", "time_window", "tier"], group_keys=False)["p_value"].apply(bh_fdr)
    return cohort, meta


def leave_one_out_pathways(effects: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, g in effects.groupby(["time_window", "pathway", "tier"]):
        if len(g) < 3:
            continue
        for omitted in g["dataset"]:
            h = g[g["dataset"] != omitted]
            m = reml_meta(h["mean_delta_z"].to_numpy(float), h["se_delta_z"].to_numpy(float))
            rows.append({"time_window": keys[0], "pathway": keys[1], "tier": keys[2], "omitted_dataset": omitted,
                         "remaining_cohorts": len(h), "pooled_delta_z": m["mu"], "ci_low": m["lower"], "ci_high": m["upper"]})
    return pd.DataFrame(rows)


def main() -> None:
    SOURCE.mkdir(exist_ok=True)
    PATHDATA.mkdir(exist_ok=True)
    gene_sets = frozen_gene_sets()
    save_pathway_lock(gene_sets)
    score_cache = PATHDATA / "sample_level_pathway_scores.parquet"
    coverage_cache = SOURCE / "04_09_pathway_gene_coverage.csv"
    expression_cache = SOURCE / "04_09_full_transcriptome_expression_audit.csv"
    if score_cache.exists() and coverage_cache.exists() and expression_cache.exists() and not FORCE_REBUILD:
        scores = pd.read_parquet(score_cache)
        coverage = pd.read_csv(coverage_cache)
        expression_audit = pd.read_csv(expression_cache).to_dict("records")
    else:
        maps = platform_maps()
        gencode = parse_gencode()
        coverage_parts = []
        score_parts = []
        expression_audit = []
        for dataset in COHORTS:
            expr = load_full_expression(dataset, maps, gencode)
            expr = expr.loc[~expr.index.duplicated(keep="first")]
            manifest_samples = set(pd.read_csv(INPUT / "metadata/longitudinal_sample_manifest.csv", dtype=str)["sample_id"])
            usable_samples = [s for s in expr.columns if s in manifest_samples]
            expr = expr[usable_samples]
            cov, dataset_scores = score_pathways(expr, gene_sets, dataset)
            coverage_parts.append(cov)
            score_parts.append(dataset_scores)
            expression_audit.append({
                "dataset": dataset, "platform": PLATFORM[dataset], "mapped_gene_count": len(expr),
                "sample_count_in_expression": len(expr.columns), "expression_min": float(np.nanmin(expr.to_numpy())),
                "expression_max": float(np.nanmax(expr.to_numpy())), "constant_gene_count": int((expr.std(axis=1) == 0).sum()),
                "nonfinite_count": int((~np.isfinite(expr.to_numpy())).sum()),
            })
            del expr
        coverage = pd.concat(coverage_parts, ignore_index=True)
        scores = pd.concat(score_parts, ignore_index=True)
    changes = paired_pathway_changes(scores)
    cohort_effects = cohort_pathway_effects(changes)
    meta_effects = meta_pathway_effects(cohort_effects)
    coupling_cohort, coupling_meta = signature_pathway_coupling(changes)
    coupling_cohort_ssgsea, coupling_meta_ssgsea = signature_pathway_coupling(
        changes, pathway_change_column="ssgsea_delta_z", method="SSGSEA_RANK_WALK"
    )
    sensitivity_effects = cohort_effects.rename(columns={
        "mean_delta_z": "mean_singscore_delta_z", "se_delta_z": "se_singscore_delta_z",
        "mean_ssgsea_delta_z": "mean_delta_z", "se_ssgsea_delta_z": "se_delta_z",
    })
    meta_effects_ssgsea = meta_pathway_effects(sensitivity_effects)
    loo = leave_one_out_pathways(cohort_effects)

    gate_support = changes.groupby(["dataset", "time_window"])["patient_id"].nunique().reset_index()
    primary_cov = coverage[coverage["tier"] == "PRESET_PRIMARY"]
    primary_support = primary_cov.groupby("pathway")["compatible"].apply(lambda x: int((x == "YES").sum()))
    gate = {
        "independent_cohorts_with_full_transcriptome": int(pd.Series([r["dataset"] for r in expression_audit]).nunique()),
        "cohorts_supporting_T24": int((gate_support["time_window"] == "T1").sum()),
        "cohorts_supporting_T48": int((gate_support["time_window"] == "T2").sum()),
        "primary_compatible_fraction": float((primary_cov["compatible"] == "YES").mean()),
        "primary_incompatible_count": int((primary_cov["compatible"] != "YES").sum()),
        "minimum_compatible_cohorts_across_primary_pathways": int(primary_support.min()),
        "primary_pathways_supported_by_at_least_four_cohorts": int((primary_support >= 4).sum()),
        "primary_pathway_count": int(len(primary_support)),
        "constant_primary_score_combinations": int(scores[scores["tier"] == "PRESET_PRIMARY"].groupby(["dataset", "pathway"])["singscore"].nunique().le(1).sum()),
    }
    gate["decision"] = "START" if (
        gate["independent_cohorts_with_full_transcriptome"] >= 4
        and max(gate["cohorts_supporting_T24"], gate["cohorts_supporting_T48"]) >= 4
        and gate["minimum_compatible_cohorts_across_primary_pathways"] >= 4
        and gate["constant_primary_score_combinations"] == 0
    ) else "STOP"

    pd.DataFrame(expression_audit).to_csv(SOURCE / "04_09_full_transcriptome_expression_audit.csv", index=False)
    coverage.to_csv(SOURCE / "04_09_pathway_gene_coverage.csv", index=False)
    gate_support.to_csv(SOURCE / "04_09_pathway_time_window_support.csv", index=False)
    (SOURCE / "04_09_pathway_gate_decision.json").write_text(json.dumps(gate, indent=2) + "\n")
    scores.to_parquet(PATHDATA / "sample_level_pathway_scores.parquet", index=False)
    changes.to_parquet(PATHDATA / "patient_level_pathway_changes.parquet", index=False)
    cohort_effects.to_csv(SOURCE / "04_11_cohort_longitudinal_pathway_changes.csv", index=False)
    meta_effects.to_csv(SOURCE / "04_11_meta_longitudinal_pathway_changes.csv", index=False)
    loo.to_csv(SOURCE / "04_11_pathway_leave_one_out.csv", index=False)
    coupling_cohort.to_csv(SOURCE / "04_12_cohort_signature_pathway_coupling.csv", index=False)
    coupling_meta.to_csv(SOURCE / "04_12_meta_signature_pathway_coupling.csv", index=False)
    meta_effects_ssgsea.to_csv(SOURCE / "04_11_meta_longitudinal_pathway_changes_ssgsea.csv", index=False)
    coupling_cohort_ssgsea.to_csv(SOURCE / "04_12_cohort_signature_pathway_coupling_ssgsea.csv", index=False)
    coupling_meta_ssgsea.to_csv(SOURCE / "04_12_meta_signature_pathway_coupling_ssgsea.csv", index=False)
    summary = {
        "gate": gate,
        "expression_audit_rows": len(expression_audit),
        "coverage_rows": len(coverage), "sample_score_rows": len(scores), "patient_change_rows": len(changes),
        "cohort_effect_rows": len(cohort_effects), "meta_effect_rows": len(meta_effects),
        "coupling_cohort_rows": len(coupling_cohort), "coupling_meta_rows": len(coupling_meta),
        "ssgsea_meta_effect_rows": len(meta_effects_ssgsea),
        "ssgsea_coupling_cohort_rows": len(coupling_cohort_ssgsea),
        "ssgsea_coupling_meta_rows": len(coupling_meta_ssgsea),
    }
    (SOURCE / "04_pathway_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
