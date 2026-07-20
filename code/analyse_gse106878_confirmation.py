#!/usr/bin/env python3
"""Post-freeze confirmation analysis for GSE106878.

This script does not alter the six-cohort primary analysis. It verifies the
updated-search eligibility criteria, reconstructs all eight source-locked scores
from the GEO series matrix, and reports paired T24 changes overall and by
randomized treatment arm using the same cohort-baseline scaling convention.
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
MATRIX = ROOT / "downloads" / "GSE106878_series_matrix.txt.gz"
PLATFORM = ROOT / "downloads" / "GPL10295_family.soft.gz"
RELEASE = ROOT.parent
SIGNATURES_PY = RELEASE / "code" / "portable_analysis" / "03_00_09_environment_lock" / "signatures.py"
COEFFICIENTS = RELEASE / "reproduction_project" / "longitudinal_stage3" / "04_source_data" / "04_03_signature_gene_coefficients.csv"
OUTDIR = ROOT / "output"
OUTDIR.mkdir(parents=True, exist_ok=True)

ALIAS_TO_CURRENT = {
    "KIAA1370": "FAM214A",
    "C9ORF103": "IDNK",
    "C9ORF95": "NMRK1",
    "NALP1": "NLRP1",
    "FCMR": "FAIM3",
}


def load_signature_module():
    spec = importlib.util.spec_from_file_location("frozen_signatures", SIGNATURES_PY)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_platform(required: set[str]) -> tuple[dict[str, list[str]], pd.DataFrame]:
    gene_to_probes: dict[str, list[str]] = defaultdict(list)
    evidence: list[dict[str, str]] = []
    inside = False
    with gzip.open(PLATFORM, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                inside = True
                header = next(handle).rstrip("\n").split("\t")
                idx = {name: i for i, name in enumerate(header)}
                continue
            if line.startswith("!platform_table_end"):
                break
            if not inside:
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) <= max(idx["ID"], idx["Symbol"]):
                continue
            probe = row[idx["ID"]].strip()
            raw_symbols = re.split(r"\s*(?:///|;|\|)\s*", row[idx["Symbol"]].strip())
            for raw in raw_symbols:
                if not raw:
                    continue
                canonical = ALIAS_TO_CURRENT.get(raw.upper(), raw.upper())
                if canonical in required:
                    gene_to_probes[canonical].append(probe)
                    evidence.append(
                        {
                            "probe_id": probe,
                            "platform_symbol": raw,
                            "canonical_symbol": canonical,
                            "mapping_basis": "historical alias" if canonical != raw.upper() else "exact symbol",
                        }
                    )
    return gene_to_probes, pd.DataFrame(evidence)


def parse_matrix(selected_probes: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata: dict[str, list[str]] = {}
    rows: list[list[str]] = []
    with gzip.open(MATRIX, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                metadata["title"] = next(csv.reader([line], delimiter="\t"))[1:]
            elif line.startswith("!Sample_geo_accession"):
                metadata["geo_accession"] = next(csv.reader([line], delimiter="\t"))[1:]
            elif line.startswith("!series_matrix_table_begin"):
                header = next(csv.reader([next(handle)], delimiter="\t"))
                for row in csv.reader(handle, delimiter="\t"):
                    if row and row[0].startswith("!series_matrix_table_end"):
                        break
                    if row and row[0] in selected_probes:
                        rows.append(row)
                break
    matrix = pd.DataFrame(rows, columns=header).set_index("ID_REF").astype(float)
    sample = pd.DataFrame(metadata)
    parsed = sample["title"].str.extract(
        r"^(?P<patient>C\d+)\s+(?P<timepoint>Pre|Post\(24h\))\s+(?P<treatment>placebo|hydrocortisone)$"
    )
    sample = pd.concat([sample, parsed], axis=1)
    if sample[["patient", "timepoint", "treatment"]].isna().any().any():
        raise ValueError("Could not parse all GSE106878 sample titles")
    return matrix, sample


def bootstrap_mean(values: np.ndarray, rng: np.random.Generator, n_boot: int = 2000) -> tuple[float, float, float]:
    draws = np.empty(n_boot, dtype=float)
    n = len(values)
    for i in range(n_boot):
        draws[i] = rng.choice(values, size=n, replace=True).mean()
    return float(draws.std(ddof=1)), float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def main() -> None:
    signatures = load_signature_module()
    coefficients = pd.read_csv(COEFFICIENTS)
    required = set(coefficients["gene"].astype(str).str.upper())
    gene_to_probes, mapping = parse_platform(required)
    missing = sorted(required - set(gene_to_probes))
    if missing:
        raise RuntimeError("Incomplete 74-gene coverage: " + ", ".join(missing))

    probe_matrix, sample = parse_matrix({p for ps in gene_to_probes.values() for p in ps})
    missing_matrix = sorted(g for g, ps in gene_to_probes.items() if not any(p in probe_matrix.index for p in ps))
    if missing_matrix:
        raise RuntimeError("Mapped probes absent from series matrix: " + ", ".join(missing_matrix))

    # Frozen rule: mean across unambiguous probes after symbol resolution.
    gene_matrix = pd.DataFrame(index=sorted(required), columns=probe_matrix.columns, dtype=float)
    for gene in gene_matrix.index:
        probes = [p for p in gene_to_probes[gene] if p in probe_matrix.index]
        gene_matrix.loc[gene] = probe_matrix.loc[probes].mean(axis=0)

    if (gene_matrix <= 0).any().any():
        raise RuntimeError("Non-positive expression prevents geometric-mean score reconstruction")

    functions = signatures.STAGE3_FUNCTIONS
    directions = signatures.SCORE_DIRECTION
    records: list[dict[str, object]] = []
    scores: dict[str, dict[str, float]] = defaultdict(dict)
    for sig, func in functions.items():
        genes = signatures.STAGE3_REQUIRED_GENES[sig]
        for gsm in gene_matrix.columns:
            expr = {gene: float(gene_matrix.at[gene, gsm]) for gene in genes}
            scores[sig][gsm] = directions[sig] * float(func(expr))

    rng = np.random.default_rng(20260720)
    for sig in functions:
        baseline_gsms = sample.loc[sample["timepoint"].eq("Pre"), "geo_accession"].tolist()
        baseline_values = np.array([scores[sig][gsm] for gsm in baseline_gsms], dtype=float)
        baseline_sd = float(baseline_values.std(ddof=1))
        if not np.isfinite(baseline_sd) or baseline_sd <= 0:
            raise RuntimeError(f"{sig}: invalid baseline SD")

        paired: list[dict[str, object]] = []
        for patient, group in sample.groupby("patient", sort=False):
            if set(group["timepoint"]) != {"Pre", "Post(24h)"} or len(group) != 2:
                raise RuntimeError(f"{patient}: pair is incomplete or duplicated")
            pre = group.loc[group["timepoint"].eq("Pre")].iloc[0]
            post = group.loc[group["timepoint"].eq("Post(24h)")].iloc[0]
            if pre["treatment"] != post["treatment"]:
                raise RuntimeError(f"{patient}: treatment label changed within pair")
            paired.append(
                {
                    "patient": patient,
                    "treatment": pre["treatment"],
                    "delta_z": (scores[sig][post["geo_accession"]] - scores[sig][pre["geo_accession"]]) / baseline_sd,
                }
            )
        pair_df = pd.DataFrame(paired)
        for stratum, frame in [("all", pair_df)] + list(pair_df.groupby("treatment", sort=True)):
            values = frame["delta_z"].to_numpy(float)
            se, lo, hi = bootstrap_mean(values, rng)
            records.append(
                {
                    "dataset": "GSE106878",
                    "time_window": "T24",
                    "signature_id": sig,
                    "stratum": stratum,
                    "n_pairs": len(values),
                    "baseline_sd": baseline_sd,
                    "mean_delta_z": float(values.mean()),
                    "bootstrap_se": se,
                    "bootstrap_ci_low": lo,
                    "bootstrap_ci_high": hi,
                    "bootstrap_replicates": 2000,
                    "seed": 20260720,
                }
            )

    pd.DataFrame(records).to_csv(OUTDIR / "GSE106878_postfreeze_confirmation_effects.csv", index=False)
    mapping.sort_values(["canonical_symbol", "probe_id"]).to_csv(
        OUTDIR / "GSE106878_platform_feature_map.csv", index=False
    )
    gene_matrix.to_csv(OUTDIR / "GSE106878_signature_gene_expression.csv")
    sample.to_csv(OUTDIR / "GSE106878_sample_registry.csv", index=False)
    pd.DataFrame(
        [
            {
                "dataset": "GSE106878",
                "decision": "INCLUDE_POST_FREEZE_CONFIRMATION",
                "baseline_cases": sample["patient"].nunique(),
                "paired_T24": sample.loc[sample["timepoint"].eq("Post(24h)"), "patient"].nunique(),
                "required_genes": len(required),
                "covered_genes": len(gene_to_probes),
                "public_matrix": "GEO series matrix",
                "platform": "GPL10295 Illumina human-6 v2.0 expression beadchip",
                "normalization": "GEO submitter: log2 transformation followed by quantile normalization",
                "role": "post-freeze confirmation only; randomized treatment arms retained and reported",
            }
        ]
    ).to_csv(OUTDIR / "GSE106878_eligibility_summary.csv", index=False)


if __name__ == "__main__":
    main()
