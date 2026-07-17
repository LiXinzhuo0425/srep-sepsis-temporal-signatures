#!/usr/bin/env python3
"""Result-blind rescue preprocessing for GSE236713 Agilent raw files.

The deposited series matrix is 75th-percentile normalized and then baseline
transformed to the global median.  That last centering step creates negative
values and destroys the positive absolute scale required by several frozen
ratio/geometric-mean signatures.  This script reconstructs the pre-baseline
scale from the public Feature Extraction files without using phenotype,
timepoint or outcome information:

1. read gProcessedSignal for non-control features;
2. compute each array's 75th percentile across non-control features;
3. scale every array to the median 75th percentile across all arrays;
4. log2 transform the positive, scaled processed signal;
5. retain only the probes locked before unblinding.

No sample label or signature score is read by this script.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import math
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd


GSM_RE = re.compile(r"(GSM\d+)")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_agilent_member(handle, required_probes: set[str]) -> tuple[float, dict[str, list[float]], int]:
    q75_values: list[float] = []
    retained: dict[str, list[float]] = {}
    feature_header: dict[str, int] | None = None
    with gzip.GzipFile(fileobj=handle, mode="rb") as zipped:
        text = io.TextIOWrapper(zipped, encoding="utf-8", errors="replace", newline="")
        for raw in text:
            fields = raw.rstrip("\r\n").split("\t")
            if fields and fields[0] == "FEATURES":
                feature_header = {name: index for index, name in enumerate(fields)}
                needed = {"ControlType", "ProbeName", "gProcessedSignal"}
                if not needed.issubset(feature_header):
                    raise RuntimeError(f"Agilent FEATURES header missing {sorted(needed - set(feature_header))}")
                continue
            if not feature_header or not fields or fields[0] != "DATA":
                continue
            try:
                control = int(fields[feature_header["ControlType"]])
                signal = float(fields[feature_header["gProcessedSignal"]])
            except (ValueError, IndexError) as exc:
                raise RuntimeError("Malformed Agilent feature row") from exc
            if control != 0 or not math.isfinite(signal) or signal <= 0:
                continue
            q75_values.append(signal)
            probe = fields[feature_header["ProbeName"]]
            if probe in required_probes:
                retained.setdefault(probe, []).append(signal)
    if not q75_values:
        raise RuntimeError("No positive non-control gProcessedSignal values")
    return float(np.quantile(np.asarray(q75_values, dtype=float), 0.75)), retained, len(q75_values)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("raw_tar", type=Path)
    parser.add_argument("feature_map", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("qc_json", type=Path)
    args = parser.parse_args()

    fmap = pd.read_csv(args.feature_map, dtype=str, keep_default_na=False)
    fmap = fmap[(fmap["dataset"] == "GSE236713") & (fmap["included_in_lock"] == "YES")]
    required_probes = set(fmap["feature_id"])
    if not required_probes:
        raise RuntimeError("No locked GSE236713 probes")

    arrays: dict[str, dict] = {}
    member_count = 0
    with tarfile.open(args.raw_tar, mode="r:") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".txt.gz"):
                continue
            match = GSM_RE.search(member.name)
            if not match:
                continue
            sample_id = match.group(1)
            extracted = archive.extractfile(member)
            if extracted is None:
                raise RuntimeError(f"Cannot extract {member.name}")
            q75, retained, noncontrol_n = read_agilent_member(extracted, required_probes)
            arrays[sample_id] = {"q75": q75, "retained": retained, "noncontrol_n": noncontrol_n}
            member_count += 1

    if member_count != 447 or len(arrays) != 447:
        raise RuntimeError(f"Expected 447 unique raw arrays, found members={member_count}, unique={len(arrays)}")
    target_q75 = float(np.median([item["q75"] for item in arrays.values()]))
    if not math.isfinite(target_q75) or target_q75 <= 0:
        raise RuntimeError("Invalid cohort target 75th percentile")

    rows: list[dict] = []
    missing: dict[str, list[str]] = {}
    for sample_id, item in sorted(arrays.items()):
        scale = target_q75 / item["q75"]
        missing_probes = sorted(required_probes - set(item["retained"]))
        if missing_probes:
            missing[sample_id] = missing_probes
        for probe, values in item["retained"].items():
            # Replicate features with the same ProbeName are averaged after
            # log transformation, matching the frozen multi-probe mean rule.
            normalized = [math.log2(value * scale) for value in values]
            rows.append({
                "sample_id": sample_id,
                "feature_id": probe,
                "expression": float(np.mean(normalized)),
                "raw_replicates": len(values),
            })
    if missing:
        example = next(iter(missing.items()))
        raise RuntimeError(f"Locked probes missing from raw data; example {example}")

    result = pd.DataFrame(rows)
    if result.duplicated(["sample_id", "feature_id"]).any():
        raise RuntimeError("Duplicate sample/probe rows after replicate aggregation")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(args.output, index=False)
    qc = {
        "status": "PASS",
        "raw_tar": str(args.raw_tar),
        "raw_tar_sha256": sha256(args.raw_tar),
        "array_count": len(arrays),
        "locked_probe_count": len(required_probes),
        "output_rows": len(result),
        "target_q75": target_q75,
        "sample_q75_min": min(item["q75"] for item in arrays.values()),
        "sample_q75_max": max(item["q75"] for item in arrays.values()),
        "expression_min": float(result["expression"].min()),
        "expression_max": float(result["expression"].max()),
        "phenotype_or_outcome_used": False,
        "method": "gProcessedSignal; non-control per-array q75; scale to cohort median q75; log2; locked probes only",
    }
    args.qc_json.write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
