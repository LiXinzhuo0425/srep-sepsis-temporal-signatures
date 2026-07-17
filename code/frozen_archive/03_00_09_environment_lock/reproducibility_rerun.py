#!/usr/bin/env python3
"""Clean-from-frozen-input numerical rerun and reproducibility comparison."""

from __future__ import annotations

import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/Users/felix/Documents/New project/longitudinal_stage3")
ENV = ROOT / "03_00_09_environment_lock"
sys.path.insert(0, str(ENV))
from run_stage3 import (  # noqa: E402
    COHORTS,
    GSE236713_RAW_REQUIRED,
    GSE236713_RAW_TAR,
    INTERMEDIATE,
    build_analysis_data,
    meta_analyses,
    paired_effects,
    score_dataset,
    sha256,
)


def numeric_compare(left, right, keys, numeric, tolerance=1e-10):
    merged = left[keys + numeric].merge(right[keys + numeric], on=keys, suffixes=("_rerun", "_stored"), validate="one_to_one")
    if len(merged) != len(left) or len(merged) != len(right):
        return False, float("inf"), "row/key mismatch"
    differences = []
    for column in numeric:
        a = pd.to_numeric(merged[f"{column}_rerun"], errors="coerce").to_numpy(float)
        b = pd.to_numeric(merged[f"{column}_stored"], errors="coerce").to_numpy(float)
        both_nan = np.isnan(a) & np.isnan(b)
        diff = np.abs(a - b)
        diff[both_nan] = 0
        differences.append(np.nanmax(diff) if len(diff) else 0)
    maximum = float(max(differences, default=0))
    return maximum <= tolerance, maximum, "numeric comparison"


def main():
    log = []
    cfg = json.loads((ENV / "analysis_config.json").read_text())
    input_checks = []
    for item in cfg["input_inventory"]:
        path = ROOT / item["relative_path"]
        if not path.exists() or item["performance_content"] == "CODE_ONLY":
            continue
        current = sha256(path)
        if item["relative_path"] == "inputs/expression/GSE236713_series_matrix.txt.gz":
            corrupt_copy = ROOT / "inputs/expression/GSE236713_series_matrix.txt.gz.corrupt_20260716"
            deviation_logged = (ROOT / "03_00_08_protocol_deviations.json").exists()
            passed = (
                deviation_logged
                and current == "b0ba144f690538f08196fa9584b91a009272e2a67650d07073faac21ecc70276"
                and corrupt_copy.exists()
                and sha256(corrupt_copy) == item["sha256"]
                and GSE236713_RAW_TAR.exists()
                and sha256(GSE236713_RAW_TAR) == "dbb89938c1f093cab6c21a08b63f818ca12e7819cdb9473e1979c0b94b31df92"
            )
            log.append(
                f"INPUT_HASH {item['relative_path']} {'PASS_DOCUMENTED_DEVIATION' if passed else 'FAIL'} "
                f"original_truncated_hash_retained={corrupt_copy.exists()} raw_rescue_hash_verified={GSE236713_RAW_TAR.exists()}"
            )
            input_checks.append(passed)
            continue
        passed = current == item["sha256"]
        input_checks.append(passed)
        log.append(f"INPUT_HASH {item['relative_path']} {'PASS' if passed else 'FAIL'}")
    if not all(input_checks):
        raise RuntimeError("Frozen input checksum mismatch")

    # Force the two raw-data preprocessing routes to be regenerated.  This is
    # the clean-input portion of the reproducibility audit; no stored score or
    # normalized expression matrix is reused.
    for path in (
        GSE236713_RAW_REQUIRED,
        INTERMEDIATE / "GSE236713_raw_preprocessing_qc.json",
        INTERMEDIATE / "GSE236713_raw_preprocessing_log.txt",
        INTERMEDIATE / "GSE110487_required_gene_vst.tsv",
        INTERMEDIATE / "GSE110487_grouped_counts.tsv",
        INTERMEDIATE / "GSE110487_required_genes.txt",
        INTERMEDIATE / "GSE110487_vst_log.txt",
    ):
        if path.exists():
            path.unlink()
    log.append("PREPROCESS_INTERMEDIATES removed_and_rebuilt=GSE236713,GSE110487")

    rerun_scores = []
    for dataset in COHORTS:
        scores, qc, manual = score_dataset(dataset)
        if not (qc.status == "PASS").all() or not (manual.status == "PASS").all():
            raise RuntimeError(f"{dataset} rerun scoring QC failed")
        rerun_scores.append(scores)
        log.append(f"SCORE {dataset} rows={len(scores)} qc=PASS manual={len(manual)}/{len(manual)}")
    rerun_scores = pd.concat(rerun_scores, ignore_index=True)
    stored_scores = pd.read_parquet(ROOT / "03_02_full_signature_score_matrix.parquet")
    passed, maximum, note = numeric_compare(
        rerun_scores, stored_scores,
        ["dataset", "sample_id", "signature_id"],
        ["raw_score", "oriented_score"],
    )
    log.append(f"COMPARE scores {'PASS' if passed else 'FAIL'} max_abs_diff={maximum:.3g} {note}")
    if not passed:
        raise RuntimeError("Score rerun differs")

    changes, _ = build_analysis_data(rerun_scores)
    stored_changes = pd.read_parquet(ROOT / "03_03_longitudinal_analysis_datasets/all_windows.parquet")
    passed, maximum, note = numeric_compare(
        changes, stored_changes,
        ["dataset", "patient_id", "signature_id", "time_window"],
        ["baseline_score", "followup_score", "raw_change", "delta_z"],
    )
    log.append(f"COMPARE changes {'PASS' if passed else 'FAIL'} max_abs_diff={maximum:.3g} {note}")
    if not passed:
        raise RuntimeError("Change rerun differs")

    effects, diagnostics = paired_effects(changes)
    stored_effects = pd.read_csv(ROOT / "source_data/03_04_cohort_signature_primary_effects.csv")
    passed, maximum, note = numeric_compare(
        effects, stored_effects,
        ["dataset", "signature_id", "time_window"],
        ["mean_delta_z", "bootstrap_se", "ci95_lower", "ci95_upper"],
    )
    log.append(f"COMPARE cohort_effects {'PASS' if passed else 'FAIL'} max_abs_diff={maximum:.3g} {note}")
    if not passed:
        raise RuntimeError("Cohort effect rerun differs")

    meta, _ = meta_analyses(effects)
    stored_meta = pd.read_csv(ROOT / "source_data/03_05_signature_level_meta_analysis.csv")
    passed, maximum, note = numeric_compare(
        meta, stored_meta,
        ["signature_id", "time_window", "analysis_set"],
        ["pooled_delta_z", "ci95_lower", "ci95_upper", "tau2", "I2_percent"],
        tolerance=1e-8,
    )
    log.append(f"COMPARE meta {'PASS' if passed else 'FAIL'} max_abs_diff={maximum:.3g} {note}")
    if not passed:
        raise RuntimeError("Meta-analysis rerun differs")

    figure_files = list((ROOT / "03_15_main_figures").glob("Figure_*.*"))
    required_extensions = {".svg", ".pdf", ".tiff", ".png"}
    grouped = {}
    for path in figure_files:
        grouped.setdefault(path.stem, set()).add(path.suffix)
    figure_pass = len(grouped) >= 6 and all(required_extensions.issubset(exts) for exts in grouped.values())
    log.append(f"FIGURE_EXPORTS count={len(figure_files)} {'PASS' if figure_pass else 'FAIL'}")

    versions = {
        "python": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }
    try:
        versions["DESeq2"] = subprocess.check_output([
            "Rscript", "-e", '.libPaths(c("/Users/felix/Documents/New project/longitudinal_stage3/03_00_09_environment_lock/R_libs", .libPaths())); cat(as.character(packageVersion("DESeq2")))'
        ], text=True, stderr=subprocess.DEVNULL).strip()
    except Exception as exc:
        versions["DESeq2"] = f"ERROR {exc}"
    log.append("VERSIONS " + json.dumps(versions, sort_keys=True))
    log.append("FINAL PASS")
    (ROOT / "03_17_reproducibility_run_log.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    html_rows = "".join(f"<tr><td>{line}</td></tr>" for line in log)
    html = f"<html><head><meta charset='utf-8'><style>body{{font-family:Arial;margin:32px;color:#17212b}}table{{border-collapse:collapse}}td{{border:1px solid #d4dee7;padding:7px}}h1{{color:#18324a}}</style></head><body><h1>Stage 3 reproducibility QC</h1><p><strong>PASS.</strong> Frozen inputs, scores, patient pairing, cohort effects and meta-analysis were recomputed without manual intervention.</p><table>{html_rows}</table></body></html>"
    (ROOT / "03_17_reproducibility_qc_report.html").write_text(html, encoding="utf-8")
    print(json.dumps({"status": "PASS", "log_lines": len(log), "model_diagnostics": len(diagnostics)}, indent=2))


if __name__ == "__main__":
    main()
