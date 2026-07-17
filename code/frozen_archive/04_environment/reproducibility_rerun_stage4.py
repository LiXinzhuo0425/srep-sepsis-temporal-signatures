#!/usr/bin/env python3
"""Clean Stage 4 rerun with semantic numerical comparison to the frozen outputs."""
from __future__ import annotations

import hashlib
import html
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path("/Users/felix/Documents/New project/longitudinal_stage3")
ENV = ROOT / "04_environment"
SOURCE = ROOT / "04_source_data"
LOG = ROOT / "04_23_reproducibility_run_log.txt"
REPORT = ROOT / "04_23_reproducibility_qc_report.html"
NODE = "/Users/felix/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node"
BUNDLED_PY = "/Users/felix/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3"

CANONICAL = [
    ROOT / "04_04_gene_level_longitudinal_dataset.parquet",
    SOURCE / "04_05_decomposition_reconstruction_qc.csv",
    SOURCE / "04_06_cohort_level_gene_contributions.csv",
    SOURCE / "04_07_gene_contribution_meta_analysis.csv",
    SOURCE / "04_08_signature_drift_architecture.csv",
    SOURCE / "04_09_pathway_gene_coverage.csv",
    ROOT / "04_pathway_data/patient_level_pathway_changes.parquet",
    SOURCE / "04_11_meta_longitudinal_pathway_changes.csv",
    SOURCE / "04_12_meta_signature_pathway_coupling.csv",
    SOURCE / "04_13_cell_source_annotation.csv",
    SOURCE / "04_15_clinical_anchor_availability_audit.csv",
    SOURCE / "04_20_integrated_signature_interpretation_matrix.csv",
]


def semantic_bytes(file: Path) -> bytes:
    if file.suffix == ".parquet":
        frame = pd.read_parquet(file)
    else:
        frame = pd.read_csv(file)
    frame = frame.reindex(sorted(frame.columns), axis=1)
    for col in frame.select_dtypes(include="number").columns:
        frame[col] = frame[col].round(12)
    frame = frame.fillna("")
    if len(frame):
        frame = frame.sort_values(list(frame.columns), key=lambda s: s.astype(str), kind="mergesort").reset_index(drop=True)
    return frame.to_csv(index=False, lineterminator="\n").encode("utf-8")


def digest(file: Path) -> str:
    return hashlib.sha256(semantic_bytes(file)).hexdigest()


def run_step(name: str, command: list[str], env: dict | None, records: list[dict], log_handle) -> None:
    started = time.time()
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{stamp}] START {name}: {' '.join(command)}"
    print(line, flush=True)
    log_handle.write(line + "\n")
    completed = subprocess.run(command, cwd=ROOT, env=env, text=True, capture_output=True)
    elapsed = time.time() - started
    log_handle.write(completed.stdout)
    log_handle.write(completed.stderr)
    log_handle.write(f"\nEND {name}: exit={completed.returncode}; elapsed={elapsed:.1f}s\n\n")
    log_handle.flush()
    records.append({"step": name, "exit_code": completed.returncode, "elapsed_seconds": round(elapsed, 1)})
    if completed.returncode != 0:
        raise RuntimeError(f"{name} failed with exit code {completed.returncode}")
    print(f"DONE {name} ({elapsed:.1f}s)", flush=True)


before = {str(p.relative_to(ROOT)): digest(p) for p in CANONICAL}
records: list[dict] = []
with LOG.open("w", encoding="utf-8") as log_handle:
    log_handle.write("Stage 4 reproducibility rerun\n")
    log_handle.write("Date: 2026-07-16\n")
    log_handle.write(f"Python: {sys.version}\n")
    log_handle.write("Random seed: 20260716\n")
    log_handle.write("Force-rebuild pathway caches: true\n\n")
    run_step("Exact gene decomposition", [sys.executable, str(ENV / "stage4_decomposition.py")], None, records, log_handle)
    force_env = os.environ.copy()
    force_env["STAGE4_FORCE_REBUILD"] = "1"
    run_step("Full-transcriptome pathway analysis", [sys.executable, str(ENV / "stage4_pathway.py")], force_env, records, log_handle)
    run_step("Cell-source annotation", [sys.executable, str(ENV / "stage4_cell_annotation.py")], None, records, log_handle)
    run_step("Clinical-anchor gate audit", [sys.executable, str(ENV / "stage4_clinical_audit.py")], None, records, log_handle)
    run_step("Integrated interpretation", [sys.executable, str(ENV / "stage4_integrate.py")], None, records, log_handle)
    run_step("Publication figure generation", [sys.executable, str(ENV / "stage4_figures.py")], None, records, log_handle)
    run_step("Independent verification", [sys.executable, str(ENV / "verify_stage4.py")], None, records, log_handle)
    run_step("Method documents", [BUNDLED_PY, str(ENV / "make_stage4_documents_en.py")], None, records, log_handle)
    run_step("Artifact-tool workbooks", [NODE, str(ENV / "build_stage4_workbooks.mjs")], None, records, log_handle)

after = {str(p.relative_to(ROOT)): digest(p) for p in CANONICAL}
comparisons = []
for name in before:
    comparisons.append({"file": name, "before": before[name], "after": after[name], "status": "PASS" if before[name] == after[name] else "FAIL"})
all_steps = all(r["exit_code"] == 0 for r in records)
all_equal = all(r["status"] == "PASS" for r in comparisons)
overall = "PASS" if all_steps and all_equal else "FAIL"

qa_csv = SOURCE / "04_23_reproducibility_comparison.csv"
pd.DataFrame(comparisons).to_csv(qa_csv, index=False)

step_rows = "".join(
    f"<tr><td>{html.escape(r['step'])}</td><td>{r['exit_code']}</td><td>{r['elapsed_seconds']:.1f}</td><td class='pass'>PASS</td></tr>"
    for r in records
)
file_rows = "".join(
    f"<tr><td>{html.escape(r['file'])}</td><td><code>{r['before'][:16]}</code></td><td><code>{r['after'][:16]}</code></td><td class='{r['status'].lower()}'>{r['status']}</td></tr>"
    for r in comparisons
)
REPORT.write_text(
    f"""<!doctype html><html><head><meta charset='utf-8'><title>Stage 4 reproducibility QC</title>
<style>body{{font-family:Arial,sans-serif;margin:38px;color:#1f2937}}h1{{color:#17365D}}.banner{{padding:14px 18px;background:{'#DFF0D8' if overall=='PASS' else '#F8D7DA'};font-weight:700;border-left:6px solid {'#2E7D32' if overall=='PASS' else '#B64342'}}}table{{border-collapse:collapse;width:100%;margin:16px 0 28px}}th{{background:#17365D;color:white;text-align:left}}th,td{{padding:8px;border:1px solid #c8d1dc;font-size:13px}}.pass{{color:#2E7D32;font-weight:700}}.fail{{color:#B64342;font-weight:700}}code{{font-size:11px}}</style></head><body>
<h1>Stage 4 reproducibility quality-control report</h1><div class='banner'>Overall status: {overall}</div>
<p>The rerun started from the frozen Stage 3 inputs, forced regeneration of cohort-level full-transcriptome matrices and pathway scores, reran every Stage 4 analytical module, and rebuilt figures, documents and workbooks. Numerical equality uses canonicalized tables with numeric values rounded to 12 decimal places.</p>
<h2>Execution steps</h2><table><tr><th>Step</th><th>Exit code</th><th>Elapsed seconds</th><th>Status</th></tr>{step_rows}</table>
<h2>Canonical numerical comparisons</h2><table><tr><th>Output</th><th>Before hash</th><th>After hash</th><th>Status</th></tr>{file_rows}</table>
<p>Random seed: 20260716. Full command output and software version are recorded in 04_23_reproducibility_run_log.txt.</p></body></html>""",
    encoding="utf-8",
)

# Move artifact-tool inspect sidecars into the QA directory without changing the workbooks.
sidecar_dir = ROOT / "04_qa/workbooks/export_inspect"
sidecar_dir.mkdir(parents=True, exist_ok=True)
for file in ROOT.glob("*.xlsx.inspect.ndjson"):
    shutil.move(str(file), str(sidecar_dir / file.name))

print(json.dumps({"overall_status": overall, "steps": len(records), "canonical_outputs": len(comparisons), "matches": sum(r['status'] == 'PASS' for r in comparisons)}, indent=2), flush=True)
if overall != "PASS":
    raise SystemExit(1)
