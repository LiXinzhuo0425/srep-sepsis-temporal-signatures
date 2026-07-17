#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


RELEASE = Path(__file__).resolve().parents[1]
PROJECT = RELEASE / "reproduction_project"
SCRIPT = PROJECT / "stage6_submission/06_02_revision/08_code/stage6_figures.py"
OUT = PROJECT / "stage6_submission/06_02_revision/06_figures/main"
REFERENCE = RELEASE / "reference_outputs/main_figures"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


env = os.environ.copy()
env["SEPSIS_SIGNATURE_PROJECT_ROOT"] = str(PROJECT)
completed = subprocess.run([sys.executable, str(SCRIPT)], cwd=PROJECT, env=env, text=True, capture_output=True)
checks = []
for ref in sorted(REFERENCE.glob("Figure_*.png")):
    produced = OUT / ref.name
    checks.append({
        "figure": ref.stem,
        "reference_exists": ref.exists(),
        "produced_exists": produced.exists(),
        "reference_sha256": sha256(ref) if ref.exists() else None,
        "produced_sha256": sha256(produced) if produced.exists() else None,
    })
for item in checks:
    item["status"] = "PASS" if item["reference_sha256"] == item["produced_sha256"] else "FAIL"
result = {
    "release_id": "srep-sepsis-temporal-signatures-v1.0.0",
    "command_exit_code": completed.returncode,
    "stdout": completed.stdout,
    "stderr": completed.stderr,
    "figure_checks": checks,
    "status": "PASS" if completed.returncode == 0 and len(checks) == 6 and all(x["status"] == "PASS" for x in checks) else "FAIL",
}
(PROJECT / "reproduction_result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["status"] == "PASS" else 1)
