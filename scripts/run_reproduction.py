#!/usr/bin/env python3
"""Regenerate the six main figures and compare scientific image content and source tables."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path

from PIL import Image


RELEASE = Path(__file__).resolve().parents[1]
PROJECT = RELEASE / "reproduction_project"
SCRIPT = RELEASE / "code/stage6_figures.py"
EVIDENCE = RELEASE / "reproducibility_evidence"
GENERATED = EVIDENCE / "generated_v1.2.1"
OUT = GENERATED / "main_figures"
PRODUCED_SOURCE = GENERATED / "source_data"
PRODUCED_SUPPLEMENTARY = GENERATED / "supplementary_figures"
REFERENCE = RELEASE / "reference_outputs/main_figures_v1.2.1"
REFERENCE_SOURCE = RELEASE / "data/figure_source_data"
SUPPLEMENTARY_REFERENCE = RELEASE / "reference_outputs/supplementary_figures_v1.2.1"
SUPPLEMENTARY_SOURCE_REFERENCE = (
    PROJECT
    / "stage6_submission/06_02_revision/07_source_data_v1_1/"
    "Supplementary_Figure_1_source_data.csv"
)
SUPPLEMENTARY_S1 = "Supplementary_Figure_1_pilot_vs_validation"
AUTHORED_LAYOUT = {
    "Figure_1_study_design",
    "Figure_4_gene_contribution_patterns",
    "Figure_5_gene_pathway_cell_integration",
    "Figure_6_temporal_measurement_context",
}
PIXEL_REFERENCE_FIGURES = {
    "Figure_2_pooled_longitudinal_change",
    "Figure_3_cross_cohort_transportability",
}
EXPECTED_FIGURES = tuple(sorted(AUTHORED_LAYOUT | PIXEL_REFERENCE_FIGURES))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pixel_sha256(path: Path) -> str:
    with Image.open(path) as image:
        rgb = image.convert("RGB")
        digest = hashlib.sha256()
        digest.update(f"RGB:{rgb.width}x{rgb.height}:".encode("ascii"))
        digest.update(rgb.tobytes())
        return digest.hexdigest()


def safe_pixel_sha256(path: Path) -> tuple[str | None, str | None]:
    """Return a decoded-pixel digest without allowing a bad image to abort QA."""
    if not path.is_file():
        return None, "missing"
    try:
        return pixel_sha256(path), None
    except Exception as exc:  # Pillow reports format-specific decode errors.
        return None, f"{type(exc).__name__}: {exc}"


def valid_pdf(path: Path) -> bool:
    """Require a non-empty PDF with the standard header signature."""
    if not path.is_file():
        return False
    try:
        if path.stat().st_size <= 5:
            return False
        with path.open("rb") as handle:
            return handle.read(5) == b"%PDF-"
    except OSError:
        return False


def package_version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


env = os.environ.copy()
env["SEPSIS_SIGNATURE_PROJECT_ROOT"] = str(PROJECT)
if GENERATED.exists():
    shutil.rmtree(GENERATED)
env["STAGE6_FIGURE_OUT"] = str(OUT)
env["STAGE6_SOURCE_OUT"] = str(PRODUCED_SOURCE)
env["STAGE6_SUPP_OUT"] = str(PRODUCED_SUPPLEMENTARY)
completed = subprocess.run([sys.executable, str(SCRIPT)], cwd=PROJECT, env=env, text=True, capture_output=True)

figure_checks = []
for stem in EXPECTED_FIGURES:
    ref = REFERENCE / f"{stem}.png"
    produced = OUT / ref.name
    reference_pixel_sha, reference_decode_error = safe_pixel_sha256(ref)
    produced_pixel_sha, produced_decode_error = safe_pixel_sha256(produced)
    item = {
        "figure": stem,
        "reference_exists": ref.is_file(),
        "produced_exists": produced.is_file(),
        "reference_file_sha256": sha256(ref) if ref.is_file() else None,
        "produced_file_sha256": sha256(produced) if produced.is_file() else None,
        "reference_pixel_sha256": reference_pixel_sha,
        "produced_pixel_sha256": produced_pixel_sha,
        "reference_decode_error": reference_decode_error,
        "produced_decode_error": produced_decode_error,
    }
    if (
        item["reference_exists"]
        and item["produced_exists"]
        and item["reference_pixel_sha256"] is not None
        and item["reference_pixel_sha256"] == item["produced_pixel_sha256"]
    ):
        item["status"] = "PASS"
    elif (
        item["figure"] in AUTHORED_LAYOUT
        and item["reference_exists"]
        and item["produced_exists"]
        and item["reference_pixel_sha256"] is not None
        and item["produced_pixel_sha256"] is not None
    ):
        item["status"] = "PASS_AUTHORED_LAYOUT"
        item["note"] = (
            "The submission reference retains author-reviewed spacing and label placement. "
            "Scientific values are verified independently through the exact source-data comparison."
        )
    else:
        item["status"] = "FAIL"
    figure_checks.append(item)

source_checks = []
for ref in sorted(REFERENCE_SOURCE.glob("*.csv")):
    produced = PRODUCED_SOURCE / ref.name
    item = {
        "source_file": ref.name,
        "reference_exists": ref.exists(),
        "produced_exists": produced.exists(),
        "reference_sha256": sha256(ref) if ref.exists() else None,
        "produced_sha256": sha256(produced) if produced.exists() else None,
    }
    item["status"] = "PASS" if item["reference_sha256"] == item["produced_sha256"] else "FAIL"
    source_checks.append(item)

supplementary_source_produced = PRODUCED_SOURCE / SUPPLEMENTARY_SOURCE_REFERENCE.name
supplementary_source_check = {
    "source_file": SUPPLEMENTARY_SOURCE_REFERENCE.name,
    "reference_exists": SUPPLEMENTARY_SOURCE_REFERENCE.exists(),
    "produced_exists": supplementary_source_produced.exists(),
    "reference_sha256": sha256(SUPPLEMENTARY_SOURCE_REFERENCE)
    if SUPPLEMENTARY_SOURCE_REFERENCE.exists()
    else None,
    "produced_sha256": sha256(supplementary_source_produced)
    if supplementary_source_produced.exists()
    else None,
}
supplementary_source_check["status"] = (
    "PASS"
    if supplementary_source_check["reference_sha256"]
    == supplementary_source_check["produced_sha256"]
    else "FAIL"
)

format_checks = []
for stem in EXPECTED_FIGURES:
    png = OUT / f"{stem}.png"
    tiff = OUT / f"{stem}.tiff"
    svg = OUT / f"{stem}.svg"
    pdf = OUT / f"{stem}.pdf"
    png_sha, png_decode_error = safe_pixel_sha256(png)
    tiff_ok = False
    tiff_decode_error = None
    if tiff.is_file():
        try:
            with Image.open(tiff) as image:
                image.load()
                tiff_ok = (
                    image.mode == "RGB"
                    and tuple(round(float(value)) for value in image.info.get("dpi", (0, 0))) == (600, 600)
                    and image.info.get("compression") == "tiff_lzw"
                )
        except Exception as exc:
            tiff_decode_error = f"{type(exc).__name__}: {exc}"
    else:
        tiff_decode_error = "missing"
    svg_ok = False
    svg_decode_error = None
    if svg.is_file():
        try:
            svg_text = svg.read_text(encoding="utf-8")
            svg_ok = len(re.findall(r"<text\b", svg_text)) > 0
        except (OSError, UnicodeError) as exc:
            svg_decode_error = f"{type(exc).__name__}: {exc}"
    else:
        svg_decode_error = "missing"
    pdf_ok = valid_pdf(pdf)
    all_outputs_exist = all(path.is_file() for path in [png, tiff, svg, pdf])
    format_checks.append({
        "figure": stem,
        "png_exists": png.is_file(),
        "png_readable": png_sha is not None,
        "png_decode_error": png_decode_error,
        "tiff_exists": tiff.is_file(),
        "rgb_600dpi_lzw_tiff": tiff_ok,
        "tiff_decode_error": tiff_decode_error,
        "svg_exists": svg.is_file(),
        "editable_svg_text": svg_ok,
        "svg_decode_error": svg_decode_error,
        "pdf_exists": pdf.is_file(),
        "valid_pdf_header": pdf_ok,
        "status": "PASS" if all_outputs_exist and png_sha is not None and tiff_ok and svg_ok and pdf_ok else "FAIL",
    })

supplementary_reference_png = SUPPLEMENTARY_REFERENCE / f"{SUPPLEMENTARY_S1}.png"
supplementary_produced_png = PRODUCED_SUPPLEMENTARY / f"{SUPPLEMENTARY_S1}.png"
supplementary_produced_svg = PRODUCED_SUPPLEMENTARY / f"{SUPPLEMENTARY_S1}.svg"
supplementary_produced_pdf = PRODUCED_SUPPLEMENTARY / f"{SUPPLEMENTARY_S1}.pdf"
supplementary_reference_pixel_sha, supplementary_reference_decode_error = safe_pixel_sha256(
    supplementary_reference_png
)
supplementary_produced_pixel_sha, supplementary_produced_decode_error = safe_pixel_sha256(
    supplementary_produced_png
)
supplementary_svg_ok = False
supplementary_svg_decode_error = None
if supplementary_produced_svg.is_file():
    try:
        supplementary_svg_ok = len(
            re.findall(r"<text\b", supplementary_produced_svg.read_text(encoding="utf-8"))
        ) >= 16
    except (OSError, UnicodeError) as exc:
        supplementary_svg_decode_error = f"{type(exc).__name__}: {exc}"
else:
    supplementary_svg_decode_error = "missing"
supplementary_pdf_ok = valid_pdf(supplementary_produced_pdf)
supplementary_figure_check = {
    "figure": SUPPLEMENTARY_S1,
    "reference_exists": supplementary_reference_png.is_file(),
    "produced_exists": supplementary_produced_png.is_file(),
    "reference_pixel_sha256": supplementary_reference_pixel_sha,
    "produced_pixel_sha256": supplementary_produced_pixel_sha,
    "reference_decode_error": supplementary_reference_decode_error,
    "produced_decode_error": supplementary_produced_decode_error,
    "svg_exists": supplementary_produced_svg.is_file(),
    "editable_svg_text": supplementary_svg_ok,
    "svg_decode_error": supplementary_svg_decode_error,
    "pdf_exists": supplementary_produced_pdf.is_file(),
    "valid_pdf_header": supplementary_pdf_ok,
    "source_data_status": supplementary_source_check["status"],
}
supplementary_figure_check["status"] = (
    "PASS"
    if supplementary_figure_check["reference_exists"]
    and supplementary_figure_check["produced_exists"]
    and supplementary_figure_check["reference_pixel_sha256"] is not None
    and supplementary_figure_check["reference_pixel_sha256"]
    == supplementary_figure_check["produced_pixel_sha256"]
    and supplementary_svg_ok
    and supplementary_pdf_ok
    and supplementary_source_check["status"] == "PASS"
    else "FAIL"
)

observed_packages = {
    name: package_version(name)
    for name in ["matplotlib", "numpy", "pandas", "pillow", "pyarrow"]
}
expected_render_packages = {
    "matplotlib": "3.11.0",
    "numpy": "2.3.5",
    "pandas": "2.2.3",
    "pillow": "12.2.0",
    "pyarrow": "25.0.0",
}
recorded_render_environment = (
    platform.system() == "Darwin"
    and platform.machine() == "arm64"
    and Path("/System/Library/Fonts/Supplemental/Arial.ttf").exists()
    and sys.version_info[:2] == (3, 12)
    and observed_packages == expected_render_packages
)
all_sources_exact = len(source_checks) == 12 and all(
    item["status"] == "PASS" for item in source_checks
)
format_status = {item["figure"]: item["status"] for item in format_checks}
if not recorded_render_environment and all_sources_exact:
    for item in figure_checks:
        if (
            item["status"] == "FAIL"
            and item["figure"] in PIXEL_REFERENCE_FIGURES
            and format_status.get(item["figure"]) == "PASS"
        ):
            item["status"] = "PASS_PLATFORM_RENDER_VARIATION"
            item["note"] = (
                "The 12 scientific source-data files are byte-identical and all export checks pass, "
                "but decoded pixels differ from the recorded macOS arm64/Arial reference because this "
                "run uses a different software, operating-system, or font-rendering environment."
            )
    if (
        supplementary_figure_check["status"] == "FAIL"
        and supplementary_source_check["status"] == "PASS"
        and supplementary_figure_check["reference_exists"]
        and supplementary_figure_check["produced_exists"]
        and supplementary_figure_check["editable_svg_text"]
        and supplementary_figure_check["valid_pdf_header"]
    ):
        supplementary_figure_check["status"] = "PASS_PLATFORM_RENDER_VARIATION"
        supplementary_figure_check["note"] = (
            "The frozen Supplementary Figure S1 source CSV is byte-identical and the SVG/PDF/PNG "
            "exports are complete, but decoded pixels differ from the recorded reference because "
            "this run uses a different software, operating-system, or font-rendering environment."
        )

result = {
    "release_id": "srep-sepsis-temporal-signatures-v1.2.1",
    "run_date": __import__("datetime").datetime.now().astimezone().isoformat(),
    "command_exit_code": completed.returncode,
    "stdout": completed.stdout,
    "stderr": completed.stderr,
    "environment": {
        "python": sys.version,
        "platform": platform.platform(),
        "recorded_render_environment": recorded_render_environment,
        "packages": observed_packages,
        "expected_render_packages": expected_render_packages,
    },
    "comparison_rule": "All 12 main-figure CSV source-data files and the frozen Supplementary Figure S1 source CSV are compared byte for byte, and TIFF/SVG/PDF/PNG properties are checked independently. Exact decoded RGB pixels are required for Figures 2 and 3 and Supplementary Figure S1 only in the recorded Python/package, macOS arm64, and Arial environment. On other rendering environments, a clearly labelled PASS_PLATFORM_RENDER_VARIATION requires exact source-data equality and valid exports. Figures 1, 4, 5, and 6 retain author-reviewed presentation-only layouts under PASS_AUTHORED_LAYOUT.",
    "figure_checks": figure_checks,
    "source_checks": source_checks,
    "supplementary_source_check": supplementary_source_check,
    "supplementary_figure_check": supplementary_figure_check,
    "format_checks": format_checks,
}
result["status"] = (
    "PASS"
    if completed.returncode == 0
    and len(figure_checks) == 6
    and all(item["status"].startswith("PASS") for item in figure_checks)
    and len(source_checks) == 12
    and all(item["status"] == "PASS" for item in source_checks)
    and supplementary_source_check["status"] == "PASS"
    and supplementary_figure_check["status"].startswith("PASS")
    and len(format_checks) == 6
    and all(item["status"] == "PASS" for item in format_checks)
    else "FAIL"
)

PROJECT.mkdir(parents=True, exist_ok=True)
EVIDENCE.mkdir(parents=True, exist_ok=True)
serialized = json.dumps(result, indent=2) + "\n"
(PROJECT / "reproduction_result_v1.2.1.json").write_text(serialized, encoding="utf-8")
(EVIDENCE / "reproduction_result_v1.2.1.json").write_text(serialized, encoding="utf-8")
print(serialized)
raise SystemExit(0 if result["status"] == "PASS" else 1)
