#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXCLUDED = {".git", "MANIFEST_SHA256.txt", "release_manifest.json"}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


files = []
for path in sorted(ROOT.rglob("*")):
    rel = path.relative_to(ROOT)
    if not path.is_file() or any(part in EXCLUDED for part in rel.parts):
        continue
    files.append(
        {
            "path": rel.as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": digest(path),
        }
    )

manifest = {
    "release_id": "srep-sepsis-temporal-signatures-v1.0.1",
    "doi": "10.5281/zenodo.21417810",
    "concept_doi": "10.5281/zenodo.21415496",
    "repository": "https://github.com/LiXinzhuo0425/srep-sepsis-temporal-signatures",
    "files": files,
}
(ROOT / "release_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
(ROOT / "MANIFEST_SHA256.txt").write_text(
    "".join(f"{item['sha256']}  {item['path']}\n" for item in files), encoding="utf-8"
)
print(json.dumps({"release_id": manifest["release_id"], "files": len(files)}, indent=2))
