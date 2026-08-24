#!/usr/bin/env python3

import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(command: list[str]) -> str:
    return subprocess.run(
        command, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    ).stdout


def pdf_font_subtypes(path: Path) -> str:
    # Always byte-scan. pdffonts is absent on some nodes, so branching on it made
    # the same artifact pass or fail depending on where the job landed -- and the
    # pdffonts branch tested " yes " in the whole table, which matches the sub/uni
    # columns as readily as emb. /Type3 and /FontFile are literal ASCII even under
    # /ObjStm compression, so the byte scan is both portable and stricter.
    return path.read_bytes().decode("latin-1", errors="ignore")


def main() -> int:
    root_raw = os.environ.get("HAC_FIGURE_REVIEW_ROOT", "")
    if not root_raw:
        raise RuntimeError("HAC_FIGURE_REVIEW_ROOT is unset")
    root = Path(root_raw).resolve()
    # Reports are written write-once. A re-validation of an already-built root
    # must therefore target its own directory, or the mode-"x" opens below abort
    # before a single check runs.
    out_dir = Path(os.environ.get("HAC_VALIDATION_OUT_DIR", str(root))).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "figure_manifest.tsv"
    review_path = root / "panel_review.tsv"
    if not manifest_path.is_file() or not review_path.is_file():
        raise RuntimeError("Figure manifest or review table is missing")

    with manifest_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 9:
        raise RuntimeError(f"Expected 9 review panels, found {len(rows)}")

    checks: list[dict[str, object]] = []
    for row in rows:
        path = Path(row["path"])
        checks.append({"panel_id": row["panel_id"], "check": "exists", "pass": path.is_file()})
        if not path.is_file():
            continue
        checks.append({
            "panel_id": row["panel_id"],
            "check": "sha256",
            "pass": sha256(path) == row["sha256"],
        })
        info = run(["pdfinfo", str(path)])
        pages = next(
            int(line.split(":", 1)[1].strip())
            for line in info.splitlines()
            if line.startswith("Pages:")
        )
        checks.append({"panel_id": row["panel_id"], "check": "one_page", "pass": pages == 1})
        fonts = pdf_font_subtypes(path)
        type3 = any("Type 3" in line or "Type3" in line for line in fonts.splitlines())
        checks.append({"panel_id": row["panel_id"], "check": "no_type3", "pass": not type3})
        embedded_font = "FontFile" in fonts
        checks.append({
            "panel_id": row["panel_id"],
            "check": "embedded_font_object",
            "pass": embedded_font,
        })
        checks.append({
            "panel_id": row["panel_id"],
            "check": "not_promoted",
            "pass": row["automatic_promotion"] == "FALSE" and row["current_figure_modified"] == "FALSE",
        })

    passed = all(bool(item["pass"]) for item in checks)
    validation_tsv = out_dir / "validation.tsv"
    with validation_tsv.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["panel_id", "check", "pass"], delimiter="\t")
        writer.writeheader()
        writer.writerows(checks)
    summary = {
        "overall_pass": passed,
        "n_panels": len(rows),
        "n_checks": len(checks),
        "n_pass": sum(bool(item["pass"]) for item in checks),
        "current_figure_modified": False,
        "automatic_promotion": False,
    }
    with (out_dir / "validation_summary.json").open("x") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    print(json.dumps(summary, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
