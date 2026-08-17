#!/usr/bin/env python3
"""Validate compact Figure 3E and refresh the cloned candidate manifests."""

import csv
import hashlib
import os
from pathlib import Path
import re
import subprocess


root = Path(os.environ["FIGURE_CANDIDATE_ROOT"]).resolve()
pdf = root / "figure3/panels/fig3e_stage_remodeling.pdf"
relative = str(pdf.relative_to(root))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


info = subprocess.run(["pdfinfo", str(pdf)], check=True, text=True,
                      capture_output=True).stdout
pages = int(next(line.split(":", 1)[1] for line in info.splitlines()
                 if line.startswith("Pages:")))
page_size = next(line for line in info.splitlines() if line.startswith("Page size:"))
match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", page_size)
# Two-column (fibrosis | NAS) compact panel: 5.00 x 2.27 in = 360.0 x 163.4 pts.
if (pages != 1 or match is None or
        abs(float(match.group(1)) - 360.0) > 1.0 or
        abs(float(match.group(2)) - 163.0) > 1.0):
    raise RuntimeError(f"Unexpected compact panel geometry: {page_size}")
gs = subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-dPDFDEBUG",
                     "-sDEVICE=nullpage", str(pdf)], capture_output=True)
debug = gs.stdout + gs.stderr
raw = pdf.read_bytes()
has_type3 = b"/Subtype /Type3" in raw or b"/FontType 3" in raw or b"/FontType 3" in debug
if gs.returncode or has_type3 or b"/FontFamily (Helvetica)" not in raw:
    raise RuntimeError("Compact Figure 3E failed PDF/font validation")

validation_path = root / "pdf_validation.tsv"
with validation_path.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
replacement = {
    "relative_path": relative,
    "pages": "1",
    "ghostscript": "pass",
    "type3_check": "pass",
    "size_bytes": str(pdf.stat().st_size),
    "sha256": sha256(pdf),
}
rows = [replacement if row["relative_path"] == relative else row for row in rows]
if not any(row["relative_path"] == relative for row in rows):
    rows.append(replacement)
with validation_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=replacement, delimiter="\t",
                            lineterminator="\n")
    writer.writeheader()
    writer.writerows(sorted(rows, key=lambda row: row["relative_path"]))

manifest_path = root / "output_manifest.tsv"
outputs = sorted(path for path in root.rglob("*")
                 if path.is_file() and path != manifest_path)
with manifest_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle,
                            fieldnames=["relative_path", "size_bytes", "sha256"],
                            delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for path in outputs:
        writer.writerow({"relative_path": str(path.relative_to(root)),
                         "size_bytes": path.stat().st_size,
                         "sha256": sha256(path)})
print("FIG3E_COMPACT_VALIDATION PASS")
