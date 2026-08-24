#!/usr/bin/env python3
"""Validate removal of the Figure 4E point-size encoding."""

from __future__ import annotations

import csv
import hashlib
import math
import re
import subprocess
from pathlib import Path

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FIG = ROOT / "figures/main/fig4_singlecell_programs"
CAND = ROOT / "figures/candidates/fig4e-uniform-points-2026-08-17-v1"
OUT = FIG / "provenance/fig4e_uniform_points_2026-08-17_v1/promotion_validation.tsv"
checks: list[tuple[str, bool, str]] = []


def add(name: str, passed: bool, observed: object) -> None:
    checks.append((name, bool(passed), str(observed)))


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


active = FIG / "panels/fig4e_bulk_tissue_state_transport.pdf"
candidate = CAND / "panels/fig4e_bulk_tissue_state_transport.pdf"
source = FIG / "source_tables/current_candidate/fig4e_bulk_tissue_state_transport.tsv"
candidate_source = CAND / "source_tables/fig4e_bulk_tissue_state_transport.tsv"
archive = FIG / "panels/_superseded_main_4e_2026-08-17_v2/fig4e_bulk_tissue_state_transport.pdf"
active_4d = FIG / "panels/fig4d_same_atlas_lineage_specificity.pdf"

add("panel_exact_candidate_bytes", active.read_bytes() == candidate.read_bytes(), digest(active))
add("source_exact_candidate_bytes", source.read_bytes() == candidate_source.read_bytes(), digest(source))
add("old_sized_panel_archived", digest(archive) == "1f45726a22f6b2cccb8ddeafe01180a55db9380d52a775716a812df74e2292de", digest(archive))
add("active_4d_unchanged", digest(active_4d) == "e13c836db28c2ba9c8a3c22bd43963025f43d4870026ff26f48163ce4b316c8e", digest(active_4d))

with source.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
add("source_rows_8", len(rows) == 8, len(rows))
add("two_programs", len({row["program_uid"] for row in rows}) == 2, len({row["program_uid"] for row in rows}))
add("four_stage_contrasts", {row["stage"] for row in rows} == {"F1", "F2", "F3", "F4"}, ";".join(sorted({row["stage"] for row in rows})))
add("effects_finite", all(math.isfinite(float(row["effect"])) for row in rows), len(rows))
fields = ("biological_unit", "multiple_testing_family", "evidence_state", "state_reason", "unresolved_alternative", "claim_boundary")
add("contract_fields_complete", all(all(row.get(field, "") for field in fields) for row in rows), len(rows))

with (FIG / "manifests/current_candidate_manifest.tsv").open(newline="") as handle:
    manifest = {row["path"]: row for row in csv.DictReader(handle, delimiter="\t")}
for relative, path in (
    ("panels/fig4e_bulk_tissue_state_transport.pdf", active),
    ("source_tables/current_candidate/fig4e_bulk_tissue_state_transport.tsv", source),
    ("panels/fig4d_same_atlas_lineage_specificity.pdf", active_4d),
):
    row = manifest[relative]
    add(f"manifest_hash:{relative}", row["sha256"] == digest(path), row["sha256"])
    add(f"manifest_size:{relative}", int(row["size_bytes"]) == path.stat().st_size, row["size_bytes"])

pdfinfo = subprocess.run(["pdfinfo", str(active)], text=True, capture_output=True, check=False)
pages = re.search(r"^Pages:\s+(\d+)", pdfinfo.stdout, flags=re.MULTILINE)
add("pdf_one_page", pdfinfo.returncode == 0 and pages is not None and pages.group(1) == "1", pages.group(1) if pages else "missing")
ghostscript = subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=nullpage", str(active)], text=True, capture_output=True, check=False)
add("ghostscript", ghostscript.returncode == 0, ghostscript.returncode)
raw = active.read_bytes()
add("no_type3", re.search(rb"/Subtype\s*/Type3\b|/FontType\s*3\b", raw) is None, "binary scan")
text = subprocess.run(["pdftotext", str(active), "-"], text=True, capture_output=True, check=False)
add("size_legend_removed", text.returncode == 0 and "BH-supported" not in text.stdout and "program weight" not in text.stdout, text.stdout.replace("\n", " ").strip())

scope = subprocess.run(["python3", str(ROOT / "scripts/manuscript/validate_resource_scope.py")], text=True, capture_output=True, check=False)
add("resource_scope", scope.returncode == 0 and "RESOURCE_SCOPE_VALIDATION\tPASS" in scope.stdout, scope.returncode)

OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["check", "status", "observed"])
    for name, passed, observed in checks:
        writer.writerow([name, "PASS" if passed else "FAIL", observed])

failed = [name for name, passed, _ in checks if not passed]
if failed:
    raise SystemExit("FAIL: " + ", ".join(failed))
print(f"FIG4E_UNIFORM_PROMOTION_VALIDATION\tPASS\t{len(checks)} checks")
