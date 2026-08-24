#!/usr/bin/env python3
"""Validate the approved Figure 4D geometry promotion without mutating files."""

from __future__ import annotations

import csv
import hashlib
import re
import subprocess
from pathlib import Path

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FIG = ROOT / "figures/main/fig4_singlecell_programs"
CAND = ROOT / "figures/candidates/fig4de-geometry-review-2026-08-17-v5"
OUT = FIG / "provenance/fig4d_lineage_geometry_2026-08-17_v1/promotion_validation.tsv"

checks: list[tuple[str, bool, str]] = []


def add(name: str, passed: bool, observed: object) -> None:
    checks.append((name, bool(passed), str(observed)))


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


active_panel = FIG / "panels/fig4d_same_atlas_lineage_specificity.pdf"
candidate_panel = CAND / "panels/fig4d_lineage_enrichment_forest.pdf"
active_source = FIG / "source_tables/current_candidate/fig4d_same_atlas_lineage_specificity.tsv"
candidate_source = CAND / "source_tables/fig4d_compact_lineage_forest.tsv"
bulk_panel = FIG / "panels/fig4e_bulk_tissue_state_transport.pdf"
archive = FIG / "panels/_superseded_main_4d_2026-08-17_v2"

add("panel_exact_candidate_bytes", active_panel.read_bytes() == candidate_panel.read_bytes(), digest(active_panel))
add("source_exact_candidate_bytes", active_source.read_bytes() == candidate_source.read_bytes(), digest(active_source))
add("bulk_gradual_line_unchanged", digest(bulk_panel) == "01305a2a2b82cd8746544359007a856f243dd733f89f951a4e1909d3d5f76f87", digest(bulk_panel))
add("old_panel_archived", digest(archive / "fig4d_same_atlas_lineage_specificity.pdf") == "e8027909eb42235993ede0852e3e363cefa459c34e35c973be1eb2c05b2eb719", digest(archive / "fig4d_same_atlas_lineage_specificity.pdf"))
add("old_source_archived", digest(archive / "fig4d_same_atlas_lineage_specificity.tsv") == "37664ca62ee918a11caa1bc9a2d86add4f891726ecb5345085b5a7f71354eea1", digest(archive / "fig4d_same_atlas_lineage_specificity.tsv"))

with active_source.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
programs = {row["program_uid"] for row in rows}
add("source_rows_10", len(rows) == 10, len(rows))
add("two_frozen_programs", len(programs) == 2, len(programs))
counts = sorted(sum(row["program_uid"] == uid for row in rows) for uid in programs)
add("five_lineages_per_program", counts == [5, 5], ";".join(map(str, counts)))
add("no_reference_rows", all(row["comparison_lineage"] != "Hepatocytes" for row in rows), sum(row["comparison_lineage"] == "Hepatocytes" for row in rows))
fields = ("biological_unit", "multiple_testing_family", "evidence_state", "state_reason", "unresolved_alternative", "claim_boundary")
add("contract_fields_complete", all(all(row.get(field, "") for field in fields) for row in rows), len(rows))
targets = [row for row in rows if row.get("target", "").upper() == "TRUE"]
add("two_highlighted_targets", len(targets) == 2, len(targets))
target_lineages = {row["comparison_lineage"] for row in targets}
add("target_lineages_exact", target_lineages == {"Fibroblasts", "Cholangiocytes"}, ";".join(sorted(target_lineages)))
add("targets_hc3_family_supported", all(float(row["hc3_qvalue"]) < 0.05 for row in targets), ";".join(row["hc3_qvalue"] for row in targets))

with (FIG / "manifests/current_candidate_manifest.tsv").open(newline="") as handle:
    manifest = {row["path"]: row for row in csv.DictReader(handle, delimiter="\t")}
for relative, path in (
    ("panels/fig4d_same_atlas_lineage_specificity.pdf", active_panel),
    ("source_tables/current_candidate/fig4d_same_atlas_lineage_specificity.tsv", active_source),
    ("panels/fig4e_bulk_tissue_state_transport.pdf", bulk_panel),
):
    row = manifest[relative]
    add(f"manifest_hash:{relative}", row["sha256"] == digest(path), row["sha256"])
    add(f"manifest_size:{relative}", int(row["size_bytes"]) == path.stat().st_size, row["size_bytes"])

pdfinfo = subprocess.run(["pdfinfo", str(active_panel)], text=True, capture_output=True, check=False)
pages = re.search(r"^Pages:\s+(\d+)", pdfinfo.stdout, flags=re.MULTILINE)
add("pdf_one_page", pdfinfo.returncode == 0 and pages is not None and pages.group(1) == "1", pages.group(1) if pages else "missing")
ghostscript = subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=nullpage", str(active_panel)], text=True, capture_output=True, check=False)
add("ghostscript", ghostscript.returncode == 0, ghostscript.returncode)
raw = active_panel.read_bytes()
add("no_type3", re.search(rb"/Subtype\s*/Type3\b|/FontType\s*3\b", raw) is None, "binary scan")

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
print(f"FIG4D_GEOMETRY_PROMOTION_VALIDATION\tPASS\t{len(checks)} checks")
