#!/usr/bin/env python3
"""Fail-closed validation for the user-approved Figure 4D/E promotion."""

from __future__ import annotations

import csv
import hashlib
import re
import subprocess
import sys
from pathlib import Path


BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FIG = BASE / "figures/main/fig4_singlecell_programs"
CANDIDATE = BASE / "figures/candidates/fig4def-multicellular-story-2026-08-17-v2"
PROVENANCE = FIG / "provenance/fig4de_multicellular_story_2026-08-17_v1"
REPORT = PROVENANCE / "promotion_validation_v2.tsv"

checks: list[tuple[str, str, str]] = []


def add(name: str, passed: bool, detail: str) -> None:
    checks.append((name, "PASS" if passed else "FAIL", detail))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def shell(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, capture_output=True, text=True, check=False)


expected_panels = {
    FIG / "panels/fig4d_same_atlas_lineage_specificity.pdf": (
        CANDIDATE / "panels/fig4e_same_atlas_lineage_specificity.pdf",
        "e8027909eb42235993ede0852e3e363cefa459c34e35c973be1eb2c05b2eb719",
    ),
    FIG / "panels/fig4e_bulk_tissue_state_transport.pdf": (
        CANDIDATE / "panels/fig4f_bulk_tissue_state_transport.pdf",
        "1f45726a22f6b2cccb8ddeafe01180a55db9380d52a775716a812df74e2292de",
    ),
}
for active, (source, expected_hash) in expected_panels.items():
    add(f"exists:{active.name}", active.is_file(), str(active))
    add(f"candidate_exists:{source.name}", source.is_file(), str(source))
    if active.is_file() and source.is_file():
        add(f"byte_reproduction:{active.name}", active.read_bytes() == source.read_bytes(),
            sha256(active))
        add(f"checksum:{active.name}", sha256(active) == expected_hash, sha256(active))
        info = shell(["pdfinfo", str(active)])
        pages = re.search(r"^Pages:\s+(\d+)", info.stdout, flags=re.MULTILINE)
        add(f"one_page:{active.name}", info.returncode == 0 and pages is not None and pages.group(1) == "1",
            pages.group(1) if pages else info.stderr.strip())
        gs = shell(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-dPDFDEBUG", "-sDEVICE=nullpage", str(active)])
        raw = active.read_bytes()
        has_type3 = (
            re.search(rb"/Subtype\s*/Type3\b", raw) is not None
            or re.search(rb"/FontType\s*3\b", raw) is not None
            or "/Subtype /Type3" in (gs.stdout + gs.stderr)
            or "/FontType 3" in (gs.stdout + gs.stderr)
        )
        add(f"ghostscript:{active.name}", gs.returncode == 0, str(gs.returncode))
        add(f"no_type3:{active.name}", not has_type3, "no Type 3 marker")

old_active = [
    FIG / "panels/fig4e_bulk_projection_pre_genetics.pdf",
    FIG / "panels/fig4f_tf_activity.pdf",
]
for path in old_active:
    add(f"old_active_removed:{path.name}", not path.exists(), str(path))
archived = [
    FIG / "panels/_superseded_main_4ef_2026-08-17_v1/fig4e_bulk_projection_pre_genetics.pdf",
    FIG / "panels/_superseded_main_4ef_2026-08-17_v1/fig4f_tf_activity.pdf",
]
for path in archived:
    add(f"old_panel_archived:{path.name}", path.is_file(), str(path))
add("ambient_panel_not_promoted",
    not (FIG / "panels/fig4d_ambient_effect_transport.pdf").exists(),
    "no active ambient panel")

active_sources = {
    FIG / "source_tables/current_candidate/fig4d_same_atlas_lineage_specificity.tsv":
        CANDIDATE / "source_tables/fig4e_same_atlas_lineage_specificity.tsv",
    FIG / "source_tables/current_candidate/fig4e_bulk_tissue_state_transport.tsv":
        CANDIDATE / "source_tables/fig4f_bulk_stage_transport.tsv",
}
contract_columns = {
    "biological_unit", "multiple_testing_family", "evidence_state",
    "state_reason", "unresolved_alternative", "claim_boundary",
}
for active, source in active_sources.items():
    add(f"source_exists:{active.name}", active.is_file(), str(active))
    if active.is_file() and source.is_file():
        add(f"source_byte_reproduction:{active.name}", active.read_bytes() == source.read_bytes(),
            sha256(active))
        with active.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
            fields = set(rows[0]) if rows else set()
        add(f"source_contract:{active.name}", contract_columns <= fields,
            ";".join(sorted(contract_columns - fields)) or "all columns present")
        add(f"source_nonempty:{active.name}",
            all(all(row.get(column, "") for column in contract_columns) for row in rows),
            f"{len(rows)} rows")

with (FIG / "source_tables/current_candidate/fig4d_same_atlas_lineage_specificity.tsv").open(newline="") as handle:
    lineage_rows = list(csv.DictReader(handle, delimiter="\t"))
nonreference = [row for row in lineage_rows if "reference" not in row["comparison"]]
add("lineage_rows", len(lineage_rows) == 12, str(len(lineage_rows)))
add("lineage_prespecified_family", len(nonreference) == 10, str(len(nonreference)))
focus = {
    (row["program_name"], row["comparison_lineage"]): row
    for row in nonreference
}
add("lineage_ecm_fibroblast_hc3",
    float(focus[("ECM/IGFBP7", "Fibroblasts")]["hc3_qvalue"]) < 0.05,
    focus[("ECM/IGFBP7", "Fibroblasts")]["hc3_qvalue"])
add("lineage_ductular_cholangiocyte_hc3",
    float(focus[("Ductular-injury/BICC1", "Cholangiocytes")]["hc3_qvalue"]) < 0.05,
    focus[("Ductular-injury/BICC1", "Cholangiocytes")]["hc3_qvalue"])

with (FIG / "source_tables/current_candidate/fig4e_bulk_tissue_state_transport.tsv").open(newline="") as handle:
    bulk_rows = list(csv.DictReader(handle, delimiter="\t"))
add("bulk_transport_rows", len(bulk_rows) == 8, str(len(bulk_rows)))
add("bulk_transport_stages", {row["stage"] for row in bulk_rows} == {"F1", "F2", "F3", "F4"},
    ";".join(sorted({row["stage"] for row in bulk_rows})))
add("bulk_transport_programs", {row["program_name"] for row in bulk_rows} == {"ECM/IGFBP7", "Ductular-injury/BICC1"},
    ";".join(sorted({row["program_name"] for row in bulk_rows})))

manifest = FIG / "manifests/current_candidate_manifest.tsv"
with manifest.open(newline="") as handle:
    manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
for row in manifest_rows:
    path = FIG / row["path"]
    if path.exists():
        add(f"manifest:{row['role']}",
            path.stat().st_size == int(row["size_bytes"]) and sha256(path) == row["sha256"],
            sha256(path))
    else:
        add(f"manifest:{row['role']}", False, f"missing {path}")

index = (FIG / "CANDIDATE_PANEL_INDEX.tsv").read_text()
add("index_has_main_4d", "fig4d_same_atlas_lineage_specificity.pdf" in index, "panel path")
add("index_has_main_4e", "fig4e_bulk_tissue_state_transport.pdf" in index, "panel path")
add("index_omits_main_4f", "\n4F\t" not in index, "no active 4F row")

scope = shell(["python3", str(BASE / "scripts/manuscript/validate_resource_scope.py")])
add("resource_scope", scope.returncode == 0, (scope.stdout + scope.stderr).strip())

PROVENANCE.mkdir(parents=True, exist_ok=True)
with REPORT.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["check", "status", "detail"])
    writer.writerows(checks)

for check, status, detail in checks:
    print(f"{status}\t{check}\t{detail}")
if any(status == "FAIL" for _, status, _ in checks):
    sys.exit(1)
print("FIG4DE_PROMOTION_VALIDATION PASS")
