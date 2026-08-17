#!/usr/bin/env python3
"""Validate the isolated Figure 4C–D program-to-CCC bridge candidate."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
import re
import subprocess
import sys


CONTRACT_FIELDS = {
    "biological_unit",
    "multiple_testing_family",
    "evidence_state",
    "state_reason",
    "unresolved_alternative",
    "claim_boundary",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: validate_fig4cd_ccc_bridge.py CANDIDATE_ROOT CURRENT_4C")
    root = Path(sys.argv[1]).resolve()
    current_4c = Path(sys.argv[2]).resolve()
    panels = root / "panels"
    sources = root / "source_tables"
    provenance = root / "provenance"
    report: list[tuple[str, bool, str]] = []

    def check(name: str, passed: bool, detail: object) -> None:
        report.append((name, bool(passed), str(detail)))

    pdfs = [
        panels / "fig4c_hotspot_stage_heatmap.pdf",
        panels / "fig4d_program_linked_communication.pdf",
    ]
    check("two_panel_pdfs", all(path.is_file() and path.stat().st_size > 0 for path in pdfs), pdfs)
    if pdfs[0].exists() and current_4c.exists():
        check("figure_4c_byte_identical", sha256(pdfs[0]) == sha256(current_4c), sha256(pdfs[0]))

    for pdf in pdfs:
        if not pdf.exists():
            continue
        info = run(["pdfinfo", str(pdf)])
        pages = [line.split(":", 1)[1].strip() for line in info.stdout.splitlines() if line.startswith("Pages:")]
        check(f"one_page:{pdf.name}", info.returncode == 0 and pages == ["1"], pages)
        gs = run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-dPDFDEBUG", "-sDEVICE=nullpage", str(pdf)])
        debug = f"{gs.stdout}\n{gs.stderr}"
        raw = pdf.read_bytes()
        has_type3 = (
            re.search(rb"/Subtype\s*/Type3\b", raw) is not None
            or re.search(rb"/FontType\s+3\b", raw) is not None
            or "/Subtype /Type3" in debug
            or "/FontType 3" in debug
        )
        check(f"ghostscript:{pdf.name}", gs.returncode == 0, gs.stderr.strip() or "passed")
        check(f"no_type3:{pdf.name}", not has_type3, "passed" if not has_type3 else "Type 3 detected")

    expected_sources = {
        "fig4d_fixed_13_stage_context.tsv": 39,
        "fig4d_program_ccc_associations.tsv": 26,
        "fig4d_program_ccc_member_overlap.tsv": 26,
        "fig4d_program_ccc_lodo.tsv": None,
        "fig4d_program_ccc_matched_donors.tsv": None,
        "fig4d_supported_bridge_partial_residuals.tsv": 60,
    }
    loaded: dict[str, list[dict[str, str]]] = {}
    for name, expected_rows in expected_sources.items():
        path = sources / name
        check(f"source_exists:{name}", path.is_file() and path.stat().st_size > 0, path.stat().st_size if path.exists() else 0)
        if not path.exists():
            continue
        rows = read_tsv(path)
        loaded[name] = rows
        if expected_rows is not None:
            check(f"source_rows:{name}", len(rows) == expected_rows, len(rows))
        fields = set(rows[0]) if rows else set()
        check(f"contract_fields:{name}", CONTRACT_FIELDS.issubset(fields), sorted(CONTRACT_FIELDS - fields))

    stage = loaded.get("fig4d_fixed_13_stage_context.tsv", [])
    check("fixed_13_stage_pairs", len({row.get("headline_label") for row in stage}) == 13, len({row.get("headline_label") for row in stage}))
    check("fixed_13_zero_supported", sum(row.get("donor_state") == "BH-supported" for row in stage) == 0, "expected 0")

    assoc = loaded.get("fig4d_program_ccc_associations.tsv", [])
    association_pairs = {(row.get("ct_pair"), row.get("lr_pair")) for row in assoc}
    check("association_2x13", len({row.get("program_uid") for row in assoc}) == 2 and len(association_pairs) == 13, "2 programs × 13 pairs")
    check("association_family_declared_26", all("26 planned tests" in row.get("multiple_testing_family", "") for row in assoc), "26 planned tests")
    supported = [
        row for row in assoc
        if row.get("qvalue") and row.get("hc3_qvalue")
        and float(row["qvalue"]) < 0.05 and float(row["hc3_qvalue"]) < 0.05
    ]
    check(
        "one_supported_program_ccc_bridge",
        len(supported) == 1
        and supported[0].get("program_short") == "Ductular / BICC1"
        and supported[0].get("lr_pair") == "CDH1__PTPRM",
        [(row.get("program_short"), row.get("lr_pair")) for row in supported],
    )
    check("hc3_se_complete_for_testable", all(row.get("hc3_se") for row in assoc if row.get("testable") == "TRUE"), "testable HC3 SE fields")

    partial = loaded.get("fig4d_supported_bridge_partial_residuals.tsv", [])
    check(
        "supported_bridge_partial_residual_donors",
        len(partial) == 60 and len({row.get("donor") for row in partial}) == 60,
        f"rows={len(partial)};donors={len({row.get('donor') for row in partial})}",
    )

    overlap = loaded.get("fig4d_program_ccc_member_overlap.tsv", [])
    observed = [row for row in overlap if row.get("overlap_gene")]
    observed_genes = {row.get("overlap_gene") for row in observed}
    observed_programs = {row.get("program_uid") for row in observed}
    check("member_overlap_two_ecm_only", len(observed) == 2 and observed_genes == {"COL4A1", "COL4A2"} and len(observed_programs) == 1, f"rows={len(observed)};genes={sorted(observed_genes)}")

    assertion_path = provenance / "render_assertions.tsv"
    assertions = read_tsv(assertion_path) if assertion_path.exists() else []
    failed_assertions = [row for row in assertions if row.get("passed", "").upper() != "TRUE"]
    check("render_assertions_pass", bool(assertions) and not failed_assertions, f"rows={len(assertions)};failed={len(failed_assertions)}")

    report_path = provenance / "validation_report.tsv"
    provenance.mkdir(parents=True, exist_ok=True)
    with report_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["check", "passed", "detail"])
        writer.writerows(report)

    failed = [row for row in report if not row[1]]
    if failed:
        for name, _, detail in failed:
            print(f"[FAIL] {name}: {detail}", file=sys.stderr)
        return 1

    manifest_path = provenance / "output_manifest.tsv"
    excluded = {manifest_path.resolve(), (root / "VALIDATED").resolve()}
    artifacts = sorted(path for path in root.rglob("*") if path.is_file() and path.resolve() not in excluded)
    with manifest_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["relative_path", "size_bytes", "sha256"])
        for path in artifacts:
            writer.writerow([path.relative_to(root), path.stat().st_size, sha256(path)])

    with (root / "VALIDATED").open("w") as handle:
        handle.write("status\tvalidated_candidate\n")
        handle.write(f"validation_report_sha256\t{sha256(report_path)}\n")
        handle.write(f"output_manifest_sha256\t{sha256(manifest_path)}\n")
        handle.write("promotion_state\tnot_promoted_requires_user_review\n")

    print(f"[validated] {root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
