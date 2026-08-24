#!/usr/bin/env python3
"""Validate release-candidate Fig. 2A, 2E, and 2F artifacts."""

from __future__ import annotations

import csv
import os
import re
import subprocess
import sys
from pathlib import Path


EXPECTED = {
    "fig2A_gwas_cascade.pdf": (3.53, 2.35),
    "Fig2B_PIP_vs_SuSiE-coloc.pdf": (1.80, 2.35),
    "Fig2E_ancestry_unique_coloc_GWS.pdf": (2.00, 2.10),
    "Fig2F_crossancestry_coloc.pdf": (2.58, 2.10),
}
PROHIBITED = ("ABF-fallback", "exploratory", "causal gene", "ancestry-specific")
REQUIRED_TERMS = {
    "fig2A_gwas_cascade.pdf": ("Multi-signal COLOC", "Single-signal", "COLOC only"),
    "Fig2B_PIP_vs_SuSiE-coloc.pdf": ("Multi-signal COLOC",),
    "Fig2E_ancestry_unique_coloc_GWS.pdf": ("multi-signal COLOC",),
    "Fig2F_crossancestry_coloc.pdf": ("multi-signal", "single-signal only"),
}
TOLERANCE_IN = 0.015


def pdf_info(path: Path) -> tuple[int, float, float]:
    text = subprocess.check_output(["pdfinfo", str(path)], text=True)
    pages = int(re.search(r"^Pages:\s+(\d+)", text, re.M).group(1))
    match = re.search(r"^Page size:\s+([\d.]+) x ([\d.]+) pts", text, re.M)
    return pages, float(match.group(1)) / 72, float(match.group(2)) / 72


def pdf_text(path: Path) -> str:
    return subprocess.check_output(["pdftotext", str(path), "-"], text=True)


def check_source_table(path: Path) -> list[str]:
    errors: list[str] = []
    required = {
        "gene", "ancestry", "gwas_ancestry", "eqtl_ancestry", "eqtl_source",
        "evidence_state", "multi_pp4", "multi_study", "multi_trait",
        "multi_lead_variant", "multi_ld_panel", "multi_ld_reliability",
        "single_pp4", "single_study", "single_trait", "single_lead_variant",
        "single_ld_panel", "single_ld_reliability", "multi_signal_available",
        "single_signal_available", "display_pp4", "display_method",
    }
    delimiter = "\t" if path.suffix == ".tsv" else ","
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter=delimiter))
    missing = required - set(rows[0] if rows else ())
    if missing:
        errors.append(f"{path.name}: missing source fields {sorted(missing)}")
        return errors
    allowed = {"multi_signal", "single_signal_only", "evaluated_no_support", "not_evaluable"}
    for row in rows:
        state = row["evidence_state"]
        if state not in allowed:
            errors.append(f"{path.name}: invalid state {state}")
        multi = float(row["multi_pp4"]) if row["multi_pp4"] not in ("", "NA") else None
        single = float(row["single_pp4"]) if row["single_pp4"] not in ("", "NA") else None
        expected = (
            "multi_signal" if multi is not None and multi > 0.5 else
            "single_signal_only" if single is not None and single > 0.5 else
            "evaluated_no_support" if multi is not None or single is not None else
            "not_evaluable"
        )
        if state != expected:
            errors.append(f"{path.name}: {row['gene']}/{row['ancestry']} is {state}, expected {expected}")
        expected_display = (
            multi if state == "multi_signal" else
            single if state == "single_signal_only" else
            max(value for value in (multi, single) if value is not None)
            if state == "evaluated_no_support" else None
        )
        display = float(row["display_pp4"]) if row["display_pp4"] not in ("", "NA") else None
        if expected_display is None and display is not None or (
            expected_display is not None and (display is None or abs(display - expected_display) > 1e-12)
        ):
            errors.append(f"{path.name}: invalid displayed PP.H4 for {row['gene']}/{row['ancestry']}")
        expected_method = {
            "multi_signal": "multi_signal",
            "single_signal_only": "single_signal",
            "evaluated_no_support": "maximum_available",
            "not_evaluable": "not_evaluable",
        }[state]
        if row["display_method"] != expected_method:
            errors.append(f"{path.name}: invalid displayed method for {row['gene']}/{row['ancestry']}")
        if row["gwas_ancestry"] != row["ancestry"] or row["eqtl_ancestry"] != "EUR":
            errors.append(f"{path.name}: ancestry provenance mismatch for {row['gene']}/{row['ancestry']}")
    return errors


def main() -> int:
    candidate = Path(os.environ.get("FIG2_CANDIDATE_DIR", sys.argv[1] if len(sys.argv) > 1 else ""))
    if not candidate.is_dir():
        print(f"ERROR: candidate directory not found: {candidate}", file=sys.stderr)
        return 2
    errors: list[str] = []
    for name, (want_w, want_h) in EXPECTED.items():
        path = candidate / name
        if not path.exists():
            errors.append(f"missing {name}")
            continue
        pages, width, height = pdf_info(path)
        if pages != 1:
            errors.append(f"{name}: {pages} pages, expected 1")
        if abs(width - want_w) > TOLERANCE_IN or abs(height - want_h) > TOLERANCE_IN:
            errors.append(f"{name}: {width:.3f}x{height:.3f} in, expected {want_w:.2f}x{want_h:.2f}")
        text = pdf_text(path)
        for term in PROHIBITED:
            if term.lower() in text.lower():
                errors.append(f"{name}: prohibited artwork label {term!r}")
        for term in REQUIRED_TERMS.get(name, ()):
            if term.lower() not in text.lower():
                errors.append(f"{name}: missing required artwork label {term!r}")
        raw = path.read_bytes()
        if b"/FontFile" not in raw:
            errors.append(f"{name}: no embedded FontFile object detected")

    source_tsv = candidate / "Fig2F_crossancestry_coloc_source.tsv"
    source_csv = candidate / "Fig2F_crossancestry_coloc_source.csv"
    for source in (source_tsv, source_csv):
        if not source.exists():
            errors.append(f"missing {source.name}")
        else:
            errors.extend(check_source_table(source))
    if source_tsv.exists() and source_csv.exists():
        with source_tsv.open(newline="") as handle:
            tsv_rows = list(csv.DictReader(handle, delimiter="\t"))
        with source_csv.open(newline="") as handle:
            csv_rows = list(csv.DictReader(handle))
        if tsv_rows != csv_rows:
            errors.append("Fig2F CSV and TSV source tables disagree")
    audit = candidate / "Fig2H_noneur_candidate_audit.tsv"
    blocked = candidate / "Fig2H_SELECTION_BLOCKED.txt"
    if not audit.exists() or not blocked.exists():
        errors.append("Fig2H release-block audit/status is incomplete")

    report = candidate / "validation_report.txt"
    if errors:
        report.write_text("FAIL\n" + "\n".join(f"- {error}" for error in errors) + "\n")
        print(report.read_text(), end="")
        return 1
    report.write_text(
        "PASS\n- one page and exact panel dimensions\n- embedded font objects detected\n"
        "- prohibited artwork labels absent\n- Fig2F thresholds and provenance rederived\n"
        "- Fig2H remains release-blocked\n"
    )
    print(report.read_text(), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
