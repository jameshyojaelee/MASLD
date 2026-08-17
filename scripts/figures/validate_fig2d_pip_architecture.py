#!/usr/bin/env python3
"""Validate the promoted working-main Figure 2D PIP-composition panel."""

from __future__ import annotations

import csv
import hashlib
import math
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PANELS = ROOT / "figures/main/fig2_genetics/panels"
PDF = PANELS / "Fig2D_pip_architecture_by_trait_directness.pdf"
SOURCE = PANELS / "Fig2D_pip_architecture_by_trait_directness_source.tsv"
MANIFEST = PANELS / "Fig2D_pip_architecture_by_trait_directness_manifest.tsv"
ARCHIVE = PANELS / "_archive/2026-08-13_pre_pip_architecture_fig2J"
EXPECTED = {
    ("tier1_direct_masld_pdff", "protein_altering_pip_mass"): 0.22368342851023765,
    ("tier1_direct_masld_pdff", "canonical_splice_pip_mass"): 0.0031211800797125142,
    ("tier1_direct_masld_pdff", "synonymous_or_utr_pip_mass"): 0.0470738652226514,
    ("tier1_direct_masld_pdff", "other_noncoding_pip_mass"): 0.7261215261873982,
    ("tier2_liver_enzyme", "protein_altering_pip_mass"): 0.03389892569783767,
    ("tier2_liver_enzyme", "canonical_splice_pip_mass"): 0.0,
    ("tier2_liver_enzyme", "synonymous_or_utr_pip_mass"): 0.03571531824327637,
    ("tier2_liver_enzyme", "other_noncoding_pip_mass"): 0.9303857560588858,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    errors: list[str] = []
    required = [
        PDF,
        SOURCE,
        MANIFEST,
        PANELS / "Fig2D_pip_architecture_by_trait_directness_contrasts.tsv",
        PANELS / "Fig2D_pip_architecture_by_trait_directness_region_grouping.tsv",
        PANELS / "Fig2D_pip_architecture_by_trait_directness_caption.txt",
        ARCHIVE / "fig2J_eqtl_power_universes.pdf",
        ARCHIVE / "fig2J_eqtl_power_universes_source.csv",
    ]
    for path in required:
        if not path.exists():
            errors.append(f"missing {path.relative_to(ROOT)}")

    if PDF.exists():
        info = subprocess.check_output(["pdfinfo", str(PDF)], text=True)
        pages = re.search(r"^Pages:\s+(\d+)", info, re.M)
        size = re.search(r"^Page size:\s+([\d.]+) x ([\d.]+) pts", info, re.M)
        if not pages or pages.group(1) != "1":
            errors.append("promoted PDF is not one page")
        if not size or abs(float(size.group(1)) / 72 - 2.10) > 0.015 or abs(
            float(size.group(2)) / 72 - 2.44
        ) > 0.015:
            errors.append("promoted PDF is not 2.10 x 2.44 inches")
        text = subprocess.check_output(["pdftotext", str(PDF), "-"], text=True)
        prohibited = (
            "Share of fine-mapping probability by variant consequence",
            "ABF-fallback",
            "exploratory",
            "causal gene",
            "ancestry-specific",
        )
        for term in prohibited:
            if term.lower() in text.lower():
                errors.append(f"prohibited artwork text: {term}")
        raw = PDF.read_bytes()
        if b"/FontFile" not in raw:
            errors.append("no embedded font object")
        result = subprocess.run(
            ["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=nullpage", str(PDF)],
            check=False,
        )
        if result.returncode:
            errors.append("Ghostscript failed")

    if SOURCE.exists():
        with SOURCE.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        observed = {
            (row["trait_scope"], row["category"]): float(row["locus_weighted_fraction"])
            for row in rows
            if row["family"] == "consequence"
        }
        if set(observed) != set(EXPECTED):
            errors.append("consequence source rows do not match the eight expected cells")
        for key, expected in EXPECTED.items():
            if key in observed and not math.isclose(observed[key], expected, abs_tol=1e-12):
                errors.append(f"unexpected value for {key}: {observed[key]}")
        for trait in ("tier1_direct_masld_pdff", "tier2_liver_enzyme"):
            total = sum(value for (scope, _), value in observed.items() if scope == trait)
            if not math.isclose(total, 1.0, abs_tol=1e-12):
                errors.append(f"consequence fractions do not sum to one for {trait}: {total}")

    if MANIFEST.exists():
        with MANIFEST.open(newline="") as handle:
            manifest = {row["field"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
        if PDF.exists() and manifest.get("promoted_pdf_sha256") != sha256(PDF):
            errors.append("promoted PDF hash disagrees with manifest")
        if SOURCE.exists() and manifest.get("source_estimates_sha256") != sha256(SOURCE):
            errors.append("source-table hash disagrees with manifest")

    size_spec = (ROOT / "figures/layout_specs/figure2_panel_sizes.tsv").read_text()
    if "Fig2D_pip_architecture_by_trait_directness.pdf\t2.10\t2.44" not in size_spec:
        errors.append("Figure 2D is absent from the panel-size index")

    if errors:
        print("FAIL")
        print("\n".join(f"- {error}" for error in errors))
        return 1
    print("PASS")
    print(f"pdf_sha256\t{sha256(PDF)}")
    print("consequence_cells\t8")
    print("trait_sums\t1.0,1.0")
    print("old_fig2j_archived\tTRUE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
