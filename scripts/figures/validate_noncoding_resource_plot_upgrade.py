#!/usr/bin/env python3
"""Independent structural and scientific validator for the noncoding plot upgrade."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path


EXPECTED_ROWS = {
    "fig3_pooled_deg_biotype_v2_source.tsv": 6,
    "fig3_stage_deg_biotype_v2_source.tsv": 24,
    "figS_lncrna_universe_composition_source.tsv": 6,
    "figS_lncrna_robustness_heatmap_v2_source.tsv": 271 * 11,
    "figS_lncrna_genomic_classes_v2_source.tsv": 5,
    "figS_lncrna_example_forest_v2_source.tsv": 11,
    "figS_program_lncrna_content_source.tsv": 117,
    "figS_program_without_lncrna_sensitivity_source.tsv": 49,
    "fig5_noncoding_object_routes_source.tsv": 5,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate(source: Path, review: Path) -> None:
    if review.exists():
        raise FileExistsError(f"refusing to overwrite {review}")
    review.mkdir(parents=True)

    manifest = rows(source / "output_manifest.tsv")
    for row in manifest:
        path = source / row["relative_path"]
        assert path.stat().st_size == int(row["size_bytes"])
        assert sha256(path) == row["sha256"]

    pdfs = sorted(source.glob("*.pdf"))
    assert len(pdfs) == 9
    for pdf in pdfs:
        info = subprocess.check_output(["pdfinfo", str(pdf)], text=True)
        assert "Pages:           1" in info
        subprocess.run(
            ["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=nullpage", str(pdf)],
            check=True,
        )
        subprocess.run(
            [
                "pdftoppm",
                "-png",
                "-r",
                "180",
                "-singlefile",
                str(pdf),
                str(review / pdf.stem),
            ],
            check=True,
        )

    tables = {name: rows(source / name) for name in EXPECTED_ROWS}
    for name, expected in EXPECTED_ROWS.items():
        assert len(tables[name]) == expected, (name, len(tables[name]), expected)

    pooled = tables["fig3_pooled_deg_biotype_v2_source.tsv"]
    assert sum(int(row["n_significant"]) for row in pooled) == 1616
    composition = tables["figS_lncrna_universe_composition_source.tsv"]
    lnc = {
        row["universe"]: float(row["percent"])
        for row in composition
        if row["display_biotype"] == "lncRNA"
    }
    assert abs(lnc["All tested genes"] - 26.7394) < 0.001
    assert abs(lnc["TREAT-positive genes"] - 26.8564) < 0.001
    sensitivity = tables["figS_program_without_lncrna_sensitivity_source.tsv"]
    assert sum(row["sensitivity_passed"].lower() == "true" for row in sensitivity) == 48
    assert (
        sum(row["stage_direction_preserved_y"].lower() == "true" for row in sensitivity)
        == 48
    )
    assert {
        row["rule_id"] for row in tables["fig5_noncoding_object_routes_source.tsv"]
    } == {
        "EXP_REGULATORY_DNA_V2",
        "EXP_LNCRNA_RNA_PRODUCT_V2",
        "EXP_LNCRNA_LOCUS_DISAMBIGUATION_V2",
        "EXP_PROTEIN_STATE_CONTEXT_V2",
        "EXP_MEASURE_MISSING_ASSAY_V2",
    }

    verdict = {
        "status": "pass",
        "source": str(source.resolve()),
        "pdf_count": len(pdfs),
        "source_table_count": len(tables),
        "lncrna_share_tested_percent": lnc["All tested genes"],
        "lncrna_share_treat_positive_percent": lnc["TREAT-positive genes"],
        "program_sensitivity_passed": 48,
        "program_sensitivity_tested": 49,
        "contains_corrected_genetics_outcomes": False,
        "source_manifest_sha256": sha256(source / "output_manifest.tsv"),
    }
    (review / "VALIDATED.json").write_text(
        json.dumps(verdict, indent=2, sort_keys=True) + "\n"
    )
    with (review / "validation.tsv").open("w") as handle:
        handle.write("check\tstatus\n")
        for check in (
            "manifest_size_and_sha256",
            "nine_one_page_pdfs",
            "ghostscript_parse",
            "png_previews",
            "source_table_cardinality",
            "validated_scientific_censuses",
            "catalog_route_family",
        ):
            handle.write(f"{check}\tPASS\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--review", type=Path, required=True)
    args = parser.parse_args()
    validate(args.source.resolve(), args.review.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
