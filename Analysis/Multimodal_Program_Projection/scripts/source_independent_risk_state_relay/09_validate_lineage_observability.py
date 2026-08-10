#!/usr/bin/env python3
"""Validate the independent donor-collapsed Plan 45 lineage reference."""

from __future__ import annotations

import csv
import json
import math
import os
from collections import Counter
from pathlib import Path

from relay_common import PROJECT_ROOT, sha256_file


CANDIDATE_ID = os.environ.get(
    "PLAN45_LINEAGE_CANDIDATE_ID",
    "source-independent-risk-state-relay-lineage-observability-2026-08-10",
).strip()
CANDIDATE_ROOT = (
    PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates" / CANDIDATE_ID
)
EXPECTED_CELL_TYPES = {
    "B_cells",
    "Basophils",
    "Cholangiocytes",
    "Circulating_NK_NKT",
    "Endothelial_cells",
    "Fibroblasts",
    "Hepatocytes",
    "Macrophages",
    "Mono+mono_derived_cells",
    "Neutrophils",
    "Plasma_cells",
    "Resident_NK",
    "T_cells",
    "cDC1s",
    "cDC2s",
    "pDCs",
}


def rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def main() -> None:
    seal_path = CANDIDATE_ROOT / "LINEAGE_REFERENCE_SEALED.json"
    long_path = CANDIDATE_ROOT / "donor_lineage_observability.tsv"
    summary_path = CANDIDATE_ROOT / "donor_lineage_gene_summary.tsv"
    audit_path = CANDIDATE_ROOT / "lineage_source_audit.tsv"
    manifest_path = CANDIDATE_ROOT / "lineage_source_manifest.tsv"
    required = [seal_path, long_path, summary_path, audit_path, manifest_path]
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing lineage-reference artifact: {path}")
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("experimental_targets_frozen") is not False:
        raise RuntimeError("Lineage reference improperly claims frozen targets")
    if seal.get("disease_or_genetic_outcomes_inspected") is not False:
        raise RuntimeError("Lineage reference indicates outcome-driven construction")
    expected_hashes = seal["output_sha256"]
    for key, path in {
        "donor_lineage_observability": long_path,
        "donor_lineage_gene_summary": summary_path,
        "lineage_source_audit": audit_path,
        "lineage_source_manifest": manifest_path,
    }.items():
        if sha256_file(path) != expected_hashes[key]:
            raise RuntimeError(f"Lineage-reference hash mismatch: {key}")

    audit = list(rows(audit_path))
    if {row["cell_type"] for row in audit} != EXPECTED_CELL_TYPES:
        raise RuntimeError("Lineage source audit differs from the frozen 16-cell-type inventory")
    if any(int(row["n_usable_donors"]) < 5 for row in audit):
        raise RuntimeError("Lineage source contains a cell type with fewer than five donors")
    if sum(int(row["n_collapsed_columns"]) for row in audit) <= 0:
        raise RuntimeError("Lineage source did not actually collapse any repeated runs")

    long_count = 0
    cells = set()
    genes = set()
    for row in rows(long_path):
        long_count += 1
        cells.add(row["cell_type"])
        genes.add(row["gene_symbol"])
        n_donors = int(row["n_donors"])
        detect = float(row["donor_detection_fraction_cpm1"])
        if n_donors < 5 or not 0 <= detect <= 1:
            raise RuntimeError(f"Invalid long lineage row: {row['gene_symbol']} / {row['cell_type']}")
    if cells != EXPECTED_CELL_TYPES or len(genes) != int(seal["n_genes"]):
        raise RuntimeError("Lineage long-table dimensions differ from the seal")
    if long_count != int(seal["n_long_rows"]):
        raise RuntimeError("Lineage long-table row count differs from the seal")

    summary_count = 0
    labels = Counter()
    for row in rows(summary_path):
        summary_count += 1
        labels[row["interpretation"]] += 1
        dominant = float(row["dominant_mean_cpm"])
        detected = float(row["dominant_detection_fraction_cpm1"])
        expression_expected = dominant >= 1 and detected >= 0.25
        if (row["lineage_expression_gate"].lower() == "true") != expression_expected:
            raise RuntimeError(f"Lineage expression gate mismatch: {row['gene_symbol']}")
        ratio_text = row["dominant_to_second_ratio"].strip()
        ratio = float(ratio_text) if ratio_text else float("nan")
        specificity_expected = expression_expected and not math.isnan(ratio) and ratio >= 2
        if (row["lineage_specificity_gate"].lower() == "true") != specificity_expected:
            raise RuntimeError(f"Lineage specificity gate mismatch: {row['gene_symbol']}")
    if summary_count != len(genes):
        raise RuntimeError("Lineage summary does not contain one row per gene")

    print(
        "Plan 45 donor-lineage validation passed: "
        f"genes={len(genes)}; celltypes={len(cells)}; long_rows={long_count}; "
        f"labels={dict(labels)}; targets not frozen"
    )


if __name__ == "__main__":
    main()
