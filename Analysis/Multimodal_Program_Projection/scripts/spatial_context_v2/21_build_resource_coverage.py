#!/usr/bin/env python3
"""Build outcome-free 117-program coverage and gene-context products."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from spatial_resource_lib import (
    COVERAGE_COLUMNS,
    DETECTION_MIN_FRACTION,
    EXPECTED_REGISTRY_SHA256,
    PROGRAM_RELEASE_ID,
    RESOURCE_RELEASE_ID,
    SpatialResourceError,
    coverage_state,
    h5ad_gene_sets,
    load_dataset_registry,
    load_frozen_programs,
    parse_bool,
    require_complete_coverage,
    sha256_file,
    write_tsv,
)


GENE_CONTEXT_COLUMNS = (
    "release_id", "program_release_id", "dataset_id", "assay_id", "gene_symbol",
    "measurement_scope", "measured", "detected", "testable", "testability_reason",
    "assay_native_moran_i", "assay_native_svg_statistic", "assay_native_zonation",
    "effect_unit", "evidence_state", "dataset_gate", "source_dependence",
    "biological_unit", "biological_unit_resolution", "provenance",
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--candidate-root", type=Path, required=True)
    args = parser.parse_args()
    project = args.project_root.resolve()
    candidate = args.candidate_root.resolve()
    config = Path(__file__).with_name("spatial_resource_datasets.tsv")
    registry = load_dataset_registry(config)
    hotspot = project / "Analysis/Multimodal_Program_Projection/candidates" / PROGRAM_RELEASE_ID / "hotspot"
    programs, weights = load_frozen_programs(hotspot)
    output = candidate / "coverage"
    if (output / "READY").exists():
        raise SpatialResourceError(f"immutable coverage candidate already exists: {output}")

    coverage_rows = []
    gene_rows = []
    gene_union = sorted({gene for program_weights in weights.values() for gene in program_weights})
    source_audit = []
    for dataset in registry:
        if not parse_bool(dataset["coverage_required"], dataset["dataset_id"]):
            continue
        relative = dataset["gene_axis_source"].strip()
        if relative:
            source = project / relative
            measured, detected = h5ad_gene_sets(source, dataset["detection_matrix"])
            source_sha256 = sha256_file(source)
        else:
            measured, detected = set(), set()
            source_sha256 = ""
        source_audit.append({
            "dataset_id": dataset["dataset_id"],
            "assay_id": dataset["assay_id"],
            "gene_universe_denominator": len(measured),
            "n_detected_genes": len(detected),
            "min_detection_fraction": DETECTION_MIN_FRACTION,
            "gene_axis_source": relative,
            "detection_matrix": dataset["detection_matrix"],
            "gene_axis_source_sha256": source_sha256,
        })
        for uid, program in programs.items():
            program_weights = weights[uid]
            genes_on_axis = set(program_weights) & measured
            detected_genes = set(program_weights) & detected
            retained = sum(program_weights[gene] for gene in detected_genes)
            status, reason = coverage_state(len(detected_genes), retained)
            coverage_rows.append({
                "release_id": RESOURCE_RELEASE_ID,
                "program_release_id": PROGRAM_RELEASE_ID,
                "registry_sha256": EXPECTED_REGISTRY_SHA256,
                "dataset_id": dataset["dataset_id"],
                "assay_id": dataset["assay_id"],
                "program_uid": uid,
                "program_label": program["module_name"],
                "cell_type": program["cell_type"],
                "membership_sha256": program["membership_sha256"],
                "n_program_genes": len(program_weights),
                "n_genes_on_axis": len(genes_on_axis),
                "n_genes_measured": len(detected_genes),
                "n_genes_detected": len(detected_genes),
                "retained_l1_weight": f"{retained:.17g}",
                "detection_fraction": f"{len(detected_genes) / len(program_weights):.17g}",
                "gene_universe_denominator": len(measured),
                "lineage_context": dataset["lineage_context"],
                "coverage_status": status,
                "testability_reason": reason,
                "dataset_gate": dataset["dataset_gate"],
                "biological_unit": dataset["biological_unit"],
                "biological_unit_resolution": dataset["biological_unit_resolution"],
                "source_dependence": dataset["source_dependence"],
            })
        for gene in gene_union:
            is_measured = gene in measured
            is_detected = gene in detected
            gene_rows.append({
                "release_id": RESOURCE_RELEASE_ID,
                "program_release_id": PROGRAM_RELEASE_ID,
                "dataset_id": dataset["dataset_id"],
                "assay_id": dataset["assay_id"],
                "gene_symbol": gene,
                "measurement_scope": "frozen_117_program_gene_union",
                "measured": is_measured,
                "detected": is_detected,
                "testable": False,
                "testability_reason": "gene_level_inference_not_predeclared" if is_measured else "gene_not_measured_or_axis_unavailable",
                "assay_native_moran_i": None,
                "assay_native_svg_statistic": None,
                "assay_native_zonation": None,
                "effect_unit": dataset["abundance_unit"],
                "evidence_state": "indeterminate" if is_measured else "untestable",
                "dataset_gate": dataset["dataset_gate"],
                "source_dependence": dataset["source_dependence"],
                "biological_unit": dataset["biological_unit"],
                "biological_unit_resolution": dataset["biological_unit_resolution"],
                "provenance": relative or "source_gate_no_gene_axis_read",
            })

    require_complete_coverage(coverage_rows, registry)
    output.mkdir(parents=True, exist_ok=True)
    write_tsv(output / "spatial_program_coverage.tsv", COVERAGE_COLUMNS, coverage_rows)
    pd.DataFrame(coverage_rows, columns=COVERAGE_COLUMNS).to_parquet(output / "spatial_program_coverage.parquet", index=False)
    pd.DataFrame(gene_rows, columns=GENE_CONTEXT_COLUMNS).to_parquet(output / "spatial_gene_context.parquet", index=False)
    write_tsv(
        output / "gene_axis_audit.tsv",
        ("dataset_id", "assay_id", "gene_universe_denominator", "n_detected_genes", "min_detection_fraction", "gene_axis_source", "detection_matrix", "gene_axis_source_sha256"),
        source_audit,
    )
    write_tsv(
        output / "READY",
        ("release_id", "status", "n_programs", "n_coverage_assays", "n_coverage_rows", "n_gene_context_rows", "coverage_sha256", "gene_context_sha256"),
        [{
            "release_id": RESOURCE_RELEASE_ID,
            "status": "complete_outcome_free_coverage_no_confirmatory_inference",
            "n_programs": len(programs),
            "n_coverage_assays": len(source_audit),
            "n_coverage_rows": len(coverage_rows),
            "n_gene_context_rows": len(gene_rows),
            "coverage_sha256": sha256_file(output / "spatial_program_coverage.parquet"),
            "gene_context_sha256": sha256_file(output / "spatial_gene_context.parquet"),
        }],
    )
    print(json.dumps({"coverage_rows": len(coverage_rows), "gene_context_rows": len(gene_rows)}, indent=2))


if __name__ == "__main__":
    main()
