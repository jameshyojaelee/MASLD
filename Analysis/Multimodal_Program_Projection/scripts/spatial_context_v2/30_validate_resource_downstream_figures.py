#!/usr/bin/env python3
"""Independently validate downstream spatial figure panels and source tables."""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from spatial_resource_lib import RESOURCE_RELEASE_ID, SpatialResourceError, sha256_file, write_tsv


HERE = Path(__file__).resolve().parent
RENDERER_PATH = HERE / "29_render_resource_downstream_figures.py"
RENDERER_SPEC = importlib.util.spec_from_file_location("downstream_renderer", RENDERER_PATH)
if RENDERER_SPEC is None or RENDERER_SPEC.loader is None:
    raise RuntimeError(f"cannot load downstream renderer: {RENDERER_PATH}")
renderer = importlib.util.module_from_spec(RENDERER_SPEC)
RENDERER_SPEC.loader.exec_module(renderer)
INTEGRATION_RELEASE_ID = renderer.INTEGRATION_RELEASE_ID
PANELS = {
    "1-spatial": ("fig1_spatial_assay_observability.pdf", "main_figure_1_candidate_inset"),
    "5-spatial": ("fig5_spatial_passport_examples.pdf", "main_figure_5_candidate_inset"),
    "S-coverage": ("figS_spatial_program_coverage.pdf", "supplement_all_117_program_coverage"),
    "S-HMSMA": ("figS_hmsma_label_blind_organization.pdf", "supplement_label_blind_array_organization"),
}
SOURCES = {
    "1-spatial": "fig1_spatial_assay_observability.tsv",
    "5-spatial": "fig5_spatial_passport_examples.tsv",
    "S-coverage": "figS_spatial_program_coverage.tsv",
    "S-HMSMA": "figS_hmsma_label_blind_organization.tsv",
}


def command_output(command: list[str]) -> str:
    return subprocess.run(command, check=True, capture_output=True, text=True).stdout


def pdf_pages(path: Path) -> int:
    for line in command_output(["pdfinfo", str(path)]).splitlines():
        if line.startswith("Pages:"):
            return int(line.split(":", 1)[1].strip())
    raise SpatialResourceError(f"pdfinfo did not report pages for {path}")


def reject_type3_fonts(path: Path) -> None:
    pdffonts = shutil.which("pdffonts")
    if pdffonts:
        output = command_output([pdffonts, str(path)])
        is_type3 = "Type 3" in output or "Type3" in output
    else:
        is_type3 = re.search(rb"/Subtype\s*/Type3\b", path.read_bytes()) is not None
    if is_type3:
        raise SpatialResourceError(f"panel contains Type 3 fonts: {path}")


def assert_exact_numeric(observed: object, expected: object, label: str) -> None:
    if not np.isclose(float(observed), float(expected), rtol=0, atol=1e-12, equal_nan=True):
        raise SpatialResourceError(f"{label} drifted: observed={observed} expected={expected}")


def validate_fig1(output: Path, upstream: Path) -> None:
    table = pd.read_csv(output / "data/fig1_spatial_assay_observability.tsv", sep="\t", dtype=str, keep_default_na=False)
    if table["dataset_id"].tolist() != renderer.COVERAGE_DATASETS:
        raise SpatialResourceError("Figure 1 datasets are missing or reordered")
    numeric = table[["n_programs", "n_observable", "n_partial", "n_untestable"]].astype(int)
    if not (numeric["n_programs"] == 117).all() or not ((numeric[["n_observable", "n_partial", "n_untestable"]].sum(axis=1)) == 117).all():
        raise SpatialResourceError("Figure 1 coverage counts do not partition 117 programs")
    registry_payload = json.loads((upstream / "registry/spatial_dataset_registry.json").read_text(encoding="utf-8"))
    registry = pd.DataFrame(registry_payload["datasets"]).set_index("dataset_id")
    for _, row in table.iterrows():
        source = registry.loc[row["dataset_id"]]
        for field in ("dataset_gate", "biological_unit_resolution", "donor_join", "phenotype_join", "histology_join", "spatial_coordinate_join"):
            if row[field] != str(source[field]):
                raise SpatialResourceError(f"Figure 1 registry field drifted: {row['dataset_id']}/{field}")
        if row["n_biological"] and row["biological_unit_resolution"] != "resolved":
            raise SpatialResourceError(f"Figure 1 assigns a biological count to unresolved source {row['dataset_id']}")


def validate_fig5(output: Path, upstream: Path) -> None:
    table = pd.read_csv(output / "data/fig5_spatial_passport_examples.tsv", sep="\t")
    if table["primary_evidence_state"].tolist() != ["supported", "indeterminate", "untestable"]:
        raise SpatialResourceError("Figure 5 must show supported, indeterminate, and untestable examples")
    effects = pd.read_parquet(upstream / "effects/spatial_program_effects.parquet")
    mapping = {
        "igfbp7_program": "hotspot_hepatocytes_f05c535ae5bbc0b9",
        "bicc1_program": "hotspot_hepatocytes_48f39dd4d817a10e",
    }
    for example, uid in mapping.items():
        observed = table[table["example_id"] == example].iloc[0]
        for prefix, dataset in (("gse", "GSE192741"), ("vu", "Vu_et_al_2025")):
            expected = effects[(effects["dataset_id"] == dataset) & (effects["program_uid"] == uid)].iloc[0]
            assert_exact_numeric(observed[f"{prefix}_estimate"], expected["estimate"], f"{example}/{dataset}/estimate")
            assert_exact_numeric(observed[f"{prefix}_qvalue"], expected["qvalue"], f"{example}/{dataset}/qvalue")
    hmsma = pd.read_csv(upstream / "hmsma_label_blind/per_array_program_organization.tsv", sep="\t")
    expected_median = hmsma[hmsma["program_uid"] == mapping["igfbp7_program"]]["residual_moran_i"].median()
    observed_median = table[table["example_id"] == "hmsma_clinical_gate"].iloc[0]["hmsma_median_moran_i"]
    assert_exact_numeric(observed_median, expected_median, "HMSMA passport median")
    if not table[table["example_id"] == "hmsma_clinical_gate"].iloc[0]["inference_unit"].startswith("35 arrays; donor key absent"):
        raise SpatialResourceError("Figure 5 HMSMA example does not preserve the unit gate")


def validate_coverage(output: Path, upstream: Path) -> None:
    observed = pd.read_csv(output / "data/figS_spatial_program_coverage.tsv", sep="\t", dtype=str)
    if len(observed) != 819 or observed["dataset_id"].nunique() != 7 or observed["program_uid"].nunique() != 117:
        raise SpatialResourceError("supplement coverage source is not the complete 117-by-7 family")
    banned = {"pvalue", "qvalue", "evidence_state", "validation_label"}.intersection(observed.columns)
    if banned:
        raise SpatialResourceError(f"outcome or validation fields leaked into coverage display: {sorted(banned)}")
    expected = pd.read_parquet(upstream / "coverage/spatial_program_coverage.parquet")
    expected = expected[expected["dataset_id"].isin(renderer.COVERAGE_DATASETS)]
    keys = ["dataset_id", "program_uid"]
    merged = observed.merge(expected[keys + ["coverage_status"]], on=keys, how="outer", suffixes=("_observed", "_expected"), indicator=True)
    if not (merged["_merge"] == "both").all() or not (merged["coverage_status_observed"] == merged["coverage_status_expected"]).all():
        raise SpatialResourceError("supplement coverage states drifted from the sealed table")


def validate_hmsma(output: Path, upstream: Path) -> None:
    observed = pd.read_csv(output / "data/figS_hmsma_label_blind_organization.tsv", sep="\t", dtype=str)
    if len(observed) != 70 or observed["array_id"].nunique() != 35 or observed["program_uid"].nunique() != 2:
        raise SpatialResourceError("HMSMA panel source is not 35 arrays by two programs")
    if set(observed["inferential_pvalue_authorized"]) != {"FALSE"}:
        raise SpatialResourceError("HMSMA panel authorizes inferential P values")
    if any(column.lower() in {"pvalue", "qvalue", "population_pvalue"} for column in observed.columns):
        raise SpatialResourceError("HMSMA label-blind panel contains inferential P-value fields")
    expected = pd.read_csv(upstream / "hmsma_label_blind/per_array_program_organization.tsv", sep="\t", dtype=str)
    keys = ["array_id", "program_uid"]
    merged = observed.merge(expected, on=keys, how="outer", suffixes=("_observed", "_expected"), indicator=True)
    if not (merged["_merge"] == "both").all():
        raise SpatialResourceError("HMSMA panel array/program family drifted from sealed output")
    numeric_fields = {
        "n_spots", "n_graph_eligible_spots", "n_directed_edges", "n_genes_measured",
        "retained_l1_weight", "residual_moran_i",
    }
    for field in expected.columns:
        if field in keys:
            continue
        left = merged[f"{field}_observed"]
        right = merged[f"{field}_expected"]
        if field in numeric_fields:
            if not np.allclose(left.astype(float), right.astype(float), rtol=0, atol=1e-12):
                raise SpatialResourceError(f"HMSMA numeric source field drifted: {field}")
        elif field == "inferential_pvalue_authorized":
            if not left.str.lower().isin({"false", "0"}).all() or not right.str.lower().isin({"false", "0"}).all():
                raise SpatialResourceError("HMSMA authorization flag drifted")
        elif not (left == right).all():
            raise SpatialResourceError(f"HMSMA categorical source field drifted: {field}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    upstream = args.upstream_root.resolve()
    output = args.output_root.resolve()
    if not (output / "BUILD_COMPLETE").is_file() or (output / "READY").exists():
        raise SpatialResourceError("downstream candidate is absent, incomplete, or already immutable")
    renderer.require_ready(upstream)
    validate_fig1(output, upstream)
    validate_fig5(output, upstream)
    validate_coverage(output, upstream)
    validate_hmsma(output, upstream)

    panel_rows, source_rows = [], []
    for callout, (filename, role) in PANELS.items():
        panel = output / "panels" / filename
        if not panel.is_file() or panel.stat().st_size == 0:
            raise SpatialResourceError(f"missing downstream panel: {panel}")
        pages = pdf_pages(panel)
        if pages != 1:
            raise SpatialResourceError(f"downstream panel has {pages} pages: {panel}")
        reject_type3_fonts(panel)
        panel_rows.append({
            "integration_release_id": INTEGRATION_RELEASE_ID,
            "spatial_release_id": RESOURCE_RELEASE_ID,
            "callout": callout,
            "relative_path": panel.relative_to(output).as_posix(),
            "role": role,
            "bytes": panel.stat().st_size,
            "sha256": sha256_file(panel),
            "pdf_pages": pages,
            "status": "candidate_not_promoted",
        })
        source = output / "data" / SOURCES[callout]
        source_rows.append({
            "integration_release_id": INTEGRATION_RELEASE_ID,
            "spatial_release_id": RESOURCE_RELEASE_ID,
            "callout": callout,
            "relative_path": source.relative_to(output).as_posix(),
            "bytes": source.stat().st_size,
            "sha256": sha256_file(source),
        })
    write_tsv(output / "panel_manifest.tsv", tuple(panel_rows[0]), panel_rows)
    write_tsv(output / "source_table_manifest.tsv", tuple(source_rows[0]), source_rows)
    write_tsv(
        output / "READY",
        ("integration_release_id", "spatial_release_id", "status", "n_panels", "n_source_tables", "panel_manifest_sha256", "source_manifest_sha256", "canonical_figure_written"),
        [{
            "integration_release_id": INTEGRATION_RELEASE_ID,
            "spatial_release_id": RESOURCE_RELEASE_ID,
            "status": "validated_downstream_figure_candidate_awaiting_adjudication",
            "n_panels": len(panel_rows),
            "n_source_tables": len(source_rows),
            "panel_manifest_sha256": sha256_file(output / "panel_manifest.tsv"),
            "source_manifest_sha256": sha256_file(output / "source_table_manifest.tsv"),
            "canonical_figure_written": "FALSE",
        }],
    )
    print("PASS four one-page downstream spatial panels and exact source tables")


if __name__ == "__main__":
    main()
