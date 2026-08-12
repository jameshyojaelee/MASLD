#!/usr/bin/env python3
"""Independently validate spatial unit, matched-null, and sensitivity diagnostics."""

from __future__ import annotations

import argparse
import importlib.util
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from spatial_resource_lib import PROGRAM_RELEASE_ID, SpatialResourceError, sha256_file, write_tsv


HERE = Path(__file__).resolve().parent
RENDERER_PATH = HERE / "33_render_unit_null_diagnostics.py"
SPEC = importlib.util.spec_from_file_location("unit_null_renderer", RENDERER_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load diagnostic renderer: {RENDERER_PATH}")
renderer = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(renderer)

PANELS = {
    "S-unit": ("figS_spatial_unit_heterogeneity.pdf", "figS_spatial_unit_heterogeneity.tsv"),
    "S-null": ("figS_spatial_matched_null_calibration.pdf", "figS_spatial_matched_null_calibration.tsv"),
    "S-sensitivity": ("figS_spatial_weight_sensitivity.pdf", "figS_spatial_weight_sensitivity.tsv"),
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


def assert_close(left: object, right: object, label: str) -> None:
    if not np.isclose(float(left), float(right), rtol=0, atol=1e-12, equal_nan=True):
        raise SpatialResourceError(f"{label} drifted: observed={left} expected={right}")


def validate_upstream_manifest(output: Path, upstream: Path) -> None:
    manifest = pd.read_csv(output / "upstream_source_manifest.tsv", sep="\t", dtype=str, keep_default_na=False)
    expected_files = {
        "READY", "spatial_section_results.tsv", "spatial_program_results.tsv",
        "spatial_null_summary.tsv", "spatial_sensitivity_audit.tsv",
    }
    if {Path(path).name for path in manifest["source_path"]} != expected_files or len(manifest) != len(expected_files):
        raise SpatialResourceError("upstream source manifest is incomplete")
    for row in manifest.itertuples(index=False):
        path = Path(row.source_path)
        if path.parent.resolve() != upstream or not path.is_file():
            raise SpatialResourceError(f"manifested upstream source is outside the sealed candidate: {path}")
        if str(path.stat().st_size) != row.bytes or sha256_file(path) != row.sha256:
            raise SpatialResourceError(f"manifested upstream source drifted: {path}")


def validate_unit_table(output: Path, upstream: Path) -> None:
    observed = pd.read_csv(output / "data/figS_spatial_unit_heterogeneity.tsv", sep="\t", dtype=str, keep_default_na=False)
    expected = pd.read_csv(upstream / "spatial_section_results.tsv", sep="\t", dtype=str, keep_default_na=False)
    if len(observed) != 30 or observed["program_id"].nunique() != 2:
        raise SpatialResourceError("unit diagnostic must contain 30 physical-unit rows across two programs")
    if set(observed[observed["dataset"] == "GSE192741"]["sample_id"]) != set(renderer.GSE_SECTION_ORDER):
        raise SpatialResourceError("unit diagnostic does not contain all five GSE192741 sections")
    if set(observed[observed["dataset"] == "Vu_et_al_2025"]["sample_id"]) != set(renderer.VU_ARRAY_ORDER):
        raise SpatialResourceError("unit diagnostic does not contain all ten Vu physical arrays")
    h35 = observed[(observed["dataset"] == "GSE192741") & (observed["aggregation_group"] == "H35")]
    if set(h35["sample_id"]) != {"JBO014", "JBO015"} or len(h35) != 4:
        raise SpatialResourceError("H35 two-section collapse is not explicit for both programs")
    if observed["inferential_pvalue_authorized"].ne("FALSE").any():
        raise SpatialResourceError("unit diagnostic authorizes a physical-unit P value")
    if any(column.lower() in {"pvalue", "qvalue", "padj"} for column in observed.columns):
        raise SpatialResourceError("physical-unit source contains inferential P-value fields")
    keys = ["program_id", "dataset", "sample_id"]
    merged = observed.merge(expected, on=keys, how="outer", suffixes=("_observed", "_expected"), indicator=True)
    if not merged["_merge"].eq("both").all():
        raise SpatialResourceError("physical-unit source rows drifted from the sealed v2 results")
    for field in ("individual", "condition"):
        if not merged[f"{field}_observed"].eq(merged[f"{field}_expected"]).all():
            raise SpatialResourceError(f"physical-unit categorical field drifted: {field}")
    for field in ("raw_moran_i", "residual_moran_i", "zonation_moran_i"):
        if not np.allclose(
            merged[f"{field}_observed"].astype(float), merged[f"{field}_expected"].astype(float), rtol=0, atol=1e-12
        ):
            raise SpatialResourceError(f"physical-unit numeric field drifted: {field}")
    numeric = observed.assign(
        residual_moran_i=observed["residual_moran_i"].astype(float),
        collapsed_residual_moran_i=observed["collapsed_residual_moran_i"].astype(float),
    )
    rederived = numeric.groupby(["dataset", "program_id", "aggregation_group"], observed=True)["residual_moran_i"].transform("mean")
    if not np.allclose(numeric["collapsed_residual_moran_i"], rederived, rtol=0, atol=1e-12):
        raise SpatialResourceError("physical-unit collapse values do not rederive from section/array values")
    gse = numeric[numeric["dataset"] == "GSE192741"]
    if gse["aggregation_group"].nunique() != 4:
        raise SpatialResourceError("GSE192741 does not collapse to exactly four donors")


def validate_null_table(output: Path, upstream: Path) -> None:
    observed = pd.read_csv(output / "data/figS_spatial_matched_null_calibration.tsv", sep="\t")
    programs = pd.read_csv(upstream / "spatial_program_results.tsv", sep="\t")
    null = pd.read_csv(upstream / "spatial_null_summary.tsv", sep="\t")
    sensitivity = pd.read_csv(upstream / "spatial_sensitivity_audit.tsv", sep="\t")
    if len(observed) != 4 or set(zip(observed["dataset"], observed["program_id"])) != {
        (dataset, uid) for dataset in renderer.DATASETS for uid in renderer.PROGRAMS
    }:
        raise SpatialResourceError("matched-null source is not the complete two-program by two-dataset family")
    if not observed["uncertainty_semantics"].eq("matched_gene_null_quantiles_not_sampling_confidence_interval").all():
        raise SpatialResourceError("matched-null dispersion is mislabeled as sampling uncertainty")
    for row in observed.itertuples(index=False):
        program = programs[(programs["dataset"] == row.dataset) & (programs["program_id"] == row.program_id)].iloc[0]
        null_row = null[(null["dataset"] == row.dataset) & (null["program_id"] == row.program_id) & (null["statistic"] == "residual_moran_i")].iloc[0]
        primary = sensitivity[
            (sensitivity["dataset"] == row.dataset)
            & (sensitivity["program_id"] == row.program_id)
            & (sensitivity["sensitivity_id"] == "primary_original_l1_weight")
        ].iloc[0]
        for field in ("residual_moran_i", "residual_pvalue", "residual_padj", "residual_null_mean", "residual_null_sd"):
            assert_close(getattr(row, field), program[field], f"{row.dataset}/{row.program_id}/{field}")
        for field in ("null_mean", "null_sd", "q010", "q050", "q500", "q950", "q990"):
            assert_close(getattr(row, field), null_row[field], f"{row.dataset}/{row.program_id}/{field}")
        centered = float(row.residual_moran_i) - float(row.null_mean)
        assert_close(row.centered_residual_moran_i, centered, f"{row.dataset}/{row.program_id}/centered")
        assert_close(row.sensitivity_primary_centered_residual_moran_i, primary["centered_residual_moran_i"], f"{row.dataset}/{row.program_id}/primary")
        assert_close(centered, primary["centered_residual_moran_i"], f"{row.dataset}/{row.program_id}/independent centered agreement")
        if not float(row.q010) <= float(row.q050) <= float(row.q500) <= float(row.q950) <= float(row.q990):
            raise SpatialResourceError(f"matched-null quantiles are not monotonic for {row.dataset}/{row.program_id}")
        expected_counts = (4, 5) if row.dataset == "GSE192741" else (10, 10)
        if (int(row.n_aggregation_units), int(row.n_technical_units)) != expected_counts:
            raise SpatialResourceError(f"unit counts drifted for {row.dataset}/{row.program_id}")
        if row.dataset == "Vu_et_al_2025":
            source_dependent = str(row.source_dependent).strip().lower()
            if source_dependent not in {"true", "1"}:
                raise SpatialResourceError("Vu matched-null row lost source dependence")


def validate_sensitivity_table(output: Path, upstream: Path) -> None:
    observed = pd.read_csv(output / "data/figS_spatial_weight_sensitivity.tsv", sep="\t", dtype=str, keep_default_na=False)
    expected = pd.read_csv(upstream / "spatial_sensitivity_audit.tsv", sep="\t", dtype=str, keep_default_na=False)
    if len(observed) != 12 or set(observed["sensitivity_id"]) != set(renderer.SENSITIVITY_ORDER):
        raise SpatialResourceError("weight-sensitivity source is not the exact 4-by-3 family")
    if observed["inferential_pvalue_authorized"].ne("FALSE").any():
        raise SpatialResourceError("weight-sensitivity source authorizes a new P value")
    if any(column.lower() in {"pvalue", "qvalue", "padj"} for column in observed.columns):
        raise SpatialResourceError("weight-sensitivity source contains inferential P-value fields")
    keys = ["release_id", "registry_version", "dataset", "program_id", "membership_sha256", "sensitivity_id"]
    merged = observed.merge(expected, on=keys, how="outer", suffixes=("_observed", "_expected"), indicator=True)
    if not merged["_merge"].eq("both").all():
        raise SpatialResourceError("weight-sensitivity family drifted from the sealed audit")
    if not np.allclose(
        merged["centered_residual_moran_i_observed"].astype(float),
        merged["centered_residual_moran_i_expected"].astype(float), rtol=0, atol=1e-12,
    ):
        raise SpatialResourceError("weight-sensitivity estimates drifted from the sealed audit")
    if not merged["direction_observed"].eq(merged["direction_expected"]).all():
        raise SpatialResourceError("weight-sensitivity field drifted: direction")
    observed_sign = merged["sign_agree_with_primary_observed"].str.strip().str.lower()
    expected_sign = merged["sign_agree_with_primary_expected"].str.strip().str.lower()
    true_values = {"true", "1"}
    false_values = {"false", "0"}
    if not observed_sign.isin(true_values | false_values).all() or not expected_sign.isin(true_values | false_values).all():
        raise SpatialResourceError("weight-sensitivity sign flag is not boolean")
    if not observed_sign.isin(true_values).eq(expected_sign.isin(true_values)).all():
        raise SpatialResourceError("weight-sensitivity field drifted: sign_agree_with_primary")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    upstream = args.upstream_root.resolve()
    output = args.output_root.resolve()
    if not (output / "BUILD_COMPLETE").is_file() or (output / "READY").exists():
        raise SpatialResourceError("diagnostic candidate is absent, incomplete, or already immutable")
    renderer.require_upstream(upstream)
    validate_upstream_manifest(output, upstream)
    validate_unit_table(output, upstream)
    validate_null_table(output, upstream)
    validate_sensitivity_table(output, upstream)

    panel_rows = []
    source_rows = []
    for callout, (panel_name, source_name) in PANELS.items():
        panel = output / "panels" / panel_name
        source = output / "data" / source_name
        if not panel.is_file() or panel.stat().st_size == 0 or not source.is_file() or source.stat().st_size == 0:
            raise SpatialResourceError(f"missing panel or exact source table for {callout}")
        pages = pdf_pages(panel)
        if pages != 1:
            raise SpatialResourceError(f"diagnostic panel has {pages} pages: {panel}")
        reject_type3_fonts(panel)
        panel_rows.append({
            "diagnostic_release_id": renderer.DIAGNOSTIC_RELEASE_ID,
            "upstream_release_id": PROGRAM_RELEASE_ID,
            "callout": callout,
            "relative_path": panel.relative_to(output).as_posix(),
            "bytes": panel.stat().st_size,
            "sha256": sha256_file(panel),
            "pdf_pages": pages,
            "status": "candidate_not_promoted",
        })
        source_rows.append({
            "diagnostic_release_id": renderer.DIAGNOSTIC_RELEASE_ID,
            "upstream_release_id": PROGRAM_RELEASE_ID,
            "callout": callout,
            "relative_path": source.relative_to(output).as_posix(),
            "bytes": source.stat().st_size,
            "sha256": sha256_file(source),
        })
    write_tsv(output / "panel_manifest.tsv", tuple(panel_rows[0]), panel_rows)
    write_tsv(output / "source_table_manifest.tsv", tuple(source_rows[0]), source_rows)
    write_tsv(
        output / "READY",
        (
            "diagnostic_release_id", "upstream_release_id", "status", "n_panels", "n_source_tables",
            "panel_manifest_sha256", "source_manifest_sha256", "canonical_figure_written",
        ),
        [{
            "diagnostic_release_id": renderer.DIAGNOSTIC_RELEASE_ID,
            "upstream_release_id": PROGRAM_RELEASE_ID,
            "status": "validated_spatial_unit_null_diagnostic_candidate_awaiting_adjudication",
            "n_panels": len(panel_rows),
            "n_source_tables": len(source_rows),
            "panel_manifest_sha256": sha256_file(output / "panel_manifest.tsv"),
            "source_manifest_sha256": sha256_file(output / "source_table_manifest.tsv"),
            "canonical_figure_written": "FALSE",
        }],
    )
    print("PASS three one-page spatial unit/null diagnostic panels and exact source tables")


if __name__ == "__main__":
    main()
