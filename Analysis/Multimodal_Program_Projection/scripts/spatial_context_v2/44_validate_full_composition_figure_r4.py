#!/usr/bin/env python3
"""Validate r4 composition figure values, rendering, and full provenance."""

from __future__ import annotations

import argparse
import importlib.util
import re
import subprocess
import sys
import tempfile
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


def import_script(name: str, module_name: str):
    path = SCRIPT_DIR / name
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import script: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# Importing the accepted r3 validator traverses 42 -> 41 -> 39 and 40.  The
# r4 producer traverses 43 -> 41 -> 39.  The release manifest records and this
# validator rechecks every member of that executable dependency graph.
R3_VALIDATOR = import_script("42_validate_full_composition_figure_r3.py", "accepted_r3_validator_loader")
R4 = import_script("43_render_full_composition_figure_r4.py", "full_composition_figure_r4")
HELPERS = R3_VALIDATOR.VALIDATOR


def raster_sha256(pdf: Path, temporary_root: Path) -> str:
    stem = temporary_root / pdf.parent.name
    subprocess.run(
        ["pdftoppm", "-png", "-r", "200", "-singlefile", str(pdf), str(stem)],
        check=True,
        capture_output=True,
        text=True,
    )
    return R4.sha256_file(stem.with_suffix(".png"))


def validate(project_root: Path, input_root: Path | None, output_root: Path | None):
    import pandas as pd

    root = project_root.resolve()
    input_path = (input_root or R4.default_input(root)).resolve()
    output_path = (output_root or R4.default_output(root)).resolve()
    if not output_path.is_dir():
        raise RuntimeError(f"r4 figure candidate is missing: {output_path}")
    if (output_path / "READY").exists() or (output_path / "validation_report.tsv").exists():
        raise RuntimeError("refusing to overwrite validated r4 artifacts")
    checks = HELPERS.Checks()
    R4.verify_input(input_path)

    release_manifest_path = output_path / "release_manifest.tsv"
    release_manifest = pd.read_csv(release_manifest_path, sep="\t", dtype=str)
    expected_dependencies = set(R4.SCRIPT_NAMES)
    dependencies = release_manifest[release_manifest["record_type"] == "executable_dependency"].copy()
    observed_dependencies = {Path(value).name for value in dependencies["relative_path"]}
    checks.require(
        observed_dependencies == expected_dependencies and len(dependencies) == len(expected_dependencies),
        "complete_executable_dependency_graph",
        f"dependencies={sorted(observed_dependencies)}",
    )
    expected_roles = {
        "validated_numeric_sensitivity_ready",
        "validated_numeric_sensitivity_comparison",
        "visually_accepted_r3_ready",
        "visually_accepted_r3_panel_reference",
        "exact_figure_source_table",
        "rendered_one_page_panel",
    }
    artifact_rows = release_manifest[release_manifest["record_type"] == "input_or_output_artifact"]
    checks.require(set(artifact_rows["role"]) == expected_roles, "complete_artifact_provenance", str(sorted(expected_roles)))

    for row in release_manifest.itertuples(index=False):
        path = (
            output_path / row.relative_path
            if row.path_scope == "candidate_root"
            else root / row.relative_path
        )
        ok = (
            path.is_file()
            and path.stat().st_size == int(row.bytes)
            and R4.sha256_file(path) == row.sha256
        )
        checks.require(ok, f"provenance_{row.role}", f"{row.path_scope}:{row.relative_path}")

    manifest = pd.read_csv(output_path / "panel_manifest.tsv", sep="\t", dtype=str)
    checks.require(len(manifest) == 1, "one_panel_manifest_row", f"rows={len(manifest)}")
    row = manifest.iloc[0]
    panel = output_path / R4.PANEL_NAME
    source_path = output_path / R4.SOURCE_NAME
    checks.require(row["figure_release_id"] == R4.FIGURE_RELEASE_ID, "figure_release", R4.FIGURE_RELEASE_ID)
    checks.require(row["canonical_written"] == "FALSE", "canonical_write_forbidden", "FALSE")
    checks.require(
        row["input_release_id"] == R4.INPUT_RELEASE_ID
        and row["input_ready_sha256"] == R4.INPUT_READY_SHA256,
        "input_ready_link",
        R4.INPUT_RELEASE_ID,
    )
    checks.require(
        row["release_manifest_sha256"] == R4.sha256_file(release_manifest_path),
        "release_manifest_checksum",
        "panel manifest links full dependency/artifact manifest",
    )
    checks.require(panel.is_file() and row["panel_sha256"] == R4.sha256_file(panel), "panel_checksum", "SHA256 rederived")
    checks.require(source_path.is_file() and row["source_sha256"] == R4.sha256_file(source_path), "source_checksum", "SHA256 rederived")
    checks.require(
        row["supersedes_visual_release"] == R4.R3_FIGURE_RELEASE_ID
        and row["layout_identity"] == "pixel_layout_specification_identical_to_visually_accepted_r3",
        "r3_supersession_contract",
        "r3 visually accepted; r4 provenance-complete",
    )

    observed = pd.read_csv(source_path, sep="\t")
    expected = R4.build_source(input_path)[list(observed.columns)]
    try:
        pd.testing.assert_frame_equal(observed, expected, check_dtype=False, rtol=1e-12, atol=1e-14)
        exact = True
        detail = "four rows independently rebuilt from numeric r2"
    except AssertionError as exc:
        exact = False
        detail = str(exc).splitlines()[0]
    checks.require(exact, "source_table_exact", detail)
    checks.require(
        len(observed) == 4 and not observed.duplicated(["dataset", "program_id"]).any(),
        "complete_family",
        "2 programs x 2 datasets",
    )
    checks.require(
        set(observed.loc[observed["dataset"] == "Vu_et_al_2025", "evidence_state"]) == {"source_dependent"},
        "vu_source_dependence",
        "both Vu rows remain source-dependent",
    )
    checks.require(
        observed["claim_boundary"].eq(
            "attenuation_is_composition_linked_context_not_mediation_causality_or_cell_intrinsic_proof"
        ).all(),
        "claim_boundary",
        "no mediation, causality, or cell-intrinsic upgrade",
    )
    checks.require(observed["matched_set_hash_identical"].astype(str).str.lower().eq("true").all(), "matched_sets", "4/4 identical")
    checks.require(HELPERS.pdf_pages(panel) == 1, "one_page_pdf", "Pages=1")
    checks.require(
        re.search(rb"/Subtype\s*/Type3\b", panel.read_bytes()) is None,
        "no_type3_fonts",
        "raw PDF object scan",
    )
    text = HELPERS.command_output(["pdftotext", str(panel), "-"])
    checks.require(
        "All 16 cell types + QC" in text and "Hepatocyte + QC adjustment" in text,
        "adjustment_labels",
        "both adjustments rendered",
    )
    checks.require(
        "source-dependent" in text and "not evidence of mediation or causality" in text,
        "visible_claim_boundaries",
        "source dependence and noncausal wording rendered",
    )
    checks.require(
        "Stromal ECM (IGFBP7)" in text and "Ductular injury (BICC1)" in text,
        "program_labels",
        "both frozen program identities visible",
    )

    r3_root = root / "Analysis/Multimodal_Program_Projection/candidates" / R4.R3_FIGURE_RELEASE_ID
    with tempfile.TemporaryDirectory(prefix="composition_r4_visual_") as temporary:
        temp_root = Path(temporary)
        r3_raster = raster_sha256(r3_root / R4.PANEL_NAME, temp_root)
        r4_raster = raster_sha256(panel, temp_root)
    checks.require(
        r3_raster == r4_raster,
        "r3_r4_200dpi_pixel_identity",
        f"raster_sha256={r4_raster}",
    )

    report = output_path / "validation_report.tsv"
    pd.DataFrame(checks.rows).to_csv(report, sep="\t", index=False)
    if not checks.passed:
        failures = [row["check_id"] for row in checks.rows if row["status"] == "FAIL"]
        raise RuntimeError("r4 figure validation failed: " + ", ".join(failures))
    pd.DataFrame(
        [{
            "figure_release_id": R4.FIGURE_RELEASE_ID,
            "status": "pass_full_composition_supplementary_figure_provenance_complete",
            "panel_sha256": R4.sha256_file(panel),
            "source_sha256": R4.sha256_file(source_path),
            "release_manifest_sha256": R4.sha256_file(release_manifest_path),
            "validation_report_sha256": R4.sha256_file(report),
            "input_ready_sha256": R4.INPUT_READY_SHA256,
            "r3_r4_200dpi_raster_sha256": r4_raster,
            "canonical_written": "FALSE",
        }]
    ).to_csv(output_path / "READY", sep="\t", index=False)
    print(f"full-composition figure r4 validated: {len(checks.rows)} checks passed", flush=True)
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, default=None)
    parser.add_argument("--output-root", type=Path, default=None)
    args = parser.parse_args()
    try:
        validate(args.project_root, args.input_root, args.output_root)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
