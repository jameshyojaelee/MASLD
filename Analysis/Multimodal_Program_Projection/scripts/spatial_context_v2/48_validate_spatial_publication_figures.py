#!/usr/bin/env python3
"""Independently validate the compact spatial publication-figure candidate."""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

from spatial_resource_lib import SpatialResourceError, sha256_file


RELEASE_ID = "spatial-publication-figures-candidate-2026-08-11-r4"
PANEL_NAMES = (
    "fig1_spatial_assay_observability.pdf",
    "fig4a_spatial_firewall.pdf",
    "fig4e_spatial_program_maps.pdf",
    "fig4f_multimodal_program_summary.pdf",
    "fig5_spatial_passport_examples.pdf",
    "figS_spatial_program_coverage.pdf",
    "figS_hmsma_label_blind_organization.pdf",
    "figS_full_composition_sensitivity.pdf",
    "figS_cell_context_association.pdf",
    "figS_spatial_unit_heterogeneity.pdf",
    "figS_spatial_matched_null_calibration.pdf",
    "figS_spatial_weight_sensitivity.pdf",
    "figS_hmsma_proxy_mask_sensitivity.pdf",
)

OLD_PANELS = {
    "fig1_spatial_assay_observability.pdf": "Analysis/Multimodal_Program_Projection/candidates/spatial-figure-integration-candidate-2026-08-11-r2/panels/fig1_spatial_assay_observability.pdf",
    "fig4a_spatial_firewall.pdf": "Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11/figures/main/fig4_validation/panels/fig4a_overview_cascade.pdf",
    "fig4f_multimodal_program_summary.pdf": "Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11/figures/main/fig4_validation/panels/fig4f_multimodal_program_summary.pdf",
    "fig5_spatial_passport_examples.pdf": "Analysis/Multimodal_Program_Projection/candidates/spatial-figure-integration-candidate-2026-08-11-r2/panels/fig5_spatial_passport_examples.pdf",
    "figS_hmsma_label_blind_organization.pdf": "Analysis/Multimodal_Program_Projection/candidates/spatial-figure-integration-candidate-2026-08-11-r2/panels/figS_hmsma_label_blind_organization.pdf",
    "figS_full_composition_sensitivity.pdf": "Analysis/Multimodal_Program_Projection/candidates/spatial-full-composition-figure-candidate-2026-08-11-r4/figS_full_composition_sensitivity.pdf",
    "figS_cell_context_association.pdf": "Analysis/Multimodal_Program_Projection/candidates/spatial-cell-context-attribution-candidate-2026-08-11-r2/panels/figS_cell_context_association.pdf",
    "figS_spatial_unit_heterogeneity.pdf": "Analysis/Multimodal_Program_Projection/candidates/spatial-unit-null-diagnostics-candidate-2026-08-11-r2/panels/figS_spatial_unit_heterogeneity.pdf",
    "figS_spatial_matched_null_calibration.pdf": "Analysis/Multimodal_Program_Projection/candidates/spatial-unit-null-diagnostics-candidate-2026-08-11-r2/panels/figS_spatial_matched_null_calibration.pdf",
    "figS_spatial_weight_sensitivity.pdf": "Analysis/Multimodal_Program_Projection/candidates/spatial-unit-null-diagnostics-candidate-2026-08-11-r2/panels/figS_spatial_weight_sensitivity.pdf",
    "figS_hmsma_proxy_mask_sensitivity.pdf": "Analysis/Multimodal_Program_Projection/candidates/hmsma-mask-sensitivity-candidate-2026-08-11-r2/panels/figS_hmsma_proxy_mask_sensitivity.pdf",
}


def command(command: list[str]) -> str:
    return subprocess.run(command, check=True, capture_output=True, text=True).stdout


def pdf_info(path: Path) -> tuple[int, float, float]:
    output = command(["pdfinfo", str(path)])
    pages = None
    width = height = None
    for line in output.splitlines():
        if line.startswith("Pages:"):
            pages = int(line.split(":", 1)[1].strip())
        if line.startswith("Page size:"):
            match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", line)
            if match:
                width, height = float(match.group(1)), float(match.group(2))
    if pages is None or width is None or height is None:
        raise SpatialResourceError(f"incomplete pdfinfo output: {path}")
    return pages, width, height


def reject_type3(path: Path) -> None:
    if shutil.which("pdffonts"):
        output = command(["pdffonts", str(path)])
        found = "Type 3" in output or "Type3" in output
    else:
        found = re.search(rb"/Subtype\s*/Type3\b", path.read_bytes()) is not None
    if found:
        raise SpatialResourceError(f"Type 3 font found: {path}")


def word_count(path: Path) -> int:
    text = command(["pdftotext", str(path), "-"])
    return len(re.findall(r"\b[\w*†↔+–-]+\b", text, flags=re.UNICODE))


def assert_frame_equal(left: pd.DataFrame, right: pd.DataFrame, keys: list[str], label: str) -> None:
    left = left.sort_values(keys).reset_index(drop=True)
    right = right.sort_values(keys).reset_index(drop=True)
    if list(left.columns) != list(right.columns) or not left.equals(right):
        raise SpatialResourceError(f"figure source drifted: {label}")


def validate_sources(project: Path, output: Path) -> None:
    downstream = project / "Analysis/Multimodal_Program_Projection/candidates/spatial-figure-integration-candidate-2026-08-11-r2/data"
    resource = project / "Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11"
    unit_root = project / "Analysis/Multimodal_Program_Projection/candidates/spatial-unit-null-diagnostics-candidate-2026-08-11-r2/data"
    data = output / "data"

    comparisons = [
        (
            data / "fig1_spatial_assay_observability.tsv",
            downstream / "fig1_spatial_assay_observability.tsv",
            ["dataset_id"],
        ),
        (
            data / "fig5_spatial_passport_examples.tsv",
            downstream / "fig5_spatial_passport_examples.tsv",
            ["example_id"],
        ),
        (
            data / "figS_hmsma_label_blind_organization.tsv",
            downstream / "figS_hmsma_label_blind_organization.tsv",
            ["array_id", "program_uid"],
        ),
        (
            data / "figS_full_composition_sensitivity.tsv",
            project / "Analysis/Multimodal_Program_Projection/candidates/spatial-full-composition-figure-candidate-2026-08-11-r4/figS_full_composition_sensitivity.tsv",
            ["program_id", "dataset"],
        ),
        (
            data / "figS_spatial_unit_heterogeneity.tsv",
            unit_root / "figS_spatial_unit_heterogeneity.tsv",
            ["program_id", "dataset", "sample_id"],
        ),
        (
            data / "figS_spatial_matched_null_calibration.tsv",
            unit_root / "figS_spatial_matched_null_calibration.tsv",
            ["program_id", "dataset"],
        ),
        (
            data / "figS_spatial_weight_sensitivity.tsv",
            unit_root / "figS_spatial_weight_sensitivity.tsv",
            ["program_id", "dataset", "sensitivity_id"],
        ),
        (
            data / "figS_hmsma_proxy_mask_sensitivity.tsv",
            project / "Analysis/Multimodal_Program_Projection/candidates/hmsma-mask-sensitivity-candidate-2026-08-11-r2/data/hmsma_proxy_mask_summary.tsv",
            ["program_uid", "min_umi"],
        ),
        (
            data / "figS_spatial_program_coverage.tsv",
            downstream / "figS_spatial_program_coverage.tsv",
            ["dataset_id", "program_uid"],
        ),
    ]
    display_only = {
        "fig1_spatial_assay_observability.tsv": {"unit_state", "clinical_state", "histology_state", "coverage_state", "unit_label", "gate_mark"},
        "fig5_spatial_passport_examples.tsv": {"spatial_message", "short_next_test"},
    }
    for observed_path, expected_path, keys in comparisons:
        observed = pd.read_csv(observed_path, sep="\t", dtype=str, keep_default_na=False)
        expected = pd.read_csv(expected_path, sep="\t", dtype=str, keep_default_na=False)
        observed = observed.drop(columns=sorted(display_only.get(observed_path.name, set())))
        assert_frame_equal(observed, expected, keys, observed_path.name)

    cell_observed = pd.read_csv(data / "figS_cell_context_association.tsv", sep="\t", dtype=str, keep_default_na=False)
    cell_expected = pd.read_csv(
        project / "Analysis/Multimodal_Program_Projection/candidates/spatial-cell-context-attribution-candidate-2026-08-11-r2/data/figS_cell_context_association.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    focus = cell_expected[
        ((cell_expected["program_uid"] == "hotspot_hepatocytes_f05c535ae5bbc0b9") & (cell_expected["cell2location_factor"] == "Fibroblasts"))
        | ((cell_expected["program_uid"] == "hotspot_hepatocytes_48f39dd4d817a10e") & (cell_expected["cell2location_factor"] == "Cholangiocytes"))
    ]
    if len(cell_observed) != 4:
        raise SpatialResourceError("focused cell-context figure source must contain exactly four rows")
    assert_frame_equal(cell_observed, focus, ["program_uid", "dataset", "cell2location_factor"], "cell context focus")

    effects_observed = pd.read_csv(data / "fig4f_multimodal_program_summary.tsv", sep="\t", dtype=str, keep_default_na=False)
    effects_expected = pd.read_parquet(resource / "effects/spatial_program_effects.parquet").astype(str)
    dataset_order = {
        "GSE192741", "Vu_et_al_2025", "HRA007511_HMSMA", "Yakubovsky2026",
        "Govaere2026_CosMx", "Govaere2026_GeoMx", "GSE287826", "GSE244832", "GSE281367", "PXD051911",
    }
    effects_expected = effects_expected[effects_expected["dataset_id"].isin(dataset_order)]
    if len(effects_observed) != 20:
        raise SpatialResourceError("Figure 4F source must contain 20 rows")
    assert_frame_equal(effects_observed, effects_expected, ["dataset_id", "program_uid"], "Figure 4F")

    firewall = pd.read_csv(data / "fig4a_spatial_firewall.tsv", sep="\t")
    if firewall[["path", "n_programs"]].to_records(index=False).tolist() != [("coverage", 117), ("inference", 2)]:
        raise SpatialResourceError("Figure 4A firewall counts drifted")

    if len(pd.read_csv(data / "figS_spatial_program_coverage.tsv", sep="\t")) != 819:
        raise SpatialResourceError("coverage source is not 117 programs by seven spatial sources")


def validate_manifests(project: Path, output: Path) -> None:
    panels = pd.read_csv(output / "panel_manifest.tsv", sep="\t", dtype=str)
    if len(panels) != len(PANEL_NAMES) or set(panels["panel"]) != {Path(x).stem for x in PANEL_NAMES}:
        raise SpatialResourceError("panel manifest is incomplete")
    if set(panels["canonical_written"]) != {"FALSE"}:
        raise SpatialResourceError("publication candidate claims a canonical write")
    for row in panels.itertuples(index=False):
        path = output / row.relative_path
        if not path.is_file() or path.stat().st_size != int(row.bytes) or sha256_file(path) != row.sha256:
            raise SpatialResourceError(f"panel manifest drift: {row.relative_path}")

    sources = pd.read_csv(output / "source_manifest.tsv", sep="\t", dtype=str)
    required_roles = {
        "sealed_resource_ready", "sealed_downstream_ready", "sealed_composition_ready",
        "sealed_cell_context_ready", "sealed_unit_null_ready", "sealed_mask_ready",
        "figure1_exact_source", "figure5_exact_source", "hmsma_label_blind_exact_source",
        "coverage_exact_source", "figure4f_exact_source", "composition_exact_source",
        "cell_context_complete_source", "unit_exact_source", "matched_null_exact_source",
        "weight_sensitivity_exact_source", "mask_summary_exact_source", "retained_figure4e_panel",
        "retained_coverage_panel", "figure4e_selection_source", "figure4e_per_spot_source",
        "producer", "validator", "execution_wrapper",
    }
    if len(sources) != len(required_roles) or set(sources["source_role"]) != required_roles:
        raise SpatialResourceError("source manifest does not seal the complete input and code dependency set")
    for row in sources.itertuples(index=False):
        path = project / row.relative_path
        if not path.is_file() or path.stat().st_size != int(row.bytes) or sha256_file(path) != row.sha256:
            raise SpatialResourceError(f"source manifest drift: {row.relative_path}")

    disposition = pd.read_csv(output / "review_disposition.tsv", sep="\t")
    if len(disposition) != len(PANEL_NAMES) or disposition["panel"].nunique() != len(PANEL_NAMES):
        raise SpatialResourceError("review disposition is incomplete")
    if disposition["decision"].value_counts().to_dict() != {"simplified": 11, "retained": 2}:
        raise SpatialResourceError("review disposition counts drifted")


def validate_retained(project: Path, output: Path) -> None:
    pairs = [
        (
            project / "Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11/figures/main/fig4_validation/panels/fig4e_spatial_program_maps.pdf",
            output / "panels/fig4e_spatial_program_maps.pdf",
        ),
        (
            project / "Analysis/Multimodal_Program_Projection/candidates/spatial-figure-integration-candidate-2026-08-11-r2/panels/figS_spatial_program_coverage.pdf",
            output / "panels/figS_spatial_program_coverage.pdf",
        ),
    ]
    for old, new in pairs:
        if sha256_file(old) != sha256_file(new):
            raise SpatialResourceError(f"retained panel is not byte-identical: {new.name}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    project = args.project_root.resolve()
    output = args.output_root.resolve()
    try:
        if output.name != RELEASE_ID or not (output / "BUILD_COMPLETE").is_file() or (output / "READY").exists():
            raise SpatialResourceError("publication candidate is absent, incomplete, or already immutable")
        validate_sources(project, output)
        validate_manifests(project, output)
        validate_retained(project, output)

        panel_audit = []
        for name in PANEL_NAMES:
            path = output / "panels" / name
            pages, width, height = pdf_info(path)
            if pages != 1:
                raise SpatialResourceError(f"panel is not one page: {name}")
            if width > 400 or height > 400:
                raise SpatialResourceError(f"panel canvas exceeds compact 5.55-inch bound: {name}: {width}x{height}")
            reject_type3(path)
            panel_audit.append(
                {
                    "figure_release_id": RELEASE_ID,
                    "panel": Path(name).stem,
                    "pages": pages,
                    "width_points": width,
                    "height_points": height,
                    "type3_fonts": "FALSE",
                }
            )
        pd.DataFrame(panel_audit).to_csv(output / "layout_audit.tsv", sep="\t", index=False)

        cognitive = []
        for name, old_relative in OLD_PANELS.items():
            old = project / old_relative
            new = output / "panels" / name
            old_words, new_words = word_count(old), word_count(new)
            cognitive.append(
                {
                    "figure_release_id": RELEASE_ID,
                    "panel": Path(name).stem,
                    "old_words": old_words,
                    "new_words": new_words,
                    "fraction_retained": new_words / old_words if old_words else float("nan"),
                    "words_removed": old_words - new_words,
                }
            )
        cognitive = pd.DataFrame(cognitive)
        total_fraction = cognitive["new_words"].sum() / cognitive["old_words"].sum()
        if total_fraction >= 0.75:
            raise SpatialResourceError(f"aggregate text reduction is inadequate: retained fraction={total_fraction:.3f}")
        focused_limits = {
            "fig5_spatial_passport_examples": 0.60,
            "figS_cell_context_association": 0.60,
            "figS_hmsma_proxy_mask_sensitivity": 0.70,
        }
        for panel, limit in focused_limits.items():
            fraction = float(cognitive.loc[cognitive["panel"] == panel, "fraction_retained"].iloc[0])
            if fraction >= limit:
                raise SpatialResourceError(f"text reduction is inadequate for {panel}: {fraction:.3f} >= {limit}")
        cognitive.to_csv(output / "cognitive_load_audit.tsv", sep="\t", index=False)

        checks = pd.DataFrame(
            [
                ("PUB-FAMILY", "pass", len(PANEL_NAMES), len(PANEL_NAMES), "complete panel family"),
                ("PUB-SOURCE", "pass", "exact", "exact", "all plotted tables rederive from sealed inputs"),
                ("PUB-PDF", "pass", "one page; no Type 3", "one page; no Type 3", "editable compact PDFs"),
                ("PUB-TEXT", "pass", total_fraction, "<0.75", "aggregate plot text reduced"),
                ("PUB-RETAIN", "pass", 2, 2, "already-clear panels retained byte-for-byte"),
                ("PUB-CANONICAL", "pass", "FALSE", "FALSE", "canonical figures untouched"),
            ],
            columns=["check_id", "status", "observed", "expected", "detail"],
        )
        checks.to_csv(output / "validation_report.tsv", sep="\t", index=False)
        ready = pd.DataFrame(
            [
                {
                    "figure_release_id": RELEASE_ID,
                    "status": "validated_publication_layout_candidate_awaiting_adjudication",
                    "n_panels": len(PANEL_NAMES),
                    "n_simplified": 11,
                    "n_retained": 2,
                    "aggregate_text_fraction_retained": total_fraction,
                    "panel_manifest_sha256": sha256_file(output / "panel_manifest.tsv"),
                    "source_manifest_sha256": sha256_file(output / "source_manifest.tsv"),
                    "validation_report_sha256": sha256_file(output / "validation_report.tsv"),
                    "layout_audit_sha256": sha256_file(output / "layout_audit.tsv"),
                    "cognitive_load_audit_sha256": sha256_file(output / "cognitive_load_audit.tsv"),
                    "canonical_written": "FALSE",
                }
            ]
        )
        ready.to_csv(output / "READY", sep="\t", index=False)
        print(f"PASS compact spatial publication figures: {len(PANEL_NAMES)} panels; text fraction {total_fraction:.3f}")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
