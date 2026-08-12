#!/usr/bin/env python3
# KEY MESSAGE: IGFBP7 spatial organization attenuates but remains supported after all-cell-type adjustment, whereas BICC1 remains indeterminate.
"""Render provenance-complete r4 of the accepted composition panel."""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path


FIGURE_RELEASE_ID = "spatial-full-composition-figure-candidate-2026-08-11-r4"
R3_FIGURE_RELEASE_ID = "spatial-full-composition-figure-candidate-2026-08-11-r3"
SCRIPT_NAMES = (
    "39_render_full_composition_figure.py",
    "40_validate_full_composition_figure.py",
    "41_render_full_composition_figure_r3.py",
    "42_validate_full_composition_figure_r3.py",
    "43_render_full_composition_figure_r4.py",
    "44_validate_full_composition_figure_r4.py",
)


def import_r3(script_dir: Path):
    path = script_dir / "41_render_full_composition_figure_r3.py"
    spec = importlib.util.spec_from_file_location("accepted_full_composition_figure_r3", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import accepted r3 renderer: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.FIGURE_RELEASE_ID = FIGURE_RELEASE_ID
    module.BASE.FIGURE_RELEASE_ID = FIGURE_RELEASE_ID
    return module


SCRIPT_DIR = Path(__file__).resolve().parent
R3 = import_r3(SCRIPT_DIR)
INPUT_RELEASE_ID = R3.INPUT_RELEASE_ID
INPUT_READY_SHA256 = R3.INPUT_READY_SHA256
PANEL_NAME = R3.PANEL_NAME
SOURCE_NAME = R3.SOURCE_NAME
sha256_file = R3.sha256_file
verify_input = R3.verify_input
build_source = R3.build_source
render = R3.render


def default_input(project_root: Path) -> Path:
    return R3.default_input(project_root)


def default_output(project_root: Path) -> Path:
    return project_root / "Analysis/Multimodal_Program_Projection/candidates" / FIGURE_RELEASE_ID


def dependency_rows(project_root: Path) -> list[dict[str, object]]:
    rows = []
    for name in SCRIPT_NAMES:
        path = SCRIPT_DIR / name
        if not path.is_file():
            raise RuntimeError(f"missing executable dependency: {path}")
        role = {
            "39_render_full_composition_figure.py": "transitive_source_builder_and_constants",
            "40_validate_full_composition_figure.py": "validation_helpers_and_base_contract",
            "41_render_full_composition_figure_r3.py": "accepted_r3_layout_renderer",
            "42_validate_full_composition_figure_r3.py": "accepted_r3_validator_dependency_loader",
            "43_render_full_composition_figure_r4.py": "r4_candidate_producer",
            "44_validate_full_composition_figure_r4.py": "r4_independent_validator",
        }[name]
        rows.append(
            {
                "record_type": "executable_dependency",
                "role": role,
                "path_scope": "project_root",
                "relative_path": path.relative_to(project_root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    return rows


def run(project_root: Path, input_root: Path | None, output_root: Path | None) -> Path:
    import pandas as pd

    root = project_root.resolve()
    input_path = (input_root or default_input(root)).resolve()
    output_path = (output_root or default_output(root)).resolve()
    candidate_parent = (root / "Analysis/Multimodal_Program_Projection/candidates").resolve()
    if output_path.parent != candidate_parent or output_path.name != FIGURE_RELEASE_ID:
        raise RuntimeError("figure output must be the named isolated r4 candidate root")
    if output_path.exists():
        raise RuntimeError(f"refusing to overwrite figure candidate: {output_path}")
    verify_input(input_path)
    staging = output_path.parent / f".{output_path.name}.incomplete.{os.environ.get('SLURM_JOB_ID', os.getpid())}"
    if staging.exists():
        raise RuntimeError(f"staging path already exists: {staging}")
    staging.mkdir()

    source = build_source(input_path)
    columns = [
        "figure_release_id", "release_id", "program_id", "program_label", "dataset",
        "dataset_label", "display_order", "primary_adjustment", "full_adjustment",
        "primary_centered_moran_i", "full_composition_centered_moran_i", "centered_change",
        "primary_residual_pvalue", "primary_residual_padj", "residual_pvalue",
        "residual_padj", "primary_robust", "robust", "evidence_state", "interpretation",
        "matched_set_sha256", "matched_set_hash_identical", "claim_boundary",
    ]
    source[columns].to_csv(staging / SOURCE_NAME, sep="\t", index=False)
    render(source, staging / PANEL_NAME)

    manifest_rows = dependency_rows(root)
    r3_root = root / "Analysis/Multimodal_Program_Projection/candidates" / R3_FIGURE_RELEASE_ID
    input_files = (
        (input_path / "READY", "validated_numeric_sensitivity_ready"),
        (input_path / "comparison_to_primary.tsv", "validated_numeric_sensitivity_comparison"),
        (r3_root / "READY", "visually_accepted_r3_ready"),
        (r3_root / PANEL_NAME, "visually_accepted_r3_panel_reference"),
        (staging / SOURCE_NAME, "exact_figure_source_table"),
        (staging / PANEL_NAME, "rendered_one_page_panel"),
    )
    for path, role in input_files:
        manifest_rows.append(
            {
                "record_type": "input_or_output_artifact",
                "role": role,
                "path_scope": "candidate_root" if path.parent == staging else "project_root",
                "relative_path": path.name if path.parent == staging else path.relative_to(root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    pd.DataFrame(manifest_rows).to_csv(staging / "release_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{
            "figure_release_id": FIGURE_RELEASE_ID,
            "panel": PANEL_NAME,
            "panel_sha256": sha256_file(staging / PANEL_NAME),
            "source_table": SOURCE_NAME,
            "source_sha256": sha256_file(staging / SOURCE_NAME),
            "release_manifest": "release_manifest.tsv",
            "release_manifest_sha256": sha256_file(staging / "release_manifest.tsv"),
            "input_release_id": INPUT_RELEASE_ID,
            "input_ready_sha256": INPUT_READY_SHA256,
            "producer_sha256": sha256_file(Path(__file__).resolve()),
            "supersedes_visual_release": "spatial-full-composition-figure-candidate-2026-08-11-r3",
            "supersession_reason": "complete_transitive_code_and_artifact_provenance",
            "layout_identity": "pixel_layout_specification_identical_to_visually_accepted_r3",
            "canonical_written": "FALSE",
        }]
    ).to_csv(staging / "panel_manifest.tsv", sep="\t", index=False)
    os.replace(staging, output_path)
    print(f"full-composition supplementary figure r4 rendered: {output_path}", flush=True)
    return output_path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, default=None)
    parser.add_argument("--output-root", type=Path, default=None)
    args = parser.parse_args()
    try:
        run(args.project_root, args.input_root, args.output_root)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
