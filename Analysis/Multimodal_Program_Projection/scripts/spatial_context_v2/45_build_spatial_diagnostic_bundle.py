#!/usr/bin/env python3
"""Index the accepted spatial diagnostic candidates in one immutable bundle."""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import tempfile
from pathlib import Path

import pandas as pd


RELEASE_ID = "spatial-diagnostic-bundle-candidate-2026-08-11"
COMPONENTS = (
    {
        "component_id": "full_composition_numeric",
        "candidate": "spatial-full-composition-sensitivity-candidate-2026-08-11-r2",
        "status": "pass_full_composition_sensitivity",
        "paper_role": "supplementary_falsification",
        "biological_unit": "GSE_donor_first;Vu_unresolved_physical_array",
        "permitted_claim": "composition_linked_spatial_organization_sensitivity",
        "prohibited_claim": "mediation_causality_or_cell_intrinsic_mechanism",
        "expected_pdfs": 0,
    },
    {
        "component_id": "full_composition_figure",
        "candidate": "spatial-full-composition-figure-candidate-2026-08-11-r4",
        "status": "pass_full_composition_supplementary_figure_provenance_complete",
        "paper_role": "supplementary_panel",
        "biological_unit": "GSE_donor_first;Vu_unresolved_physical_array",
        "permitted_claim": "visualize_prespecified_composition_attenuation",
        "prohibited_claim": "mediation_causality_or_cell_intrinsic_mechanism",
        "expected_pdfs": 1,
    },
    {
        "component_id": "unit_null_weight_diagnostics",
        "candidate": "spatial-unit-null-diagnostics-candidate-2026-08-11-r2",
        "status": "validated_spatial_unit_null_diagnostic_candidate_awaiting_adjudication",
        "paper_role": "supplementary_diagnostics",
        "biological_unit": "GSE_donor_first;Vu_unresolved_physical_array",
        "permitted_claim": "unit_heterogeneity_matched_null_and_weight_sensitivity",
        "prohibited_claim": "section_or_array_level_population_inference",
        "expected_pdfs": 3,
    },
    {
        "component_id": "cell_context_attribution",
        "candidate": "spatial-cell-context-attribution-candidate-2026-08-11-r2",
        "status": "validated_descriptive_cell_context_candidate_awaiting_adjudication",
        "paper_role": "supplementary_descriptive_context",
        "biological_unit": "GSE_donor_first;Vu_unresolved_physical_array",
        "permitted_claim": "shared_RNA_matrix_spatial_context_covariation",
        "prohibited_claim": "lineage_origin_independent_validation_or_mechanism",
        "expected_pdfs": 1,
    },
    {
        "component_id": "hmsma_proxy_mask_sensitivity",
        "candidate": "hmsma-mask-sensitivity-candidate-2026-08-11-r2",
        "status": "validated_hmsma_proxy_mask_sensitivity_awaiting_adjudication",
        "paper_role": "supplementary_descriptive_sensitivity",
        "biological_unit": "35_unresolved_physical_arrays",
        "permitted_claim": "stability_under_stricter_UMI_proxy_masks",
        "prohibited_claim": "donor_disease_histology_or_population_inference",
        "expected_pdfs": 1,
    },
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def read_ready_status(path: Path) -> str:
    ready = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if len(ready) != 1 or "status" not in ready:
        raise RuntimeError(f"invalid component READY: {path}")
    return ready.iloc[0]["status"]


def code_paths(script_root: Path, project: Path) -> list[Path]:
    names = [
        "22_build_hmsma_label_blind.py",
        "31_run_full_composition_sensitivity.py",
        "32_validate_full_composition_sensitivity.py",
        "33_render_unit_null_diagnostics.py",
        "34_validate_unit_null_diagnostics.py",
        "35_build_cell_context_attribution.py",
        "36_validate_cell_context_attribution.py",
        "37_build_hmsma_mask_sensitivity.py",
        "38_validate_hmsma_mask_sensitivity.py",
        "39_render_full_composition_figure.py",
        "40_validate_full_composition_figure.py",
        "41_render_full_composition_figure_r3.py",
        "42_validate_full_composition_figure_r3.py",
        "43_render_full_composition_figure_r4.py",
        "44_validate_full_composition_figure_r4.py",
        "45_build_spatial_diagnostic_bundle.py",
        "46_validate_spatial_diagnostic_bundle.py",
        "spatial_resource_lib.py",
        "visium_rerun_lib.py",
        "run_full_composition_sensitivity.sbatch",
        "run_full_composition_figure_r4.sbatch",
        "run_unit_null_diagnostics.sbatch",
        "run_cell_context_attribution.sbatch",
        "run_hmsma_mask_sensitivity.sbatch",
        "run_spatial_diagnostic_bundle.sbatch",
        "tests/test_full_composition_sensitivity.py",
    ]
    paths = [script_root / name for name in names]
    paths.append(project / "Analysis/Multimodal_Program_Projection/scripts/04_spatial_projection.py")
    return paths


def build(project: Path, output: Path) -> None:
    project = project.resolve()
    output = output.resolve()
    candidate_parent = (project / "Analysis/Multimodal_Program_Projection/candidates").resolve()
    if output.parent != candidate_parent or output.name != RELEASE_ID:
        raise RuntimeError("diagnostic bundle must be the named direct candidate child")
    if output.exists():
        raise RuntimeError(f"refusing to overwrite diagnostic bundle: {output}")
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
    try:
        registry_rows = []
        artifact_rows = []
        for component in COMPONENTS:
            root = candidate_parent / component["candidate"]
            ready = root / "READY"
            observed_status = read_ready_status(ready)
            if observed_status != component["status"]:
                raise RuntimeError(
                    f"component status drift for {component['component_id']}: {observed_status}"
                )
            pdfs = sorted(root.rglob("*.pdf"))
            if len(pdfs) != component["expected_pdfs"]:
                raise RuntimeError(f"unexpected PDF family for {component['component_id']}")
            registry_rows.append({
                "bundle_release_id": RELEASE_ID,
                "component_id": component["component_id"],
                "candidate_release_id": component["candidate"],
                "relative_root": root.relative_to(project).as_posix(),
                "ready_status": observed_status,
                "ready_sha256": sha256_file(ready),
                "paper_role": component["paper_role"],
                "biological_unit": component["biological_unit"],
                "permitted_claim": component["permitted_claim"],
                "prohibited_claim": component["prohibited_claim"],
                "n_pdfs": len(pdfs),
                "canonical_write_allowed": "FALSE",
            })
            for path in sorted(item for item in root.rglob("*") if item.is_file()):
                artifact_rows.append({
                    "bundle_release_id": RELEASE_ID,
                    "component_id": component["component_id"],
                    "relative_path": path.relative_to(project).as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                })
        pd.DataFrame(registry_rows).to_csv(staging / "diagnostic_registry.tsv", sep="\t", index=False)
        pd.DataFrame(artifact_rows).to_csv(staging / "artifact_manifest.tsv", sep="\t", index=False)

        script_root = Path(__file__).resolve().parent
        code_rows = []
        for path in code_paths(script_root, project):
            if not path.is_file():
                raise RuntimeError(f"missing diagnostic code dependency: {path}")
            code_rows.append({
                "bundle_release_id": RELEASE_ID,
                "relative_path": path.relative_to(project).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            })
        pd.DataFrame(code_rows).to_csv(staging / "code_manifest.tsv", sep="\t", index=False)
        pd.DataFrame([{
            "bundle_release_id": RELEASE_ID,
            "status": "built_pending_independent_validation",
            "n_components": len(registry_rows),
            "n_pdfs": sum(row["n_pdfs"] for row in registry_rows),
            "n_artifacts": len(artifact_rows),
            "n_code_files": len(code_rows),
            "diagnostic_registry_sha256": sha256_file(staging / "diagnostic_registry.tsv"),
            "artifact_manifest_sha256": sha256_file(staging / "artifact_manifest.tsv"),
            "code_manifest_sha256": sha256_file(staging / "code_manifest.tsv"),
            "canonical_write_allowed": "FALSE",
        }]).to_csv(staging / "BUILD_COMPLETE", sep="\t", index=False)
        os.replace(staging, output)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    print(f"built spatial diagnostic bundle: {output}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    build(args.project_root, args.output_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
