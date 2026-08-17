#!/usr/bin/env python3
"""Independently validate the result-led spatial impact figure family."""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from spatial_resource_lib import SpatialResourceError, sha256_file


RELEASE_ID = "spatial-impact-figures-candidate-2026-08-12"
R4_RELEASE_ID = "spatial-publication-figures-candidate-2026-08-11-r4"
PROGRAM_RELEASE_ID = "program-context-v2-candidate-2026-08-07"
PROGRAMS = [
    "hotspot_hepatocytes_f05c535ae5bbc0b9",
    "hotspot_hepatocytes_48f39dd4d817a10e",
]
DATASETS = ["GSE192741", "Vu_et_al_2025"]


def command(args: list[str]) -> str:
    return subprocess.run(args, check=True, capture_output=True, text=True).stdout


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


def text(path: Path) -> str:
    return command(["pdftotext", str(path), "-"])


def close(observed: object, expected: object, label: str) -> None:
    if not np.isclose(float(observed), float(expected), rtol=0, atol=1e-12, equal_nan=True):
        raise SpatialResourceError(f"{label} drifted: observed={observed}, expected={expected}")


def validate_profile(project: Path, output: Path) -> pd.DataFrame:
    candidates = project / "Analysis/Multimodal_Program_Projection/candidates"
    profile = pd.read_csv(output / "data/spatial_transportability_profile.tsv", sep="\t")
    figure = pd.read_csv(output / "data/fig4f_matched_null_discrimination.tsv", sep="\t")
    if len(profile) != 4 or set(profile["program_id"]) != set(PROGRAMS) or set(profile["dataset"]) != set(DATASETS):
        raise SpatialResourceError("transportability profile is not the complete two-program by two-source family")
    keys = ["program_id", "dataset"]
    if not profile.sort_values(keys).reset_index(drop=True).equals(figure.sort_values(keys).reset_index(drop=True)):
        raise SpatialResourceError("Figure 4F source does not exactly equal the impact profile")

    registry = pd.read_csv(candidates / PROGRAM_RELEASE_ID / "hotspot/program_registry_v2.tsv", sep="\t")
    registry = registry[registry["program_uid"].isin(PROGRAMS)].set_index("program_uid")
    for row in profile.itertuples(index=False):
        source = registry.loc[row.program_id]
        if str(row.membership_sha256) != str(source.membership_sha256):
            raise SpatialResourceError(f"membership hash drift: {row.program_id}")
        for field in ("primary_qvalue", "primary_hc3_qvalue", "stability_median"):
            close(getattr(row, field), source[field], f"{row.program_id}/{field}")
        if str(row.robust_display).lower() not in {"true", "1"} or str(source.robust_display).lower() not in {"true", "1"}:
            raise SpatialResourceError(f"program is not frozen robust-display: {row.program_id}")
        if str(row.external_outcomes_read).upper() != "FALSE":
            raise SpatialResourceError("spatial outcomes leaked into frozen program selection")

    null = pd.read_csv(candidates / "spatial-unit-null-diagnostics-candidate-2026-08-11-r2/data/figS_spatial_matched_null_calibration.tsv", sep="\t")
    null = null.set_index(keys)
    for row in profile.itertuples(index=False):
        source = null.loc[(row.program_id, row.dataset)]
        mapping = {
            "residual_moran_i": "residual_moran_i",
            "primary_qvalue_spatial": "residual_padj",
            "null_mean": "null_mean",
            "q050": "q050",
            "q500": "q500",
            "q950": "q950",
            "centered_residual_moran_i": "centered_residual_moran_i",
        }
        for observed_field, source_field in mapping.items():
            close(getattr(row, observed_field), source[source_field], f"{row.program_id}/{row.dataset}/{observed_field}")
        if int(row.n_null) != 9999 or str(row.matched_set_sha256) != str(source.matched_set_sha256):
            raise SpatialResourceError(f"matched-null contract drift: {row.program_id}/{row.dataset}")
        if row.uncertainty_semantics != "matched_gene_null_quantiles_not_sampling_confidence_interval":
            raise SpatialResourceError("null quantiles were relabeled as sampling uncertainty")

    unit = pd.read_csv(candidates / "spatial-unit-null-diagnostics-candidate-2026-08-11-r2/data/figS_spatial_unit_heterogeneity.tsv", sep="\t")
    expected_counts = {}
    for (program_id, dataset), part in unit.groupby(keys, observed=True):
        if dataset == "GSE192741":
            values = part.groupby("aggregation_group", observed=True)["residual_moran_i"].mean()
        else:
            values = part.set_index("sample_id")["residual_moran_i"]
        expected_counts[(program_id, dataset)] = (len(values), int((values > 0).sum()))
    for row in profile.itertuples(index=False):
        expected = expected_counts[(row.program_id, row.dataset)]
        if (int(row.n_reporting_units), int(row.n_positive_units)) != expected:
            raise SpatialResourceError(f"reporting-unit sign count drift: {row.program_id}/{row.dataset}")
        close(row.positive_fraction, expected[1] / expected[0], f"{row.program_id}/{row.dataset}/positive_fraction")

    expected_signs = {
        (PROGRAMS[0], "GSE192741"): (4, 4),
        (PROGRAMS[0], "Vu_et_al_2025"): (10, 10),
        (PROGRAMS[1], "GSE192741"): (4, 4),
        (PROGRAMS[1], "Vu_et_al_2025"): (10, 9),
    }
    if expected_counts != expected_signs:
        raise SpatialResourceError(f"unexpected unit-sign profile: {expected_counts}")
    if not (profile["residual_moran_i"] > 0).all():
        raise SpatialResourceError("not all four aggregate Moran statistics are positive")

    effects = pd.read_parquet(candidates / "spatial-resource-candidate-2026-08-11/effects/spatial_program_effects.parquet")
    effects = effects[effects["program_uid"].isin(PROGRAMS) & effects["dataset_id"].isin(DATASETS)].set_index(["program_uid", "dataset_id"])
    for row in profile.itertuples(index=False):
        source = effects.loc[(row.program_id, row.dataset)]
        if row.primary_evidence_state != source.evidence_state:
            raise SpatialResourceError(f"accepted evidence state drift: {row.program_id}/{row.dataset}")

    full = pd.read_csv(candidates / "spatial-full-composition-sensitivity-candidate-2026-08-11-r2/full_composition_program_results.tsv", sep="\t").set_index(keys)
    reuse = pd.read_csv(candidates / "spatial-full-composition-sensitivity-candidate-2026-08-11-r2/matched_set_reuse_audit.tsv", sep="\t")
    if len(reuse) != 4 or not reuse["identical"].astype(bool).all():
        raise SpatialResourceError("composition sensitivity did not preserve all four matched sets")
    for row in profile.itertuples(index=False):
        source = full.loc[(row.program_id, row.dataset)]
        close(row.full_residual_moran_i, source.residual_moran_i, f"{row.program_id}/{row.dataset}/full_residual")
        close(row.full_null_mean, source.residual_null_mean, f"{row.program_id}/{row.dataset}/full_null")
        close(row.full_centered_residual_moran_i, source.residual_moran_i - source.residual_null_mean, f"{row.program_id}/{row.dataset}/full_centered")
        close(row.full_qvalue, source.residual_padj, f"{row.program_id}/{row.dataset}/full_q")
        if bool(row.full_robust) != bool(source.robust) or bool(row.full_within_source_support) != bool(source.within_source_support):
            raise SpatialResourceError(f"full-composition robustness drift: {row.program_id}/{row.dataset}")

    ig = profile[profile["program_id"] == PROGRAMS[0]]
    bi = profile[profile["program_id"] == PROGRAMS[1]]
    if not (ig["within_source_call"] == "supported").all() or not (bi["within_source_call"] == "indeterminate").all():
        raise SpatialResourceError("impact contrast does not match the accepted asymmetric result")
    if not ig["full_robust"].astype(bool).all() or bi["full_robust"].astype(bool).any():
        raise SpatialResourceError("composition sensitivity does not preserve the asymmetric result")
    return profile


def validate_manifests(project: Path, output: Path) -> None:
    panel_manifest = pd.read_csv(output / "panel_manifest.tsv", sep="\t", dtype=str)
    if len(panel_manifest) != 14 or panel_manifest["panel"].nunique() != 14:
        raise SpatialResourceError("impact panel manifest must contain 14 unique panels")
    if set(panel_manifest["canonical_written"]) != {"FALSE"}:
        raise SpatialResourceError("impact candidate claims a canonical write")
    for row in panel_manifest.itertuples(index=False):
        path = output / row.relative_path
        if not path.is_file() or path.stat().st_size != int(row.bytes) or sha256_file(path) != row.sha256:
            raise SpatialResourceError(f"impact panel manifest drift: {row.relative_path}")

    required_roles = {
        "accepted_r4_publication_family", "frozen_program_registry", "physical_unit_results",
        "matched_null_results", "accepted_spatial_evidence_states", "all_cell_composition_sensitivity",
        "all_cell_composition_ready", "matched_set_reuse_audit", "producer", "validator", "execution_wrapper",
    }
    sources = pd.read_csv(output / "source_manifest.tsv", sep="\t", dtype=str)
    if len(sources) != len(required_roles) or set(sources["source_role"]) != required_roles:
        raise SpatialResourceError("impact source manifest is incomplete")
    for row in sources.itertuples(index=False):
        path = project / row.relative_path
        if not path.is_file() or path.stat().st_size != int(row.bytes) or sha256_file(path) != row.sha256:
            raise SpatialResourceError(f"impact source manifest drift: {row.relative_path}")


def validate_retained(project: Path, output: Path) -> None:
    r4 = project / "Analysis/Multimodal_Program_Projection/candidates" / R4_RELEASE_ID
    for source in sorted((r4 / "panels").glob("*.pdf")):
        if source.name == "fig4a_spatial_firewall.pdf":
            continue
        target_name = "figS_multimodal_program_source_matrix.pdf" if source.name == "fig4f_multimodal_program_summary.pdf" else source.name
        target = output / "panels" / target_name
        if not target.is_file() or sha256_file(source) != sha256_file(target):
            raise SpatialResourceError(f"retained r4 panel drift: {target_name}")


def validate_claim_contract(output: Path) -> None:
    contract = pd.read_csv(output / "impact_claim_contract.tsv", sep="\t")
    if contract["claim_class"].value_counts().to_dict() != {"allowed": 4, "prohibited": 5}:
        raise SpatialResourceError("impact claim contract is incomplete")
    joined = "\n".join(text(path) for path in sorted((output / "panels").glob("*.pdf"))).lower()
    for prohibited in ("bicc1 is not spatial", "cell autonomous", "causal", "donor-level replication"):
        if prohibited in joined:
            raise SpatialResourceError(f"prohibited impact claim leaked into artwork: {prohibited}")
    logic_text = text(output / "panels/fig4a_spatial_calibration_logic.pdf").lower()
    result_text = text(output / "panels/fig4f_matched_null_discrimination.pdf").lower()
    for required in ("positive", "moran", "matched-gene", "igfbp7", "bicc1", "supported", "indeterminate"):
        if required not in logic_text:
            raise SpatialResourceError(f"Figure 4A omits required result logic: {required}")
    for required in ("matched genes", "observed program", "4/4 donors", "10/10 arrays", "9/10 arrays", "supported", "indeterminate"):
        if required not in result_text:
            raise SpatialResourceError(f"Figure 4F omits required discrimination evidence: {required}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    project = args.project_root.resolve()
    output = args.output_root.resolve()
    try:
        if output.name != RELEASE_ID or not (output / "BUILD_COMPLETE").is_file() or (output / "READY").exists():
            raise SpatialResourceError("impact candidate is absent, incomplete, or already immutable")
        profile = validate_profile(project, output)
        validate_manifests(project, output)
        validate_retained(project, output)
        validate_claim_contract(output)

        layout_rows = []
        for panel in sorted((output / "panels").glob("*.pdf")):
            pages, width, height = pdf_info(panel)
            if pages != 1 or width > 400 or height > 400:
                raise SpatialResourceError(f"impact panel is not a compact one-page PDF: {panel.name}: {pages}, {width}x{height}")
            reject_type3(panel)
            layout_rows.append(
                {
                    "impact_release_id": RELEASE_ID,
                    "panel": panel.stem,
                    "pages": pages,
                    "width_points": width,
                    "height_points": height,
                    "type3_fonts": "FALSE",
                }
            )
        pd.DataFrame(layout_rows).to_csv(output / "layout_audit.tsv", sep="\t", index=False)

        build = pd.read_csv(output / "BUILD_COMPLETE", sep="\t", dtype=str)
        if len(build) != 1 or build.iloc[0]["new_inferential_family"] != "FALSE" or build.iloc[0]["canonical_written"] != "FALSE":
            raise SpatialResourceError("impact build broadened inference or wrote canonical outputs")

        checks = pd.DataFrame(
            [
                ("IMPACT-FAMILY", "pass", 14, 14, "complete compact publication family"),
                ("IMPACT-FROZEN", "pass", 2, 2, "exact frozen robust-display programs"),
                ("IMPACT-POSITIVE", "pass", 4, 4, "all aggregate Moran statistics are positive"),
                ("IMPACT-NULL", "pass", "2 supported / 2 indeterminate", "2 supported / 2 indeterminate", "matched-gene null discriminates programs"),
                ("IMPACT-UNITS", "pass", "4/4, 10/10, 4/4, 9/10", "4/4, 10/10, 4/4, 9/10", "positive unit signs independently rederived"),
                ("IMPACT-COMPOSITION", "pass", "IGFBP7 persists; BICC1 indeterminate", "IGFBP7 persists; BICC1 indeterminate", "all-cell sensitivity remains supplementary support"),
                ("IMPACT-INFERENCE", "pass", "FALSE", "FALSE", "no new inferential family or score"),
                ("IMPACT-CANONICAL", "pass", "FALSE", "FALSE", "canonical figures untouched"),
            ],
            columns=["check_id", "status", "observed", "expected", "detail"],
        )
        checks.to_csv(output / "validation_report.tsv", sep="\t", index=False)
        pd.DataFrame(
            [
                {
                    "impact_release_id": RELEASE_ID,
                    "status": "validated_result_led_spatial_candidate_awaiting_adjudication",
                    "n_panels": len(layout_rows),
                    "n_profile_rows": len(profile),
                    "n_new_main_panels": 2,
                    "n_retained_r4_panels": 12,
                    "main_claim": "positive_spatial_autocorrelation_requires_matched_null_calibration",
                    "panel_manifest_sha256": sha256_file(output / "panel_manifest.tsv"),
                    "source_manifest_sha256": sha256_file(output / "source_manifest.tsv"),
                    "validation_report_sha256": sha256_file(output / "validation_report.tsv"),
                    "layout_audit_sha256": sha256_file(output / "layout_audit.tsv"),
                    "claim_contract_sha256": sha256_file(output / "impact_claim_contract.tsv"),
                    "new_inferential_family": "FALSE",
                    "canonical_written": "FALSE",
                }
            ]
        ).to_csv(output / "READY", sep="\t", index=False)
        print("PASS result-led spatial impact figures: 14 panels; matched-null contrast rederived", flush=True)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
