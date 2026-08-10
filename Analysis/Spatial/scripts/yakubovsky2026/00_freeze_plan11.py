#!/usr/bin/env python3
"""Freeze Plan 11 inputs/specification before binary lipid labels are opened.

``--registry-only --check-only`` is deliberately safe before the Plan 10 gate:
it validates only the already sealed Plan 20 registry.  Production ``--write``
also validates the one-row Plan 10 gate and hashes the row-level join file, but
does not parse that row-level file or read expression/lipid values.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from yakubovsky_common import (
    BASE_SEED,
    BLOCK_WIDTH_MULTIPLIER,
    MAX_FAILED_BOOTSTRAP_FRACTION,
    FREEZE_STAGE_ARTIFACTS,
    MIN_BLOCK_SPOTS,
    MIN_DONOR_SPOTS,
    MIN_FINAL_BLOCKS,
    MIN_PROGRAM_GENES,
    MIN_RETAINED_L1,
    N_BOOTSTRAP,
    PRIMARY_MODEL,
    PRIMARY_SPLINE_DF,
    PRODUCTION_CODE_FILENAMES,
    RELEASE_ID,
    ROBUST_RULE,
    SPLINE_SENSITIVITY_DF,
    ContractError,
    atomic_write_frame,
    bool_value,
    default_paths,
    sha256_file,
    validate_registry_contract,
    validate_analysis_freeze,
    validate_gene_axis_manifest,
    validate_source_gate_summary,
    validate_stage_seal,
    write_stage_seal,
)


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def validate_v1_preservation(path: Path, base: Path) -> tuple[int, str]:
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {
        "relative_path",
        "baseline_bytes",
        "baseline_sha256",
        "final_bytes",
        "final_sha256",
        "unchanged",
    }
    if not required.issubset(table.columns):
        raise ContractError(f"Plan 20 v1 preservation table lacks {sorted(required)}")
    if len(table) != 44:
        raise ContractError(f"Expected 44 protected v1 files, found {len(table)}")
    for row in table.itertuples(index=False):
        if row.unchanged != "TRUE" or row.baseline_sha256 != row.final_sha256:
            raise ContractError(f"Plan 20 reports v1 drift for {row.relative_path}")
        target = base / row.relative_path
        if not target.is_file():
            raise ContractError(f"Protected v1 file is missing: {target}")
        if target.stat().st_size != int(row.final_bytes):
            raise ContractError(f"Protected v1 byte-size drift: {target}")
        if sha256_file(target) != row.final_sha256:
            raise ContractError(f"Protected v1 hash drift: {target}")
    return len(table), sha256_file(path)


def validate_gene_axis_verification(
    path: Path,
    *,
    script_dir: Path,
    v_mat_sha256: str,
    manifest_sha256: str,
    common_sha256: str,
    n_donors: int,
    n_genes: int,
) -> str:
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {
        "release_id",
        "status",
        "n_donors",
        "n_genes",
        "common_gene_axis_sha256",
        "v_mat_sha256",
        "manifest_sha256",
        "producer",
        "producer_sha256",
        "matrix_values_read",
        "lipid_fields_read",
    }
    if len(table) != 1 or not required.issubset(table.columns):
        raise ContractError("Gene-axis verification artifact is not one complete row")
    row = table.iloc[0]
    verifier = (script_dir / "07_verify_gene_axis.py").resolve()
    expected = {
        "release_id": RELEASE_ID,
        "status": "exact_all_donor_gene_axes_verified",
        "n_donors": str(n_donors),
        "n_genes": str(n_genes),
        "common_gene_axis_sha256": common_sha256,
        "v_mat_sha256": v_mat_sha256,
        "manifest_sha256": manifest_sha256,
        "producer": str(verifier),
        "producer_sha256": sha256_file(verifier),
    }
    for field, value in expected.items():
        if row[field] != value:
            raise ContractError(
                f"Gene-axis verification artifact drift for {field}: "
                f"{row[field]!r} != {value!r}"
            )
    if bool_value(row["matrix_values_read"]) or bool_value(row["lipid_fields_read"]):
        raise ContractError("Gene-axis verification artifact is not outcome-blind")
    return sha256_file(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--registry-only", action="store_true")
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.check_only == args.write:
        raise SystemExit("Choose exactly one of --check-only or --write")
    if args.write and args.registry_only:
        raise SystemExit("--registry-only is available only with --check-only")

    paths = default_paths(args.base)
    output = (args.output_dir or paths["candidate"]).resolve()
    gene_axis_verification = output / "gene_axis_verification.tsv"
    registry = validate_registry_contract(paths)
    result: dict[str, object] = {
        "release_id": RELEASE_ID,
        "mode": "registry_only" if args.registry_only else "full_source_gate",
        "write_requested": args.write,
        "external_programs": len(registry.family),
        "external_modules": sorted(
            pd.to_numeric(registry.family["module"]).astype(int).tolist()
        ),
        "registry_sha256": registry.registry_sha256,
        "membership_sha256": registry.membership_sha256,
        "hotspot_validation_sha256": registry.validation_sha256,
        "hotspot_release_manifest_sha256": registry.release_manifest_sha256,
        "hotspot_ready_sha256": registry.ready_sha256,
        "external_outcomes_read": False,
    }
    if args.registry_only:
        print(json.dumps(result, indent=2, sort_keys=True))
        return

    source_gate = validate_source_gate_summary(paths["source_gate"])
    required_source = [
        paths["source_join"],
        paths["source_files"],
        paths["source_git"],
        paths["v_mat"],
        paths["gencode"],
        paths["source_importer"],
        paths["source_background"],
        paths["v1_preservation"],
        paths["gene_axis_manifest"],
        gene_axis_verification,
    ]
    missing = [str(path) for path in required_source if not path.is_file()]
    if missing:
        raise ContractError("Missing production input(s): " + ", ".join(missing))
    n_v1, v1_table_hash = validate_v1_preservation(
        paths["v1_preservation"], paths["base"]
    )
    v_mat_sha256 = sha256_file(paths["v_mat"])
    gene_axis = validate_gene_axis_manifest(
        paths["gene_axis_manifest"], v_mat_sha256
    )
    script_dir = Path(__file__).resolve().parent
    gene_axis_verification_sha256 = validate_gene_axis_verification(
        gene_axis_verification,
        script_dir=script_dir,
        v_mat_sha256=v_mat_sha256,
        manifest_sha256=gene_axis.manifest_sha256,
        common_sha256=gene_axis.common_sha256,
        n_donors=len(gene_axis.rows),
        n_genes=gene_axis.n_genes,
    )
    result.update(
        {
            "source_gate": source_gate.row["terminal_verdict"],
            "passing_binary_lipid_donors": list(source_gate.passing_donors),
            "source_gate_sha256": source_gate.source_gate_sha256,
            "source_join_sha256": sha256_file(paths["source_join"]),
            "source_files_sha256": sha256_file(paths["source_files"]),
            "source_git_manifest_sha256": sha256_file(paths["source_git"]),
            "v_mat_sha256": v_mat_sha256,
            "gene_axis_manifest_sha256": gene_axis.manifest_sha256,
            "common_gene_axis_sha256": gene_axis.common_sha256,
            "gene_axis_verification_sha256": gene_axis_verification_sha256,
            "gencode_sha256": sha256_file(paths["gencode"]),
            "protected_v1_files": n_v1,
            "v1_preservation_table_sha256": v1_table_hash,
            "binary_labels_opened": False,
            "continuous_lipid_fields_read": False,
        }
    )
    if args.check_only:
        freeze_path = output / "analysis_freeze_manifest.tsv"
        if freeze_path.is_file():
            validate_analysis_freeze(output)
            validate_stage_seal(output, "freeze", FREEZE_STAGE_ARTIFACTS)
            result["pre_outcome_producer_seal"] = "rederived"
            result["analysis_freeze_manifest_sha256"] = sha256_file(freeze_path)
        else:
            result["pre_outcome_producer_seal"] = "not_yet_written"
        print(json.dumps(result, indent=2, sort_keys=True))
        return

    output.mkdir(parents=True, exist_ok=True)
    protected_outputs = [
        output / "analysis_freeze_manifest.tsv",
        output / "analysis_specification.tsv",
        output / "registry_manifest.tsv",
        output / "freeze_execution_manifest.tsv",
        output / "freeze_stage_seal.tsv",
    ]
    existing = [str(path) for path in protected_outputs if path.exists()]
    if existing:
        raise ContractError(
            "Refusing to overwrite frozen Plan 11 pre-outcome artifact(s): "
            + ", ".join(existing)
        )

    frozen_at = utc_now()
    manifest_rows = []
    file_inputs = {
        "plan10_source_gate": paths["source_gate"],
        "plan10_authoritative_join": paths["source_join"],
        "plan10_source_files": paths["source_files"],
        "plan10_source_git_manifest": paths["source_git"],
        "yakubovsky_v_mat": paths["v_mat"],
        "gencode_v49": paths["gencode"],
        "yakubovsky_source_importer": paths["source_importer"],
        "yakubovsky_background_estimator": paths["source_background"],
        "yakubovsky_gene_axis_equivalence": paths["gene_axis_manifest"],
        "yakubovsky_gene_axis_verification": gene_axis_verification,
        "plan20_registry": paths["registry"],
        "plan20_membership": paths["membership"],
        "plan20_external_family": paths["external"],
        "plan20_tested_universe": paths["tested_universe"],
        "plan20_validation": paths["hotspot_validation"],
        "plan20_release_manifest": paths["hotspot_release_manifest"],
        "plan20_ready": paths["hotspot_ready"],
        "plan20_v1_preservation": paths["v1_preservation"],
    }
    for filename in PRODUCTION_CODE_FILENAMES:
        file_inputs[f"producer_code::{filename}"] = script_dir / filename
    for role, path in file_inputs.items():
        manifest_rows.append(
            {
                "release_id": RELEASE_ID,
                "record_type": "file_input",
                "name": role,
                "value": str(path),
                "sha256": sha256_file(path),
                "frozen_at_utc": frozen_at,
                "lipid_outcome_read": False,
            }
        )
    scalar_inputs = {
        "source_gate_verdict": source_gate.row["terminal_verdict"],
        "passing_binary_lipid_donors": ";".join(source_gate.passing_donors),
        "test_family_program_uids": ";".join(registry.family["program_uid"]),
        "expected_directions": ";".join(
            f"{row.program_uid}:{row.expected_direction}"
            for row in registry.family.itertuples(index=False)
        ),
        "primary_model": PRIMARY_MODEL,
        "primary_spline_df": str(PRIMARY_SPLINE_DF),
        "spline_sensitivity_df": ";".join(map(str, SPLINE_SENSITIVITY_DF)),
        "model_variants": (
            "primary;equal_weight;leave_top_gene;zonation_df3;zonation_df5;"
            "zonation_landmark_overadjusted"
        ),
        "zonation_landmark_interpretation": (
            "LM_pp_ind_and_LM_pc_ind_are_zonation_landmarks_not_hepatocyte_"
            "identity_purity_or_eligibility"
        ),
        "n_bootstrap": str(N_BOOTSTRAP),
        "base_seed": str(BASE_SEED),
        "minimum_program_genes": str(MIN_PROGRAM_GENES),
        "minimum_retained_l1_weight": str(MIN_RETAINED_L1),
        "minimum_donor_spots": str(MIN_DONOR_SPOTS),
        "minimum_final_blocks": str(MIN_FINAL_BLOCKS),
        "minimum_block_spots": str(MIN_BLOCK_SPOTS),
        "block_width_nearest_neighbor_multiplier": str(BLOCK_WIDTH_MULTIPLIER),
        "maximum_failed_bootstrap_fraction": str(MAX_FAILED_BOOTSTRAP_FRACTION),
        "multiplicity": "BH_once_across_complete_testable_frozen_family",
        "effect_unit": "source_lipid_zone_minus_non_lipid_zone_score",
        "robust_rule": ROBUST_RULE,
        "prohibited_interpretation": (
            "continuous_gradient;ordered_trend;high_vs_low_severity;dose_response;"
            "lipid_free_complement"
        ),
    }
    for name, value in scalar_inputs.items():
        manifest_rows.append(
            {
                "release_id": RELEASE_ID,
                "record_type": "frozen_specification",
                "name": name,
                "value": value,
                "sha256": "",
                "frozen_at_utc": frozen_at,
                "lipid_outcome_read": False,
            }
        )
    atomic_write_frame(output / "analysis_freeze_manifest.tsv", pd.DataFrame(manifest_rows))

    specification = pd.DataFrame(
        [
            {"parameter": name, "value": value, "frozen_at_utc": frozen_at}
            for name, value in scalar_inputs.items()
        ]
    )
    atomic_write_frame(output / "analysis_specification.tsv", specification)

    registry_manifest = registry.family.copy()
    registry_manifest.insert(0, "candidate_release_id", RELEASE_ID)
    registry_manifest["registry_sha256_rederived"] = registry.registry_sha256
    registry_manifest["membership_table_sha256_rederived"] = registry.membership_sha256
    registry_manifest["frozen_before_binary_labels_opened"] = True
    registry_manifest["frozen_at_utc"] = frozen_at
    atomic_write_frame(output / "registry_manifest.tsv", registry_manifest)

    environment = pd.DataFrame(
        [
            {
                "release_id": RELEASE_ID,
                "stage": "outcome_blind_freeze",
                "python": sys.version.replace("\n", " "),
                "platform": platform.platform(),
                "conda_prefix": os.environ.get("CONDA_PREFIX", ""),
                "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
                "completed_utc": utc_now(),
                "binary_labels_opened": False,
                "continuous_lipid_fields_read": False,
            }
        ]
    )
    atomic_write_frame(output / "freeze_execution_manifest.tsv", environment)
    write_stage_seal(output, "freeze", FREEZE_STAGE_ARTIFACTS)
    result["analysis_freeze_manifest_sha256"] = sha256_file(
        output / "analysis_freeze_manifest.tsv"
    )
    validate_analysis_freeze(output)
    result["pre_outcome_producer_seal"] = "written_and_rederived"
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
