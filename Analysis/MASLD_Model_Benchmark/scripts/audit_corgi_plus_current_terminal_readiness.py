#!/usr/bin/env python3
"""Re-audit current Corgi+ release readiness without opening biological outcomes."""

from __future__ import annotations

import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path
import tomllib
from typing import Any


class CorgiPlusCurrentReadinessError(RuntimeError):
    """Raised when current primary sources or frozen authorities differ."""


CHECKPOINT_SHA256 = "01855614bfeffa72d5bcaaa28729c022ccfbc8f5fdcc0add8e04ba4533e94ccd"
CHECKPOINT_MD5 = "f7dbe2ca301fe64814bebcd4fc3377bd"
CHECKPOINT_SIZE = 2_384_030_158
REGULAR_CHECKPOINT_SHA256 = "cd51539f01de1e66a3a4f9b89e772bdeae07df2f0178d5b692f2368db2f4d2a4"

AUTHORITIES = {
    "corgi_plus_checkpoint": (
        "config/artifacts/models/corgi_plus/checkpoints.json",
        "936718cd13ffd621a4d77e2a340e9ede5f328208da4f72055a41f8c1b87af50f",
    ),
    "corgi_plus_exposure": (
        "config/artifacts/models/corgi_plus/exposure_audit.json",
        "60c176303db6beb14d1b20519ada8d7aace8e2cdba94703a85c5e4428a864148",
    ),
    "corgi_plus_crosswalk": (
        "config/artifacts/models/corgi_plus/development_crosswalk.json",
        "af1e5f6c53fa3cfc729faa14be549f259f8e2f08e0e0378d587865a48c28ae8c",
    ),
    "corgi_plus_prior_disposition": (
        "config/artifacts/models/corgi_plus/release_disposition_20260824.json",
        "43210eabff99797a7c026e52f5ec19ead336424ba0ff13a237e71203232a173f",
    ),
    "regular_corgi_checkpoint": (
        "config/artifacts/models/corgi/checkpoints.json",
        "1df7633507cd716f85ab03ab7ab6a171e3e56c0e3e1f29f81b04f565134a089e",
    ),
    "regular_corgi_smoke": (
        "config/artifacts/models/corgi/bounded_prediction_smoke_disposition_20260824.json",
        "92543e0ee28e0cdbdd9e2bb13a32b0b4f2960a9320dc417aa881d4d1f8f63635",
    ),
    "context_borzoi_interface": (
        "config/artifacts/models/context_borzoi/component_interface_readiness_20260825.json",
        "8195a56a52468077d089942b1cd7dec93ce1b728216a94208d07c0c13e231f70",
    ),
    "context_borzoi_fixture": (
        "config/artifacts/models/context_borzoi/rna_context_adapter_fixture_20260825.json",
        "1b92694d16ba9c48ab0d1063aa7f362d60c8bfd95337095ffdbf59d8538be15c",
    ),
    "resource_atlas": (
        "config/datasets/resource_atlas_current.toml",
        "ce711352dbeb1c95fb997d615308bbfcf7f08d23aed67622a584f99785ae2baa",
    ),
    "bulk_five_cohort": (
        "config/datasets/bulk_five_cohort.toml",
        "566c96580a5895a66c0d941588a0fde2fe142616806fa8d234932e391cd99990",
    ),
    "gse244832": (
        "config/datasets/gse244832.toml",
        "8356e8badbfeb6838519e772f4707330dede879ccb77f0d93ac8456d0dcfe302",
    ),
    "gse296875": (
        "config/datasets/gse296875.toml",
        "9465b012bbf46611641555a5acae4d1f033d053c900bf944dcec02979c735b6f",
    ),
    "gse281367": (
        "config/datasets/gse281367.toml",
        "55986a4c4ee43d08abd5517cc4f95bc95a21c796597d70a703acb0e11d33fdca",
    ),
    "rna_conditioned_capabilities": (
        "config/evaluation/rna_conditioned_atac_capabilities.toml",
        "c779f06dedf69a15ec44280323a7701e0938b75bfa2660d11dd2a6f5de37868d",
    ),
}

FROZEN_ARTIFACTS = {
    "official_acquisition": (
        "executions/upstream-acquisitions/corgi-d84c000-zenodo18630048-21064463",
        "73aad9ab9ad3e080a2ccaa0cee1a1d727d126a454f0e69f74f343a900b168c0f",
    ),
    "checkpoint_schema": (
        "executions/corgi-official-checkpoint-schema-21064532",
        "312bb3c90ae6da3eff20f43d920b03fbf6a5b05b2b480bbc98d24aaad71d868e",
    ),
    "l40s_runtime": (
        "executions/corgi-native-runtime-l40s-r8-21065022",
        "09d8811b612cba4f75361feeba29bd0676efde83f0e4785a37990c2d2f52a0a2",
    ),
    "regular_prediction_smoke": (
        "executions/corgi-hepatocyte-tile-smoke-bundle-v3-21069656",
        "bf53c7a23a1d023442a3de64150fd98799a5be3991e95e6a2f4057f6998d4eda",
    ),
}


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file() or path.stat().st_size > 5_000_000:
        raise CorgiPlusCurrentReadinessError(f"unsafe or unexpected JSON: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CorgiPlusCurrentReadinessError(f"invalid JSON: {path}") from error
    if not isinstance(value, dict):
        raise CorgiPlusCurrentReadinessError(f"expected JSON object: {path}")
    return value


def load_toml(path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file() or path.stat().st_size > 2_000_000:
        raise CorgiPlusCurrentReadinessError(f"unsafe or unexpected TOML: {path}")
    try:
        value = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise CorgiPlusCurrentReadinessError(f"invalid TOML: {path}") from error
    if not isinstance(value, dict):
        raise CorgiPlusCurrentReadinessError(f"expected TOML table: {path}")
    return value


def verify_authorities(root: Path) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for name, (relative, expected) in AUTHORITIES.items():
        observed = digest(root / relative)
        if observed != expected:
            raise CorgiPlusCurrentReadinessError(
                f"authority differs: {relative}: {observed}"
            )
        result[name] = {"path": relative, "sha256": expected}
    return result


def verify_frozen_manifests(root: Path) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for name, (relative, expected) in FROZEN_ARTIFACTS.items():
        artifact = root / relative
        if not (artifact / "COMPLETE").is_file():
            raise CorgiPlusCurrentReadinessError(f"artifact is incomplete: {relative}")
        manifest = artifact / "ARTIFACTS.json"
        if digest(manifest) != expected:
            raise CorgiPlusCurrentReadinessError(f"artifact differs: {relative}")
        result[name] = load_json(manifest).get("metadata", {})
    return result


def class_methods_raise_not_implemented(
    source: str, class_name: str, method_names: tuple[str, ...]
) -> bool:
    tree = ast.parse(source)
    class_node = next(
        (
            node
            for node in tree.body
            if isinstance(node, ast.ClassDef) and node.name == class_name
        ),
        None,
    )
    if class_node is None:
        return False
    methods = {
        node.name: node
        for node in class_node.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    for name in method_names:
        method = methods.get(name)
        if method is None:
            return False
        raises = [node for node in ast.walk(method) if isinstance(node, ast.Raise)]
        if not any(
            isinstance(node.exc, ast.Call)
            and isinstance(node.exc.func, ast.Name)
            and node.exc.func.id == "NotImplementedError"
            for node in raises
        ):
            return False
    return True


def zenodo_terms_and_files(
    record: dict[str, Any], expected_id: int
) -> dict[str, Any]:
    if record.get("id") != expected_id:
        raise CorgiPlusCurrentReadinessError("Zenodo record identity differs")
    metadata = record.get("metadata", {})
    license_value = metadata.get("license")
    license_id = (
        license_value.get("id") if isinstance(license_value, dict) else license_value
    )
    files = {
        str(row.get("key")): {
            "size_bytes": row.get("size"),
            "checksum": row.get("checksum"),
        }
        for row in record.get("files", [])
        if isinstance(row, dict)
    }
    checkpoint = files.get("corgiplus_model.pt")
    if (
        str(license_id).lower() != "cc-by-4.0"
        or checkpoint is None
        or checkpoint["size_bytes"] != CHECKPOINT_SIZE
        or checkpoint["checksum"] != f"md5:{CHECKPOINT_MD5}"
    ):
        raise CorgiPlusCurrentReadinessError(
            f"Zenodo Corgi+ record differs: {expected_id}"
        )
    return {
        "record_id": expected_id,
        "created": record.get("created"),
        "updated": record.get("updated"),
        "license": str(license_id),
        "files": files,
    }


def fixture_decisions() -> dict[str, dict[str, Any]]:
    return {
        "gse296875": {
            "topology": "same_nucleus_RNA_ATAC",
            "global_trans_context": "derivable_from_gene_counts_inside_training_folds",
            "local_stranded_total_RNA": "not_native_to_3prime_single_nucleus_gene_expression",
            "native_22_track_mean_baseline": "unavailable",
            "corgi_plus_family_native_input_complete": False,
        },
        "gse244832": {
            "topology": "same_sample_different_aliquot_RNA_ATAC",
            "false_cell_pairing_forbidden": True,
            "native_22_track_mean_baseline": "unavailable",
            "corgi_plus_family_native_input_complete": False,
        },
        "resource_atlas_current": {
            "topology": "single_cell_RNA_multi_source",
            "available_input": "gene_by_cell_UMI_counts",
            "local_stranded_total_RNA": "structurally_missing",
            "corgi_plus_family_native_input_complete": False,
        },
        "bulk_five_cohort": {
            "topology": "bulk_RNA_same_study_unpaired",
            "cell_state_specific_paired_ATAC": "absent",
            "native_22_track_mean_baseline": "unavailable",
            "corgi_plus_family_native_input_complete": False,
        },
        "gse281367": {
            "topology": "ATAC_only",
            "RNA": "structurally_missing",
            "corgi_plus_family_native_input_complete": False,
        },
    }


def audit(args: argparse.Namespace) -> dict[str, Any]:
    root = args.project_root.resolve(strict=True)
    source = args.source_cache.resolve(strict=True)
    if args.output.exists() or args.output.is_symlink():
        raise CorgiPlusCurrentReadinessError("output already exists")
    authority_receipts = verify_authorities(root)
    frozen = verify_frozen_manifests(root)

    plus = load_json(root / AUTHORITIES["corgi_plus_checkpoint"][0])
    prior = load_json(root / AUTHORITIES["corgi_plus_prior_disposition"][0])
    regular = load_json(root / AUTHORITIES["regular_corgi_smoke"][0])
    context_borzoi = load_json(root / AUTHORITIES["context_borzoi_interface"][0])
    if (
        plus.get("checkpoint", {}).get("sha256") != CHECKPOINT_SHA256
        or plus.get("checkpoint", {}).get("size_bytes") != CHECKPOINT_SIZE
        or plus.get("architecture_contract", {}).get("checkpoint_state_key_count")
        != 190
        or plus.get("architecture_contract", {}).get("checkpoint_state_numel")
        != 198_664_365
        or prior.get("status") != "TERMINAL_BLOCKED_CURRENT_UPSTREAM_RELEASE"
    ):
        raise CorgiPlusCurrentReadinessError("frozen Corgi+ identity differs")
    if (
        frozen["official_acquisition"].get("checkpoint_deserialized")
        or frozen["checkpoint_schema"].get("safe_weights_only_load") is not True
        or frozen["l40s_runtime"].get("plus_strict_restore") is not True
        or frozen["l40s_runtime"].get("corgi_plus_released_input_bundle_complete")
        is not False
    ):
        raise CorgiPlusCurrentReadinessError("frozen compatibility evidence differs")

    corgi_repo = load_json(source / "corgi.repo.json")
    corgi_commit = load_json(source / "corgi.commit.json")
    corgi_tags = json.loads((source / "corgi.tags.json").read_text(encoding="utf-8"))
    reproduction_repo = load_json(source / "reproduction.repo.json")
    reproduction_commit = load_json(source / "reproduction.commit.json")
    if not isinstance(corgi_tags, list):
        raise CorgiPlusCurrentReadinessError("GitHub tags response differs")
    current_code_revision = corgi_commit.get("sha")
    current_reproduction_revision = reproduction_commit.get("sha")
    if (
        corgi_repo.get("full_name") != "ekinda/corgi"
        or reproduction_repo.get("full_name") != "ekinda/corgi-reproduction"
        or not isinstance(current_code_revision, str)
        or len(current_code_revision) != 40
        or not isinstance(current_reproduction_revision, str)
        or len(current_reproduction_revision) != 40
    ):
        raise CorgiPlusCurrentReadinessError("GitHub repository identity differs")

    code_dir = source / "corgi-current"
    license_text = (code_dir / "LICENSE").read_text(encoding="utf-8")
    pyproject_text = (code_dir / "pyproject.toml").read_text(encoding="utf-8")
    predict_text = (code_dir / "corgi/predict.py").read_text(encoding="utf-8")
    config_plus_text = (code_dir / "corgi/config_corgiplus.py").read_text(
        encoding="utf-8"
    )
    trainer_plus_text = (code_dir / "corgi/trainer_corgiplus.py").read_text(
        encoding="utf-8"
    )
    reproduction_text = (
        source / "reproduction-current/code/corgiplus_delta_bench_tta.py"
    ).read_text(encoding="utf-8")
    placeholder = class_methods_raise_not_implemented(
        predict_text,
        "corgiplus_pretrained",
        ("predict", "predict_regions", "predict_regions_with_bigwig"),
    )
    if (
        "Apache License" not in license_text
        or 'license = {text = "MIT"}' not in pyproject_text
        or not placeholder
        or "corgiplusplus_rna" not in trainer_plus_text
        or "/project/" not in config_plus_text
        or "/project/" not in reproduction_text
    ):
        raise CorgiPlusCurrentReadinessError(
            "current source no longer matches the terminal blocker surface"
        )

    zenodo = {
        str(record_id): zenodo_terms_and_files(
            load_json(source / f"zenodo-{record_id}.json"), record_id
        )
        for record_id in (18_630_048, 20_686_797)
    }
    paper_record_files = zenodo["20686797"]["files"]
    if "source_data.zip" not in paper_record_files:
        raise CorgiPlusCurrentReadinessError("paper source-data archive disappeared")

    datasets = {
        name: load_toml(root / AUTHORITIES[name][0])
        for name in (
            "resource_atlas",
            "bulk_five_cohort",
            "gse244832",
            "gse296875",
            "gse281367",
        )
    }
    if (
        datasets["resource_atlas"].get("modalities") != ["single_cell_rna"]
        or datasets["bulk_five_cohort"].get("modalities") != ["bulk_rna"]
        or datasets["gse244832"].get("pairing_levels")
        != ["same_sample_different_aliquot"]
        or datasets["gse296875"].get("pairing_levels") != ["same_nucleus"]
        or datasets["gse281367"].get("modalities") != ["single_nucleus_atac"]
    ):
        raise CorgiPlusCurrentReadinessError("dataset topology authority differs")

    fixtures = fixture_decisions()
    if any(
        row["corgi_plus_family_native_input_complete"] for row in fixtures.values()
    ):
        raise CorgiPlusCurrentReadinessError("a fixture was incorrectly activated")
    if (
        regular.get("model_identity", {}).get("checkpoint_sha256")
        != REGULAR_CHECKPOINT_SHA256
        or regular.get("status")
        != "BOUNDED_OUTCOME_FREE_PREDICTION_SMOKE_COMPLETE"
        or regular.get("tournament_disposition", {}).get("native_lane_admitted")
        or context_borzoi.get("conditional_slot", {}).get("architecture_built")
        or context_borzoi.get("execution_disposition", {}).get("training_authorized")
    ):
        raise CorgiPlusCurrentReadinessError("substitute-lane disposition differs")

    result = {
        "schema_version": "masld-bench-corgi-plus-current-terminal-readiness-v1",
        "audit_status": "passed_terminal_blocked",
        "current_primary_sources": {
            "corgi": {
                "repository": "ekinda/corgi",
                "default_branch": corgi_repo.get("default_branch"),
                "current_head": current_code_revision,
                "prior_audit_head": "d84c00082da3d04eb05169e290e53f1c7da2b9c5",
                "tags": [row.get("name") for row in corgi_tags if isinstance(row, dict)],
                "high_level_corgi_plus_methods_are_placeholders": placeholder,
                "absolute_project_paths_remain": True,
            },
            "reproduction": {
                "repository": "ekinda/corgi-reproduction",
                "default_branch": reproduction_repo.get("default_branch"),
                "current_head": current_reproduction_revision,
                "prior_audit_head": "2367c42cf36630b7cf561eccc2e781abc22b36be",
                "absolute_project_paths_remain": True,
            },
            "zenodo": zenodo,
        },
        "checkpoint": {
            "sha256": CHECKPOINT_SHA256,
            "md5": CHECKPOINT_MD5,
            "size_bytes": CHECKPOINT_SIZE,
            "state_keys": 190,
            "state_numel": 198_664_365,
            "auxiliary_width": 26,
            "safe_schema_inventory_passed": True,
            "strict_restore_passed": True,
            "synthetic_full_window_forward_passed": True,
            "shortened_backward_passed": True,
            "checkpoint_producing_source_revision": "UNRESOLVED",
            "biological_numeric_parity": False,
            "checkpoint_resume": False,
        },
        "licenses": {
            "repository_root": "Apache-2.0",
            "pyproject_metadata": "MIT",
            "code_conflict": "UNRESOLVED",
            "checkpoint_records": "CC-BY-4.0",
            "derivative_weight_redistribution": "REQUIRES_UPSTREAM_OR_INSTITUTIONAL_DETERMINATION",
            "open_champion_eligible": False,
        },
        "native_input_contract": {
            "DNA_bp": 524_288,
            "global_trans_regulators": 2_891,
            "local_auxiliary_bins": 8_192,
            "local_auxiliary_resolution_bp": 64,
            "local_auxiliary_channels": 26,
            "channel_families": {
                "stranded_local_total_RNA": 2,
                "training_context_mean_baseline_in_native_output_order": 22,
                "local_RNA_delta_from_corresponding_baseline_RNA": 2,
            },
            "channel_order_checkpoint_bound": False,
            "native_input_bundle_released": False,
            "gene_count_matrix_is_valid_local_RNA": False,
            "missing_local_RNA_may_be_zero_filled": False,
        },
        "native_output_contract": {
            "channels": 22,
            "central_bins": 6_144,
            "resolution_bp": 64,
            "activation": "Conv1d_then_Softplus",
            "RNA_channels_given_as_input_are_scored_as_outcomes": False,
        },
        "reference_contract": {
            "native_build_label": "hg38",
            "checkpoint_bound_modified_FASTA": "UNRESOLVED",
            "benchmark_reference": "GRCh38.p14_plus_GENCODE_v49",
            "local_RNA_baseline_output_crop_shift_and_reverse_complement_parity": "UNRESOLVED",
            "biological_execution_gate": "blocked",
        },
        "project_fixture_fit": fixtures,
        "terminal_disposition": {
            "status": "TERMINAL_BLOCKED_CURRENT_UPSTREAM_RELEASE",
            "rna_conditioned_atac_execution_authorized": False,
            "model_training_or_adaptation_authorized": False,
            "production_or_GPU_item_authorized": False,
            "comparative_claim_eligible": False,
            "reevaluation_trigger": "new_immutable_upstream_release_or_clarification_closing_source_input_loader_reference_numeric_and_terms_gates",
        },
        "nearest_executable_lane": {
            "model_id": "corgi_regular",
            "identity": "restricted_development_comparator_not_Corgi_plus",
            "checkpoint_sha256": REGULAR_CHECKPOINT_SHA256,
            "outcome_free_GSE296875_prediction_smoke_complete": True,
            "resubmission_needed": False,
            "development_evaluation_authorized_by_this_audit": False,
            "open_champion_eligible": False,
            "relabel_as_Corgi_plus": False,
            "why_nearest": "same_long_range_DNA_global_FiLM_transformer_and_22_output_family_without_missing_local_auxiliary_bundle",
        },
        "context_borzoi_boundary": {
            "status": "interface_only_build_blocked",
            "architecture_built": False,
            "component_selection_locked": False,
            "training_authorized": False,
            "nearest_executable_now": False,
        },
        "next_action": "retain_Corgi_plus_as_terminal_blocked_use_existing_regular_Corgi_outcome_free_comparator_and_do_not_build_context_Borzoi_until_complementarity_selection_passes",
        "authority_hashes": authority_receipts,
        "frozen_artifact_manifests": {
            name: {"path": relative, "artifacts_sha256": expected}
            for name, (relative, expected) in FROZEN_ARTIFACTS.items()
        },
        "checkpoint_deserialized_by_this_audit": False,
        "model_forward_executed_by_this_audit": False,
        "biological_data_read": False,
        "development_outcomes_read": False,
        "sealed_features_read": False,
        "sealed_outcomes_read": False,
        "gpu_job_submitted": False,
    }
    args.output.mkdir(parents=True, exist_ok=False, mode=0o750)
    (args.output / "audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--source-cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    print(json.dumps(audit(parse_args()), sort_keys=True))


if __name__ == "__main__":
    main()
