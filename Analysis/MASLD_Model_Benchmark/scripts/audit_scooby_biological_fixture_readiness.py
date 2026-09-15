#!/usr/bin/env python3
"""Resolve released Scooby biological-fixture readiness without reading outcomes."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import tomllib
from typing import Any


class ScoobyBiologicalReadinessError(RuntimeError):
    """Raised when a frozen authority differs from the fail-closed requirement."""


CODE_REVISION = "dbfe6cd150d018d46a6d24ffa0406d7f2e89b471"
CODE_LICENSE_SHA256 = "0a09080986ba1b09bc670290c172fc831a1f73dff1f0c6ab23039b79f80d8ed3"

CHECKPOINTS = {
    "scooby_onek1k": {
        "repository": "lauradmartens/onek1k-scooby",
        "revision": "dd7a2b7be54c906b43de4bc0159b169caaf057c5",
        "sha256": "2ae335c651182040462d965bb879118be40a6cb174550a3686fb10b04a74f7cb",
        "size_bytes": 765_092_032,
        "context_dim": 10,
        "tracks": ["rna_plus", "rna_minus"],
        "context_encoder": "scPoli",
        "context_topology": "RNA_only_OneK1K_PBMC",
        "pinned_card_license": "mit",
    },
    "scooby_epicardioids": {
        "repository": "lauradmartens/epicardioids-scooby",
        "revision": "cefc3a626f5080baaf83878754ea924df6bd7e8d",
        "sha256": "f945e1937a46ddd5cb31249ef71e721e7a7b1381aa207fb2de449f23ac75a9e8",
        "size_bytes": 772_988_628,
        "context_dim": 50,
        "tracks": ["rna_plus", "rna_minus", "atac_insertion"],
        "context_encoder": "scGLUE",
        "context_topology": "same_study_unpaired_RNA_ATAC_pseudomatches",
        "pinned_card_license": "mit",
    },
    "scooby_neurips": {
        "repository": "johahi/neurips-scooby",
        "revision": "1591d15696c65f14053e4495e52c74812a86fd8d",
        "sha256": "8f2e2de3e86016378116c90accf94572e75963d5445beb6cabd585c8372170dd",
        "size_bytes": 772_970_196,
        "context_dim": 14,
        "tracks": ["rna_plus", "rna_minus", "atac_insertion"],
        "context_encoder": "custom_Poisson_MultiVI",
        "context_topology": "same_nucleus_RNA_ATAC",
        "pinned_card_license": None,
    },
}

CURRENT_AUTHORITIES = {
    "scooby_epicardioids_checkpoint": (
        "config/artifacts/models/scooby_epicardioids/checkpoints.json",
        "bc125a6f3e146dcc891f6f3327bd6c86be506c67b41e1327696225f589a19520",
    ),
    "scooby_epicardioids_crosswalk": (
        "config/artifacts/models/scooby_epicardioids/development_crosswalk.json",
        "32a6e8c2eef0de1004719cfc1b41c93609ad9adc75e7d2911502e34456994c52",
    ),
    "scooby_epicardioids_exposure": (
        "config/artifacts/models/scooby_epicardioids/exposure_audit.json",
        "8f3cad157e701cb8e11c34ae7e4d8cf46a23b7b9426b1fd0fb138d880b200881",
    ),
    "scooby_neurips_checkpoint": (
        "config/artifacts/models/scooby_neurips/checkpoints.json",
        "d76f024290f2b4d7c7ef8da2d2e915a3f87f8e337d3413949d230347f585d412",
    ),
    "scooby_onek1k_checkpoint": (
        "config/artifacts/models/scooby_onek1k/checkpoints.json",
        "a223a282bf7030ce2bc8d754854d434ddd6d46f30052393fae4df0c7e290313e",
    ),
    "canonical_borzoi_checkpoint": (
        "config/artifacts/models/borzoi/checkpoints.json",
        "1c5d7213dab05be8bb819d61ff4d03b2bcd7ca307954ab328de739c476c9e257",
    ),
    "gse296875_dataset": (
        "config/datasets/gse296875.toml",
        "9465b012bbf46611641555a5acae4d1f033d053c900bf944dcec02979c735b6f",
    ),
    "gse281367_dataset": (
        "config/datasets/gse281367.toml",
        "55986a4c4ee43d08abd5517cc4f95bc95a21c796597d70a703acb0e11d33fdca",
    ),
    "observed_multiome_capabilities": (
        "config/evaluation/observed_multiome_task_capabilities.toml",
        "6c2e92a78af315a0412a949de208daea5ec28d30fdd2e421a7762d3ffe6252b8",
    ),
    "synthetic_probe_config": (
        "config/scooby_epicardioids_full_backbone_synthetic_probe.json",
        "b9331c5e8243b7bc336d4e209d90323f63518a45b21e6654ea57cda96a48229e",
    ),
    "synthetic_probe_queue_item": (
        "config/campaigns/gpu_bundle_queue/062-model-probe-284.json",
        "403c363b00c74b7bfc46b636e17afcc7f4bdf26d92e2b620b96d9c641d267a8a",
    ),
    "synthetic_probe_wrapper": (
        "slurm/probe_scooby_full_backbone_synthetic_l40s.sbatch",
        "f94b874d19dc62093d5851b50eee1467b09c44a6227785ec3e64633d5adbc4b6",
    ),
}

FROZEN_ARTIFACTS = {
    "released_lane": (
        "executions/model-audit-201-21085741",
        "01d5dabca02430ac18fd5cf26cc9ed64c1ea1fc37baeb49506512ff2fa4d6198",
        "audit/released_lane_audit.json",
    ),
    "runtime_schema": (
        "executions/model-check-279-21097163",
        "835244d33a0fcc2808efc76551eb0f0f771e08ec17e0acf12599cd7ceb9365ea",
        "audit/receipt.json",
    ),
    "checkpoint_acquisition": (
        "executions/model-data-281-21097353",
        "657f20fb7fe110fdd4f690b874d450f73fc6fc9b3e40c0276f0bade73bc9047a",
        "audit/receipt.json",
    ),
    "decoder_probe": (
        "executions/model-check-282-21097437",
        "17aad708f568157785905fbed57dc78fc6bf0af01f7aacaa493108340615a88a",
        "audit/receipt.json",
    ),
    "native_context": (
        "executions/model-data-283-21097537",
        "4f7995b1810ab213b17adfe2ab1f6260e447446864ac1bbc07174285bfe1fe89",
        "context/receipt.json",
    ),
    "synthetic_queue_preflight": (
        "executions/model-check-285-21097650",
        "45d0918a33feb557bde3baf205d3c541791cd12cfc4575beb0af29317da3d4bb",
        "ARTIFACTS.json",
    ),
    "gse296875_mask_plan": (
        "executions/model-check-296-21099109",
        "9df7d83623d609a6fb894b921bb914168489df0b7d172c5addf794f418df84e2",
        "receipt.json",
    ),
    "gse296875_materialization_contract": (
        "executions/model-check-297-21099157",
        "99044a88a1ccf3e8b2c2a8563b80d1f05905b3e236ddec21dce0df0b6d6f4a3b",
        "receipt.json",
    ),
    "gse296875_multivi_prepare": (
        "executions/multivi-smoke-prepared-21065093",
        "6c10206033545fe7ed1ecbd9f7d64e107e2b3742de06db55a838bbdb46a1cdaa",
        "prepared/prepare_receipt.json",
    ),
    "gse281367_activation": (
        "executions/gse281367-activation-readiness-21099679",
        "e36b86815bdb9f66c372c933c2ce3052ba5c2b475138f45f2bfa2a26805f9f01",
        "task_eligibility.json",
    ),
    "gse281367_production": (
        "executions/gse281367-atac-transport-production-readiness-21099780",
        "e7769dc43d99ed286bb7e8c8b848da8a28020d56cbb1b73ebfa775a78a074801",
        "production_decision.json",
    ),
}


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file() or path.stat().st_size > 2_000_000:
        raise ScoobyBiologicalReadinessError(f"unsafe or unexpected JSON: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScoobyBiologicalReadinessError(f"invalid JSON: {path}") from error
    if not isinstance(value, dict):
        raise ScoobyBiologicalReadinessError(f"expected JSON object: {path}")
    return value


def load_toml(path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file() or path.stat().st_size > 2_000_000:
        raise ScoobyBiologicalReadinessError(f"unsafe or unexpected TOML: {path}")
    try:
        value = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise ScoobyBiologicalReadinessError(f"invalid TOML: {path}") from error
    if not isinstance(value, dict):
        raise ScoobyBiologicalReadinessError(f"expected TOML table: {path}")
    return value


def verify_current_authorities(root: Path) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for name, (relative, expected) in CURRENT_AUTHORITIES.items():
        path = root / relative
        observed = digest(path)
        if observed != expected:
            raise ScoobyBiologicalReadinessError(
                f"current authority differs: {relative}: {observed}"
            )
        result[name] = {"path": relative, "sha256": expected}
    return result


def read_frozen_member(
    root: Path, relative_root: str, expected_manifest: str, relative_member: str
) -> dict[str, Any]:
    artifact = root / relative_root
    manifest_path = artifact / "ARTIFACTS.json"
    if digest(manifest_path) != expected_manifest or not (artifact / "COMPLETE").is_file():
        raise ScoobyBiologicalReadinessError(
            f"frozen artifact identity differs: {relative_root}"
        )
    manifest = load_json(manifest_path)
    if relative_member == "ARTIFACTS.json":
        return manifest
    rows = {
        str(row.get("path")): row
        for row in manifest.get("artifacts", [])
        if isinstance(row, dict)
    }
    row = rows.get(relative_member)
    member = artifact / relative_member
    if row is None or digest(member) != row.get("sha256"):
        raise ScoobyBiologicalReadinessError(
            f"frozen member differs: {relative_root}/{relative_member}"
        )
    if member.stat().st_size != row.get("size_bytes"):
        raise ScoobyBiologicalReadinessError(
            f"frozen member size differs: {relative_root}/{relative_member}"
        )
    return load_json(member)


def hf_terms(metadata: dict[str, Any], expected_revision: str | None) -> dict[str, Any]:
    revision = metadata.get("sha")
    if expected_revision is not None and revision != expected_revision:
        raise ScoobyBiologicalReadinessError(
            f"Hugging Face revision differs: {revision!r} != {expected_revision!r}"
        )
    card = metadata.get("cardData") or {}
    siblings = metadata.get("siblings") or []
    if not isinstance(card, dict) or not isinstance(siblings, list):
        raise ScoobyBiologicalReadinessError("Hugging Face metadata structure differs")
    license_files = sorted(
        str(row.get("rfilename"))
        for row in siblings
        if isinstance(row, dict)
        and str(row.get("rfilename", "")).lower()
        in {"license", "license.md", "license.txt", "copying"}
    )
    declared = card.get("license")
    return {
        "revision": revision,
        "declared_license": declared,
        "license_files": license_files,
        "terms_declared": bool(declared or license_files),
    }


def parse_sbatch_header(source: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in source.splitlines():
        match = re.fullmatch(r"#SBATCH --([^=]+)=(.+)", line)
        if match:
            result[match.group(1)] = match.group(2)
    return result


def biological_decisions() -> dict[str, dict[str, Any]]:
    return {
        "released_checkpoint_on_gse296875": {
            "status": "blocked",
            "reason": "no_coordinate_compatible_released_query_encoder_and_inherited_weight_terms_unresolved",
            "family_native_zero_shot": False,
            "biological_execution_authorized": False,
        },
        "project_fitted_scooby_adaptation_on_gse296875": {
            "status": "conditional_not_authorized",
            "reason": "same_nucleus_topology_is_suitable_but_a_new_latent_requires_a_fold_fitted_bridge_or_decoder_adaptation_and_cleared_derivative_terms",
            "identity": "new_project_adaptation_not_released_checkpoint_inference",
            "biological_execution_authorized": False,
        },
        "released_checkpoint_on_gse281367": {
            "status": "not_applicable",
            "reason": "RNA_is_structurally_missing_and_raw_ATAC_is_not_a_checkpoint_context_vector",
            "family_native_zero_shot": False,
            "biological_execution_authorized": False,
        },
        "gse281367_atac_transport_evaluation": {
            "status": "blocked_for_production_exchange",
            "reason": "ATAC_only_development_evaluator_requires_label_agnostic_prediction_and_mask_contracts",
            "biological_execution_authorized": False,
        },
        "synthetic_epicardioids_full_backbone": {
            "status": "admitted_existing_central_dispatch_item",
            "reason": "architecture_compatibility_only_with_synthetic_sequence_and_context",
            "biological_execution_authorized": False,
        },
    }


def audit(args: argparse.Namespace) -> dict[str, Any]:
    root = args.project_root.resolve(strict=True)
    if args.output.exists() or args.output.is_symlink():
        raise ScoobyBiologicalReadinessError("output already exists")
    current = verify_current_authorities(root)
    frozen = {
        name: read_frozen_member(root, relative, manifest, member)
        for name, (relative, manifest, member) in FROZEN_ARTIFACTS.items()
    }

    model_configs = {
        model_id: load_json(root / CURRENT_AUTHORITIES[f"{model_id}_checkpoint"][0])
        for model_id in ("scooby_onek1k", "scooby_epicardioids", "scooby_neurips")
    }
    for model_id, expected in CHECKPOINTS.items():
        config = model_configs[model_id]
        checkpoint = config.get("checkpoint", {})
        architecture = config.get("architecture_contract", {})
        if (
            checkpoint.get("hugging_face_repository") != expected["repository"]
            or checkpoint.get("revision") != expected["revision"]
            or checkpoint.get("sha256") != expected["sha256"]
            or checkpoint.get("size_bytes") != expected["size_bytes"]
            or architecture.get("input_length_bp") != 524_288
            or architecture.get("output_central_bins") != 6_144
        ):
            raise ScoobyBiologicalReadinessError(
                f"checkpoint or architecture identity differs: {model_id}"
            )

    source_cache = args.source_cache.resolve(strict=True)
    pinned_terms: dict[str, dict[str, Any]] = {}
    current_terms: dict[str, dict[str, Any]] = {}
    for model_id, expected in CHECKPOINTS.items():
        pinned_terms[model_id] = hf_terms(
            load_json(source_cache / f"{model_id}.pinned.json"), expected["revision"]
        )
        if pinned_terms[model_id]["declared_license"] != expected["pinned_card_license"]:
            raise ScoobyBiologicalReadinessError(
                f"pinned checkpoint-card terms differ: {model_id}"
            )
        current_terms[model_id] = hf_terms(
            load_json(source_cache / f"{model_id}.current.json"), None
        )
    base_converter = hf_terms(load_json(source_cache / "borzoi_converter.current.json"), None)
    code_license = source_cache / "scooby.LICENSE"
    if digest(code_license) != CODE_LICENSE_SHA256:
        raise ScoobyBiologicalReadinessError("Scooby code license differs")

    released = frozen["released_lane"]
    if (
        released.get("status") != "pass"
        or released.get("checkpoint_forward_allowed")
        or released.get("open_champion_eligible")
        or released.get("released_checkpoint_count") != 5
    ):
        raise ScoobyBiologicalReadinessError("released-lane disposition differs")
    acquisition = frozen["checkpoint_acquisition"]
    decoder = frozen["decoder_probe"]
    native_context = frozen["native_context"]
    queue_preflight = frozen["synthetic_queue_preflight"].get("metadata", {})
    if (
        acquisition.get("checkpoint", {}).get("sha256")
        != CHECKPOINTS["scooby_epicardioids"]["sha256"]
        or not acquisition.get("checkpoint", {}).get("full_payload_hash_verified")
        or decoder.get("status")
        != "pass_checkpoint_decoder_executable_full_sequence_and_native_context_still_blocked"
        or native_context.get("status")
        != "pass_native_context_rows_released_query_encoder_still_absent"
        or queue_preflight.get("status") != "pass_ready_for_central_dispatcher"
    ):
        raise ScoobyBiologicalReadinessError("Epicardioids execution evidence differs")

    gse296875 = load_toml(root / CURRENT_AUTHORITIES["gse296875_dataset"][0])
    gse281367 = load_toml(root / CURRENT_AUTHORITIES["gse281367_dataset"][0])
    capabilities = load_toml(
        root / CURRENT_AUTHORITIES["observed_multiome_capabilities"][0]
    )
    multivi = frozen["gse296875_multivi_prepare"]
    if (
        gse296875.get("expected_biological_units") != 39
        or gse296875.get("pairing_levels") != ["same_nucleus"]
        or multivi.get("pairing_topology") != "same_nucleus"
        or multivi.get("outer_unit") != "donor"
        or not all(not row.get("query_atac_exported") for row in multivi["folds"].values())
    ):
        raise ScoobyBiologicalReadinessError("GSE296875 donor-safe topology differs")
    if (
        gse281367.get("modalities") != ["single_nucleus_atac"]
        or gse281367.get("expected_biological_units") != 12
        or frozen["gse281367_activation"].get("rna_conditioning_active")
        or frozen["gse281367_production"].get("family_native_tournament_ready")
        or capabilities.get("execution_authorized")
    ):
        raise ScoobyBiologicalReadinessError("observed-multiome task state differs")

    queue = load_json(root / CURRENT_AUTHORITIES["synthetic_probe_queue_item"][0])
    probe_config = load_json(root / CURRENT_AUTHORITIES["synthetic_probe_config"][0])
    wrapper_text = (
        root / CURRENT_AUTHORITIES["synthetic_probe_wrapper"][0]
    ).read_text(encoding="utf-8")
    header = parse_sbatch_header(wrapper_text)
    if (
        not queue.get("enabled")
        or queue.get("bundle_id") != "model-probe-284"
        or queue.get("wrapper_sha256")
        != CURRENT_AUTHORITIES["synthetic_probe_wrapper"][1]
        or header.get("job-name") != "model-probe-284"
        or header.get("partition") != "gpu"
        or header.get("qos") != "nslab"
        or header.get("gres") != "gpu:l40s:1"
        or probe_config.get("task_scope", {}).get("biological_task_executed")
    ):
        raise ScoobyBiologicalReadinessError("central-dispatch queue item differs")

    decisions = biological_decisions()
    result = {
        "schema_version": "masld-bench-scooby-biological-fixture-readiness-v1",
        "audit_status": "passed_fail_closed",
        "exact_code": {
            "repository": "gagneurlab/scooby",
            "revision": CODE_REVISION,
            "license": "MIT",
            "license_sha256": CODE_LICENSE_SHA256,
            "checkpoint_bound_revision": "UNRESOLVED",
        },
        "released_checkpoints": {
            model_id: {
                **expected,
                "pinned_terms": pinned_terms[model_id],
                "current_repository_terms_snapshot": current_terms[model_id],
                "coordinate_compatible_liver_query_encoder_released": False,
                "biological_smoke_authorized": False,
            }
            for model_id, expected in CHECKPOINTS.items()
        },
        "terms": {
            "scooby_code": "MIT",
            "epicardioids_weight_card": "MIT",
            "onek1k_weight_card": "MIT",
            "neurips_weight_card": "UNDECLARED",
            "borzoi_converter_current_snapshot": base_converter,
            "canonical_calico_borzoi_weight_terms": "UNDECLARED",
            "merged_derivative_redistribution": "UNRESOLVED",
            "open_champion_eligible": False,
            "finding": "downstream_or_converter_card_labels_do_not_resolve_canonical_inherited_weight_authority",
        },
        "native_contract": {
            "sequence_input_bp": 524_288,
            "sequence_encoding": "one_hot_hg38",
            "checkpoint_bound_base_to_channel_mapping": "UNRESOLVED",
            "output_central_bins": 6_144,
            "output_resolution_bp": 32,
            "output_span_bp": 196_608,
            "native_annotation": "GENCODE_v32",
            "native_reference_bundles": [
                "refdata-cellranger-arc-GRCh38-2020-A-2.0.0",
                "refdata-gex-GRCh38-2020-A",
            ],
            "benchmark_reference": "GRCh38.p14_plus_GENCODE_v49",
            "exact_native_FASTA_identity": "UNRESOLVED",
            "post_release_sequence_window_parity": "UNRESOLVED",
            "reference_and_v32_to_v49_biological_execution_gate": "blocked",
            "rna_inverse_transform": "Borzoi_like_squashed_scale_clip_soft_5",
            "atac_inverse_transform": "multiply_by_20_where_ATAC_head_exists",
        },
        "dataset_fit": {
            "gse296875": {
                "donors": 39,
                "topology": "same_nucleus_RNA_ATAC",
                "role": "project_exposed_development",
                "donor_safe_fold_authority_available": True,
                "supports_new_fold_fitted_context_model_in_principle": True,
                "supports_frozen_released_scooby_context": False,
                "reason": "a_new_scGLUE_or_MultiVI_latent_has_an_arbitrary_coordinate_system_relative_to_the_frozen_decoder",
            },
            "gse281367": {
                "donors": 12,
                "topology": "ATAC_only_same_study_unpaired",
                "rna_state": "structurally_missing",
                "role": "project_exposed_ATAC_only_development_transport",
                "supports_frozen_released_scooby_context": False,
                "production_transport_ready": False,
            },
        },
        "biological_fixture_decisions": decisions,
        "central_dispatch": {
            "existing_item": "config/campaigns/gpu_bundle_queue/062-model-probe-284.json",
            "enabled": True,
            "duplicate_item_created": False,
            "manual_gpu_submission": False,
            "job_name": header["job-name"],
            "partition": header["partition"],
            "qos": header["qos"],
            "gpu": header["gres"],
            "cpus": int(header["cpus-per-task"]),
            "memory": header["mem"],
            "time": header["time"],
            "scope": "synthetic_architecture_compatibility_only",
        },
        "next_gate": {
            "released_biological_smoke": "recover_checkpoint_coordinate_compatible_query_encoder_then_clear_reference_and_inherited_terms",
            "project_adaptation": "clear_derivative_terms_then_prespecify_training_fold_only_context_bridge_or_decoder_adaptation",
            "gse281367_transport": "freeze_label_agnostic_prediction_mask_and_scoring_exchange",
            "gpu": "let_the_existing_central_dispatcher_run_model_probe_284_no_duplicate_submission",
        },
        "frozen_authority_manifests": {
            name: {"path": relative, "artifacts_sha256": manifest}
            for name, (relative, manifest, _member) in FROZEN_ARTIFACTS.items()
        },
        "current_authorities": current,
        "checkpoint_deserialized": False,
        "model_forward_executed": False,
        "biological_matrices_read": False,
        "project_outcomes_read": False,
        "sealed_features_read": False,
        "sealed_outcomes_read": False,
        "gpu_job_submitted": False,
        "terminal_disposition": "synthetic_probe_already_queued_released_biological_smoke_blocked_new_gse296875_adaptation_conditional_gse281367_not_a_native_query",
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
    result = audit(parse_args())
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
