#!/usr/bin/env python3
"""Bind every released scooby checkpoint to a fail-closed task disposition."""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from typing import Any

from masld_bench.artifacts import ArtifactError, verify_frozen_tree


class ScoobyReleasedLaneError(RuntimeError):
    """Raised when a frozen release or current authority differs."""


FROZEN_ARTIFACTS = {
    "cross_review": "27668bff0943a8baac6489f073e000c2ba181e729a6f3360ea1dc0a3cfd6aa7c",
    "archive_inventory": "8cb983369ff43a93fbd2bc96bd6683613ce7d3cdf45efd6a4efe1c21715d2aa7",
    "bounded_headers": "569736f2d15e3c00326c25313a4cc33c352d9aeb6da296d86aad86e582e3ed9c",
}

REGISTERED = {
    "scooby_onek1k": {
        "checkpoint_sha256": "2ae335c651182040462d965bb879118be40a6cb174550a3686fb10b04a74f7cb",
        "outputs": ["rna_plus", "rna_minus"],
        "observed_atac_context": False,
    },
    "scooby_epicardioids": {
        "checkpoint_sha256": "f945e1937a46ddd5cb31249ef71e721e7a7b1381aa207fb2de449f23ac75a9e8",
        "outputs": ["rna_plus", "rna_minus", "atac_insertion"],
        "observed_atac_context": True,
    },
    "scooby_neurips": {
        "checkpoint_sha256": "8f2e2de3e86016378116c90accf94572e75963d5445beb6cabd585c8372170dd",
        "outputs": ["rna_plus", "rna_minus", "atac_insertion"],
        "observed_atac_context": True,
    },
}

VARIANTS = {
    "scooby_neurips_rna": {
        "revision": "e2a0f23c9231c84fe3d887a2ec81692d505c6882",
        "checkpoint_sha256": "7bfdb8beb63be285e563b962a46cf398669935e43401ac158010f5ada537273e",
        "size": 765094080,
        "parameters": 191267821,
        "tracks": ["rna_plus", "rna_minus"],
    },
    "scooby_neurips_flash": {
        "revision": "06bbc89694730faf040c2f070d244fef144e5174",
        "checkpoint_sha256": "18cca385e950dbb9baa4dd6535e0ce632ee64aa2618ba2faf5f3fe7a4bbc5c11",
        "size": 754774916,
        "parameters": 188688926,
        "tracks": ["rna_plus", "rna_minus", "atac_insertion"],
    },
}


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load(path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file() or path.stat().st_size > 1_000_000:
        raise ScoobyReleasedLaneError(f"authority path differs: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScoobyReleasedLaneError(f"authority JSON differs: {path}") from error
    if not isinstance(value, dict):
        raise ScoobyReleasedLaneError(f"authority JSON is not an object: {path}")
    return value


def require_artifact(root: Path, expected: str) -> dict[str, Any]:
    try:
        verified = verify_frozen_tree(root)
    except ArtifactError as error:
        raise ScoobyReleasedLaneError(f"frozen artifact differs: {root}") from error
    if digest(root / "ARTIFACTS.json") != expected:
        raise ScoobyReleasedLaneError(f"frozen artifact manifest differs: {root}")
    return dict(verified)


def validate_variants(authority: dict[str, Any]) -> dict[str, dict[str, Any]]:
    if (
        authority.get("schema_version")
        != "masld-bench-scooby-released-variants-v1"
        or authority.get("checkpoint_payloads_downloaded")
        or authority.get("full_payload_sha256_verified_locally")
        or set(authority.get("variants", {})) != set(VARIANTS)
    ):
        raise ScoobyReleasedLaneError("released variant authority differs")
    result: dict[str, dict[str, Any]] = {}
    for model_id, expected in VARIANTS.items():
        observed = authority["variants"][model_id]
        if (
            observed.get("revision") != expected["revision"]
            or observed.get("declared_weight_sha256")
            != expected["checkpoint_sha256"]
            or observed.get("file_size_bytes") != expected["size"]
            or observed.get("header_parameter_count") != expected["parameters"]
            or observed.get("decoder_tracks") != expected["tracks"]
            or observed.get("rna_conditioned_atac_eligible")
            or observed.get("sealed_masld_champion_eligible")
            or not str(observed.get("terminal_disposition", "")).startswith(
                "blocked_"
            )
        ):
            raise ScoobyReleasedLaneError(f"released variant differs: {model_id}")
        result[model_id] = {
            "checkpoint_sha256": expected["checkpoint_sha256"],
            "revision": expected["revision"],
            "decoder_tracks": expected["tracks"],
            "registered": False,
            "checkpoint_payload_downloaded": False,
            "checkpoint_deserialized": False,
            "checkpoint_forward_allowed": False,
            "rna_conditioned_atac_eligible": False,
            "sealed_masld_champion_eligible": False,
            "terminal_disposition": observed["terminal_disposition"],
        }
    return result


def load_cross_review_module(project_root: Path) -> Any:
    path = project_root / "scripts" / "scooby_cross_review.py"
    spec = importlib.util.spec_from_file_location("scooby_current_cross_review", path)
    if spec is None or spec.loader is None:
        raise ScoobyReleasedLaneError("current cross-review module is unavailable")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def audit(args: argparse.Namespace) -> dict[str, Any]:
    root = args.project_root.resolve(strict=True)
    if args.output.exists() or args.output.is_symlink():
        raise ScoobyReleasedLaneError("output already exists")
    frozen = {
        "cross_review": require_artifact(
            args.cross_review, FROZEN_ARTIFACTS["cross_review"]
        ),
        "archive_inventory": require_artifact(
            args.archive_inventory, FROZEN_ARTIFACTS["archive_inventory"]
        ),
        "bounded_headers": require_artifact(
            args.bounded_headers, FROZEN_ARTIFACTS["bounded_headers"]
        ),
    }
    cross = load(args.cross_review / "scooby_cross_review_receipt.json")
    archive = load(
        args.archive_inventory / "scooby_released_archive_inventory_receipt.json"
    )
    headers = load(
        args.bounded_headers / "scooby_bounded_header_admission_receipt.json"
    )
    if (
        cross.get("status") != "pass"
        or cross.get("checkpoint_forward_allowed")
        or cross.get("rna_conditioned_atac_models")
        or cross.get("sealed_inference_allowed")
        or archive.get("archive_payload_downloaded")
        or archive.get("checkpoint_bytes_downloaded")
        or headers.get("checkpoint_payload_downloaded")
        or headers.get("checkpoint_deserialized")
        or headers.get("model_forward_executed")
    ):
        raise ScoobyReleasedLaneError("frozen release disposition differs")

    args.output.mkdir(parents=True, exist_ok=False, mode=0o750)
    current_dir = args.output / "current_registry_review"
    current = load_cross_review_module(root).review(root, current_dir)
    if (
        current.get("checkpoint_forward_allowed")
        or current.get("rna_conditioned_atac_models")
        or current.get("sealed_inference_allowed")
        or set(current.get("models", {})) != set(REGISTERED)
    ):
        raise ScoobyReleasedLaneError("current registry disposition differs")
    registered: dict[str, dict[str, Any]] = {}
    for model_id, expected in REGISTERED.items():
        observed = current["models"][model_id]
        if (
            observed.get("decoder_tracks") != expected["outputs"]
            or observed.get("observed_atac_in_checkpoint_context")
            is not expected["observed_atac_context"]
            or observed.get("rna_conditioned_atac_eligible")
            or observed.get("checkpoint_execution_allowed")
        ):
            raise ScoobyReleasedLaneError(f"registered checkpoint differs: {model_id}")
        registered[model_id] = {
            "checkpoint_sha256": expected["checkpoint_sha256"],
            "decoder_tracks": expected["outputs"],
            "observed_atac_in_context": expected["observed_atac_context"],
            "registered": True,
            "checkpoint_payload_downloaded": False,
            "checkpoint_deserialized": False,
            "checkpoint_forward_allowed": False,
            "rna_conditioned_atac_eligible": False,
            "sealed_masld_champion_eligible": False,
            "terminal_disposition": observed["terminal_disposition"],
        }

    variant_path = (
        root
        / "config"
        / "artifacts"
        / "models"
        / "scooby_neurips"
        / "released_variants.json"
    )
    variants = validate_variants(load(variant_path))
    roster = {**registered, **variants}
    if len(roster) != 5:
        raise ScoobyReleasedLaneError("released checkpoint roster differs")

    result = {
        "schema_version": "masld-bench-scooby-released-lane-audit-v1",
        "status": "pass",
        "frozen_artifact_manifest_sha256": FROZEN_ARTIFACTS,
        "frozen_artifact_metadata": {
            name: value["metadata"] for name, value in frozen.items()
        },
        "current_authority_sha256": {
            **current["authority_files"],
            "released_variants": {
                "path": str(variant_path.relative_to(root)),
                "sha256": digest(variant_path),
            },
        },
        "released_checkpoint_count": len(roster),
        "registered_checkpoint_count": len(registered),
        "unregistered_released_variant_count": len(variants),
        "models": roster,
        "rna_conditioned_atac_models": [],
        "observed_multiome_development_only_models": [
            "scooby_epicardioids",
            "scooby_neurips",
            "scooby_neurips_flash",
        ],
        "checkpoint_forward_allowed": [],
        "sealed_inference_allowed": [],
        "open_champion_eligible": [],
        "project_data_read": False,
        "sealed_features_read": False,
        "sealed_outcomes_read": False,
        "checkpoint_payloads_downloaded": False,
        "checkpoint_deserialized": False,
        "model_forward_executed": False,
        "validation_job": {
            "eligible_now": False,
            "conditional_header": {
                "job_name": "model-training-201",
                "partition": "gpu",
                "qos": "nslab",
                "gres": "gpu:l40s:1",
                "cpus_per_task": 8,
                "memory": "128G",
                "time": "08:00:00",
            },
            "activation_gate": "Resolve model-specific terms, checkpoint-producing source/runtime, exact context encoder, native reference, safe strict restoration, and topology-valid fixed fixture first.",
        },
        "production_job": {
            "eligible_now": False,
            "header": None,
            "reason": "No released scooby checkpoint currently has a topology-valid liver query encoder, admitted runtime/reference contract, or eligible external MASLD evaluation.",
        },
        "terminal_disposition": "all_five_released_scooby_checkpoint_objects_have_fail_closed_terminal_dispositions; no_checkpoint_forward_or_production_job_is_eligible",
    }
    (args.output / "released_lane_audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--cross-review", type=Path, required=True)
    parser.add_argument("--archive-inventory", type=Path, required=True)
    parser.add_argument("--bounded-headers", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    result = audit(arguments)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
