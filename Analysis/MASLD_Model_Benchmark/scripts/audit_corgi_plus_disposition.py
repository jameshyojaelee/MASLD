#!/usr/bin/env python3
"""Validate the frozen current-release Corgi+ terminal disposition."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "masld-bench-corgi-plus-release-disposition-v1"
BLOCKED_STATUS = "TERMINAL_BLOCKED_CURRENT_UPSTREAM_RELEASE"


class CorgiPlusDispositionError(ValueError):
    """Raised when the frozen Corgi+ disposition or an authority drifts."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _resolve_member(root: Path, relative: str) -> Path:
    member = (root / relative).resolve(strict=True)
    try:
        member.relative_to(root.resolve(strict=True))
    except ValueError as error:
        raise CorgiPlusDispositionError("authority path escapes benchmark root") from error
    if not member.is_file() or member.is_symlink():
        raise CorgiPlusDispositionError(f"authority is not a regular file: {relative}")
    return member


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise CorgiPlusDispositionError(f"JSON authority is not an object: {path}")
    return value


def audit_disposition(root: Path, authority_path: Path) -> dict[str, Any]:
    root = root.resolve(strict=True)
    authority_path = authority_path.resolve(strict=True)
    authority = _load_json(authority_path)
    if authority.get("schema_version") != SCHEMA_VERSION:
        raise CorgiPlusDispositionError("Corgi+ disposition schema differs")
    if authority.get("status") != BLOCKED_STATUS:
        raise CorgiPlusDispositionError("Corgi+ current-release status differs")

    authority_hashes = authority.get("authority_hashes")
    if not isinstance(authority_hashes, dict) or len(authority_hashes) != 7:
        raise CorgiPlusDispositionError("Corgi authority hash family differs")
    verified: dict[str, str] = {}
    loaded: dict[str, dict[str, Any]] = {}
    for relative, expected in sorted(authority_hashes.items()):
        if not isinstance(relative, str) or not isinstance(expected, str):
            raise CorgiPlusDispositionError("Corgi authority hash record differs")
        member = _resolve_member(root, relative)
        observed = _digest(member)
        if observed != expected:
            raise CorgiPlusDispositionError(f"Corgi authority drifted: {relative}")
        verified[relative] = observed
        loaded[relative] = _load_json(member)

    plus = loaded["config/artifacts/models/corgi_plus/checkpoints.json"]
    regular = loaded["config/artifacts/models/corgi/checkpoints.json"]
    regular_release = loaded["config/artifacts/models/corgi/release.json"]
    model_identity = authority.get("model_identity", {})
    comparator = authority.get("nearest_legitimate_comparator", {})
    released = authority.get("released_surface", {})
    tournament = authority.get("tournament_disposition", {})

    if (
        plus.get("checkpoint", {}).get("sha256")
        != model_identity.get("checkpoint_sha256")
        or plus.get("checkpoint", {}).get("size_bytes")
        != model_identity.get("checkpoint_size_bytes")
        or plus.get("architecture_contract", {}).get("checkpoint_state_key_count")
        != model_identity.get("checkpoint_state_key_count")
        or plus.get("architecture_contract", {}).get("checkpoint_state_numel")
        != model_identity.get("checkpoint_state_numel")
    ):
        raise CorgiPlusDispositionError("Corgi+ checkpoint identity differs")
    if (
        plus.get("loader_audit", {}).get("status")
        != "NO_COMPLETE_RELEASED_INFERENCE_PATH"
        or plus.get("verified_artifacts", {}).get("released_input_bundle_complete")
        is not False
        or released.get("inference", {}).get("portable_complete_loader_available")
        is not False
        or released.get("checkpoint", {}).get("biological_input_or_numeric_parity_established")
        is not False
    ):
        raise CorgiPlusDispositionError("Corgi+ released execution surface differs")
    if (
        released.get("training", {}).get("source_defined_auxiliary_width") != 26
        or len(released.get("training", {}).get("source_defined_candidate_auxiliary_order", []))
        != 5
        or released.get("training", {}).get("candidate_order_checkpoint_bound") is not False
        or released.get("source_data_zip_bounded_inventory", {}).get("csv_contents_opened")
        is not False
    ):
        raise CorgiPlusDispositionError("Corgi+ input or source-data audit differs")
    if any(
        tournament.get(field) is not False
        for field in (
            "rna_conditioned_atac_execution_authorized",
            "variant_to_regulation_task_registered",
            "model_training_authorized",
            "model_adaptation_authorized",
            "production_sbatch_authorized",
            "comparative_claim_eligible",
            "conditional_component_eligible",
        )
    ):
        raise CorgiPlusDispositionError("Corgi+ tournament gate is not fail-closed")
    if len(authority.get("blocking_gates", [])) != 7 or any(
        gate.get("state") != "blocked" for gate in authority.get("blocking_gates", [])
    ):
        raise CorgiPlusDispositionError("Corgi+ blocking-gate family differs")

    regular_checkpoint = regular.get("checkpoint", {})
    release_weight = regular_release.get("weight", {})
    if (
        comparator.get("model_id") != "corgi_regular"
        or comparator.get("checkpoint_sha256") != regular_checkpoint.get("sha256")
        or comparator.get("checkpoint_sha256") != release_weight.get("sha256")
        or comparator.get("checkpoint_size_bytes") != regular_checkpoint.get("size_bytes")
        or comparator.get("code_commit") != regular_release.get("code_commit")
        or comparator.get("not_a_corgi_plus_substitute") is not True
        or comparator.get("relabeling_as_corgi_plus_prohibited") is not True
    ):
        raise CorgiPlusDispositionError("nearest regular-Corgi comparator differs")

    return {
        "schema_version": "masld-bench-corgi-plus-disposition-audit-receipt-v1",
        "status": "pass_frozen_terminal_blocker",
        "authority_sha256": _digest(authority_path),
        "verified_authority_hashes": verified,
        "corgi_plus_checkpoint_sha256": model_identity["checkpoint_sha256"],
        "corgi_plus_executable_for_registered_task": False,
        "corgi_plus_production_sbatch_authorized": False,
        "nearest_legitimate_comparator": "corgi_regular",
        "nearest_comparator_checkpoint_sha256": comparator["checkpoint_sha256"],
        "outcomes_opened": False,
        "sealed_features_opened": False,
    }


def main() -> None:
    default_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=default_root)
    parser.add_argument(
        "--authority",
        type=Path,
        default=default_root
        / "config/artifacts/models/corgi_plus/release_disposition_20260824.json",
    )
    parser.add_argument(
        "--require-executable",
        action="store_true",
        help="Fail after auditing because the frozen current release is not executable.",
    )
    arguments = parser.parse_args()
    receipt = audit_disposition(arguments.root, arguments.authority)
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if arguments.require_executable:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
