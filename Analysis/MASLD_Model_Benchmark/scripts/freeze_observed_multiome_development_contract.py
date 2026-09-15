#!/usr/bin/env python3
"""Validate and freeze the standalone observed-multiome development requirements."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import tomllib
from typing import Any, Mapping

from masld_bench.artifacts import freeze_tree, write_json_exclusive, write_text_exclusive
from masld_bench.registry import Registry, load_task_spec


class ObservedMultiomeContractError(RuntimeError):
    """Raised when the standalone task opens a forbidden global or claim check."""


TASK_PATH = Path("config/evaluation/observed_multiome_task.toml")
CAPABILITY_PATH = Path(
    "config/evaluation/observed_multiome_task_capabilities.toml"
)
GATE_PATH = Path("config/evaluation/observed_multiome_development_gate.toml")
TOURNAMENT_PATH = Path("config/evaluation/family_native_tournament.toml")
GLOBAL_GATE_PATH = Path("config/evaluation/promotion_gates.toml")
EXPECTED_DATASETS = {"gse244832", "gse281367", "gse296875"}
EXPECTED_TOPOLOGY = {
    "gse296875": {"same_nucleus"},
    "gse244832": {"same_sample_different_aliquot"},
    "gse281367": {"same_study_unpaired"},
}
REQUIRED_FUTURE_BASELINES = {
    "masked_modality",
    "observed_atac_glm",
    "observed_atac_only",
    "rna_only",
    "shuffled_modality",
}


def _load_toml(path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise ObservedMultiomeContractError(f"authority is not a regular file: {path}")
    with path.open("rb") as handle:
        value = tomllib.load(handle)
    if not isinstance(value, dict):
        raise ObservedMultiomeContractError(f"authority is not a table: {path}")
    return value


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def validate_gate(gate: Mapping[str, Any]) -> None:
    expected_false = {
        "external_seal_bound",
        "model_promotion_allowed",
        "champion_claim_allowed",
        "rna_conditioned_atac_claim_allowed",
        "universal_claim_allowed",
        "clinical_claim_allowed",
        "selection_lock_allowed",
    }
    if (
        gate.get("schema_version") != "masld-bench-development-gate-v1"
        or gate.get("gate_id") != "observed_multiome_development_v1"
        or gate.get("task_id") != "observed_multiome"
        or gate.get("promotion_mode") != "development_only_no_promotion"
        or gate.get("global_census_revision_required_before_execution") is not True
        or any(gate.get(field) is not False for field in expected_false)
    ):
        raise ObservedMultiomeContractError("development or claim gate opened")
    if (
        gate.get("bootstrap_replicates") != 10000
        or gate.get("resampling_units") != ["donor", "genomic_block"]
        or gate.get("minimum_improved_lineages") != 4
        or len(gate.get("major_lineages", [])) != 5
    ):
        raise ObservedMultiomeContractError("development uncertainty gate drifted")


def validate_capabilities(
    capability: Mapping[str, Any], *, task_baselines: tuple[str, ...], model_ids: set[str]
) -> None:
    if (
        capability.get("schema_version")
        != "masld-bench-standalone-observed-multiome-capabilities-v1"
        or capability.get("registry_id")
        != "observed_multiome_development_capabilities"
        or capability.get("task_id") != "observed_multiome"
    ):
        raise ObservedMultiomeContractError("capability identity drifted")
    required_true = {
        "false_cell_pairing_forbidden",
        "scored_mask_absent_from_query_atac",
        "matched_baseline_observed_evidence_required",
        "native_output_families_ranked_separately",
        "admission_blocking",
    }
    required_false = {
        "missing_modality_encoded_as_zero",
        "rna_conditioned_atac_claim_allowed",
        "external_champion_claim_allowed",
        "global_registry_integration_allowed",
        "execution_authorized",
    }
    if any(capability.get(field) is not True for field in required_true) or any(
        capability.get(field) is not False for field in required_false
    ):
        raise ObservedMultiomeContractError("capability safety gate opened")
    if set(capability.get("development_dataset_ids", [])) != EXPECTED_DATASETS:
        raise ObservedMultiomeContractError("development dataset roster drifted")
    topology = {
        "same_nucleus_dataset_ids": {"gse296875"},
        "same_donor_different_aliquot_dataset_ids": {"gse244832"},
        "atac_only_dataset_ids": {"gse281367"},
    }
    if any(set(capability.get(field, [])) != expected for field, expected in topology.items()):
        raise ObservedMultiomeContractError("dataset topology roster drifted")
    candidate_fields = (
        "observed_atac_conditioned_track_candidates",
        "sequence_plus_rna_track_candidates",
        "joint_representation_candidates",
        "registered_baseline_models",
    )
    named_models = {
        model_id for field in candidate_fields for model_id in capability.get(field, [])
    }
    unknown = named_models.difference(model_ids)
    if unknown:
        raise ObservedMultiomeContractError(
            "standalone contract names unknown global models: " + ", ".join(sorted(unknown))
        )
    if tuple(capability.get("registered_baseline_models", [])) != task_baselines:
        raise ObservedMultiomeContractError("standalone baseline roster drifted")
    future = set(capability.get("required_future_baseline_ids", []))
    if future != REQUIRED_FUTURE_BASELINES or future.intersection(model_ids):
        raise ObservedMultiomeContractError("future baseline gate drifted")
    if capability.get("profile_fixture_passed_models") or capability.get(
        "primary_eligible_models"
    ):
        raise ObservedMultiomeContractError("a model became eligible without a census revision")


def validate_contract(root: Path) -> dict[str, Any]:
    root = root.resolve(strict=True)
    registry = Registry.load(root / "config")
    task_path = root / TASK_PATH
    capability_path = root / CAPABILITY_PATH
    gate_path = root / GATE_PATH
    task = load_task_spec(task_path)
    capability = _load_toml(capability_path)
    gate = _load_toml(gate_path)

    if task.task_id != "observed_multiome" or task.status.value != "candidate":
        raise ObservedMultiomeContractError("standalone TaskSpec identity drifted")
    if set(task.dataset_ids) != EXPECTED_DATASETS or task.datasets_sealed:
        raise ObservedMultiomeContractError("standalone TaskSpec dataset role drifted")
    if task.promotion_gate_id != gate.get("gate_id") or task.promotion_gate_config_sha256 != _digest(gate_path):
        raise ObservedMultiomeContractError("TaskSpec development-gate binding drifted")
    if "observed_multiome" in registry.tasks:
        raise ObservedMultiomeContractError("standalone task entered the frozen global registry")
    globally_supported = [
        model.model_id
        for model in registry.models.values()
        if "observed_multiome" in model.supported_tasks
    ]
    if globally_supported:
        raise ObservedMultiomeContractError(
            "historical model manifests claim standalone task support: "
            + ", ".join(sorted(globally_supported))
        )
    for dataset_id, expected_pairing in EXPECTED_TOPOLOGY.items():
        observed = {
            state.value for state in registry.get_dataset(dataset_id).pairing_levels
        }
        if observed != expected_pairing:
            raise ObservedMultiomeContractError(
                f"dataset topology differs for {dataset_id}: {sorted(observed)}"
            )
    validate_gate(gate)
    validate_capabilities(
        capability,
        task_baselines=task.baseline_model_ids,
        model_ids=set(registry.models),
    )

    tournament = _load_toml(root / TOURNAMENT_PATH)
    observed_task = tournament.get("task", {}).get("observed_multiome", {})
    if (
        "observed_multiome" in tournament.get("active_task_ids", [])
        or "observed_multiome" not in tournament.get("task_ids_to_register", [])
        or observed_task.get("status") != "planned_registration"
    ):
        raise ObservedMultiomeContractError("family-native global registration gate opened")
    if "observed_multiome" in _load_toml(root / GLOBAL_GATE_PATH).get("tasks", {}):
        raise ObservedMultiomeContractError("global promotion gate was changed")

    authorities = {
        path.as_posix(): _digest(root / path)
        for path in (
            TASK_PATH,
            CAPABILITY_PATH,
            GATE_PATH,
            TOURNAMENT_PATH,
            GLOBAL_GATE_PATH,
            Path("scripts/freeze_observed_multiome_development_contract.py"),
        )
    }
    return {
        "schema_version": "masld-bench-observed-multiome-development-freeze-v1",
        "status": "passed_standalone_development_contract_execution_blocked",
        "task_id": task.task_id,
        "global_registry_snapshot_sha256": registry.snapshot_sha256,
        "dataset_ids": sorted(EXPECTED_DATASETS),
        "global_model_count": len(registry.models),
        "primary_eligible_models": [],
        "execution_authorized": False,
        "external_champion_claim_allowed": False,
        "authorities": authorities,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = validate_contract(args.root)
    args.output.mkdir(parents=False, exist_ok=False)
    write_json_exclusive(args.output / "validation_receipt.json", receipt)
    source_lines = [
        f"{digest}  {(args.root.resolve() / relative).resolve()}"
        for relative, digest in sorted(receipt["authorities"].items())
    ]
    write_text_exclusive(args.output / "source.sha256", "\n".join(source_lines) + "\n")
    manifest = freeze_tree(
        args.output,
        metadata={
            "artifact_class": "observed_multiome_standalone_development_contract",
            "task_id": "observed_multiome",
            "execution_authorized": False,
            "external_champion_claim_allowed": False,
        },
    )
    print(json.dumps({"output": str(args.output), "artifacts_sha256": manifest}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
