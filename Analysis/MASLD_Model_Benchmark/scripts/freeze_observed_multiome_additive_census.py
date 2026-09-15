#!/usr/bin/env python3
"""Freeze an additive observed-multiome census without mutating the base registry."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import tomllib
from typing import Any

from masld_bench.artifacts import freeze_tree, write_json_exclusive
from masld_bench.registry import Registry


SCHEMA = "masld-bench-observed-multiome-additive-census-v1"
REVISION_ID = "observed_multiome_additive_census_20260825"
EXPECTED_NEW_BASELINES = {
    "masked_modality",
    "observed_atac_glm",
    "observed_atac_only",
    "rna_only",
    "shuffled_modality",
}
EXPECTED_CANDIDATES = {
    "epiagent",
    "epibert",
    "epcotv2",
    "get",
    "lsi",
    "multivi",
    "peakvi",
    "scooby_epicardioids",
    "scooby_neurips",
    "seurat_wnn",
}
EXPECTED_EXISTING_BASELINES = {
    "assay_native_pseudobulk",
    "context_only",
    "lsi",
    "mean_track",
    "peakvi",
    "sequence_only",
    "seurat_wnn",
    "shuffled_context",
    "trans_only",
}
ALLOWED_OUTPUT_FAMILIES = {"functional_track", "joint_representation"}
ALLOWED_MISSINGNESS = {
    "not_applicable_for_removed_pathway",
    "not_applicable_for_removed_pathways",
    "preserve_original_explicit_missing_state",
    "structurally_missing_when_query_atac_is_unavailable",
}


class ObservedMultiomeCensusError(ValueError):
    """Raised when the additive census or its frozen base authorities drift."""


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ObservedMultiomeCensusError(f"JSON object required: {path}")
    return value


def _load_toml(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        value = tomllib.load(handle)
    if not isinstance(value, dict):
        raise ObservedMultiomeCensusError(f"TOML object required: {path}")
    return value


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _ids(records: Any, label: str) -> set[str]:
    if not isinstance(records, list) or not records:
        raise ObservedMultiomeCensusError(f"{label} must be a non-empty list")
    identifiers = [record.get("model_id") for record in records if isinstance(record, dict)]
    if len(identifiers) != len(records) or any(not isinstance(item, str) for item in identifiers):
        raise ObservedMultiomeCensusError(f"{label} contains an invalid model_id")
    if len(set(identifiers)) != len(identifiers):
        raise ObservedMultiomeCensusError(f"{label} contains duplicate model IDs")
    return set(identifiers)


def validate_revision(root: Path, config: dict[str, Any]) -> dict[str, Any]:
    if config.get("schema_version") != SCHEMA or config.get("revision_id") != REVISION_ID:
        raise ObservedMultiomeCensusError("additive census identity differs")

    registry = Registry.load(root / "config", validate_references=False)
    base = config.get("base_registry")
    if not isinstance(base, dict):
        raise ObservedMultiomeCensusError("base registry binding is absent")
    observed_base = {
        "model_count": len(registry.models),
        "task_count": len(registry.tasks),
        "snapshot_sha256": registry.snapshot_sha256,
    }
    if base != observed_base:
        raise ObservedMultiomeCensusError(
            f"base registry drifted: expected={base!r}; observed={observed_base!r}"
        )

    contract = config.get("task_contract")
    if not isinstance(contract, dict) or contract.get("task_id") != "observed_multiome":
        raise ObservedMultiomeCensusError("task binding differs")
    authorities: dict[str, dict[str, Any]] = {}
    for prefix in ("task", "gate", "capability"):
        relative = contract.get(f"{prefix}_path")
        expected = contract.get(f"{prefix}_sha256")
        if not isinstance(relative, str) or Path(relative).is_absolute():
            raise ObservedMultiomeCensusError(f"invalid {prefix} authority path")
        path = (root / relative).resolve(strict=True)
        path.relative_to(root.resolve(strict=True))
        if _digest(path) != expected:
            raise ObservedMultiomeCensusError(f"{prefix} authority drifted")
        authorities[prefix] = _load_toml(path)

    task = authorities["task"]
    gate = authorities["gate"]
    capability = authorities["capability"]
    if (
        task.get("task_id") != "observed_multiome"
        or gate.get("task_id") != "observed_multiome"
        or capability.get("task_id") != "observed_multiome"
    ):
        raise ObservedMultiomeCensusError("authority task IDs differ")
    if task.get("datasets_train") != ["gse296875"] or set(task.get("datasets_development", [])) != {
        "gse244832",
        "gse281367",
    }:
        raise ObservedMultiomeCensusError("development dataset roster differs")
    if task.get("datasets_sealed") != [] or gate.get("external_seal_bound") is not False:
        raise ObservedMultiomeCensusError("an external seal was introduced")
    if gate.get("promotion_mode") != "development_only_no_promotion":
        raise ObservedMultiomeCensusError("development-only gate differs")

    candidates = config.get("candidate_models")
    candidate_ids = _ids(candidates, "candidate_models")
    if candidate_ids != EXPECTED_CANDIDATES:
        raise ObservedMultiomeCensusError("candidate roster differs")
    missing_models = sorted(candidate_ids.difference(registry.models))
    if missing_models:
        raise ObservedMultiomeCensusError(f"unknown candidate models: {missing_models}")
    for record in candidates:
        if record.get("output_family") not in ALLOWED_OUTPUT_FAMILIES:
            raise ObservedMultiomeCensusError("candidate output family differs")
        modalities = record.get("query_modalities")
        if not isinstance(modalities, list) or not modalities:
            raise ObservedMultiomeCensusError("candidate query modalities are absent")
        if "single_nucleus_atac" in modalities:
            raise ObservedMultiomeCensusError("observed ATAC must be explicit in the census")

    existing_baselines = config.get("existing_mandatory_baselines")
    if not isinstance(existing_baselines, list) or set(existing_baselines) != EXPECTED_EXISTING_BASELINES:
        raise ObservedMultiomeCensusError("existing baseline roster differs")
    if not EXPECTED_EXISTING_BASELINES.issubset(registry.models):
        raise ObservedMultiomeCensusError("an existing baseline is absent from the base registry")

    new_baselines = config.get("new_mandatory_baselines")
    new_baseline_ids = _ids(new_baselines, "new_mandatory_baselines")
    if new_baseline_ids != EXPECTED_NEW_BASELINES:
        raise ObservedMultiomeCensusError("new modality-ablation baseline roster differs")
    if new_baseline_ids.intersection(registry.models):
        raise ObservedMultiomeCensusError("an additive baseline collides with the base registry")
    for record in new_baselines:
        if record.get("fit_scope") != "outer_training_only":
            raise ObservedMultiomeCensusError("baseline fit scope differs")
        if record.get("missingness_output") not in ALLOWED_MISSINGNESS:
            raise ObservedMultiomeCensusError("baseline missingness state differs")
        intervention = record.get("intervention")
        if not isinstance(intervention, str) or not intervention.strip():
            raise ObservedMultiomeCensusError("baseline intervention is absent")

    execution = config.get("execution_contract")
    expected_true = {
        "scored_mask_removed_before_query_transform",
        "receptive_field_buffer_removed_before_query_transform",
        "native_output_families_ranked_separately",
    }
    if not isinstance(execution, dict) or any(execution.get(key) is not True for key in expected_true):
        raise ObservedMultiomeCensusError("masked-query execution contract differs")
    if (
        execution.get("missing_as_zero") is not False
        or execution.get("query_fitted_preprocessing") is not False
        or execution.get("external_champion_claim_allowed") is not False
        or execution.get("execution_authorized_by_census_alone") is not False
        or execution.get("same_nucleus_only_models") != ["multivi", "seurat_wnn"]
        or execution.get("different_aliquot_resolution") != "donor_by_lineage_pseudobulk_only"
        or execution.get("atac_only_rna_state") != "structurally_missing"
    ):
        raise ObservedMultiomeCensusError("topology or claim firewall differs")

    firewall = config.get("outcome_firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise ObservedMultiomeCensusError("outcome firewall is not closed")
    if set(capability.get("required_future_baseline_ids", [])) != new_baseline_ids:
        raise ObservedMultiomeCensusError("capability-to-additive-baseline binding differs")

    return {
        "schema_version": "masld-bench-observed-multiome-additive-census-receipt-v1",
        "revision_id": REVISION_ID,
        "base_registry": observed_base,
        "task_id": "observed_multiome",
        "candidate_model_count": len(candidate_ids),
        "existing_baseline_count": len(existing_baselines),
        "new_baseline_count": len(new_baseline_ids),
        "output_families": sorted(ALLOWED_OUTPUT_FAMILIES),
        "development_dataset_ids": ["gse244832", "gse281367", "gse296875"],
        "census_requirement_satisfied": True,
        "candidate_fixture_requirement_satisfied": False,
        "external_seal_bound": False,
        "execution_authorized": False,
        "outcome_firewall": firewall,
    }


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    root = args.root.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    config_path.relative_to(root)
    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    config = _load_json(config_path)
    receipt = validate_revision(root, config)
    write_json_exclusive(output / "receipt.json", receipt)
    _write_tsv(
        output / "candidate_models.tsv",
        ["model_id", "role", "output_family", "query_modalities"],
        [
            {
                **record,
                "query_modalities": ";".join(record["query_modalities"]),
            }
            for record in config["candidate_models"]
        ],
    )
    _write_tsv(
        output / "new_baselines.tsv",
        ["model_id", "role", "intervention", "missingness_output", "fit_scope"],
        config["new_mandatory_baselines"],
    )
    source_paths = [
        config_path,
        Path(__file__).resolve(strict=True),
        root / "tests/unit/test_observed_multiome_additive_census.py",
    ]
    (output / "source.sha256").write_text(
        "".join(f"{_digest(path)}  {path}\n" for path in source_paths),
        encoding="utf-8",
    )
    digest = freeze_tree(
        output,
        metadata={
            "artifact_class": "observed_multiome_additive_census_revision",
            "revision_id": REVISION_ID,
            "execution_authorized": False,
            "sealed_outcomes_accessed": False,
        },
    )
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
