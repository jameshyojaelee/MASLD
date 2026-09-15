#!/usr/bin/env python3
"""Audit the model/evaluator materialization separation without reading counts."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any, Mapping

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-observed-multiome-materialization-contract-v1"
CONTRACT_ID = "gse296875_observed_multiome_materialization_20260825"
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


class ObservedMultiomeMaterializationContractError(ValueError):
    """Raised when the model/evaluator separation or mask geometry drifts."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ObservedMultiomeMaterializationContractError("JSON object required")
    return value


def _read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ObservedMultiomeMaterializationContractError(f"TSV lacks header: {path}")
        return tuple(reader.fieldnames), [dict(row) for row in reader]


def _bound_tree(root: Path, record: Mapping[str, Any], label: str) -> Path:
    path = (root / str(record.get("path", ""))).resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise ObservedMultiomeMaterializationContractError(f"{label} artifact drifted")
    return path


def validate(root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("contract_id") != CONTRACT_ID
        or config.get("dataset_id") != "gse296875"
        or config.get("dataset_view_id") != "gse296875_rna_atac_smoke_1000_v1"
        or config.get("stage") != "smoke"
    ):
        raise ObservedMultiomeMaterializationContractError("contract identity differs")
    parents = config.get("parents")
    if not isinstance(parents, dict) or set(parents) != {
        "mask_plan",
        "target_axis_contract",
        "factorized_fixture",
    }:
        raise ObservedMultiomeMaterializationContractError("parent roster differs")
    mask_root = _bound_tree(root, parents["mask_plan"], "mask plan")
    _bound_tree(root, parents["target_axis_contract"], "target-axis contract")
    _bound_tree(root, parents["factorized_fixture"], "factorized fixture")

    unit = config.get("biological_unit")
    if unit != {
        "axis": "donor_by_lineage",
        "units": 195,
        "donors": 39,
        "lineages": 5,
        "cells_or_nuclei_as_replicates": False,
        "aggregation": "sum_raw_integer_counts_within_donor_by_lineage",
    }:
        raise ObservedMultiomeMaterializationContractError("biological unit differs")
    views = config.get("views")
    if views != {
        "genomic_folds": 5,
        "targets_per_fold": 1000,
        "observed_atac_inputs_per_fold": 2000,
        "model_input_artifacts": 5,
        "evaluator_outcome_artifacts": 5,
        "crossed_donor_by_genomic_fit_surfaces": 25,
    }:
        raise ObservedMultiomeMaterializationContractError("view census differs")

    model = config.get("model_input_role")
    if (
        not isinstance(model, dict)
        or model.get("one_artifact_per_genomic_fold") is not True
        or model.get("target_atac_count_values") != "forbidden"
        or model.get("target_metadata_without_values") != "required"
        or model.get("normalization_or_feature_selection")
        != "forbidden_during_materialization"
        or model.get("missingness_masks") != "required"
        or model.get("raw_donor_ids") != "forbidden"
    ):
        raise ObservedMultiomeMaterializationContractError("model-input role differs")
    evaluator = config.get("evaluator_role")
    if (
        not isinstance(evaluator, dict)
        or evaluator.get("one_artifact_per_genomic_fold") is not True
        or evaluator.get("allowed_values") != ["raw_target_atac_counts"]
        or evaluator.get("rna_count_values") != "forbidden"
        or evaluator.get("observed_atac_input_count_values") != "forbidden"
        or evaluator.get("model_fit_or_predict_binding") != "forbidden"
        or evaluator.get("label_join_or_metric") != "forbidden_during_materialization"
        or evaluator.get("raw_donor_ids") != "forbidden"
    ):
        raise ObservedMultiomeMaterializationContractError("evaluator role differs")
    sparse = config.get("selective_sparse_read")
    if (
        not isinstance(sparse, dict)
        or sparse.get("csr_indptr_and_indices_may_be_read_for_column_membership") is not True
        or sparse.get("model_input_data_values")
        != "read_only_positions_whose_columns_are_selected_input_peaks"
        or sparse.get("evaluator_data_values")
        != "read_only_positions_whose_columns_are_selected_target_peaks"
        or sparse.get("same_fold_model_input_and_target_value_positions_disjoint") is not True
        or sparse.get("loading_full_atac_csr_data_array_in_model_materializer") is not False
    ):
        raise ObservedMultiomeMaterializationContractError("selective sparse-read rule differs")

    binding = config.get("run_binding")
    if (
        not isinstance(binding, dict)
        or set(binding.get("fit_and_predict_may_bind", []))
        != {"one_matching_model_input_artifact", "immutable_environment_and_reference_artifacts"}
        or set(binding.get("fit_and_predict_may_not_bind", []))
        != {
            "any_evaluator_outcome_artifact",
            "another_genomic_fold_model_input_artifact",
            "raw_multimodal_h5",
        }
        or set(binding.get("evaluator_may_bind", []))
        != {"one_matching_evaluator_outcome_artifact", "frozen_prediction_bundle"}
        or binding.get("prediction_hash_committed_before_evaluator_join") is not True
    ):
        raise ObservedMultiomeMaterializationContractError("run binding differs")
    firewall = config.get("firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise ObservedMultiomeMaterializationContractError("contract-audit firewall is open")

    unit_fields, unit_rows = _read_tsv(mask_root / "donor_lineage_units.tsv")
    if unit_fields != ("donor_hash", "lineage", "donor_fold", "nuclei") or len(unit_rows) != 195:
        raise ObservedMultiomeMaterializationContractError("donor-lineage mask differs")
    if any(SHA256_PATTERN.fullmatch(row["donor_hash"]) is None for row in unit_rows):
        raise ObservedMultiomeMaterializationContractError("donor hashes differ")
    donor_lineages: dict[str, set[str]] = defaultdict(set)
    donor_folds: dict[str, str] = {}
    for row in unit_rows:
        donor_lineages[row["donor_hash"]].add(row["lineage"])
        if donor_folds.setdefault(row["donor_hash"], row["donor_fold"]) != row["donor_fold"]:
            raise ObservedMultiomeMaterializationContractError("one donor crosses folds")
        if int(row["nuclei"]) < 1:
            raise ObservedMultiomeMaterializationContractError("empty donor-lineage unit")
    if (
        len(donor_lineages) != 39
        or {len(value) for value in donor_lineages.values()} != {5}
        or Counter(map(int, donor_folds.values())) != Counter({0: 11, 1: 9, 2: 7, 3: 4, 4: 8})
    ):
        raise ObservedMultiomeMaterializationContractError("donor-fold census differs")

    target_fields, target_rows = _read_tsv(mask_root / "target_peaks.tsv")
    input_fields, input_rows = _read_tsv(mask_root / "observed_atac_input_peaks.tsv")
    if target_fields != (
        "target_hash",
        "chromosome",
        "bed_start_0based",
        "bed_end_half_open",
        "genomic_fold",
    ) or input_fields != (
        "input_hash",
        "chromosome",
        "bed_start_0based",
        "bed_end_half_open",
        "source_genomic_fold",
        "held_genomic_fold",
    ):
        raise ObservedMultiomeMaterializationContractError("peak mask schema differs")
    if len(target_rows) != 5000 or len(input_rows) != 10000:
        raise ObservedMultiomeMaterializationContractError("peak mask cardinality differs")
    target_chromosomes: dict[int, set[str]] = defaultdict(set)
    input_chromosomes: dict[int, set[str]] = defaultdict(set)
    for row in target_rows:
        fold = int(row["genomic_fold"])
        target_chromosomes[fold].add(row["chromosome"])
        if SHA256_PATTERN.fullmatch(row["target_hash"]) is None:
            raise ObservedMultiomeMaterializationContractError("target hash differs")
    for row in input_rows:
        fold = int(row["held_genomic_fold"])
        input_chromosomes[fold].add(row["chromosome"])
        if SHA256_PATTERN.fullmatch(row["input_hash"]) is None:
            raise ObservedMultiomeMaterializationContractError("input hash differs")
    if Counter(map(int, (row["genomic_fold"] for row in target_rows))) != Counter({fold: 1000 for fold in range(5)}):
        raise ObservedMultiomeMaterializationContractError("target fold budget differs")
    if Counter(map(int, (row["held_genomic_fold"] for row in input_rows))) != Counter({fold: 2000 for fold in range(5)}):
        raise ObservedMultiomeMaterializationContractError("input fold budget differs")
    if any(not target_chromosomes[fold].isdisjoint(input_chromosomes[fold]) for fold in range(5)):
        raise ObservedMultiomeMaterializationContractError("held chromosome enters model ATAC input")

    return {
        "schema_version": "masld-bench-observed-multiome-materialization-contract-receipt-v1",
        "contract_id": CONTRACT_ID,
        "dataset_id": "gse296875",
        "stage": "smoke",
        "donors": 39,
        "donor_lineage_units": 195,
        "genomic_folds": 5,
        "model_input_artifacts_planned": 5,
        "evaluator_outcome_artifacts_planned": 5,
        "target_peaks": 5000,
        "observed_atac_input_peak_views": 10000,
        "same_fold_held_chromosomes_excluded": True,
        "selective_sparse_value_reads_required": True,
        "full_atac_data_load_in_model_materializer_allowed": False,
        "raw_donor_ids_present_in_mask_plan": False,
        "biological_count_values_read": False,
        "development_outcomes_read": False,
        "sealed_outcomes_read": False,
        "benchmark_metrics_calculated": False,
        "model_fit": False,
        "production_ranking_authorized": False,
        "champion_claim_allowed": False,
        "next_jobs": config["next_jobs"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    config_path.relative_to(root)
    receipt = validate(root, _load_json(config_path))
    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [
        config_path,
        Path(__file__).resolve(strict=True),
        root / "tests/unit/test_gse296875_observed_multiome_materialization_contract.py",
    ]
    (output / "source.sha256").write_text(
        "".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8"
    )
    digest = freeze_tree(
        output,
        metadata={
            "artifact_class": "gse296875_observed_multiome_materialization_contract",
            "contract_id": CONTRACT_ID,
            "biological_count_values_read": False,
            "sealed_outcomes_accessed": False,
            "model_fit": False,
        },
    )
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
