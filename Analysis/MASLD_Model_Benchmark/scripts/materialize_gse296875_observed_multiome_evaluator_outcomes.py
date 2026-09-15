#!/usr/bin/env python3
"""Build five isolated evaluator-only GSE296875 target-count output files."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import (
    freeze_tree,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
)
from masld_bench.observed_multiome_materialization import (
    decode_strings,
    donor_lineage_aggregation,
    salted_hash,
    selected_value_positions,
    selectively_read_csr_columns,
    write_csr,
    write_string_dataset,
)


SCHEMA = "masld-bench-observed-multiome-evaluator-materialization-v1"
CAMPAIGN_ID = "gse296875_observed_multiome_evaluator_outcomes_20260825"


class EvaluatorMaterializationError(ValueError):
    """Raised when evaluator outcomes are not isolated from model inputs."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise EvaluatorMaterializationError("JSON object required")
    return value


def _tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise EvaluatorMaterializationError(f"{path.name} schema differs")
        return [dict(row) for row in reader]


def _tree(root: Path, record: Mapping[str, Any], label: str) -> Path:
    lexical = reject_symlink_components(root / str(record.get("path", "")), label=label)
    path = lexical.resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise EvaluatorMaterializationError(f"{label} drifted")
    return path


def validate_config(root: Path, config: Mapping[str, Any]) -> dict[str, Path]:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("campaign_id") != CAMPAIGN_ID
        or config.get("dataset_id") != "gse296875"
        or config.get("stage") != "smoke"
    ):
        raise EvaluatorMaterializationError("campaign identity differs")
    parents = config.get("parents")
    if not isinstance(parents, dict) or set(parents) != {"materialization_contract", "mask_plan", "verified_model_input_control"}:
        raise EvaluatorMaterializationError("parent roster differs")
    resolved = {key: _tree(root, record, key) for key, record in parents.items()}
    verification = _json(resolved["verified_model_input_control"] / "receipt.json")
    verified = parents["verified_model_input_control"]
    if (
        verified.get("biological_values_read") is not False
        or verification.get("promotion_gate_passed") is not True
        or verification.get("model_input_artifacts_sha256") != verified.get("model_input_artifacts_sha256")
    ):
        raise EvaluatorMaterializationError("verified model-input control differs")
    h5_record = config.get("input_h5")
    if not isinstance(h5_record, dict):
        raise EvaluatorMaterializationError("input HDF5 record is absent")
    lexical = reject_symlink_components(root / str(h5_record.get("path", "")), label="input HDF5")
    h5_path = lexical.resolve(strict=True)
    h5_path.relative_to(root)
    if not h5_path.is_file() or _digest(h5_path) != h5_record.get("sha256"):
        raise EvaluatorMaterializationError("input HDF5 drifted")
    resolved["input_h5"] = h5_path
    if config.get("expected") != {"nuclei": 1000, "donors": 39, "donor_lineage_units": 195, "atac_features": 306706, "genomic_folds": 5, "input_peaks_per_fold": 2000, "target_peaks_per_fold": 1000}:
        raise EvaluatorMaterializationError("expected census differs")
    if config.get("hash_namespaces") != {"donor": "gse296875:observed_multiome:donor:v1", "row": "gse296875:observed_multiome:row:v1", "peak": "gse296875:observed_multiome:peak:v1"}:
        raise EvaluatorMaterializationError("hash namespaces differ")
    artifact = config.get("evaluator_artifact")
    if (
        not isinstance(artifact, dict)
        or artifact.get("one_separately_frozen_child_per_genomic_fold") is not True
        or artifact.get("target_atac_values") != "raw_integer_donor_by_lineage_sum_selected_target_peaks_only"
        or artifact.get("rna_values") != "forbidden"
        or artifact.get("observed_atac_input_values") != "forbidden"
        or artifact.get("normalization") != "none"
        or artifact.get("feature_selection") != "none"
        or artifact.get("raw_donor_ids_exported") is not False
        or artifact.get("model_input_biological_values_read") is not False
    ):
        raise EvaluatorMaterializationError("evaluator artifact contract differs")
    selective = config.get("selective_read")
    if not isinstance(selective, dict) or selective != {"atac_indices_and_indptr_read": True, "atac_data_values_read_by_monotone_target_position_index": True, "full_atac_data_array_read": False, "observed_atac_input_value_positions_computed_without_value_read": True, "same_fold_selected_input_and_target_value_positions_disjoint": True}:
        raise EvaluatorMaterializationError("selective read contract differs")
    firewall = config.get("firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise EvaluatorMaterializationError("evaluator firewall is open")
    return resolved


def _resolve_columns(
    peak_ids: Sequence[str],
    chromosomes: Sequence[str],
    starts: Sequence[int],
    ends: Sequence[int],
    rows: Sequence[Mapping[str, str]],
    *,
    namespace: str,
    fold: int,
    kind: str,
) -> list[int]:
    field = "input_hash" if kind == "input" else "target_hash"
    by_hash: dict[str, tuple[int, str, int, int]] = {}
    for index, (peak_id, chromosome, start, end) in enumerate(zip(peak_ids, chromosomes, starts, ends, strict=True)):
        identifier = salted_hash(namespace, "input", fold, peak_id) if kind == "input" else salted_hash(namespace, "target", peak_id)
        by_hash[identifier] = (index, chromosome, int(start), int(end))
    columns: list[int] = []
    for row in rows:
        metadata = by_hash.get(row[field])
        if metadata is None or metadata[1:] != (row["chromosome"], int(row["bed_start_0based"]), int(row["bed_end_half_open"])):
            raise EvaluatorMaterializationError(f"{kind} identity differs")
        columns.append(metadata[0])
    if len(columns) != len(set(columns)):
        raise EvaluatorMaterializationError(f"{kind} columns are duplicated")
    return columns


def _write_fold_h5(path: Path, fold: int, planned: Sequence[Mapping[str, str]], row_hashes: Sequence[str], target_counts: Any, target_rows: Sequence[Mapping[str, str]]) -> dict[str, Any]:
    import h5py

    with h5py.File(path, "x") as handle:
        handle.attrs["schema_version"] = "masld-bench-observed-multiome-evaluator-outcome-h5-v1"
        handle.attrs["dataset_id"] = "gse296875"
        handle.attrs["stage"] = "smoke"
        handle.attrs["held_genomic_fold"] = fold
        handle.attrs["biological_unit"] = "donor_by_lineage"
        handle.attrs["evaluator_only"] = True
        rows = handle.create_group("rows")
        write_string_dataset(rows, "row_hash", row_hashes)
        write_string_dataset(rows, "unit_hash", [row["donor_hash"] for row in planned])
        write_string_dataset(rows, "lineage", [row["lineage"] for row in planned])
        rows.create_dataset("donor_fold", data=np.asarray([int(row["donor_fold"]) for row in planned], dtype=np.int8))
        rows.create_dataset("nuclei", data=np.asarray([int(row["nuclei"]) for row in planned], dtype=np.int16))
        target = handle.create_group("target_atac")
        identity = write_csr(target.create_group("counts_csr"), target_counts)
        write_string_dataset(target, "target_hash", [row["target_hash"] for row in target_rows])
        write_string_dataset(target, "chromosome", [row["chromosome"] for row in target_rows])
        target.create_dataset("bed_start_0based", data=np.asarray([int(row["bed_start_0based"]) for row in target_rows], dtype=np.int64))
        target.create_dataset("bed_end_half_open", data=np.asarray([int(row["bed_end_half_open"]) for row in target_rows], dtype=np.int64))
        target.attrs["scale"] = "raw_integer_donor_by_lineage_sum"
    return identity


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    config_path.relative_to(root)
    config = _json(config_path)
    resolved = validate_config(root, config)
    output = reject_symlink_components(args.output, label="evaluator output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    masks = resolved["mask_plan"]
    planned = _tsv(masks / "donor_lineage_units.tsv", ("donor_hash", "lineage", "donor_fold", "nuclei"))
    inputs = _tsv(masks / "observed_atac_input_peaks.tsv", ("input_hash", "chromosome", "bed_start_0based", "bed_end_half_open", "source_genomic_fold", "held_genomic_fold"))
    targets = _tsv(masks / "target_peaks.tsv", ("target_hash", "chromosome", "bed_start_0based", "bed_end_half_open", "genomic_fold"))
    expected = config["expected"]
    namespaces = config["hash_namespaces"]
    children: list[dict[str, Any]] = []
    import h5py

    with h5py.File(resolved["input_h5"], "r") as handle:
        donors = decode_strings(handle["obs/donor_id"][:], label="donor")
        lineages = decode_strings(handle["obs/broad_label"][:], label="lineage")
        cell_units = [(salted_hash(namespaces["donor"], donor), lineage) for donor, lineage in zip(donors, lineages, strict=True)]
        aggregation, _ = donor_lineage_aggregation(cell_unit_keys=cell_units, planned_units=planned)
        row_hashes = [salted_hash(namespaces["row"], row["donor_hash"], row["lineage"]) for row in planned]
        if len(donors) != expected["nuclei"] or len(set(donors)) != expected["donors"] or len(planned) != expected["donor_lineage_units"]:
            raise EvaluatorMaterializationError("donor or row census differs")
        peak_ids = decode_strings(handle["atac/peak_id"][:], label="peak ID")
        chromosomes = decode_strings(handle["atac/chromosome"][:], label="chromosome")
        starts = [int(value) for value in handle["atac/bed_start_0based"][:]]
        ends = [int(value) for value in handle["atac/bed_end_half_open"][:]]
        atac_group = handle["atac/counts_csr"]
        shape = tuple(int(value) for value in atac_group["shape"][:])
        indices = np.asarray(atac_group["indices"][:], dtype=np.int64)
        if shape != (expected["nuclei"], expected["atac_features"]):
            raise EvaluatorMaterializationError("ATAC shape differs")
        for fold in range(expected["genomic_folds"]):
            input_rows = [row for row in inputs if int(row["held_genomic_fold"]) == fold]
            target_rows = [row for row in targets if int(row["genomic_fold"]) == fold]
            input_columns = _resolve_columns(peak_ids, chromosomes, starts, ends, input_rows, namespace=namespaces["peak"], fold=fold, kind="input")
            target_columns = _resolve_columns(peak_ids, chromosomes, starts, ends, target_rows, namespace=namespaces["peak"], fold=fold, kind="target")
            if len(input_columns) != expected["input_peaks_per_fold"] or len(target_columns) != expected["target_peaks_per_fold"] or set(input_columns).intersection(target_columns):
                raise EvaluatorMaterializationError("fold target budget or separation differs")
            input_positions = selected_value_positions(indices, input_columns, total_columns=expected["atac_features"])
            target_positions_expected = selected_value_positions(indices, target_columns, total_columns=expected["atac_features"])
            if np.intersect1d(input_positions, target_positions_expected).size:
                raise EvaluatorMaterializationError("observed input value position overlaps target")
            target_cell, target_positions_read = selectively_read_csr_columns(atac_group, target_columns)
            if not np.array_equal(target_positions_read, target_positions_expected):
                raise EvaluatorMaterializationError("target selective read differs")
            target_counts = (aggregation @ target_cell).tocsr()
            child = output / f"genomic_fold_{fold}"
            child.mkdir()
            identity = _write_fold_h5(child / "evaluator_outcomes.h5", fold, planned, row_hashes, target_counts, target_rows)
            receipt = {"schema_version": "masld-bench-observed-multiome-evaluator-outcome-receipt-v1", "dataset_id": "gse296875", "stage": "smoke", "held_genomic_fold": fold, "donor_lineage_units": len(planned), "target_peaks": len(target_columns), "target_identity": identity, "target_sparse_value_positions_read": int(target_positions_read.size), "observed_input_sparse_value_positions": int(input_positions.size), "same_fold_input_target_value_position_overlap": 0, "rna_count_values_read": False, "observed_atac_input_count_values_read": False, "full_atac_data_array_read": False, "model_input_biological_artifact_read": False, "metric_calculated": False, "prediction_bundle_read": False, "model_fit": False, "sealed_outcomes_read": False}
            write_json_exclusive(child / "receipt.json", receipt)
            child_hash = freeze_tree(child, metadata={"artifact_class": "gse296875_observed_multiome_evaluator_outcome", "held_genomic_fold": fold, "evaluator_only": True, "sealed_outcomes_accessed": False})
            children.append({"held_genomic_fold": fold, "path": child.relative_to(output).as_posix(), "artifacts_sha256": child_hash})
    parent = {"schema_version": "masld-bench-observed-multiome-evaluator-campaign-receipt-v1", "campaign_id": CAMPAIGN_ID, "dataset_id": "gse296875", "stage": "smoke", "child_artifacts": children, "child_artifact_count": len(children), "donors": expected["donors"], "donor_lineage_units": expected["donor_lineage_units"], "rna_count_values_read": False, "observed_atac_input_count_values_read": False, "model_input_biological_artifact_read": False, "benchmark_metric_calculated": False, "prediction_bundle_read": False, "model_fit": False, "sealed_outcomes_read": False, "ranking_authorized": False, "champion_claim_allowed": False, "next_gate": "independent_evaluator_artifact_verification_then_prediction_only_smoke"}
    write_json_exclusive(output / "receipt.json", parent)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "src/masld_bench/observed_multiome_materialization.py", root / "tests/unit/test_observed_multiome_evaluator_materialization.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_evaluator_campaign", "campaign_id": CAMPAIGN_ID, "child_artifact_count": 5, "evaluator_only": True, "sealed_outcomes_accessed": False, "model_fit": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
