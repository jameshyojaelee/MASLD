#!/usr/bin/env python3
"""Build 25 held-donor-safe supervised ATAC training-label output files."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy import sparse

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
    validate_csr_structure,
    write_csr,
    write_string_dataset,
)


SCHEMA = "masld-bench-observed-multiome-training-label-materialization-v1"
CAMPAIGN_ID = "gse296875_observed_multiome_training_labels_20260825"


class TrainingLabelMaterializationError(ValueError):
    """Raised when a training label contains a held donor or genomic value."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise TrainingLabelMaterializationError("JSON object required")
    return value


def _tree(root: Path, record: Mapping[str, Any], label: str) -> Path:
    lexical = reject_symlink_components(root / str(record.get("path", "")), label=label)
    path = lexical.resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise TrainingLabelMaterializationError(f"{label} drifted")
    return path


def _tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise TrainingLabelMaterializationError(f"{path.name} schema differs")
        return [dict(row) for row in reader]


def validate_config(root: Path, config: Mapping[str, Any]) -> dict[str, Path]:
    if config.get("schema_version") != SCHEMA or config.get("campaign_id") != CAMPAIGN_ID or config.get("dataset_id") != "gse296875" or config.get("stage") != "smoke":
        raise TrainingLabelMaterializationError("campaign identity differs")
    parents = config.get("parents")
    if not isinstance(parents, dict) or set(parents) != {"supervised_training_contract", "training_target_mask_plan", "base_mask_plan"}:
        raise TrainingLabelMaterializationError("parent roster differs")
    resolved = {key: _tree(root, record, key) for key, record in parents.items()}
    contract = _json(resolved["supervised_training_contract"] / "receipt.json")
    mask = _json(resolved["training_target_mask_plan"] / "receipt.json")
    if contract.get("next_gate") != "identifier_only_training_target_mask_plan" or mask.get("next_gate") != "supervised_training_label_materialization" or mask.get("atac_count_values_read") is not False:
        raise TrainingLabelMaterializationError("parent disposition differs")
    h5 = config.get("input_h5")
    if not isinstance(h5, dict):
        raise TrainingLabelMaterializationError("input HDF5 is absent")
    lexical = reject_symlink_components(root / str(h5.get("path", "")), label="input HDF5")
    h5_path = lexical.resolve(strict=True)
    h5_path.relative_to(root)
    if not h5_path.is_file() or _digest(h5_path) != h5.get("sha256"):
        raise TrainingLabelMaterializationError("input HDF5 drifted")
    resolved["input_h5"] = h5_path
    if config.get("expected") != {"nuclei": 1000, "donors": 39, "lineages": 5, "donor_lineage_units": 195, "donor_folds": 5, "genomic_folds": 5, "children": 25, "training_targets_per_held_genomic_fold": 1000, "donor_fold_census": [11, 9, 7, 4, 8], "training_units_by_held_donor_fold": [140, 150, 160, 175, 155], "atac_features": 306706}:
        raise TrainingLabelMaterializationError("expected census differs")
    if config.get("hash_namespaces") != {"donor": "gse296875:observed_multiome:donor:v1", "row": "gse296875:observed_multiome:row:v1", "training_target": "gse296875:observed_multiome:training_target:v1"}:
        raise TrainingLabelMaterializationError("hash namespaces differ")
    artifact = config.get("artifact")
    if not isinstance(artifact, dict) or artifact != {"one_separately_frozen_child_per_held_genomic_by_held_donor_fold": True, "target_values": "raw_integer_outer_training_donor_by_lineage_sum", "held_donor_rows": "forbidden", "rna_values": "forbidden", "observed_atac_input_values": "forbidden", "evaluator_values": "forbidden", "raw_donor_ids_exported": False, "normalization": "none", "feature_selection": "none"}:
        raise TrainingLabelMaterializationError("artifact contract differs")
    if config.get("selective_read") != {"atac_indices_and_indptr_read": True, "atac_data_values_read_only_for_training_target_columns_and_outer_training_cell_rows": True, "full_atac_data_array_read": False, "held_donor_training_target_value_positions_read_per_surface": 0}:
        raise TrainingLabelMaterializationError("selective read contract differs")
    firewall = config.get("firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise TrainingLabelMaterializationError("materialization firewall is open")
    return resolved


def selectively_read_csr_rows_columns(group: Any, selected_rows: Sequence[int], selected_columns: Sequence[int]) -> tuple[Any, np.ndarray, np.ndarray]:
    shape = tuple(int(value) for value in group["shape"][:])
    indices = np.asarray(group["indices"][:], dtype=np.int64)
    indptr = np.asarray(group["indptr"][:], dtype=np.int64)
    rows, columns = validate_csr_structure(indices, indptr, shape)
    row_selection = np.asarray(selected_rows, dtype=np.int64)
    column_selection = np.asarray(selected_columns, dtype=np.int64)
    if (
        row_selection.ndim != 1
        or column_selection.ndim != 1
        or row_selection.size < 1
        or column_selection.size < 1
        or len(set(map(int, row_selection))) != row_selection.size
        or len(set(map(int, column_selection))) != column_selection.size
        or np.any(row_selection < 0)
        or np.any(row_selection >= rows)
        or np.any(column_selection < 0)
        or np.any(column_selection >= columns)
        or np.any(row_selection[1:] <= row_selection[:-1])
    ):
        raise TrainingLabelMaterializationError("selected sparse rows or columns differ")
    column_positions = selected_value_positions(indices, column_selection, total_columns=columns)
    source_rows = np.searchsorted(indptr[1:], column_positions, side="right")
    row_membership = np.zeros(rows, dtype=bool)
    row_membership[row_selection] = True
    retained = row_membership[source_rows]
    positions = column_positions[retained]
    retained_source_rows = source_rows[retained]
    row_remap = np.full(rows, -1, dtype=np.int64)
    row_remap[row_selection] = np.arange(row_selection.size, dtype=np.int64)
    column_remap = np.full(columns, -1, dtype=np.int64)
    column_remap[column_selection] = np.arange(column_selection.size, dtype=np.int64)
    data = np.asarray(group["data"][positions])
    matrix = sparse.coo_matrix(
        (data, (row_remap[retained_source_rows], column_remap[indices[positions]])),
        shape=(row_selection.size, column_selection.size),
    ).tocsr()
    if matrix.nnz != positions.size or np.any(matrix.data < 0):
        raise TrainingLabelMaterializationError("selective training CSR values differ")
    return matrix, positions, retained_source_rows


def _resolve_target_columns(peak_ids: Sequence[str], chromosomes: Sequence[str], starts: Sequence[int], ends: Sequence[int], target_rows: Sequence[Mapping[str, str]], *, held_genomic_fold: int, namespace: str) -> list[int]:
    by_hash = {
        salted_hash(namespace, held_genomic_fold, peak_id): (index, chromosome, int(start), int(end))
        for index, (peak_id, chromosome, start, end) in enumerate(zip(peak_ids, chromosomes, starts, ends, strict=True))
    }
    columns = []
    for row in target_rows:
        metadata = by_hash.get(row["training_target_hash"])
        if metadata is None or metadata[1:] != (row["chromosome"], int(row["bed_start_0based"]), int(row["bed_end_half_open"])):
            raise TrainingLabelMaterializationError("training target identity differs")
        columns.append(metadata[0])
    if len(columns) != len(set(columns)):
        raise TrainingLabelMaterializationError("training target columns are duplicated")
    return columns


def _write_child(path: Path, held_genomic_fold: int, held_donor_fold: int, planned: Sequence[Mapping[str, str]], row_hashes: Sequence[str], counts: Any, targets: Sequence[Mapping[str, str]]) -> dict[str, Any]:
    import h5py

    with h5py.File(path, "x") as handle:
        handle.attrs["schema_version"] = "masld-bench-observed-multiome-training-label-h5-v1"
        handle.attrs["dataset_id"] = "gse296875"
        handle.attrs["stage"] = "smoke"
        handle.attrs["held_genomic_fold"] = held_genomic_fold
        handle.attrs["held_donor_fold"] = held_donor_fold
        handle.attrs["biological_unit"] = "donor_by_lineage"
        handle.attrs["training_only"] = True
        rows = handle.create_group("rows")
        write_string_dataset(rows, "row_hash", row_hashes)
        write_string_dataset(rows, "unit_hash", [row["donor_hash"] for row in planned])
        write_string_dataset(rows, "lineage", [row["lineage"] for row in planned])
        rows.create_dataset("donor_fold", data=np.asarray([int(row["donor_fold"]) for row in planned], dtype=np.int8))
        rows.create_dataset("nuclei", data=np.asarray([int(row["nuclei"]) for row in planned], dtype=np.int16))
        target = handle.create_group("training_target_atac")
        identity = write_csr(target.create_group("counts_csr"), counts)
        write_string_dataset(target, "training_target_hash", [row["training_target_hash"] for row in targets])
        write_string_dataset(target, "chromosome", [row["chromosome"] for row in targets])
        target.create_dataset("bed_start_0based", data=np.asarray([int(row["bed_start_0based"]) for row in targets], dtype=np.int64))
        target.create_dataset("bed_end_half_open", data=np.asarray([int(row["bed_end_half_open"]) for row in targets], dtype=np.int64))
        target.create_dataset("source_genomic_fold", data=np.asarray([int(row["source_genomic_fold"]) for row in targets], dtype=np.int8))
        target.attrs["scale"] = "raw_integer_outer_training_donor_by_lineage_sum"
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
    expected = config["expected"]
    namespaces = config["hash_namespaces"]
    planned_all = _tsv(resolved["base_mask_plan"] / "donor_lineage_units.tsv", ("donor_hash", "lineage", "donor_fold", "nuclei"))
    target_all = _tsv(resolved["training_target_mask_plan"] / "training_target_peaks.tsv", ("training_target_hash", "chromosome", "bed_start_0based", "bed_end_half_open", "source_genomic_fold", "held_genomic_fold", "nearest_observed_input_distance_bp"))
    output = reject_symlink_components(args.output, label="training label output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    children: list[dict[str, Any]] = []
    import h5py

    with h5py.File(resolved["input_h5"], "r") as handle:
        donors = decode_strings(handle["obs/donor_id"][:], label="donor")
        lineages = decode_strings(handle["obs/broad_label"][:], label="lineage")
        cell_units = [(salted_hash(namespaces["donor"], donor), lineage) for donor, lineage in zip(donors, lineages, strict=True)]
        donor_fold_by_hash: dict[str, int] = {}
        for row in planned_all:
            prior = donor_fold_by_hash.setdefault(row["donor_hash"], int(row["donor_fold"]))
            if prior != int(row["donor_fold"]):
                raise TrainingLabelMaterializationError("donor fold differs across lineages")
        peak_ids = decode_strings(handle["atac/peak_id"][:], label="peak ID")
        chromosomes = decode_strings(handle["atac/chromosome"][:], label="chromosome")
        starts = [int(value) for value in handle["atac/bed_start_0based"][:]]
        ends = [int(value) for value in handle["atac/bed_end_half_open"][:]]
        atac = handle["atac/counts_csr"]
        if tuple(int(value) for value in atac["shape"][:]) != (expected["nuclei"], expected["atac_features"]):
            raise TrainingLabelMaterializationError("ATAC matrix shape differs")
        for held_genomic_fold in range(expected["genomic_folds"]):
            targets = [row for row in target_all if int(row["held_genomic_fold"]) == held_genomic_fold]
            columns = _resolve_target_columns(peak_ids, chromosomes, starts, ends, targets, held_genomic_fold=held_genomic_fold, namespace=namespaces["training_target"])
            if len(columns) != expected["training_targets_per_held_genomic_fold"] or any(int(row["source_genomic_fold"]) == held_genomic_fold for row in targets):
                raise TrainingLabelMaterializationError("training target fold differs")
            for held_donor_fold in range(expected["donor_folds"]):
                selected_cells = [index for index, (donor_hash, _) in enumerate(cell_units) if donor_fold_by_hash[donor_hash] != held_donor_fold]
                selected_units = [row for row in planned_all if int(row["donor_fold"]) != held_donor_fold]
                selected_cell_units = [cell_units[index] for index in selected_cells]
                aggregation, _ = donor_lineage_aggregation(cell_unit_keys=selected_cell_units, planned_units=selected_units)
                cell_counts, positions, source_rows = selectively_read_csr_rows_columns(atac, selected_cells, columns)
                if any(donor_fold_by_hash[cell_units[int(index)][0]] == held_donor_fold for index in source_rows):
                    raise TrainingLabelMaterializationError("held donor target value position was read")
                counts = (aggregation @ cell_counts).tocsr()
                if len(selected_units) != expected["training_units_by_held_donor_fold"][held_donor_fold]:
                    raise TrainingLabelMaterializationError("outer training unit census differs")
                row_hashes = [salted_hash(namespaces["row"], row["donor_hash"], row["lineage"]) for row in selected_units]
                child = output / f"genomic_{held_genomic_fold}_donor_{held_donor_fold}"
                child.mkdir()
                identity = _write_child(child / "training_labels.h5", held_genomic_fold, held_donor_fold, selected_units, row_hashes, counts, targets)
                receipt = {"schema_version": "masld-bench-observed-multiome-training-label-receipt-v1", "dataset_id": "gse296875", "stage": "smoke", "held_genomic_fold": held_genomic_fold, "held_donor_fold": held_donor_fold, "outer_training_donors": expected["donors"] - expected["donor_fold_census"][held_donor_fold], "outer_training_donor_lineage_units": len(selected_units), "training_targets": len(targets), "target_identity": identity, "training_sparse_value_positions_read": int(positions.size), "held_donor_training_target_value_positions_read": 0, "held_donor_rows_written": False, "rna_count_values_read": False, "observed_atac_input_count_values_read": False, "evaluator_artifact_read": False, "model_input_biological_artifact_read": False, "full_atac_data_array_read": False, "metric_calculated": False, "prediction_bundle_read": False, "model_fit": False, "sealed_outcomes_read": False}
                write_json_exclusive(child / "receipt.json", receipt)
                child_hash = freeze_tree(child, metadata={"artifact_class": "gse296875_observed_multiome_supervised_training_labels", "held_genomic_fold": held_genomic_fold, "held_donor_fold": held_donor_fold, "held_donor_rows_present": False, "sealed_outcomes_accessed": False})
                children.append({"held_genomic_fold": held_genomic_fold, "held_donor_fold": held_donor_fold, "path": child.relative_to(output).as_posix(), "artifacts_sha256": child_hash})
    if len(children) != expected["children"]:
        raise TrainingLabelMaterializationError("training label child census differs")
    parent = {"schema_version": "masld-bench-observed-multiome-training-label-campaign-receipt-v1", "campaign_id": CAMPAIGN_ID, "dataset_id": "gse296875", "stage": "smoke", "child_artifacts": children, "child_artifact_count": len(children), "held_donor_training_target_value_positions_read_per_surface": 0, "held_donor_rows_written": False, "rna_count_values_read": False, "observed_atac_input_count_values_read": False, "evaluator_artifact_read": False, "model_input_biological_artifact_read": False, "full_atac_data_array_read": False, "development_metric_calculated": False, "prediction_bundle_read": False, "model_fit": False, "sealed_outcomes_read": False, "ranking_authorized": False, "champion_claim_allowed": False, "next_gate": "independent_training_label_verification"}
    write_json_exclusive(output / "receipt.json", parent)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "src/masld_bench/observed_multiome_materialization.py", root / "tests/unit/test_observed_multiome_training_label_materialization.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_training_label_campaign", "campaign_id": CAMPAIGN_ID, "child_artifact_count": 25, "held_donor_rows_present": False, "sealed_outcomes_accessed": False, "model_fit": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
