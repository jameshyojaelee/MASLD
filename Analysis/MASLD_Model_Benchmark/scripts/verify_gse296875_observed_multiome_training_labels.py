#!/usr/bin/env python3
"""Verify 25 held-donor-safe supervised training-label output files."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import h5py
import numpy as np

from masld_bench.artifacts import (
    freeze_tree,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
)


SCHEMA = "masld-bench-observed-multiome-training-label-verification-v1"
VERIFICATION_ID = "gse296875_observed_multiome_training_label_verification_20260825"


class TrainingLabelVerificationError(ValueError):
    """Raised when training-label outer-fold isolation differs."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise TrainingLabelVerificationError("JSON object required")
    return value


def _strings(values: Sequence[Any]) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def _tree(root: Path, record: Mapping[str, Any], label: str) -> Path:
    lexical = reject_symlink_components(root / str(record.get("path", "")), label=label)
    path = lexical.resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise TrainingLabelVerificationError(f"{label} drifted")
    return path


def _tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise TrainingLabelVerificationError(f"{path.name} schema differs")
        return [dict(row) for row in reader]


def _csr(group: Any, shape: tuple[int, int]) -> dict[str, Any]:
    if set(group) != {"data", "indices", "indptr", "shape"}:
        raise TrainingLabelVerificationError("training CSR schema differs")
    observed_shape = tuple(int(value) for value in group["shape"][:])
    data = np.asarray(group["data"][:])
    indices = np.asarray(group["indices"][:], dtype=np.int64)
    indptr = np.asarray(group["indptr"][:], dtype=np.int64)
    if (
        observed_shape != shape
        or data.ndim != 1
        or indices.shape != data.shape
        or indptr.shape != (shape[0] + 1,)
        or int(indptr[0]) != 0
        or int(indptr[-1]) != data.size
        or np.any(indptr[1:] < indptr[:-1])
        or np.any(indices < 0)
        or np.any(indices >= shape[1])
        or np.any(data < 0)
        or not np.equal(data, np.floor(data)).all()
    ):
        raise TrainingLabelVerificationError("training CSR values differ")
    return {
        "shape": list(shape),
        "nnz": int(data.size),
        "sum": int(data.sum()),
        "data_sha256": sha256(data.tobytes()).hexdigest(),
        "indices_sha256": sha256(indices.tobytes()).hexdigest(),
        "indptr_sha256": sha256(indptr.tobytes()).hexdigest(),
    }


def validate_config(root: Path, config: Mapping[str, Any]) -> dict[str, Path]:
    if config.get("schema_version") != SCHEMA or config.get("verification_id") != VERIFICATION_ID:
        raise TrainingLabelVerificationError("verification identity differs")
    model = config.get("model_input_artifact")
    if (
        not isinstance(model, dict)
        or model.get("allowed_reads") != ["row_identifiers", "observed_input_identifiers_and_coordinates", "held_target_identifiers_and_coordinates"]
        or model.get("biological_count_values_read") is not False
    ):
        raise TrainingLabelVerificationError("model input read contract differs")
    if config.get("expected") != {"children": 25, "donor_folds": 5, "genomic_folds": 5, "training_targets": 1000, "training_units_by_held_donor_fold": [140, 150, 160, 175, 155], "minimum_observed_input_distance_bp": 524288}:
        raise TrainingLabelVerificationError("expected census differs")
    firewall = config.get("firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise TrainingLabelVerificationError("verification firewall is open")
    return {
        "labels": _tree(root, config["training_label_artifact"], "training labels"),
        "model": _tree(root, model, "model inputs"),
        "training_masks": _tree(root, config["training_target_mask_plan"], "training masks"),
        "base_masks": _tree(root, config["base_mask_plan"], "base masks"),
    }


def verify(root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    resolved = validate_config(root, config)
    expected = config["expected"]
    parent = _json(resolved["labels"] / "receipt.json")
    model_parent = _json(resolved["model"] / "receipt.json")
    children = parent.get("child_artifacts")
    model_children = model_parent.get("child_artifacts")
    if not isinstance(children, list) or len(children) != expected["children"] or not isinstance(model_children, list) or len(model_children) != expected["genomic_folds"]:
        raise TrainingLabelVerificationError("child roster differs")
    model_by_fold = {int(record["held_genomic_fold"]): record for record in model_children}
    target_masks = _tsv(resolved["training_masks"] / "training_target_peaks.tsv", ("training_target_hash", "chromosome", "bed_start_0based", "bed_end_half_open", "source_genomic_fold", "held_genomic_fold", "nearest_observed_input_distance_bp"))
    surfaces: list[dict[str, Any]] = []
    seen: set[tuple[int, int]] = set()
    for record in children:
        genomic_fold = int(record["held_genomic_fold"])
        donor_fold = int(record["held_donor_fold"])
        surface = (genomic_fold, donor_fold)
        if surface in seen or genomic_fold not in model_by_fold:
            raise TrainingLabelVerificationError("duplicate or unbound outer surface")
        seen.add(surface)
        child = resolved["labels"] / str(record["path"])
        model_child = resolved["model"] / str(model_by_fold[genomic_fold]["path"])
        verify_frozen_tree(child)
        verify_frozen_tree(model_child)
        if _digest(child / "ARTIFACTS.json") != record.get("artifacts_sha256") or _digest(model_child / "ARTIFACTS.json") != model_by_fold[genomic_fold].get("artifacts_sha256"):
            raise TrainingLabelVerificationError("surface artifact drifted")
        receipt = _json(child / "receipt.json")
        with h5py.File(child / "training_labels.h5", "r") as labels, h5py.File(model_child / "model_input.h5", "r") as model:
            if (
                set(labels) != {"rows", "training_target_atac"}
                or labels.attrs.get("schema_version") != "masld-bench-observed-multiome-training-label-h5-v1"
                or not bool(labels.attrs.get("training_only", False))
                or int(labels.attrs.get("held_genomic_fold", -1)) != genomic_fold
                or int(labels.attrs.get("held_donor_fold", -1)) != donor_fold
            ):
                raise TrainingLabelVerificationError("training HDF5 root differs")
            rows = labels["rows"]
            if set(rows) != {"row_hash", "unit_hash", "lineage", "donor_fold", "nuclei"}:
                raise TrainingLabelVerificationError("training row schema differs")
            row_hashes = _strings(rows["row_hash"][:])
            unit_hashes = _strings(rows["unit_hash"][:])
            lineages = _strings(rows["lineage"][:])
            donor_folds = [int(value) for value in rows["donor_fold"][:]]
            nuclei = [int(value) for value in rows["nuclei"][:]]
            model_row_hashes = _strings(model["rows/row_hash"][:])
            model_unit_hashes = _strings(model["rows/unit_hash"][:])
            model_lineages = _strings(model["rows/lineage"][:])
            model_donor_folds = [int(value) for value in model["rows/donor_fold"][:]]
            model_nuclei = [int(value) for value in model["rows/nuclei"][:]]
            keep = [index for index, value in enumerate(model_donor_folds) if value != donor_fold]
            if (
                len(row_hashes) != expected["training_units_by_held_donor_fold"][donor_fold]
                or donor_fold in donor_folds
                or row_hashes != [model_row_hashes[index] for index in keep]
                or unit_hashes != [model_unit_hashes[index] for index in keep]
                or lineages != [model_lineages[index] for index in keep]
                or donor_folds != [model_donor_folds[index] for index in keep]
                or nuclei != [model_nuclei[index] for index in keep]
            ):
                raise TrainingLabelVerificationError("outer-training row identity differs")
            targets = [row for row in target_masks if int(row["held_genomic_fold"]) == genomic_fold]
            target = labels["training_target_atac"]
            if (
                set(target) != {"counts_csr", "training_target_hash", "chromosome", "bed_start_0based", "bed_end_half_open", "source_genomic_fold"}
                or target.attrs.get("scale") != "raw_integer_outer_training_donor_by_lineage_sum"
                or _strings(target["training_target_hash"][:]) != [row["training_target_hash"] for row in targets]
                or _strings(target["chromosome"][:]) != [row["chromosome"] for row in targets]
                or [int(value) for value in target["bed_start_0based"][:]] != [int(row["bed_start_0based"]) for row in targets]
                or [int(value) for value in target["bed_end_half_open"][:]] != [int(row["bed_end_half_open"]) for row in targets]
                or [int(value) for value in target["source_genomic_fold"][:]] != [int(row["source_genomic_fold"]) for row in targets]
                or any(int(row["source_genomic_fold"]) == genomic_fold for row in targets)
                or any(int(row["nearest_observed_input_distance_bp"]) != -1 and int(row["nearest_observed_input_distance_bp"]) < expected["minimum_observed_input_distance_bp"] for row in targets)
            ):
                raise TrainingLabelVerificationError("training target identity or buffer differs")
            training_coordinates = {(row["chromosome"], int(row["bed_start_0based"]), int(row["bed_end_half_open"])) for row in targets}
            observed_coordinates = set(zip(_strings(model["observed_atac_input/chromosome"][:]), [int(value) for value in model["observed_atac_input/bed_start_0based"][:]], [int(value) for value in model["observed_atac_input/bed_end_half_open"][:]], strict=True))
            held_coordinates = set(zip(_strings(model["targets_without_values/chromosome"][:]), [int(value) for value in model["targets_without_values/bed_start_0based"][:]], [int(value) for value in model["targets_without_values/bed_end_half_open"][:]], strict=True))
            if training_coordinates.intersection(observed_coordinates) or training_coordinates.intersection(held_coordinates):
                raise TrainingLabelVerificationError("training targets overlap model inputs or held targets")
            identity = _csr(target["counts_csr"], (len(row_hashes), expected["training_targets"]))
        if (
            receipt.get("target_identity") != identity
            or receipt.get("held_donor_training_target_value_positions_read") != 0
            or receipt.get("held_donor_rows_written") is not False
            or receipt.get("metric_calculated") is not False
            or receipt.get("model_fit") is not False
        ):
            raise TrainingLabelVerificationError("surface receipt differs")
        surfaces.append({"held_genomic_fold": genomic_fold, "held_donor_fold": donor_fold, "shape": identity["shape"], "held_donor_rows_present": False, "nonnegative_integer_training_counts": True})
    if seen != {(genomic, donor) for genomic in range(expected["genomic_folds"]) for donor in range(expected["donor_folds"])}:
        raise TrainingLabelVerificationError("outer surface rectangle differs")
    return {
        "schema_version": "masld-bench-observed-multiome-training-label-verification-receipt-v1",
        "verification_id": VERIFICATION_ID,
        "training_label_artifacts_sha256": config["training_label_artifact"]["artifacts_sha256"],
        "model_input_artifacts_sha256": config["model_input_artifact"]["artifacts_sha256"],
        "surfaces_verified": sorted(surfaces, key=lambda row: (row["held_genomic_fold"], row["held_donor_fold"])),
        "surface_count": len(surfaces),
        "held_donor_rows_present": False,
        "model_input_biological_count_values_read": False,
        "evaluator_artifact_read": False,
        "training_count_values_read_by_verifier": True,
        "benchmark_metric_calculated": False,
        "prediction_bundle_read": False,
        "model_fit": False,
        "sealed_outcomes_read": False,
        "promotion_gate_passed": True
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
    output = reject_symlink_components(args.output, label="training label verification output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    receipt = verify(root, _json(config_path))
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "tests/unit/test_verify_observed_multiome_training_labels.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_training_label_verification", "promotion_gate_passed": True, "held_donor_rows_present": False, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
