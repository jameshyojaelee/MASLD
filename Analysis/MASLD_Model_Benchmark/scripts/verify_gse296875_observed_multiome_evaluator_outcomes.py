#!/usr/bin/env python3
"""Verify isolated GSE296875 evaluator outcomes without scoring predictions."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
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


SCHEMA = "masld-bench-observed-multiome-evaluator-verification-v1"
VERIFICATION_ID = "gse296875_observed_multiome_evaluator_verification_20260825"


class EvaluatorVerificationError(ValueError):
    """Raised when evaluator output files or their isolation differ."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise EvaluatorVerificationError("JSON object required")
    return value


def _strings(values: Sequence[Any]) -> list[str]:
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def _tree(root: Path, record: Mapping[str, Any], label: str) -> Path:
    lexical = reject_symlink_components(root / str(record.get("path", "")), label=label)
    path = lexical.resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise EvaluatorVerificationError(f"{label} drifted")
    return path


def _tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise EvaluatorVerificationError(f"{path.name} schema differs")
        return [dict(row) for row in reader]


def _csr(group: Any, shape: tuple[int, int]) -> dict[str, Any]:
    if set(group) != {"data", "indices", "indptr", "shape"}:
        raise EvaluatorVerificationError("target CSR schema differs")
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
        raise EvaluatorVerificationError("target CSR values differ")
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
        raise EvaluatorVerificationError("verification identity differs")
    model = config.get("model_input_artifact")
    if (
        not isinstance(model, dict)
        or model.get("allowed_reads") != ["row_identifiers", "target_identifiers_and_coordinates"]
        or model.get("biological_count_values_read") is not False
    ):
        raise EvaluatorVerificationError("model input read contract differs")
    firewall = config.get("firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise EvaluatorVerificationError("verification firewall is open")
    return {
        "evaluator": _tree(root, config["evaluator_artifact"], "evaluator artifact"),
        "model": _tree(root, model, "model input artifact"),
        "masks": _tree(root, config["mask_plan_artifact"], "mask plan artifact"),
    }


def verify(root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    resolved = validate_config(root, config)
    expected = config["expected"]
    evaluator_parent = _json(resolved["evaluator"] / "receipt.json")
    model_parent = _json(resolved["model"] / "receipt.json")
    evaluator_children = evaluator_parent.get("child_artifacts")
    model_children = model_parent.get("child_artifacts")
    if not isinstance(evaluator_children, list) or not isinstance(model_children, list) or len(evaluator_children) != expected["children"] or len(model_children) != expected["children"]:
        raise EvaluatorVerificationError("child roster differs")
    target_masks = _tsv(resolved["masks"] / "target_peaks.tsv", ("target_hash", "chromosome", "bed_start_0based", "bed_end_half_open", "genomic_fold"))
    planned_units = _tsv(resolved["masks"] / "donor_lineage_units.tsv", ("donor_hash", "lineage", "donor_fold", "nuclei"))
    model_by_fold = {int(record["held_genomic_fold"]): record for record in model_children}
    folds: list[dict[str, Any]] = []
    for record in evaluator_children:
        fold = int(record["held_genomic_fold"])
        if fold not in model_by_fold:
            raise EvaluatorVerificationError("model sibling fold is absent")
        evaluator_child = resolved["evaluator"] / str(record["path"])
        model_record = model_by_fold[fold]
        model_child = resolved["model"] / str(model_record["path"])
        verify_frozen_tree(evaluator_child)
        verify_frozen_tree(model_child)
        if _digest(evaluator_child / "ARTIFACTS.json") != record.get("artifacts_sha256") or _digest(model_child / "ARTIFACTS.json") != model_record.get("artifacts_sha256"):
            raise EvaluatorVerificationError("child artifact drifted")
        receipt = _json(evaluator_child / "receipt.json")
        with h5py.File(evaluator_child / "evaluator_outcomes.h5", "r") as outcome, h5py.File(model_child / "model_input.h5", "r") as model:
            if (
                set(outcome) != {"rows", "target_atac"}
                or outcome.attrs.get("schema_version") != "masld-bench-observed-multiome-evaluator-outcome-h5-v1"
                or not bool(outcome.attrs.get("evaluator_only", False))
                or int(outcome.attrs.get("held_genomic_fold", -1)) != fold
            ):
                raise EvaluatorVerificationError("evaluator HDF5 root differs")
            rows = outcome["rows"]
            if set(rows) != {"row_hash", "unit_hash", "lineage", "donor_fold", "nuclei"}:
                raise EvaluatorVerificationError("evaluator row schema differs")
            row_hashes = _strings(rows["row_hash"][:])
            unit_hashes = _strings(rows["unit_hash"][:])
            lineages = _strings(rows["lineage"][:])
            donor_folds = [int(value) for value in rows["donor_fold"][:]]
            nuclei = [int(value) for value in rows["nuclei"][:]]
            if (
                row_hashes != _strings(model["rows/row_hash"][:])
                or unit_hashes != _strings(model["rows/unit_hash"][:])
                or lineages != _strings(model["rows/lineage"][:])
                or donor_folds != [int(value) for value in model["rows/donor_fold"][:]]
                or nuclei != [int(value) for value in model["rows/nuclei"][:]]
                or unit_hashes != [row["donor_hash"] for row in planned_units]
                or lineages != [row["lineage"] for row in planned_units]
                or donor_folds != [int(row["donor_fold"]) for row in planned_units]
                or nuclei != [int(row["nuclei"]) for row in planned_units]
            ):
                raise EvaluatorVerificationError("evaluator/model row identity differs")
            unit_lineages: dict[str, set[str]] = defaultdict(set)
            unit_folds: dict[str, set[int]] = defaultdict(set)
            for unit, lineage, donor_fold in zip(unit_hashes, lineages, donor_folds, strict=True):
                unit_lineages[unit].add(lineage)
                unit_folds[unit].add(donor_fold)
            census = Counter(next(iter(values)) for values in unit_folds.values())
            if len(unit_folds) != expected["donors"] or any(len(values) != expected["lineages"] for values in unit_lineages.values()) or any(len(values) != 1 for values in unit_folds.values()) or [census[index] for index in range(5)] != expected["donor_fold_census"] or sum(nuclei) != expected["nuclei"]:
                raise EvaluatorVerificationError("evaluator donor grouping differs")
            target_rows = [row for row in target_masks if int(row["genomic_fold"]) == fold]
            target = outcome["target_atac"]
            model_targets = model["targets_without_values"]
            if (
                set(target) != {"counts_csr", "target_hash", "chromosome", "bed_start_0based", "bed_end_half_open"}
                or target.attrs.get("scale") != "raw_integer_donor_by_lineage_sum"
                or _strings(target["target_hash"][:]) != [row["target_hash"] for row in target_rows]
                or _strings(target["target_hash"][:]) != _strings(model_targets["target_hash"][:])
                or _strings(target["chromosome"][:]) != _strings(model_targets["chromosome"][:])
                or [int(value) for value in target["bed_start_0based"][:]] != [int(value) for value in model_targets["bed_start_0based"][:]]
                or [int(value) for value in target["bed_end_half_open"][:]] != [int(value) for value in model_targets["bed_end_half_open"][:]]
            ):
                raise EvaluatorVerificationError("evaluator target identity differs")
            identity = _csr(target["counts_csr"], (expected["donor_lineage_units"], expected["target_peaks_per_fold"]))
        if (
            receipt.get("held_genomic_fold") != fold
            or receipt.get("target_identity") != identity
            or receipt.get("rna_count_values_read") is not False
            or receipt.get("observed_atac_input_count_values_read") is not False
            or receipt.get("prediction_bundle_read") is not False
            or receipt.get("metric_calculated") is not False
        ):
            raise EvaluatorVerificationError("evaluator receipt differs")
        folds.append({"held_genomic_fold": fold, "shape": identity["shape"], "nonnegative_integer_target_counts": True})
    if {row["held_genomic_fold"] for row in folds} != set(range(expected["children"])):
        raise EvaluatorVerificationError("held genomic folds differ")
    return {
        "schema_version": "masld-bench-observed-multiome-evaluator-verification-receipt-v1",
        "verification_id": VERIFICATION_ID,
        "evaluator_artifacts_sha256": config["evaluator_artifact"]["artifacts_sha256"],
        "model_input_artifacts_sha256": config["model_input_artifact"]["artifacts_sha256"],
        "folds_verified": sorted(folds, key=lambda row: row["held_genomic_fold"]),
        "model_input_biological_count_values_read": False,
        "rna_count_values_read": False,
        "observed_atac_input_count_values_read": False,
        "target_atac_count_values_read_by_evaluator_verifier": True,
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
    output = reject_symlink_components(args.output, label="evaluator verification output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    receipt = verify(root, _json(config_path))
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "tests/unit/test_verify_observed_multiome_evaluator_outcomes.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_evaluator_verification", "evaluator_only": True, "promotion_gate_passed": True, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
