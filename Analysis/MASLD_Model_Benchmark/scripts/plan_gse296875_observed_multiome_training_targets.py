#!/usr/bin/env python3
"""Select count-blind supervised targets outside observed-ATAC masks."""

from __future__ import annotations

import argparse
from bisect import bisect_left
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from masld_bench.artifacts import (
    freeze_tree,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
)


SCHEMA = "masld-bench-observed-multiome-training-target-mask-plan-v1"
PLAN_ID = "gse296875_observed_multiome_training_target_masks_20260825"


class TrainingTargetMaskPlanError(ValueError):
    """Raised when a training target can leak observed or held ATAC values."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise TrainingTargetMaskPlanError("JSON object required")
    return value


def _tree(root: Path, record: Mapping[str, Any], label: str) -> Path:
    lexical = reject_symlink_components(root / str(record.get("path", "")), label=label)
    path = lexical.resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise TrainingTargetMaskPlanError(f"{label} drifted")
    return path


def _tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise TrainingTargetMaskPlanError(f"{path.name} schema differs")
        return [dict(row) for row in reader]


def _write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def validate_config(root: Path, config: Mapping[str, Any]) -> dict[str, Path]:
    if config.get("schema_version") != SCHEMA or config.get("plan_id") != PLAN_ID or config.get("dataset_id") != "gse296875" or config.get("stage") != "smoke":
        raise TrainingTargetMaskPlanError("plan identity differs")
    parents = config.get("parents")
    if not isinstance(parents, dict) or set(parents) != {"supervised_training_contract", "base_mask_plan"}:
        raise TrainingTargetMaskPlanError("parent roster differs")
    resolved = {key: _tree(root, record, key) for key, record in parents.items()}
    contract = _json(resolved["supervised_training_contract"] / "receipt.json")
    if contract.get("next_gate") != "identifier_only_training_target_mask_plan" or contract.get("supervised_fit_authorized") is not False:
        raise TrainingTargetMaskPlanError("supervised contract disposition differs")
    h5 = config.get("input_h5")
    allowed = ["atac/peak_id", "atac/chromosome", "atac/bed_start_0based", "atac/bed_end_half_open"]
    if not isinstance(h5, dict) or h5.get("allowed_datasets") != allowed:
        raise TrainingTargetMaskPlanError("HDF5 allowlist differs")
    lexical = reject_symlink_components(root / str(h5.get("path", "")), label="input HDF5")
    h5_path = lexical.resolve(strict=True)
    h5_path.relative_to(root)
    if not h5_path.is_file() or _digest(h5_path) != h5.get("sha256"):
        raise TrainingTargetMaskPlanError("input HDF5 drifted")
    resolved["input_h5"] = h5_path
    if config.get("selection") != {"genomic_folds": 5, "targets_per_nonheld_source_genomic_fold": 250, "targets_per_held_genomic_fold": 1000, "minimum_distance_from_matching_observed_input_bp": 524288, "exclude_all_base_evaluator_targets": True, "exclude_matching_observed_input_peaks": True, "identifier_hash_namespace": "gse296875:observed_multiome:training_target:v1", "selection_seed": 20260825, "selection_uses_count_values": False, "insufficient_candidates": "fail_closed"}:
        raise TrainingTargetMaskPlanError("selection contract differs")
    firewall = config.get("firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise TrainingTargetMaskPlanError("plan firewall is open")
    return resolved


def read_peak_axes(path: Path) -> list[dict[str, Any]]:
    import h5py

    with h5py.File(path, "r") as handle:
        if (
            handle.attrs.get("schema_version") != "masld-bench-multimodal-h5-v1"
            or handle.attrs.get("view_id") != "gse296875_rna_atac_smoke_1000_v1"
            or handle.attrs.get("bed_coordinate_system") != "0_based_half_open"
        ):
            raise TrainingTargetMaskPlanError("input HDF5 attributes differ")
        peak_ids = [value.decode() if isinstance(value, bytes) else str(value) for value in handle["atac/peak_id"][:]]
        chromosomes = [value.decode() if isinstance(value, bytes) else str(value) for value in handle["atac/chromosome"][:]]
        starts = [int(value) for value in handle["atac/bed_start_0based"][:]]
        ends = [int(value) for value in handle["atac/bed_end_half_open"][:]]
    if not (len(peak_ids) == len(chromosomes) == len(starts) == len(ends) == len(set(peak_ids))):
        raise TrainingTargetMaskPlanError("peak axes differ")
    if any(start < 0 or end <= start for start, end in zip(starts, ends, strict=True)):
        raise TrainingTargetMaskPlanError("peak coordinates differ")
    return [{"peak_id": peak_id, "chromosome": chromosome, "bed_start_0based": start, "bed_end_half_open": end} for peak_id, chromosome, start, end in zip(peak_ids, chromosomes, starts, ends, strict=True)]


def _expanded_intervals(inputs: Sequence[Mapping[str, Any]], buffer_bp: int) -> dict[str, tuple[list[int], list[tuple[int, int]]]]:
    grouped: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in inputs:
        grouped[str(row["chromosome"])].append((max(0, int(row["bed_start_0based"]) - buffer_bp), int(row["bed_end_half_open"]) + buffer_bp))
    result: dict[str, tuple[list[int], list[tuple[int, int]]]] = {}
    for chromosome, intervals in grouped.items():
        merged: list[tuple[int, int]] = []
        for start, end in sorted(intervals):
            if not merged or start > merged[-1][1]:
                merged.append((start, end))
            else:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        result[chromosome] = ([start for start, _ in merged], merged)
    return result


def _outside_expanded(row: Mapping[str, Any], expanded: Mapping[str, tuple[list[int], list[tuple[int, int]]]]) -> bool:
    indexed = expanded.get(str(row["chromosome"]))
    if indexed is None:
        return True
    starts, intervals = indexed
    position = bisect_left(starts, int(row["bed_end_half_open"]))
    if position == 0:
        return True
    left, right = intervals[position - 1]
    return int(row["bed_start_0based"]) >= right or int(row["bed_end_half_open"]) <= left


def _nearest_distance(row: Mapping[str, Any], inputs: Sequence[Mapping[str, Any]]) -> int:
    distances = []
    start = int(row["bed_start_0based"])
    end = int(row["bed_end_half_open"])
    for other in inputs:
        if other["chromosome"] != row["chromosome"]:
            continue
        other_start = int(other["bed_start_0based"])
        other_end = int(other["bed_end_half_open"])
        distances.append(max(0, other_start - end, start - other_end))
    return min(distances) if distances else -1


def select_training_targets(
    peaks: Sequence[Mapping[str, Any]],
    genomic_assignment: Mapping[str, int],
    inputs: Sequence[Mapping[str, Any]],
    evaluator_targets: Sequence[Mapping[str, Any]],
    *,
    folds: int,
    per_source_fold: int,
    buffer_bp: int,
    namespace: str,
    seed: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    candidates = [dict(row, source_genomic_fold=int(genomic_assignment[str(row["chromosome"])])) for row in peaks]
    excluded_targets = {(row["chromosome"], int(row["bed_start_0based"]), int(row["bed_end_half_open"])) for row in evaluator_targets}
    selected: list[dict[str, Any]] = []
    census: list[dict[str, Any]] = []
    for held_fold in range(folds):
        matching_inputs = [row for row in inputs if int(row["held_genomic_fold"]) == held_fold]
        input_coordinates = {(row["chromosome"], int(row["bed_start_0based"]), int(row["bed_end_half_open"])) for row in matching_inputs}
        expanded = _expanded_intervals(matching_inputs, buffer_bp)
        for source_fold in range(folds):
            if source_fold == held_fold:
                continue
            eligible = [
                row
                for row in candidates
                if row["source_genomic_fold"] == source_fold
                and (row["chromosome"], row["bed_start_0based"], row["bed_end_half_open"]) not in input_coordinates
                and (row["chromosome"], row["bed_start_0based"], row["bed_end_half_open"]) not in excluded_targets
                and _outside_expanded(row, expanded)
            ]
            eligible.sort(key=lambda row: (sha256(f"{seed}\0{namespace}\0{held_fold}\0{row['peak_id']}".encode()).digest(), row["peak_id"]))
            if len(eligible) < per_source_fold:
                raise TrainingTargetMaskPlanError(f"insufficient buffered candidates for held={held_fold}, source={source_fold}: {len(eligible)}")
            census.append({"held_genomic_fold": held_fold, "source_genomic_fold": source_fold, "eligible_candidates": len(eligible), "selected_targets": per_source_fold})
            for row in eligible[:per_source_fold]:
                distance = _nearest_distance(row, matching_inputs)
                if distance != -1 and distance < buffer_bp:
                    raise TrainingTargetMaskPlanError("selected target violates observed-input buffer")
                selected.append({"training_target_hash": sha256("\0".join((namespace, str(held_fold), row["peak_id"])).encode()).hexdigest(), "chromosome": row["chromosome"], "bed_start_0based": row["bed_start_0based"], "bed_end_half_open": row["bed_end_half_open"], "source_genomic_fold": source_fold, "held_genomic_fold": held_fold, "nearest_observed_input_distance_bp": distance})
    return selected, census


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
    base = resolved["base_mask_plan"]
    genomic_rows = _tsv(base / "genomic_folds.tsv", ("chromosome", "genomic_fold", "eligible_peaks"))
    assignment = {row["chromosome"]: int(row["genomic_fold"]) for row in genomic_rows}
    inputs = _tsv(base / "observed_atac_input_peaks.tsv", ("input_hash", "chromosome", "bed_start_0based", "bed_end_half_open", "source_genomic_fold", "held_genomic_fold"))
    targets = _tsv(base / "target_peaks.tsv", ("target_hash", "chromosome", "bed_start_0based", "bed_end_half_open", "genomic_fold"))
    selection = config["selection"]
    selected, census = select_training_targets(read_peak_axes(resolved["input_h5"]), assignment, inputs, targets, folds=selection["genomic_folds"], per_source_fold=selection["targets_per_nonheld_source_genomic_fold"], buffer_bp=selection["minimum_distance_from_matching_observed_input_bp"], namespace=selection["identifier_hash_namespace"], seed=selection["selection_seed"])
    expected_total = selection["genomic_folds"] * selection["targets_per_held_genomic_fold"]
    if len(selected) != expected_total or len({row["training_target_hash"] for row in selected}) != expected_total:
        raise TrainingTargetMaskPlanError("training target census differs")
    output = reject_symlink_components(args.output, label="training target plan output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    fields = ("training_target_hash", "chromosome", "bed_start_0based", "bed_end_half_open", "source_genomic_fold", "held_genomic_fold", "nearest_observed_input_distance_bp")
    _write_tsv(output / "training_target_peaks.tsv", fields, selected)
    _write_tsv(output / "candidate_census.tsv", ("held_genomic_fold", "source_genomic_fold", "eligible_candidates", "selected_targets"), census)
    finite_distances = [row["nearest_observed_input_distance_bp"] for row in selected if row["nearest_observed_input_distance_bp"] >= 0]
    fold_counts = Counter(row["held_genomic_fold"] for row in selected)
    receipt = {"schema_version": "masld-bench-observed-multiome-training-target-mask-plan-receipt-v1", "plan_id": PLAN_ID, "dataset_id": "gse296875", "stage": "smoke", "training_targets": len(selected), "training_targets_per_held_genomic_fold": {str(fold): fold_counts[fold] for fold in range(selection["genomic_folds"])}, "targets_per_nonheld_source_genomic_fold": selection["targets_per_nonheld_source_genomic_fold"], "minimum_observed_input_distance_bp": selection["minimum_distance_from_matching_observed_input_bp"], "minimum_realized_finite_distance_bp": min(finite_distances) if finite_distances else None, "selection_uses_count_values": False, "rna_count_values_read": False, "atac_count_values_read": False, "raw_donor_ids_read": False, "development_metric_calculated": False, "model_fit": False, "sealed_outcomes_read": False, "next_gate": "supervised_training_label_materialization"}
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "tests/unit/test_observed_multiome_training_target_mask_plan.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_training_target_mask_plan", "plan_id": PLAN_ID, "training_targets": len(selected), "biological_count_values_read": False, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
