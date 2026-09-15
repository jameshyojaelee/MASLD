#!/usr/bin/env python3
"""Independently rederive every outcome-free target sequence feature."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
from itertools import product
import json
import math
import os
from pathlib import Path
from typing import Any, Mapping, Sequence

import h5py
import numpy as np

from masld_bench.artifacts import freeze_tree, reject_symlink_components, verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-observed-multiome-target-sequence-feature-verification-v1"
CAMPAIGN_ID = "gse296875_observed_multiome_target_sequence_feature_verification_20260825"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


class TargetSequenceVerificationError(ValueError):
    """Raised when the independently derived target feature layer differs."""


def _digest(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise TargetSequenceVerificationError("JSON object required")
    return value


def _tree(root: Path, record: Mapping[str, Any], label: str) -> Path:
    lexical = reject_symlink_components(root / str(record.get("path", "")), label=label)
    path = lexical.resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise TargetSequenceVerificationError(f"{label} drifted")
    return path


def _tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise TargetSequenceVerificationError(f"{path.name} schema differs")
        return [dict(row) for row in reader]


def validate_config(root: Path, config: Mapping[str, Any]) -> dict[str, Path]:
    if config.get("schema_version") != SCHEMA or config.get("campaign_id") != CAMPAIGN_ID or config.get("dataset_id") != "gse296875" or config.get("stage") != "smoke":
        raise TargetSequenceVerificationError("verification identity differs")
    if config.get("expected") != {"roles": ["training", "held"], "genomic_folds": 5, "children": 10, "targets_per_child": 1000, "feature_count": 85, "window_bp": 1024, "recompute_fraction": 1.0, "absolute_tolerance": 1e-7}:
        raise TargetSequenceVerificationError("verification census differs")
    firewall = config.get("firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise TargetSequenceVerificationError("verification firewall is open")
    parents = config.get("parents")
    expected = {"target_sequence_features", "training_target_mask_plan", "base_mask_plan", "sequence_reference"}
    if not isinstance(parents, dict) or set(parents) != expected:
        raise TargetSequenceVerificationError("verification parent roster differs")
    resolved = {key: _tree(root, parents[key], key) for key in parents}
    reference = parents["sequence_reference"]
    manifest = _json(resolved["sequence_reference"] / "ARTIFACTS.json")
    identities = {row["path"]: row["sha256"] for row in manifest["artifacts"]}
    if identities.get(reference.get("fasta")) != reference.get("fasta_sha256") or identities.get(reference.get("fai")) != reference.get("fai_sha256"):
        raise TargetSequenceVerificationError("verification reference differs")
    resolved["fasta"] = resolved["sequence_reference"] / str(reference["fasta"])
    resolved["fai"] = resolved["sequence_reference"] / str(reference["fai"])
    return resolved


class PreadFasta:
    def __init__(self, fasta: Path, fai: Path) -> None:
        self.fd = os.open(fasta, os.O_RDONLY)
        self.entries: dict[str, tuple[int, int, int, int]] = {}
        for line in fai.read_text(encoding="utf-8").splitlines():
            fields = line.split("\t")
            if len(fields) < 5:
                self.close()
                raise TargetSequenceVerificationError("FAI schema differs")
            self.entries[fields[0]] = tuple(map(int, fields[1:5]))

    def close(self) -> None:
        if self.fd >= 0:
            os.close(self.fd)
            self.fd = -1

    def __enter__(self) -> "PreadFasta":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def fetch(self, chromosome: str, start: int, end: int) -> str:
        if chromosome not in self.entries or end <= start:
            raise TargetSequenceVerificationError("reference interval differs")
        length, offset, line_bases, line_width = self.entries[chromosome]
        left = max(0, -start)
        right = max(0, end - length)
        position = max(0, start)
        stop = min(length, end)
        pieces: list[bytes] = []
        while position < stop:
            within = position % line_bases
            take = min(stop - position, line_bases - within)
            payload = os.pread(self.fd, take, offset + position // line_bases * line_width + within)
            if len(payload) != take:
                raise TargetSequenceVerificationError("short independent FASTA read")
            pieces.append(payload)
            position += take
        sequence = "N" * left + b"".join(pieces).decode("ascii").upper() + "N" * right
        if len(sequence) != end - start or any(base not in "ACGTN" for base in sequence):
            raise TargetSequenceVerificationError("independent FASTA window differs")
        return sequence


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def _names() -> list[str]:
    values: list[str] = []
    for k in (1, 2, 3):
        labels = sorted({min(word, _reverse_complement(word)) for word in map("".join, product("ACGTN", repeat=k))})
        values.extend(f"canonical_k{k}:{label}" for label in labels)
    values.extend(["mono_entropy_bits", "log1p_peak_width"])
    return values


def _features(sequence: str, width: int) -> np.ndarray:
    values: list[float] = []
    for k in (1, 2, 3):
        labels = sorted({min(word, _reverse_complement(word)) for word in map("".join, product("ACGTN", repeat=k))})
        counts = dict.fromkeys(labels, 0)
        for index in range(len(sequence) - k + 1):
            word = sequence[index:index + k]
            counts[min(word, _reverse_complement(word))] += 1
        denominator = len(sequence) - k + 1
        values.extend(counts[label] / denominator for label in labels)
    fractions = [sequence.count(base) / len(sequence) for base in "ACGTN"]
    values.append(-sum(value * math.log2(value) for value in fractions if value))
    values.append(math.log1p(width))
    return np.asarray(values, dtype=np.float32)


def _strings(dataset: h5py.Dataset) -> list[str]:
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in dataset[:]]


def run(root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    resolved = validate_config(root, config)
    receipt = _json(resolved["target_sequence_features"] / "receipt.json")
    if receipt.get("child_count") != 10 or receipt.get("biological_matrix_read") is not False:
        raise TargetSequenceVerificationError("target feature parent receipt differs")
    training = _tsv(resolved["training_target_mask_plan"] / "training_target_peaks.tsv", ("training_target_hash", "chromosome", "bed_start_0based", "bed_end_half_open", "source_genomic_fold", "held_genomic_fold", "nearest_observed_input_distance_bp"))
    held = _tsv(resolved["base_mask_plan"] / "target_peaks.tsv", ("target_hash", "chromosome", "bed_start_0based", "bed_end_half_open", "genomic_fold"))
    expected_names = _names()
    checked = 0
    children = []
    child_index = {(row["role"], int(row["held_genomic_fold"])): row for row in receipt["children"]}
    with PreadFasta(resolved["fasta"], resolved["fai"]) as reference:
        for role in config["expected"]["roles"]:
            source = training if role == "training" else held
            fold_field = "held_genomic_fold" if role == "training" else "genomic_fold"
            identifier = "training_target_hash" if role == "training" else "target_hash"
            for fold in range(config["expected"]["genomic_folds"]):
                expected_rows = [row for row in source if int(row[fold_field]) == fold]
                record = child_index.get((role, fold))
                if len(expected_rows) != 1000 or record is None:
                    raise TargetSequenceVerificationError("target feature child census differs")
                child = (resolved["target_sequence_features"] / record["path"]).resolve(strict=True)
                child.relative_to(resolved["target_sequence_features"])
                verify_frozen_tree(child)
                if _digest(child / "ARTIFACTS.json") != record["artifacts_sha256"]:
                    raise TargetSequenceVerificationError("target feature child drifted")
                with h5py.File(child / "target_features.h5", "r") as handle:
                    expected_keys = {identifier, "chromosome", "bed_start_0based", "bed_end_half_open", "window_start_0based", "window_end_half_open", "feature_name", "features"}
                    if set(handle.keys()) != expected_keys or handle["features"].shape != (1000, 85):
                        raise TargetSequenceVerificationError("target feature H5 surface differs")
                    if _strings(handle[identifier]) != [row[identifier] for row in expected_rows] or _strings(handle["chromosome"]) != [row["chromosome"] for row in expected_rows] or _strings(handle["feature_name"]) != expected_names:
                        raise TargetSequenceVerificationError("target feature identifiers differ")
                    matrix = np.asarray(handle["features"][:], dtype=np.float32)
                    starts = np.asarray(handle["window_start_0based"][:], dtype=np.int64)
                    ends = np.asarray(handle["window_end_half_open"][:], dtype=np.int64)
                derived = []
                for index, row in enumerate(expected_rows):
                    start = int(row["bed_start_0based"])
                    end = int(row["bed_end_half_open"])
                    center = (start + end) // 2
                    window_start = center - 512
                    if starts[index] != window_start or ends[index] != window_start + 1024:
                        raise TargetSequenceVerificationError("target feature coordinate differs")
                    sequence = reference.fetch(row["chromosome"], window_start, window_start + 1024)
                    forward = _features(sequence, end - start)
                    reverse = _features(_reverse_complement(sequence), end - start)
                    if not np.allclose(forward, reverse, rtol=0, atol=1e-7):
                        raise TargetSequenceVerificationError("independent reverse-complement check failed")
                    derived.append(forward)
                expected_matrix = np.vstack(derived)
                if not np.allclose(matrix, expected_matrix, rtol=0, atol=float(config["expected"]["absolute_tolerance"])):
                    raise TargetSequenceVerificationError("independent target feature recomputation differs")
                child_receipt = _json(child / "receipt.json")
                if sha256(matrix.tobytes()).hexdigest() != child_receipt.get("matrix_sha256"):
                    raise TargetSequenceVerificationError("target feature matrix receipt differs")
                checked += matrix.shape[0]
                children.append({"role": role, "held_genomic_fold": fold, "targets_verified": matrix.shape[0], "features_verified": matrix.shape[1], "maximum_absolute_difference": float(np.max(np.abs(matrix - expected_matrix)))})
    return {"schema_version": "masld-bench-observed-multiome-target-sequence-feature-verification-receipt-v1", "campaign_id": CAMPAIGN_ID, "dataset_id": "gse296875", "stage": "smoke", "children_verified": len(children), "targets_verified": checked, "features_per_target": 85, "recompute_fraction": 1.0, "children": children, "reference_build": "GRCh38.p14", "biological_matrix_read": False, "training_label_artifact_read": False, "evaluator_artifact_read": False, "development_metric_calculated": False, "model_fit": False, "prediction_bundle_read": False, "sealed_outcomes_read": False, "promotion_gate_passed": True, "next_gate": "prediction_only_broker_fixture"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    config_path.relative_to(root)
    output = reject_symlink_components(args.output, label="target feature verification output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    receipt = run(root, _json(config_path))
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "tests/unit/test_verify_observed_multiome_target_sequence_features.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_target_sequence_feature_verification", "campaign_id": CAMPAIGN_ID, "targets_verified": receipt["targets_verified"], "promotion_gate_passed": True, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
