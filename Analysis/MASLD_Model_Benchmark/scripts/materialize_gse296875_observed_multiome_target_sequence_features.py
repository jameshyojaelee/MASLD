#!/usr/bin/env python3
"""Materialize outcome-free, target-generalizing sequence features."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
from typing import Any, BinaryIO, Mapping, Sequence

import h5py
import numpy as np

from masld_bench.artifacts import (
    freeze_tree,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
)
from masld_bench.observed_multiome_materialization import write_string_dataset


SCHEMA = "masld-bench-observed-multiome-target-sequence-features-v1"
CAMPAIGN_ID = "gse296875_observed_multiome_target_sequence_features_20260825"
ALPHABET = "ACGTN"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


class TargetSequenceFeatureError(ValueError):
    """Raised when target features drift from reference or outcome separation."""


def _digest(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise TargetSequenceFeatureError("JSON object required")
    return value


def _tree(root: Path, record: Mapping[str, Any], label: str) -> Path:
    lexical = reject_symlink_components(root / str(record.get("path", "")), label=label)
    path = lexical.resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise TargetSequenceFeatureError(f"{label} drifted")
    return path


def _tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise TargetSequenceFeatureError(f"{path.name} schema differs")
        return [dict(row) for row in reader]


def validate_config(root: Path, config: Mapping[str, Any]) -> dict[str, Path]:
    if config.get("schema_version") != SCHEMA or config.get("campaign_id") != CAMPAIGN_ID or config.get("dataset_id") != "gse296875" or config.get("stage") != "smoke":
        raise TargetSequenceFeatureError("campaign identity differs")
    feature = config.get("feature_contract")
    if feature != {"window_bp": 1024, "center": "floor_of_peak_midpoint", "boundary_padding": "N", "alphabet": "ACGTN", "canonical_reverse_complement_kmers": [1, 2, 3], "additional_features": ["mono_entropy_bits", "log1p_peak_width"], "feature_count": 85, "feature_dtype": "float32", "normalization": "kmer_frequency_within_window", "fit_on_outcomes": False, "target_specific_fitted_parameters": False}:
        raise TargetSequenceFeatureError("feature contract differs")
    if config.get("expected") != {"genomic_folds": 5, "roles": ["training", "held"], "children": 10, "training_targets_per_fold": 1000, "held_targets_per_fold": 1000}:
        raise TargetSequenceFeatureError("expected census differs")
    firewall = config.get("firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise TargetSequenceFeatureError("sequence feature firewall is open")
    parents = config.get("parents")
    if not isinstance(parents, dict) or set(parents) != {"training_target_mask_plan", "base_mask_plan", "sequence_reference"}:
        raise TargetSequenceFeatureError("parent roster differs")
    resolved = {key: _tree(root, record, key) for key, record in parents.items()}
    reference = parents["sequence_reference"]
    manifest = _json(resolved["sequence_reference"] / "ARTIFACTS.json")
    identities = {row["path"]: row["sha256"] for row in manifest["artifacts"]}
    contract_source = Path(str(reference.get("contract_source_fasta", ""))).resolve(strict=True)
    derivative = (root / str(reference.get("primary_analysis_derivative", ""))).resolve(strict=True)
    derivative.relative_to(root)
    if (
        reference.get("build") != "GRCh38.p14"
        or identities.get(reference.get("fasta")) != reference.get("fasta_sha256")
        or identities.get(reference.get("fai")) != reference.get("fai_sha256")
        or reference.get("contract_source_fasta_sha256") != "9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215"
        or _digest(contract_source) != reference.get("contract_source_fasta_sha256")
        or _digest(derivative) != reference.get("primary_analysis_derivative_sha256")
        or manifest.get("metadata", {}).get("source_fasta_sha256") != reference.get("primary_analysis_derivative_sha256")
        or reference.get("materialization") != "canonical_primary_contig_subset_then_bgzip_decompression_without_sequence_change"
    ):
        raise TargetSequenceFeatureError("sequence reference identity differs")
    try:
        for compressed in (contract_source, derivative):
            with gzip.open(compressed, "rb") as handle:
                while handle.read(8 * 1024 * 1024):
                    pass
    except (OSError, EOFError) as error:
        raise TargetSequenceFeatureError("contract source FASTA gzip integrity failed") from error
    resolved["fasta"] = resolved["sequence_reference"] / reference["fasta"]
    resolved["fai"] = resolved["sequence_reference"] / reference["fai"]
    resolved["contract_source_fasta"] = contract_source
    resolved["primary_analysis_derivative"] = derivative
    return resolved


class IndexedFasta:
    def __init__(self, fasta: Path, fai: Path) -> None:
        self.fasta = fasta
        self._handle: BinaryIO = fasta.open("rb")
        self.entries: dict[str, tuple[int, int, int, int]] = {}
        for line in fai.read_text(encoding="utf-8").splitlines():
            fields = line.split("\t")
            if len(fields) < 5:
                raise TargetSequenceFeatureError("FAI schema differs")
            self.entries[fields[0]] = tuple(map(int, fields[1:5]))
        if not self.entries:
            self._handle.close()
            raise TargetSequenceFeatureError("FAI is empty")

    def close(self) -> None:
        self._handle.close()

    def __enter__(self) -> "IndexedFasta":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def fetch(self, chromosome: str, start: int, end: int, *, pad: str = "N") -> str:
        if end <= start or len(pad) != 1:
            raise TargetSequenceFeatureError("FASTA interval differs")
        entry = self.entries.get(chromosome)
        if entry is None:
            raise TargetSequenceFeatureError(f"contig is absent: {chromosome}")
        length, offset, line_bases, line_width = entry
        left_pad = max(0, -start)
        right_pad = max(0, end - length)
        clipped_start = max(0, start)
        clipped_end = min(length, end)
        pieces: list[bytes] = []
        position = clipped_start
        while position < clipped_end:
            within_line = position % line_bases
            take = min(clipped_end - position, line_bases - within_line)
            byte_offset = offset + (position // line_bases) * line_width + within_line
            self._handle.seek(byte_offset)
            payload = self._handle.read(take)
            if len(payload) != take:
                raise TargetSequenceFeatureError("short FASTA read")
            pieces.append(payload)
            position += take
        sequence = pad * left_pad + b"".join(pieces).decode("ascii").upper() + pad * right_pad
        if len(sequence) != end - start:
            raise TargetSequenceFeatureError("FASTA window length differs")
        return "".join(base if base in ALPHABET else "N" for base in sequence)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def feature_names() -> list[str]:
    names = []
    for k in (1, 2, 3):
        canonical = sorted({min(kmer, reverse_complement(kmer)) for kmer in map("".join, product(ALPHABET, repeat=k))})
        names.extend(f"canonical_k{k}:{value}" for value in canonical)
    names.extend(["mono_entropy_bits", "log1p_peak_width"])
    if len(names) != 85 or len(set(names)) != len(names):
        raise TargetSequenceFeatureError("feature name census differs")
    return names


def sequence_features(sequence: str, peak_width: int) -> np.ndarray:
    if not sequence or peak_width < 1 or any(base not in ALPHABET for base in sequence):
        raise TargetSequenceFeatureError("sequence feature input differs")
    values: list[float] = []
    for k in (1, 2, 3):
        labels = sorted({min(kmer, reverse_complement(kmer)) for kmer in map("".join, product(ALPHABET, repeat=k))})
        counts = {label: 0 for label in labels}
        total = len(sequence) - k + 1
        for index in range(total):
            kmer = sequence[index : index + k]
            counts[min(kmer, reverse_complement(kmer))] += 1
        values.extend(counts[label] / total for label in labels)
    mono = [sequence.count(base) / len(sequence) for base in ALPHABET]
    entropy = -sum(value * math.log2(value) for value in mono if value > 0)
    values.extend([entropy, math.log1p(peak_width)])
    result = np.asarray(values, dtype=np.float32)
    if result.shape != (85,) or not np.isfinite(result).all():
        raise TargetSequenceFeatureError("sequence feature vector differs")
    return result


def _write_child(path: Path, role: str, fold: int, rows: Sequence[Mapping[str, str]], matrix: np.ndarray, names: Sequence[str], window_starts: Sequence[int], window_ends: Sequence[int]) -> None:
    identifier = "training_target_hash" if role == "training" else "target_hash"
    with h5py.File(path, "x") as handle:
        handle.attrs["schema_version"] = "masld-bench-observed-multiome-target-sequence-features-h5-v1"
        handle.attrs["dataset_id"] = "gse296875"
        handle.attrs["stage"] = "smoke"
        handle.attrs["role"] = role
        handle.attrs["held_genomic_fold"] = fold
        handle.attrs["reference_build"] = "GRCh38.p14"
        handle.attrs["outcome_free"] = True
        write_string_dataset(handle, identifier, [row[identifier] for row in rows])
        write_string_dataset(handle, "chromosome", [row["chromosome"] for row in rows])
        handle.create_dataset("bed_start_0based", data=np.asarray([int(row["bed_start_0based"]) for row in rows], dtype=np.int64))
        handle.create_dataset("bed_end_half_open", data=np.asarray([int(row["bed_end_half_open"]) for row in rows], dtype=np.int64))
        handle.create_dataset("window_start_0based", data=np.asarray(window_starts, dtype=np.int64))
        handle.create_dataset("window_end_half_open", data=np.asarray(window_ends, dtype=np.int64))
        write_string_dataset(handle, "feature_name", names)
        handle.create_dataset("features", data=matrix, compression="gzip", shuffle=True)


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
    training = _tsv(resolved["training_target_mask_plan"] / "training_target_peaks.tsv", ("training_target_hash", "chromosome", "bed_start_0based", "bed_end_half_open", "source_genomic_fold", "held_genomic_fold", "nearest_observed_input_distance_bp"))
    held = _tsv(resolved["base_mask_plan"] / "target_peaks.tsv", ("target_hash", "chromosome", "bed_start_0based", "bed_end_half_open", "genomic_fold"))
    names = feature_names()
    output = reject_symlink_components(args.output, label="target sequence feature output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    children: list[dict[str, Any]] = []
    window_bp = config["feature_contract"]["window_bp"]
    with IndexedFasta(resolved["fasta"], resolved["fai"]) as indexed:
        for role in config["expected"]["roles"]:
            source = training if role == "training" else held
            fold_field = "held_genomic_fold" if role == "training" else "genomic_fold"
            for fold in range(config["expected"]["genomic_folds"]):
                rows = [row for row in source if int(row[fold_field]) == fold]
                target_expected = config["expected"][f"{role}_targets_per_fold"]
                if len(rows) != target_expected:
                    raise TargetSequenceFeatureError("target role/fold census differs")
                vectors = []
                window_starts = []
                window_ends = []
                for row in rows:
                    start = int(row["bed_start_0based"])
                    end = int(row["bed_end_half_open"])
                    center = (start + end) // 2
                    window_start = center - window_bp // 2
                    window_end = window_start + window_bp
                    sequence = indexed.fetch(row["chromosome"], window_start, window_end)
                    vectors.append(sequence_features(sequence, end - start))
                    window_starts.append(window_start)
                    window_ends.append(window_end)
                matrix = np.vstack(vectors).astype(np.float32, copy=False)
                child = output / f"{role}_genomic_{fold}"
                child.mkdir()
                _write_child(child / "target_features.h5", role, fold, rows, matrix, names, window_starts, window_ends)
                receipt = {"schema_version": "masld-bench-observed-multiome-target-sequence-feature-receipt-v1", "dataset_id": "gse296875", "stage": "smoke", "role": role, "held_genomic_fold": fold, "targets": len(rows), "features": len(names), "matrix_shape": list(matrix.shape), "matrix_sha256": sha256(matrix.tobytes()).hexdigest(), "reference_build": "GRCh38.p14", "window_bp": window_bp, "reverse_complement_invariant": True, "biological_matrix_read": False, "training_label_artifact_read": False, "evaluator_artifact_read": False, "development_metric_calculated": False, "model_fit": False, "sealed_outcomes_read": False}
                write_json_exclusive(child / "receipt.json", receipt)
                child_hash = freeze_tree(child, metadata={"artifact_class": "gse296875_observed_multiome_target_sequence_features", "role": role, "held_genomic_fold": fold, "outcome_free": True, "sealed_outcomes_accessed": False})
                children.append({"role": role, "held_genomic_fold": fold, "path": child.relative_to(output).as_posix(), "artifacts_sha256": child_hash})
    parent = {"schema_version": "masld-bench-observed-multiome-target-sequence-feature-campaign-receipt-v1", "campaign_id": CAMPAIGN_ID, "dataset_id": "gse296875", "stage": "smoke", "children": children, "child_count": len(children), "feature_count": len(names), "feature_names_sha256": sha256("\n".join(names).encode()).hexdigest(), "reference_build": "GRCh38.p14", "biological_matrix_read": False, "training_label_artifact_read": False, "evaluator_artifact_read": False, "development_metric_calculated": False, "model_fit": False, "prediction_bundle_read": False, "sealed_outcomes_read": False, "next_gate": "prediction_only_broker_fixture"}
    write_json_exclusive(output / "receipt.json", parent)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "tests/unit/test_observed_multiome_target_sequence_features.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_target_sequence_feature_campaign", "campaign_id": CAMPAIGN_ID, "child_count": 10, "outcome_free": True, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
