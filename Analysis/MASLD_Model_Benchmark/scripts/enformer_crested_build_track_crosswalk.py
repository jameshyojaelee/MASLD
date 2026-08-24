#!/usr/bin/env python3
"""Bind CREsted Enformer output positions to the canonical 5,313-track manifest."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


LABEL_SHA256 = "83c67554af8900075ee542b1cbf247bf9905fc1b1bc59d5cff11ce64ef58dffb"
TARGET_SHA256 = "d90233175a9fec389e5ac04cc4aac8cfc08c4aa86217e2e5c10221de78023580"


class TrackCrosswalkError(ValueError):
    """Raised when the converted and canonical output orders differ."""


def _digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _read_targets(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "index",
            "genome",
            "identifier",
            "file",
            "clip",
            "scale",
            "sum_stat",
            "description",
        }
        if set(reader.fieldnames or []) != required:
            raise TrackCrosswalkError("canonical target fields differ")
        rows = [dict(row) for row in reader]
    if len(rows) != 5313:
        raise TrackCrosswalkError("canonical target row count differs")
    if [row["index"] for row in rows] != [str(index) for index in range(5313)]:
        raise TrackCrosswalkError("canonical target index order differs")
    if any(row["genome"] != "0" for row in rows):
        raise TrackCrosswalkError("canonical target genome code differs")
    return rows


def _assay_family(description: str) -> str:
    prefix = description.partition(":")[0].upper()
    if prefix in {"DNASE", "ATAC"}:
        return "accessibility"
    if prefix == "CAGE":
        return "cage"
    if prefix == "CHIP":
        return "chip"
    return prefix.lower() or "unclassified"


def build(labels: Path, targets: Path, output: Path) -> dict[str, Any]:
    if output.exists() or any(path.is_symlink() for path in (labels, targets)):
        raise TrackCrosswalkError("input or output contract is invalid")
    if _digest(labels) != LABEL_SHA256 or _digest(targets) != TARGET_SHA256:
        raise TrackCrosswalkError("input checksum differs")
    label_rows = labels.read_text(encoding="utf-8").splitlines()
    target_rows = _read_targets(targets)
    if label_rows != [row["description"] for row in target_rows]:
        mismatch = next(
            (
                index
                for index, (label, row) in enumerate(zip(label_rows, target_rows))
                if label != row["description"]
            ),
            min(len(label_rows), len(target_rows)),
        )
        raise TrackCrosswalkError(f"converted output order differs at position {mismatch}")

    output.mkdir(parents=True, mode=0o750)
    crosswalk_path = output / "track_crosswalk.tsv"
    fields = (
        "position",
        "canonical_index",
        "identifier",
        "description",
        "assay_family",
        "clip",
        "scale",
        "sum_stat",
    )
    assay_counts: Counter[str] = Counter()
    with crosswalk_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for position, row in enumerate(target_rows):
            assay = _assay_family(row["description"])
            assay_counts[assay] += 1
            writer.writerow(
                {
                    "position": position,
                    "canonical_index": row["index"],
                    "identifier": row["identifier"],
                    "description": row["description"],
                    "assay_family": assay,
                    "clip": row["clip"],
                    "scale": row["scale"],
                    "sum_stat": row["sum_stat"],
                }
            )

    descriptions = Counter(label_rows)
    receipt = {
        "schema_version": "masld-bench-enformer-crested-track-crosswalk-v1",
        "status": "pass",
        "model_id": "enformer_crested_restricted_port",
        "converted_output_rows": len(label_rows),
        "canonical_target_rows": len(target_rows),
        "position_by_position_description_identity": True,
        "canonical_target_identifier_bound_by_position": True,
        "description_only_track_selection_allowed": False,
        "unique_descriptions": len(descriptions),
        "duplicated_description_values": sum(
            count > 1 for count in descriptions.values()
        ),
        "maximum_description_multiplicity": max(descriptions.values()),
        "assay_family_counts": dict(sorted(assay_counts.items())),
        "track_crosswalk_sha256": _digest(crosswalk_path),
        "native_prediction_contract": {
            "sequence_input_bp": 196608,
            "output_bins": 896,
            "output_resolution_bp": 128,
            "output_tracks": 5313,
            "allele_delta_sign": "ALT_minus_REF",
            "native_SAD": "sum_over_output_bins_of_ALT_minus_REF_per_track",
            "native_SAR": "log2(sum_ALT_plus_1)-log2(sum_REF_plus_1)_per_track",
            "reverse_complement": "separately_named_adaptation_not_silent_native_inference",
        },
        "allowed_task_roles": [
            "sequence_native_regulatory_restricted_static_track_comparator",
            "variant_to_regulation_restricted_sequence_component",
        ],
        "forbidden_task_roles": [
            "RNA_conditioned_ATAC",
            "observed_multiome",
            "standalone_signed_eQTL_or_ieQTL_predictor",
            "donor_or_MASLD_context_model",
        ],
        "mandatory_baselines": [
            "zero_allele_delta",
            "mean_track",
            "gc_kmer_control",
            "deltaSVM",
            "gkm_SVM",
            "sequence_CNN_control",
            "sequence_transformer_control",
            "assay_native_linear_head",
        ],
        "development_head_rule": "track_selection_and_any_fitted_combiner_are_fit_inside_donor_and_genomic_or_LD_training_folds_only",
        "native_parity_established": False,
        "project_data_read": False,
        "outcomes_read": False,
        "sealed_features_or_labels_read": False,
        "open_champion_eligible": False,
        "terminal_disposition": "positional_track_crosswalk_passed_restricted_runtime_restore_and_native_parity_pending",
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    result = build(arguments.labels, arguments.targets, arguments.output)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
