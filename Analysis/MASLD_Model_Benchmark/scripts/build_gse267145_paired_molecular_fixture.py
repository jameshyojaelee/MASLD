#!/usr/bin/env python3
"""Materialize the outcome-free paired GSE267145 RNA and H3K27ac fixture."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import shlex
from typing import Any, Mapping, Sequence

import numpy as np


EXPECTED_PARTICIPANTS = 99
EXPECTED_RNA_SOURCE_FEATURES = 43_285
EXPECTED_RNA_ALLOWED_FEATURES = 42_163
EXPECTED_H3_FEATURES = 96_460
EXPECTED_FOLD_COUNTS = {0: 21, 1: 21, 2: 21, 3: 19, 4: 17}
PAIRING = "same_sample_different_aliquot"
RNA_JOIN_STATE = "observed_unit_semantics_unresolved"
H3_JOIN_STATE = "observed_integer_counts"
JOIN_FIELDS_USED = (
    "participant_id",
    "rna_gsm",
    "h3k27ac_gsm",
    "pairing",
    "outer_fold",
    "rna_measurement_state",
    "h3k27ac_measurement_state",
)


class PairedMolecularFixtureError(RuntimeError):
    """Raised when an input differs from the outcome-free fixture requirement."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_artifact(root: Path, expected_sha256: str) -> None:
    observed = sha256_file(root / "ARTIFACTS.json")
    if observed != expected_sha256:
        raise PairedMolecularFixtureError(
            f"input ARTIFACTS SHA-256 differs for {root}: {observed}"
        )


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise PairedMolecularFixtureError("refusing to write an empty table")
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=tuple(rows[0]),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_participant_roster(
    path: Path,
    *,
    expected_participants: int = EXPECTED_PARTICIPANTS,
    expected_fold_counts: Mapping[int, int] | None = EXPECTED_FOLD_COUNTS,
) -> list[dict[str, Any]]:
    """Read only identity, assay, pairing, and already-frozen fold columns."""

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise PairedMolecularFixtureError("participant join is empty") from error
        if len(header) != len(set(header)) or any(
            field not in header for field in JOIN_FIELDS_USED
        ):
            raise PairedMolecularFixtureError("participant join fields differ")
        index = {field: header.index(field) for field in JOIN_FIELDS_USED}
        rows: list[dict[str, Any]] = []
        for line_number, values in enumerate(reader, start=2):
            if len(values) != len(header):
                raise PairedMolecularFixtureError(
                    f"participant join width differs at line {line_number}"
                )
            selected = {field: values[index[field]] for field in JOIN_FIELDS_USED}
            try:
                fold = int(selected["outer_fold"])
            except ValueError as error:
                raise PairedMolecularFixtureError("outer fold is not an integer") from error
            if (
                selected["pairing"] != PAIRING
                or selected["rna_measurement_state"] != RNA_JOIN_STATE
                or selected["h3k27ac_measurement_state"] != H3_JOIN_STATE
                or fold not in range(5)
            ):
                raise PairedMolecularFixtureError(
                    "pairing, observation state, or outer fold differs"
                )
            rows.append(
                {
                    "participant_id": selected["participant_id"],
                    "rna_gsm": selected["rna_gsm"],
                    "h3k27ac_gsm": selected["h3k27ac_gsm"],
                    "pairing": selected["pairing"],
                    "outer_fold": fold,
                }
            )
    if len(rows) != expected_participants:
        raise PairedMolecularFixtureError("participant count differs")
    for field in ("participant_id", "rna_gsm", "h3k27ac_gsm"):
        values = [str(row[field]) for row in rows]
        if any(not value for value in values) or len(set(values)) != len(values):
            raise PairedMolecularFixtureError(f"{field} is empty or duplicated")
    if expected_fold_counts is not None:
        observed = Counter(int(row["outer_fold"]) for row in rows)
        if dict(sorted(observed.items())) != dict(sorted(expected_fold_counts.items())):
            raise PairedMolecularFixtureError("frozen outer-fold census differs")
    return rows


def read_rna_crosswalk(
    path: Path,
    *,
    expected_source_features: int = EXPECTED_RNA_SOURCE_FEATURES,
    expected_allowed_features: int = EXPECTED_RNA_ALLOWED_FEATURES,
) -> list[dict[str, Any]]:
    used = (
        "matrix_gene_id",
        "stable_gene_id",
        "mapping_state",
        "allowed_project_input",
    )
    rows: list[dict[str, Any]] = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise PairedMolecularFixtureError("RNA crosswalk is empty") from error
        if len(header) != len(set(header)) or any(field not in header for field in used):
            raise PairedMolecularFixtureError("RNA crosswalk fields differ")
        index = {field: header.index(field) for field in used}
        for source_index, values in enumerate(reader):
            if len(values) != len(header):
                raise PairedMolecularFixtureError("RNA crosswalk width differs")
            selected = {field: values[index[field]] for field in used}
            allowed = selected["allowed_project_input"] == "true"
            if selected["allowed_project_input"] not in {"true", "false"}:
                raise PairedMolecularFixtureError("RNA crosswalk admission flag differs")
            if allowed and selected["mapping_state"] != "stable_id_exact_v98_and_v49":
                raise PairedMolecularFixtureError("an admitted RNA gene is not exact")
            rows.append(
                {
                    "source_feature_index": source_index,
                    "matrix_gene_id": selected["matrix_gene_id"],
                    "stable_gene_id": selected["stable_gene_id"],
                    "allowed": allowed,
                }
            )
    if len(rows) != expected_source_features:
        raise PairedMolecularFixtureError("RNA source feature count differs")
    if sum(bool(row["allowed"]) for row in rows) != expected_allowed_features:
        raise PairedMolecularFixtureError("RNA admitted feature count differs")
    for field in ("matrix_gene_id", "stable_gene_id"):
        values = [str(row[field]) for row in rows]
        if len(values) != len(set(values)) or any(not value for value in values):
            raise PairedMolecularFixtureError(f"RNA {field} is empty or duplicated")
    return rows


def materialize_rna(
    path: Path,
    roster: Sequence[Mapping[str, Any]],
    crosswalk: Sequence[Mapping[str, Any]],
    output: Path,
) -> dict[str, Any]:
    participants = [str(row["participant_id"]) for row in roster]
    allowed = [row for row in crosswalk if bool(row["allowed"])]
    matrix = np.empty((len(participants), len(allowed)), dtype=np.float64)
    fractional_values = 0
    minimum = math.inf
    maximum = -math.inf
    allowed_index = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise PairedMolecularFixtureError("RNA matrix is empty") from error
        if not header or header[0] != "ensembl_gene_id" or len(header) != len(set(header)):
            raise PairedMolecularFixtureError("RNA matrix header differs")
        try:
            source_columns = [header.index(participant) for participant in participants]
        except ValueError as error:
            raise PairedMolecularFixtureError("RNA participant axis is incomplete") from error
        if [header[column] for column in source_columns] != participants:
            raise PairedMolecularFixtureError("RNA participant alignment differs")
        observed_rows = 0
        for source_index, values in enumerate(reader):
            if source_index >= len(crosswalk):
                raise PairedMolecularFixtureError("RNA matrix has extra features")
            if len(values) != len(header):
                raise PairedMolecularFixtureError("RNA matrix width differs")
            mapping = crosswalk[source_index]
            if values[0] != mapping["matrix_gene_id"]:
                raise PairedMolecularFixtureError("RNA matrix and crosswalk order differ")
            if bool(mapping["allowed"]):
                try:
                    parsed = np.fromiter(
                        (float(values[column]) for column in source_columns),
                        dtype=np.float64,
                        count=len(source_columns),
                    )
                except ValueError as error:
                    raise PairedMolecularFixtureError("RNA value is not numeric") from error
                if not np.isfinite(parsed).all() or np.any(parsed < 0):
                    raise PairedMolecularFixtureError("RNA values are not finite nonnegative")
                fractional_values += int(np.count_nonzero(parsed != np.floor(parsed)))
                minimum = min(minimum, float(parsed.min()))
                maximum = max(maximum, float(parsed.max()))
                matrix[:, allowed_index] = parsed
                allowed_index += 1
            observed_rows += 1
    if observed_rows != len(crosswalk) or allowed_index != len(allowed):
        raise PairedMolecularFixtureError("RNA matrix feature census differs")
    if fractional_values == 0:
        raise PairedMolecularFixtureError(
            "RNA participant matrix does not retain fractional measurements"
        )
    np.save(output / "rna_values.npy", matrix, allow_pickle=False)
    axis = [
        {
            "rna_feature_index": index,
            "source_feature_index": row["source_feature_index"],
            "matrix_gene_id": row["matrix_gene_id"],
            "stable_gene_id": row["stable_gene_id"],
        }
        for index, row in enumerate(allowed)
    ]
    write_tsv(output / "rna_feature_axis.tsv", axis)
    return {
        "shape": [len(participants), len(allowed)],
        "dtype": str(matrix.dtype),
        "minimum": minimum,
        "maximum": maximum,
        "fractional_values": fractional_values,
        "source_features": len(crosswalk),
        "admitted_features": len(allowed),
        "masked_by_frozen_crosswalk": len(crosswalk) - len(allowed),
        "normalization_applied": False,
        "variance_filter_applied": False,
    }


def _participant_for_h3_column(label: str, participants: Sequence[str]) -> str:
    matches = [participant for participant in participants if label.startswith(participant + "_")]
    if len(matches) != 1:
        raise PairedMolecularFixtureError("H3 source column does not map uniquely")
    return matches[0]


def materialize_h3(
    path: Path,
    roster: Sequence[Mapping[str, Any]],
    output: Path,
    *,
    expected_features: int = EXPECTED_H3_FEATURES,
) -> dict[str, Any]:
    participants = [str(row["participant_id"]) for row in roster]
    matrix = np.empty((len(participants), expected_features), dtype=np.uint32)
    feature_keys: list[str] = []
    feature_key_set: set[str] = set()
    minimum: int | None = None
    maximum: int | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        header_line = handle.readline()
        if not header_line:
            raise PairedMolecularFixtureError("H3 matrix is empty")
        try:
            labels = shlex.split(header_line, comments=False, posix=True)
        except ValueError as error:
            raise PairedMolecularFixtureError("H3 header quoting differs") from error
        mapped = [_participant_for_h3_column(label, participants) for label in labels]
        if len(mapped) != len(participants) or len(set(mapped)) != len(participants):
            raise PairedMolecularFixtureError("H3 participant axis differs")
        source_by_participant = {participant: index for index, participant in enumerate(mapped)}
        source_columns = [source_by_participant[participant] + 1 for participant in participants]
        for feature_index, line in enumerate(handle):
            if feature_index >= expected_features:
                raise PairedMolecularFixtureError("H3 matrix has extra features")
            values = line.rstrip("\r\n").split()
            if len(values) != len(labels) + 1:
                raise PairedMolecularFixtureError("H3 matrix width differs")
            token = values[0]
            if len(token) < 3 or token[0] != '"' or token[-1] != '"':
                raise PairedMolecularFixtureError("H3 opaque feature quoting differs")
            feature_key = token[1:-1]
            if not feature_key or feature_key in feature_key_set:
                raise PairedMolecularFixtureError("H3 opaque feature key is empty or duplicated")
            selected = [values[column] for column in source_columns]
            if any(re.fullmatch(r"[0-9]+", value) is None for value in selected):
                raise PairedMolecularFixtureError("H3 values are not nonnegative integers")
            parsed64 = np.fromiter(
                (int(value) for value in selected), dtype=np.uint64, count=len(selected)
            )
            if np.any(parsed64 > np.iinfo(np.uint32).max):
                raise PairedMolecularFixtureError("H3 value exceeds uint32")
            parsed = parsed64.astype(np.uint32)
            local_minimum = int(parsed.min())
            local_maximum = int(parsed.max())
            minimum = local_minimum if minimum is None else min(minimum, local_minimum)
            maximum = local_maximum if maximum is None else max(maximum, local_maximum)
            matrix[:, feature_index] = parsed
            feature_keys.append(feature_key)
            feature_key_set.add(feature_key)
    if len(feature_keys) != expected_features:
        raise PairedMolecularFixtureError("H3 matrix feature census differs")
    np.save(output / "h3k27ac_counts.npy", matrix, allow_pickle=False)
    write_tsv(
        output / "h3k27ac_feature_axis.tsv",
        [
            {"h3k27ac_feature_index": index, "opaque_source_feature_key": value}
            for index, value in enumerate(feature_keys)
        ],
    )
    return {
        "shape": [len(participants), len(feature_keys)],
        "dtype": str(matrix.dtype),
        "minimum": minimum,
        "maximum": maximum,
        "integer_semantics_verified": True,
        "opaque_feature_keys_preserved": True,
        "coordinate_semantics_inferred": False,
        "normalization_applied": False,
        "variance_filter_applied": False,
    }


def build_fixture(
    *,
    matrix_root: Path,
    join_path: Path,
    crosswalk_path: Path,
    output: Path,
    expected_participants: int = EXPECTED_PARTICIPANTS,
    expected_rna_source_features: int = EXPECTED_RNA_SOURCE_FEATURES,
    expected_rna_allowed_features: int = EXPECTED_RNA_ALLOWED_FEATURES,
    expected_h3_features: int = EXPECTED_H3_FEATURES,
    expected_fold_counts: Mapping[int, int] | None = EXPECTED_FOLD_COUNTS,
) -> dict[str, Any]:
    if output.exists():
        raise PairedMolecularFixtureError(f"refusing to overwrite output: {output}")
    molecular = output / "molecular"
    folds = output / "folds"
    molecular.mkdir(parents=True)
    folds.mkdir()
    roster = read_participant_roster(
        join_path,
        expected_participants=expected_participants,
        expected_fold_counts=expected_fold_counts,
    )
    crosswalk = read_rna_crosswalk(
        crosswalk_path,
        expected_source_features=expected_rna_source_features,
        expected_allowed_features=expected_rna_allowed_features,
    )
    rna = materialize_rna(
        matrix_root / "raw" / "GSE269412_RNA.txt.gz", roster, crosswalk, molecular
    )
    h3 = materialize_h3(
        matrix_root / "raw" / "GSE267119_H3K27ac.txt.gz",
        roster,
        molecular,
        expected_features=expected_h3_features,
    )
    participant_axis = [
        {
            "participant_index": index,
            "participant_id": row["participant_id"],
            "rna_source_sample_accession": row["rna_gsm"],
            "h3k27ac_source_sample_accession": row["h3k27ac_gsm"],
            "pairing": row["pairing"],
            "rna_observation_state": "observed",
            "h3k27ac_observation_state": "observed",
        }
        for index, row in enumerate(roster)
    ]
    write_tsv(molecular / "participant_axis.tsv", participant_axis)
    np.save(
        molecular / "rna_observed_mask.npy",
        np.ones(expected_participants, dtype=np.bool_),
        allow_pickle=False,
    )
    np.save(
        molecular / "h3k27ac_observed_mask.npy",
        np.ones(expected_participants, dtype=np.bool_),
        allow_pickle=False,
    )
    fold_rows = [
        {"participant_id": row["participant_id"], "outer_fold": row["outer_fold"]}
        for row in roster
    ]
    write_tsv(folds / "participant_outer_folds.tsv", fold_rows)
    fold_counts = dict(sorted(Counter(row["outer_fold"] for row in roster).items()))
    fold_receipt = {
        "schema_version": "masld-bench-gse267145-molecular-folds-v1",
        "status": "passed",
        "participants": len(roster),
        "biological_unit": "participant",
        "outer_fold_counts": {str(key): value for key, value in fold_counts.items()},
        "folds_copied_from_authoritative_join": True,
        "folds_recomputed": False,
        "outcome_values_included": False,
    }
    write_json(folds / "receipt.json", fold_receipt)
    receipt = {
        "schema_version": "masld-bench-gse267145-paired-molecular-fixture-v1",
        "status": "passed",
        "cohort_family_id": "gse267145_znf469_human_liver",
        "participants": len(roster),
        "biological_unit": "participant",
        "pairing": PAIRING,
        "participant_alignment": "authoritative_join_order",
        "source_sample_accessions_preserved": True,
        "rna": rna,
        "h3k27ac": h3,
        "missingness": {
            "rna_observed": len(roster),
            "h3k27ac_observed": len(roster),
            "both_observed": len(roster),
            "explicit_boolean_masks": True,
            "missing_encoded_as_zero": False,
        },
        "feature_firewall": {
            "outcome_values_accessed": False,
            "outcome_columns_emitted": False,
            "outcome_bearing_source_labels_emitted": False,
            "outcomes_used_for_feature_construction": False,
            "published_outcome_selected_features_used": False,
            "global_normalization_used": False,
            "global_variance_filter_used": False,
        },
        "model_fitted": False,
        "metrics_calculated": False,
    }
    write_json(molecular / "receipt.json", receipt)
    return receipt


def validate_upstream_contracts(arguments: argparse.Namespace) -> None:
    for root, expected in (
        (arguments.activation_campaign, arguments.activation_campaign_artifacts_sha256),
        (arguments.activation_view, arguments.activation_view_artifacts_sha256),
        (arguments.matrix_audit, arguments.matrix_artifacts_sha256),
        (arguments.join, arguments.join_artifacts_sha256),
        (arguments.measurement, arguments.measurement_artifacts_sha256),
        (arguments.reference, arguments.reference_artifacts_sha256),
        (arguments.qc_rights, arguments.qc_rights_artifacts_sha256),
        (arguments.coordinate_activation, arguments.coordinate_artifacts_sha256),
    ):
        verify_artifact(root, expected)
    measurement = json.loads(
        (arguments.measurement / "measurement_evidence.json").read_text(encoding="utf-8")
    )
    reference = json.loads(
        (arguments.reference / "crosswalk_summary.json").read_text(encoding="utf-8")
    )
    coordinate = json.loads(
        (arguments.coordinate_activation / "activation_contract.json").read_text(
            encoding="utf-8"
        )
    )
    activation = json.loads(
        (arguments.activation_view / "activation_contract.json").read_text(
            encoding="utf-8"
        )
    )
    if (
        measurement.get("not_raw_integer_molecule_counts") is not True
        or measurement.get("deposited_measurement")
        != "nonnegative_fractional_gene_expression_estimates"
        or reference.get("rna_matrix_genes") != EXPECTED_RNA_SOURCE_FEATURES
        or reference.get("gene_mapping_states", {}).get("stable_id_exact_v98_and_v49")
        != EXPECTED_RNA_ALLOWED_FEATURES
        or reference.get("h3k27ac", {}).get("regions") != EXPECTED_H3_FEATURES
        or coordinate.get("coordinate_semantics_resolved") is not False
        or coordinate.get("sequence_extraction_allowed") is not False
        or activation.get("task_id") != "paired_bulk_histology_state"
        or activation.get("participants") != EXPECTED_PARTICIPANTS
        or activation.get("feature_contract", {}).get("complete_deposited_assay_axes_only")
        is not True
        or activation.get("feature_contract", {}).get(
            "published_outcome_selected_features_allowed"
        )
        is not False
    ):
        raise PairedMolecularFixtureError("upstream molecular boundary differs")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        "activation-campaign",
        "activation-view",
        "matrix-audit",
        "join",
        "measurement",
        "reference",
        "qc-rights",
        "coordinate-activation",
    ):
        parser.add_argument(f"--{name}", required=True, type=Path)
    for name in (
        "activation-campaign-artifacts-sha256",
        "activation-view-artifacts-sha256",
        "matrix-artifacts-sha256",
        "join-artifacts-sha256",
        "measurement-artifacts-sha256",
        "reference-artifacts-sha256",
        "qc-rights-artifacts-sha256",
        "coordinate-artifacts-sha256",
    ):
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    arguments = parse_args()
    validate_upstream_contracts(arguments)
    receipt = build_fixture(
        matrix_root=arguments.matrix_audit,
        join_path=arguments.join / "participant_join.tsv",
        crosswalk_path=arguments.reference / "rna_gene_crosswalk.tsv",
        output=arguments.output,
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
