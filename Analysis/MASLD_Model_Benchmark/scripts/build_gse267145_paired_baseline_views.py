#!/usr/bin/env python3
"""Build fold-scoped model views and an evaluator-only GSE267145 H3 view."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
from typing import Any, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file


EXPECTED_PARTICIPANTS = 99
EXPECTED_RNA_FEATURES = 42_163
EXPECTED_H3_FEATURES = 96_460
EXPECTED_FOLD_COUNTS = {0: 21, 1: 21, 2: 21, 3: 19, 4: 17}
PAIRING = "same_sample_different_aliquot"
FORBIDDEN_FIELDS = {
    "stage",
    "stage3",
    "stage5",
    "fibrosis",
    "steatosis",
    "ballooning",
    "lobular_inflammation",
    "lobular_necrosis",
    "sex",
    "age",
    "nas",
    "nash_crn_component_sum",
}


class PairedBaselineViewError(ValueError):
    """Raised when model/evaluator separation cannot be proven."""


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise PairedBaselineViewError(f"TSV lacks a header: {path}")
        return list(reader.fieldnames), list(reader)


def write_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def build_views(
    *,
    molecular: Path,
    folds: Path,
    output: Path,
    task_spec_sha256: str,
    promotion_gate_sha256: str,
    expected_participants: int = EXPECTED_PARTICIPANTS,
    expected_rna_features: int = EXPECTED_RNA_FEATURES,
    expected_h3_features: int = EXPECTED_H3_FEATURES,
    expected_fold_counts: Mapping[int, int] = EXPECTED_FOLD_COUNTS,
) -> dict[str, Any]:
    if output.exists():
        raise PairedBaselineViewError("refusing to overwrite paired baseline views")
    participant_fields, participant_rows = read_tsv(molecular / "participant_axis.tsv")
    fold_fields, fold_rows = read_tsv(folds / "participant_outer_folds.tsv")
    rna_axis_fields, rna_axis = read_tsv(molecular / "rna_feature_axis.tsv")
    h3_axis_fields, h3_axis = read_tsv(molecular / "h3k27ac_feature_axis.tsv")
    if (
        fold_fields != ["participant_id", "outer_fold"]
        or set(participant_fields)
        != {
            "participant_index",
            "participant_id",
            "rna_source_sample_accession",
            "h3k27ac_source_sample_accession",
            "pairing",
            "rna_observation_state",
            "h3k27ac_observation_state",
        }
        or rna_axis_fields
        != [
            "rna_feature_index",
            "source_feature_index",
            "matrix_gene_id",
            "stable_gene_id",
        ]
        or h3_axis_fields
        != ["h3k27ac_feature_index", "opaque_source_feature_key"]
        or FORBIDDEN_FIELDS & set(participant_fields)
        or FORBIDDEN_FIELDS & set(fold_fields)
    ):
        raise PairedBaselineViewError("paired fixture schemas differ")
    participants = [row["participant_id"] for row in participant_rows]
    fold_participants = [row["participant_id"] for row in fold_rows]
    outer = np.asarray([int(row["outer_fold"]) for row in fold_rows], dtype=np.int8)
    if (
        len(participants) != expected_participants
        or len(set(participants)) != expected_participants
        or participants != fold_participants
        or {fold: int(np.sum(outer == fold)) for fold in range(5)}
        != dict(expected_fold_counts)
        or any(row["pairing"] != PAIRING for row in participant_rows)
        or any(row["rna_observation_state"] != "observed" for row in participant_rows)
        or any(row["h3k27ac_observation_state"] != "observed" for row in participant_rows)
        or len(rna_axis) != expected_rna_features
        or len(h3_axis) != expected_h3_features
        or [int(row["rna_feature_index"]) for row in rna_axis]
        != list(range(expected_rna_features))
        or [int(row["h3k27ac_feature_index"]) for row in h3_axis]
        != list(range(expected_h3_features))
        or len({row["opaque_source_feature_key"] for row in h3_axis})
        != expected_h3_features
    ):
        raise PairedBaselineViewError("participant, fold, or feature census differs")

    rna = np.load(molecular / "rna_values.npy", mmap_mode="r", allow_pickle=False)
    h3 = np.load(
        molecular / "h3k27ac_counts.npy", mmap_mode="r", allow_pickle=False
    )
    rna_mask = np.load(molecular / "rna_observed_mask.npy", allow_pickle=False)
    h3_mask = np.load(molecular / "h3k27ac_observed_mask.npy", allow_pickle=False)
    if (
        rna.shape != (expected_participants, expected_rna_features)
        or rna.dtype != np.float64
        or h3.shape != (expected_participants, expected_h3_features)
        or h3.dtype != np.uint32
        or rna_mask.shape != (expected_participants,)
        or h3_mask.shape != (expected_participants,)
        or rna_mask.dtype != np.bool_
        or h3_mask.dtype != np.bool_
        or not rna_mask.all()
        or not h3_mask.all()
    ):
        raise PairedBaselineViewError("paired molecular arrays or masks differ")

    model_root = output / "model"
    evaluator_root = output / "evaluator"
    model_root.mkdir(parents=True)
    evaluator_root.mkdir()
    rna_axis_sha256 = sha256_file(molecular / "rna_feature_axis.tsv")
    h3_axis_sha256 = sha256_file(molecular / "h3k27ac_feature_axis.tsv")
    fold_receipts: list[dict[str, Any]] = []
    for outer_fold in range(5):
        fold_root = model_root / f"outer_{outer_fold}"
        fold_root.mkdir()
        training_indices = np.flatnonzero(outer != outer_fold)
        query_indices = np.flatnonzero(outer == outer_fold)
        np.save(
            fold_root / "training_rna_values.npy",
            np.asarray(rna[training_indices], dtype=np.float64),
            allow_pickle=False,
        )
        np.save(
            fold_root / "query_rna_values.npy",
            np.asarray(rna[query_indices], dtype=np.float64),
            allow_pickle=False,
        )
        np.save(
            fold_root / "training_h3k27ac_counts.npy",
            np.asarray(h3[training_indices], dtype=np.uint32),
            allow_pickle=False,
        )
        training_rows = [
            {
                "participant_index": int(index),
                "participant_id": participants[index],
                "outer_fold": int(outer[index]),
            }
            for index in training_indices
        ]
        query_rows = [
            {
                "participant_index": int(index),
                "participant_id": participants[index],
                "outer_fold": int(outer[index]),
            }
            for index in query_indices
        ]
        write_tsv(
            fold_root / "training_participants.tsv",
            ["participant_index", "participant_id", "outer_fold"],
            training_rows,
        )
        write_tsv(
            fold_root / "query_participants.tsv",
            ["participant_index", "participant_id", "outer_fold"],
            query_rows,
        )
        receipt = {
            "schema_version": "masld-bench-gse267145-paired-model-view-v1",
            "status": "passed",
            "outer_fold": outer_fold,
            "training_participants": len(training_indices),
            "query_participants": len(query_indices),
            "rna_features": expected_rna_features,
            "h3k27ac_opaque_features": expected_h3_features,
            "pairing": PAIRING,
            "training_rna_included": True,
            "query_rna_included": True,
            "training_h3k27ac_included": True,
            "query_h3k27ac_included": False,
            "held_h3_library_totals_included": False,
            "outcome_or_participant_covariates_included": False,
            "masked_retired_rna_genes": 1_122,
            "rna_feature_axis_sha256": rna_axis_sha256,
            "h3_feature_axis_sha256": h3_axis_sha256,
            "h3_feature_identity": "opaque_source_feature_key",
            "coordinate_semantics_inferred": False,
            "task_spec_sha256": task_spec_sha256,
            "promotion_gate_sha256": promotion_gate_sha256,
        }
        write_json_exclusive(fold_root / "receipt.json", receipt)
        freeze_tree(
            fold_root,
            {
                "artifact_class": "gse267145_paired_model_fold_view",
                "outer_fold": outer_fold,
                "training_participants": len(training_indices),
                "query_participants": len(query_indices),
                "query_h3k27ac_included": False,
                "outcome_or_participant_covariates_included": False,
                "status": "passed",
            },
        )
        fold_receipts.append(receipt)

    np.save(
        evaluator_root / "observed_h3k27ac_counts.npy",
        np.asarray(h3, dtype=np.uint32),
        allow_pickle=False,
    )
    evaluator_rows = [
        {
            "participant_index": index,
            "participant_id": participant,
            "outer_fold": int(outer[index]),
            "h3k27ac_observation_state": "observed",
        }
        for index, participant in enumerate(participants)
    ]
    write_tsv(
        evaluator_root / "participant_axis.tsv",
        [
            "participant_index",
            "participant_id",
            "outer_fold",
            "h3k27ac_observation_state",
        ],
        evaluator_rows,
    )
    shutil.copyfile(
        molecular / "h3k27ac_feature_axis.tsv",
        evaluator_root / "h3k27ac_feature_axis.tsv",
    )
    evaluator_receipt = {
        "schema_version": "masld-bench-gse267145-paired-evaluator-view-v1",
        "status": "passed",
        "participants": expected_participants,
        "h3k27ac_opaque_features": expected_h3_features,
        "pairing": PAIRING,
        "observed_h3k27ac_included": True,
        "rna_values_included": False,
        "h3_feature_axis_sha256": h3_axis_sha256,
        "h3_feature_identity": "opaque_source_feature_key",
        "coordinate_semantics_inferred": False,
        "h3_library_totals_evaluator_only": True,
        "histology_or_participant_covariates_included": False,
        "task_spec_sha256": task_spec_sha256,
        "promotion_gate_sha256": promotion_gate_sha256,
    }
    write_json_exclusive(evaluator_root / "receipt.json", evaluator_receipt)
    freeze_tree(
        evaluator_root,
        {
            "artifact_class": "gse267145_paired_h3_evaluator_view",
            "participants": expected_participants,
            "h3k27ac_opaque_features": expected_h3_features,
            "rna_values_included": False,
            "histology_or_participant_covariates_included": False,
            "status": "passed",
        },
    )
    summary = {
        "schema_version": "masld-bench-gse267145-paired-baseline-views-v1",
        "status": "passed",
        "participants": expected_participants,
        "outer_folds": 5,
        "pairing": PAIRING,
        "rna_features": expected_rna_features,
        "h3k27ac_opaque_features": expected_h3_features,
        "masked_retired_rna_genes": 1_122,
        "model_fold_views": fold_receipts,
        "held_participant_h3_is_evaluator_only": True,
        "outcome_or_participant_covariates_included": False,
        "coordinate_semantics_inferred": False,
        "artifact_release_state": "internal_only_source_data_redistribution_prohibited",
        "model_fitted": False,
        "metrics_calculated": False,
    }
    write_json_exclusive(output / "receipt.json", summary)
    freeze_tree(
        model_root,
        {
            "artifact_class": "gse267145_paired_model_views",
            "outer_folds": 5,
            "held_participant_h3_included": False,
            "outcome_or_participant_covariates_included": False,
            "status": "passed",
        },
    )
    freeze_tree(
        output,
        {
            "artifact_class": "gse267145_paired_baseline_views",
            "participants": expected_participants,
            "outer_folds": 5,
            "held_participant_h3_is_evaluator_only": True,
            "model_fitted": False,
            "metrics_calculated": False,
            "status": "passed",
        },
    )
    return summary


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--molecular", required=True, type=Path)
    value.add_argument("--molecular-artifacts-sha256", required=True)
    value.add_argument("--folds", required=True, type=Path)
    value.add_argument("--folds-artifacts-sha256", required=True)
    value.add_argument("--task-spec", required=True, type=Path)
    value.add_argument("--task-spec-sha256", required=True)
    value.add_argument("--promotion-gate", required=True, type=Path)
    value.add_argument("--promotion-gate-sha256", required=True)
    value.add_argument("--output", required=True, type=Path)
    return value


def main() -> int:
    arguments = parser().parse_args()
    for root, expected in (
        (arguments.molecular, arguments.molecular_artifacts_sha256),
        (arguments.folds, arguments.folds_artifacts_sha256),
    ):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise PairedBaselineViewError("input ARTIFACTS SHA-256 differs")
        verify_frozen_tree(root)
    if sha256_file(arguments.task_spec) != arguments.task_spec_sha256:
        raise PairedBaselineViewError("TaskSpec SHA-256 differs")
    if sha256_file(arguments.promotion_gate) != arguments.promotion_gate_sha256:
        raise PairedBaselineViewError("promotion gate SHA-256 differs")
    gate = json.loads(arguments.promotion_gate.read_text(encoding="utf-8"))
    if gate.get("champion_eligible") is not False or gate.get("claim_mode") != "development_only":
        raise PairedBaselineViewError("promotion gate is not development-only")
    result = build_views(
        molecular=arguments.molecular,
        folds=arguments.folds,
        output=arguments.output,
        task_spec_sha256=arguments.task_spec_sha256,
        promotion_gate_sha256=arguments.promotion_gate_sha256,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
