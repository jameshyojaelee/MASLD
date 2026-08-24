#!/usr/bin/env python3
"""Evaluate frozen scGLUE embeddings with evaluator-only same-nucleus pairs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


MODELS = ("scglue", "paired_scglue", "linear_cca")


class SCGLUERetrievalEvaluationError(ValueError):
    """Raised when blinded cross-modal retrieval evaluation differs."""


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SCGLUERetrievalEvaluationError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), [dict(row) for row in reader]


def write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def normalize_rows(values: Any) -> Any:
    import numpy as np

    matrix = np.asarray(values, dtype=np.float64)
    if matrix.ndim != 2:
        raise SCGLUERetrievalEvaluationError("embedding rows are invalid")
    norms = np.linalg.norm(matrix, axis=1)
    if np.any(~np.isfinite(matrix)) or np.any(norms <= 0):
        raise SCGLUERetrievalEvaluationError("embedding rows are invalid")
    return matrix / norms[:, None]


def directional_ranks(
    query_ids: Sequence[str],
    target_ids: Sequence[str],
    query_embedding: Any,
    target_embedding: Any,
    correct_target: Mapping[str, str],
) -> dict[str, int]:
    import numpy as np

    query = normalize_rows(query_embedding)
    target = normalize_rows(target_embedding)
    if query.shape[1] != target.shape[1] or len(set(target_ids)) != len(target_ids):
        raise SCGLUERetrievalEvaluationError("retrieval embedding axes differ")
    target_index = {identifier: index for index, identifier in enumerate(target_ids)}
    scores = query @ target.T
    result: dict[str, int] = {}
    for index, identifier in enumerate(query_ids):
        expected = correct_target.get(identifier)
        if expected not in target_index:
            raise SCGLUERetrievalEvaluationError("correct target is absent")
        expected_index = target_index[expected]
        order = np.argsort(-scores[index], kind="stable")
        positions = np.flatnonzero(order == expected_index)
        if len(positions) != 1:
            raise SCGLUERetrievalEvaluationError("correct target rank differs")
        result[identifier] = int(positions[0]) + 1
    if set(result) != set(query_ids):
        raise SCGLUERetrievalEvaluationError("query rank universe differs")
    return result


def load_embeddings(path: Path) -> tuple[list[str], Any, list[str], Any]:
    import numpy as np

    with np.load(path, allow_pickle=False) as values:
        if set(values.files) != {
            "rna_ids",
            "rna_embedding",
            "atac_ids",
            "atac_embedding",
        }:
            raise SCGLUERetrievalEvaluationError("embedding artifact keys differ")
        rna_ids = [str(value) for value in values["rna_ids"]]
        atac_ids = [str(value) for value in values["atac_ids"]]
        rna = values["rna_embedding"]
        atac = values["atac_embedding"]
    if (
        len(rna_ids) != rna.shape[0]
        or len(atac_ids) != atac.shape[0]
        or len(set(rna_ids)) != len(rna_ids)
        or len(set(atac_ids)) != len(atac_ids)
        or set(rna_ids) & set(atac_ids)
    ):
        raise SCGLUERetrievalEvaluationError("embedding identity axes differ")
    return rna_ids, rna, atac_ids, atac


def evaluate(predictions: Path, evaluator_only: Path, output: Path) -> dict[str, Any]:
    if output.exists():
        raise SCGLUERetrievalEvaluationError("output exists")
    donor_values: dict[tuple[str, str], list[tuple[int, int]]] = {}
    total_pairs = 0
    for fold in range(5):
        fields, pairs = read_tsv(evaluator_only / f"fold_{fold}/hidden_pairs.tsv")
        if fields != ("pair_hash", "donor_hash", "rna_query_id", "atac_query_id"):
            raise SCGLUERetrievalEvaluationError("hidden pair schema differs")
        if len({row["pair_hash"] for row in pairs}) != len(pairs):
            raise SCGLUERetrievalEvaluationError("hidden pair is duplicated")
        total_pairs += len(pairs)
        rna_to_atac = {row["rna_query_id"]: row["atac_query_id"] for row in pairs}
        atac_to_rna = {row["atac_query_id"]: row["rna_query_id"] for row in pairs}
        donor_by_rna = {row["rna_query_id"]: row["donor_hash"] for row in pairs}
        donor_by_atac = {row["atac_query_id"]: row["donor_hash"] for row in pairs}
        for model_id in MODELS:
            rna_ids, rna, atac_ids, atac = load_embeddings(
                predictions / f"fold_{fold}/{model_id}/query_embeddings.npz"
            )
            if set(rna_ids) != set(rna_to_atac) or set(atac_ids) != set(atac_to_rna):
                raise SCGLUERetrievalEvaluationError("embedding and hidden-pair IDs differ")
            forward = directional_ranks(rna_ids, atac_ids, rna, atac, rna_to_atac)
            reverse = directional_ranks(atac_ids, rna_ids, atac, rna, atac_to_rna)
            for identifier, rank in forward.items():
                donor_values.setdefault((model_id, donor_by_rna[identifier]), []).append(
                    (rank, len(atac_ids))
                )
            for identifier, rank in reverse.items():
                donor_values.setdefault((model_id, donor_by_atac[identifier]), []).append(
                    (rank, len(rna_ids))
                )
    donor_rows: list[dict[str, Any]] = []
    for (model_id, donor_hash), ranks in sorted(donor_values.items()):
        values = [rank for rank, _size in ranks]
        donor_rows.append(
            {
                "model_id": model_id,
                "donor_hash": donor_hash,
                "directional_queries": len(values),
                "mean_reciprocal_rank": format(
                    sum(1.0 / rank for rank in values) / len(values), ".17g"
                ),
                "top1_accuracy": format(sum(rank == 1 for rank in values) / len(values), ".17g"),
                "top5_accuracy": format(sum(rank <= 5 for rank in values) / len(values), ".17g"),
                "mean_null_top1": format(
                    sum(1.0 / size for _rank, size in ranks) / len(ranks), ".17g"
                ),
            }
        )
    summaries: dict[str, Any] = {}
    for model_id in MODELS:
        rows = [row for row in donor_rows if row["model_id"] == model_id]
        if len(rows) != 39:
            raise SCGLUERetrievalEvaluationError("model donor census differs")
        summaries[model_id] = {
            "donor_macro_mean_reciprocal_rank": sum(
                float(row["mean_reciprocal_rank"]) for row in rows
            )
            / len(rows),
            "donor_macro_top1_accuracy": sum(float(row["top1_accuracy"]) for row in rows)
            / len(rows),
            "donor_macro_top5_accuracy": sum(float(row["top5_accuracy"]) for row in rows)
            / len(rows),
            "donor_macro_null_top1": sum(float(row["mean_null_top1"]) for row in rows)
            / len(rows),
        }
    baseline = summaries["linear_cca"]["donor_macro_mean_reciprocal_rank"]
    gains = {
        model_id: summaries[model_id]["donor_macro_mean_reciprocal_rank"] - baseline
        for model_id in ("scglue", "paired_scglue")
    }
    output.mkdir(parents=True, mode=0o750)
    write_tsv(
        output / "donor_retrieval_metrics.tsv",
        (
            "model_id",
            "donor_hash",
            "directional_queries",
            "mean_reciprocal_rank",
            "top1_accuracy",
            "top5_accuracy",
            "mean_null_top1",
        ),
        donor_rows,
    )
    result = {
        "schema_version": "masld-bench-scglue-retrieval-evaluation-v1",
        "status": "pass",
        "dataset_id": "gse296875",
        "task_id": "same_nucleus_cross_modal_retrieval",
        "model_roster": list(MODELS),
        "outer_unit": "donor",
        "n_donors": 39,
        "n_same_nucleus_pairs": total_pairs,
        "directions": ["rna_to_atac", "atac_to_rna"],
        "summaries": summaries,
        "mrr_gain_over_linear_cca": gains,
        "smoke_gate": {
            "minimum_absolute_mrr_gain": 0.02,
            "scglue_passed": gains["scglue"] >= 0.02,
            "paired_scglue_passed": gains["paired_scglue"] >= 0.02,
        },
        "hidden_pair_map_available_to_models": False,
        "query_modality_ids_disjoint": True,
        "query_atac_row_order_permuted": True,
        "cells_used_as_biological_replicates": False,
        "donor_macro_aggregation": True,
        "test_outcomes_read": False,
        "rna_to_atac_prediction_claim": False,
        "sealed_rna_conditioned_atac_eligible": False,
        "smoke_only": True,
        "champion_claim_allowed": False,
    }
    (output / "evaluation.json").write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, sort_keys=True))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--evaluator-only", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    evaluate(arguments.predictions, arguments.evaluator_only, arguments.output)


if __name__ == "__main__":
    main()
