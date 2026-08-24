#!/usr/bin/env python3
"""Evaluate Seurat WNN bridge embeddings with evaluator-only nucleus pairs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from evaluate_scglue_retrieval_smoke import (
    directional_ranks,
    load_embeddings,
    read_tsv,
)


MODELS = ("seurat_wnn_bridge", "linear_cca")


class SeuratWNNRetrievalEvaluationError(ValueError):
    """Raised when WNN retrieval evaluation differs."""


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


def load_wnn_embedding(path: Path) -> tuple[list[str], Any]:
    import numpy as np

    fields, rows = read_tsv(path)
    expected = ("query_id",) + tuple(f"bridge-lap-{index}" for index in range(1, 17))
    if fields != expected or not rows:
        raise SeuratWNNRetrievalEvaluationError("WNN embedding schema differs")
    identifiers = [row["query_id"] for row in rows]
    if len(set(identifiers)) != len(identifiers):
        raise SeuratWNNRetrievalEvaluationError("WNN embedding ID is duplicated")
    values = np.asarray(
        [[float(row[field]) for field in fields[1:]] for row in rows], dtype=np.float64
    )
    if values.shape != (len(identifiers), 16) or np.any(~np.isfinite(values)):
        raise SeuratWNNRetrievalEvaluationError("WNN embedding values differ")
    return identifiers, values


def evaluate(
    wnn_predictions: Path,
    baseline_predictions: Path,
    evaluator_only: Path,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise SeuratWNNRetrievalEvaluationError("output exists")
    donor_values: dict[tuple[str, str], list[tuple[int, int]]] = {}
    total_pairs = 0
    for fold in range(5):
        fields, pairs = read_tsv(evaluator_only / f"fold_{fold}/hidden_pairs.tsv")
        if fields != ("pair_hash", "donor_hash", "rna_query_id", "atac_query_id"):
            raise SeuratWNNRetrievalEvaluationError("hidden pair schema differs")
        total_pairs += len(pairs)
        rna_to_atac = {row["rna_query_id"]: row["atac_query_id"] for row in pairs}
        atac_to_rna = {row["atac_query_id"]: row["rna_query_id"] for row in pairs}
        donor_by_rna = {row["rna_query_id"]: row["donor_hash"] for row in pairs}
        donor_by_atac = {row["atac_query_id"]: row["donor_hash"] for row in pairs}
        wnn_root = wnn_predictions / f"fold_{fold}/seurat_wnn_bridge"
        wnn_rna_ids, wnn_rna = load_wnn_embedding(
            wnn_root / "query_rna_embeddings.tsv"
        )
        wnn_atac_ids, wnn_atac = load_wnn_embedding(
            wnn_root / "query_atac_embeddings.tsv"
        )
        baseline_rna_ids, baseline_rna, baseline_atac_ids, baseline_atac = load_embeddings(
            baseline_predictions / f"fold_{fold}/linear_cca/query_embeddings.npz"
        )
        embeddings = {
            "seurat_wnn_bridge": (wnn_rna_ids, wnn_rna, wnn_atac_ids, wnn_atac),
            "linear_cca": (
                baseline_rna_ids,
                baseline_rna,
                baseline_atac_ids,
                baseline_atac,
            ),
        }
        for model_id, (rna_ids, rna, atac_ids, atac) in embeddings.items():
            if set(rna_ids) != set(rna_to_atac) or set(atac_ids) != set(atac_to_rna):
                raise SeuratWNNRetrievalEvaluationError(
                    "embedding and hidden-pair IDs differ"
                )
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
                "top1_accuracy": format(
                    sum(rank == 1 for rank in values) / len(values), ".17g"
                ),
                "top5_accuracy": format(
                    sum(rank <= 5 for rank in values) / len(values), ".17g"
                ),
                "mean_null_top1": format(
                    sum(1.0 / size for _rank, size in ranks) / len(ranks), ".17g"
                ),
            }
        )
    summaries: dict[str, Any] = {}
    for model_id in MODELS:
        rows = [row for row in donor_rows if row["model_id"] == model_id]
        if len(rows) != 39:
            raise SeuratWNNRetrievalEvaluationError("model donor census differs")
        summaries[model_id] = {
            "donor_macro_mean_reciprocal_rank": sum(
                float(row["mean_reciprocal_rank"]) for row in rows
            )
            / len(rows),
            "donor_macro_top1_accuracy": sum(
                float(row["top1_accuracy"]) for row in rows
            )
            / len(rows),
            "donor_macro_top5_accuracy": sum(
                float(row["top5_accuracy"]) for row in rows
            )
            / len(rows),
            "donor_macro_null_top1": sum(
                float(row["mean_null_top1"]) for row in rows
            )
            / len(rows),
        }
    gain = (
        summaries["seurat_wnn_bridge"]["donor_macro_mean_reciprocal_rank"]
        - summaries["linear_cca"]["donor_macro_mean_reciprocal_rank"]
    )
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
        "schema_version": "masld-bench-seurat-wnn-bridge-evaluation-v1",
        "status": "pass",
        "dataset_id": "gse296875",
        "task_id": "same_nucleus_cross_modal_retrieval",
        "model_roster": list(MODELS),
        "outer_unit": "donor",
        "n_donors": 39,
        "n_same_nucleus_pairs": total_pairs,
        "directions": ["rna_to_atac", "atac_to_rna"],
        "summaries": summaries,
        "mrr_gain_over_linear_cca": gain,
        "smoke_gate": {
            "minimum_absolute_mrr_gain": 0.02,
            "seurat_wnn_bridge_passed": gain >= 0.02,
        },
        "hidden_pair_map_available_to_model": False,
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
    parser.add_argument("--wnn-predictions", type=Path, required=True)
    parser.add_argument("--baseline-predictions", type=Path, required=True)
    parser.add_argument("--evaluator-only", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    evaluate(
        arguments.wnn_predictions,
        arguments.baseline_predictions,
        arguments.evaluator_only,
        arguments.output,
    )


if __name__ == "__main__":
    main()
