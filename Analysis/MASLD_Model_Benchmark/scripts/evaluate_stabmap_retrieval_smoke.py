#!/usr/bin/env python3
"""Evaluate frozen StabMap folds with evaluator-only same-nucleus pairs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
from typing import Any, Iterable, Mapping, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))
from evaluate_scglue_retrieval_smoke import directional_ranks, load_embeddings, read_tsv


MODELS = ("stabmap", "linear_cca")


class StabMapRetrievalEvaluationError(ValueError):
    """Raised when StabMap retrieval evaluation differs."""


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


def load_stabmap_embedding(path: Path) -> tuple[list[str], Any]:
    import numpy as np

    fields, rows = read_tsv(path)
    expected = ("query_id",) + tuple(
        f"reference_rna_PC{index}" for index in range(1, 17)
    )
    if fields != expected or not rows:
        raise StabMapRetrievalEvaluationError("StabMap embedding schema differs")
    identifiers = [row["query_id"] for row in rows]
    if len(set(identifiers)) != len(identifiers):
        raise StabMapRetrievalEvaluationError("StabMap embedding ID is duplicated")
    values = np.asarray(
        [[float(row[field]) for field in fields[1:]] for row in rows], dtype=np.float64
    )
    if values.shape != (len(identifiers), 16) or np.any(~np.isfinite(values)):
        raise StabMapRetrievalEvaluationError("StabMap embedding values differ")
    return identifiers, values


def parse_fold_roots(values: Sequence[str]) -> dict[int, Path]:
    roots: dict[int, Path] = {}
    for value in values:
        fold_text, separator, path_text = value.partition("=")
        if not separator or not fold_text.isdigit():
            raise StabMapRetrievalEvaluationError("fold-root syntax differs")
        fold = int(fold_text)
        if fold in roots or fold not in range(5):
            raise StabMapRetrievalEvaluationError("fold-root census differs")
        roots[fold] = Path(path_text)
    if set(roots) != set(range(5)):
        raise StabMapRetrievalEvaluationError("five fold roots are required")
    return roots


def evaluate(
    fold_roots: Mapping[int, Path],
    baseline_predictions: Path,
    evaluator_only: Path,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise StabMapRetrievalEvaluationError("output exists")
    donor_values: dict[tuple[str, str], list[tuple[int, int]]] = {}
    total_pairs = 0
    for fold in range(5):
        fields, pairs = read_tsv(evaluator_only / f"fold_{fold}/hidden_pairs.tsv")
        if fields != ("pair_hash", "donor_hash", "rna_query_id", "atac_query_id"):
            raise StabMapRetrievalEvaluationError("hidden pair schema differs")
        if len({row["pair_hash"] for row in pairs}) != len(pairs):
            raise StabMapRetrievalEvaluationError("hidden pair is duplicated")
        total_pairs += len(pairs)
        rna_to_atac = {row["rna_query_id"]: row["atac_query_id"] for row in pairs}
        atac_to_rna = {row["atac_query_id"]: row["rna_query_id"] for row in pairs}
        donor_by_rna = {row["rna_query_id"]: row["donor_hash"] for row in pairs}
        donor_by_atac = {row["atac_query_id"]: row["donor_hash"] for row in pairs}
        candidate = fold_roots[fold] / f"predictions/fold_{fold}/stabmap"
        rna_ids, rna = load_stabmap_embedding(candidate / "query_rna_embeddings.tsv")
        atac_ids, atac = load_stabmap_embedding(candidate / "query_atac_embeddings.tsv")
        baseline_rna_ids, baseline_rna, baseline_atac_ids, baseline_atac = load_embeddings(
            baseline_predictions / f"fold_{fold}/linear_cca/query_embeddings.npz"
        )
        embeddings = {
            "stabmap": (rna_ids, rna, atac_ids, atac),
            "linear_cca": (
                baseline_rna_ids,
                baseline_rna,
                baseline_atac_ids,
                baseline_atac,
            ),
        }
        for model_id, (model_rna_ids, model_rna, model_atac_ids, model_atac) in embeddings.items():
            if set(model_rna_ids) != set(rna_to_atac) or set(model_atac_ids) != set(atac_to_rna):
                raise StabMapRetrievalEvaluationError("embedding and hidden-pair IDs differ")
            forward = directional_ranks(
                model_rna_ids, model_atac_ids, model_rna, model_atac, rna_to_atac
            )
            reverse = directional_ranks(
                model_atac_ids, model_rna_ids, model_atac, model_rna, atac_to_rna
            )
            for identifier, rank in forward.items():
                donor_values.setdefault((model_id, donor_by_rna[identifier]), []).append(
                    (rank, len(model_atac_ids))
                )
            for identifier, rank in reverse.items():
                donor_values.setdefault((model_id, donor_by_atac[identifier]), []).append(
                    (rank, len(model_rna_ids))
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
            raise StabMapRetrievalEvaluationError("model donor census differs")
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
        summaries["stabmap"]["donor_macro_mean_reciprocal_rank"]
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
        "schema_version": "masld-bench-stabmap-retrieval-evaluation-v1",
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
            "stabmap_passed": gain >= 0.02,
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
    parser.add_argument("--fold-root", action="append", required=True)
    parser.add_argument("--baseline-predictions", type=Path, required=True)
    parser.add_argument("--evaluator-only", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    evaluate(
        parse_fold_roots(arguments.fold_root),
        arguments.baseline_predictions,
        arguments.evaluator_only,
        arguments.output,
    )


if __name__ == "__main__":
    main()
