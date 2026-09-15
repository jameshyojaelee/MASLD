#!/usr/bin/env python3
"""Apply exact Atlas outer-fit TF-Sapiens heads to label-free query embeddings."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import torch

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import canonical_sha256, sha256_file


VIEW_ID = "gse296875_rna_cell_state_external_development_7500_v1"
SPLIT_SHA256 = "10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680"
ROSTER = (
    "cholangiocyte",
    "endothelial",
    "hepatocyte",
    "immune",
    "mesenchymal_stromal",
)
SEEDS = (1103, 1201, 1301)
FOLDS = tuple(range(5))
HEADS = ("linear", "two_layer_mlp")


class TranscriptFormerTransferError(ValueError):
    """Raised when source-fit or query prediction requirement differs."""


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or ()), list(reader)


def _write_tsv(path: Path, rows: Iterable[Mapping[str, Any]]) -> None:
    fields = ("row_id", "predicted_class", *(f"probability::{label}" for label in ROSTER))
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _head_logits(state_path: Path, features: np.ndarray, head_id: str) -> np.ndarray:
    with np.load(state_path, allow_pickle=False) as values:
        expected = (
            {"feature_mean", "feature_scale", "state::weight", "state::bias"}
            if head_id == "linear"
            else {
                "feature_mean", "feature_scale", "state::0.weight", "state::0.bias",
                "state::3.weight", "state::3.bias",
            }
        )
        if set(values.files) != expected:
            raise TranscriptFormerTransferError(f"{head_id} head state schema differs")
        mean = np.asarray(values["feature_mean"], dtype=np.float32)
        scale = np.asarray(values["feature_scale"], dtype=np.float32)
        tensors = {key: torch.from_numpy(np.asarray(values[key])) for key in values.files if key.startswith("state::")}
    if mean.shape != (2_048,) or scale.shape != (2_048,) or np.any(scale <= 0):
        raise TranscriptFormerTransferError("head standardization state differs")
    x = torch.from_numpy(((features - mean) / scale).astype(np.float32, copy=False))
    with torch.inference_mode():
        if head_id == "linear":
            return torch.nn.functional.linear(
                x, tensors["state::weight"], tensors["state::bias"]
            ).numpy()
        hidden = torch.nn.functional.gelu(
            torch.nn.functional.linear(x, tensors["state::0.weight"], tensors["state::0.bias"])
        )
        return torch.nn.functional.linear(
            hidden, tensors["state::3.weight"], tensors["state::3.bias"]
        ).numpy()


def _softmax(logits: np.ndarray, temperature: float) -> np.ndarray:
    values = np.asarray(logits, dtype=np.float64) / float(temperature)
    values -= values.max(axis=1, keepdims=True)
    values = np.exp(values)
    values /= values.sum(axis=1, keepdims=True)
    return values


def run(args: argparse.Namespace) -> dict[str, Any]:
    if args.output.exists():
        raise TranscriptFormerTransferError("refusing to overwrite TF-Sapiens transfer")
    query_root = args.query_embeddings.resolve(strict=True)
    source_embeddings = args.source_embeddings.resolve(strict=True)
    heads = args.source_heads.resolve(strict=True)
    split = args.split.resolve(strict=True)
    for root, expected in (
        (query_root, args.expected_query_embeddings_sha256),
        (source_embeddings, args.expected_source_embeddings_sha256),
        (heads, args.expected_source_heads_sha256),
        (split, SPLIT_SHA256),
    ):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise TranscriptFormerTransferError(f"input ARTIFACTS SHA-256 differs: {root}")
        verify_frozen_tree(root)
    query_meta = verify_frozen_tree(query_root)["metadata"]
    source_meta = verify_frozen_tree(source_embeddings)["metadata"]
    head_meta = verify_frozen_tree(heads)["metadata"]
    if (
        query_meta.get("artifact_class") != "gse296875_transcriptformer_outcome_blind_embeddings"
        or query_meta.get("view_id") != VIEW_ID
        or query_meta.get("rows") != 7_500
        or query_meta.get("labels_read") is not False
        or query_meta.get("donor_ids_read") is not False
        or query_meta.get("target_adaptation_performed") is not False
        or source_meta.get("artifact_class") != "cell_foundation_frozen_screen_embeddings"
        or source_meta.get("model_id") != "transcriptformer_tf_sapiens"
        or source_meta.get("development_rows") != 50_000
        or head_meta.get("artifact_class") != "cell_foundation_common_head_shard_bundle"
        or head_meta.get("logical_shards") != 30
    ):
        raise TranscriptFormerTransferError("TF-Sapiens input authority differs")

    with np.load(query_root / "embeddings/common_embeddings.npz", allow_pickle=False) as bundle:
        if set(bundle.files) != {"embeddings", "outer_folds", "row_ids"}:
            raise TranscriptFormerTransferError("query embedding fields differ")
        query = np.asarray(bundle["embeddings"], dtype=np.float32)
        query_rows = np.asarray(bundle["row_ids"]).astype(str)
        query_folds = np.asarray(bundle["outer_folds"])
    with np.load(source_embeddings / "embeddings/common_embeddings.npz", allow_pickle=False) as bundle:
        source = np.asarray(bundle["embeddings"], dtype=np.float32)
        source_rows = np.asarray(bundle["row_ids"]).astype(str)
        source_activation_outer = np.asarray(bundle["outer_folds"], dtype=np.int8)
    split_fields, split_rows = _read_tsv(split / "row_outer_folds.tsv")
    if split_fields != ["row_id", "donor_id", "dataset", "outer_fold"]:
        raise TranscriptFormerTransferError("study-held split schema differs")
    source_outer = np.asarray([int(row["outer_fold"]) for row in split_rows], dtype=np.int8)
    if (
        query.shape != (7_500, 2_048)
        or len(set(query_rows)) != 7_500
        or np.any(query_folds != -1)
        or source.shape != (50_000, 2_048)
        or len(set(source_rows)) != 50_000
        or [row["row_id"] for row in split_rows] != source_rows.tolist()
        or set(source_outer.tolist()) != set(FOLDS)
    ):
        raise TranscriptFormerTransferError("embedding tensor or row contract differs")
    activation_split_mismatch_rows = int(np.sum(source_activation_outer != source_outer))
    if activation_split_mismatch_rows != 42_584:
        raise TranscriptFormerTransferError("activation-versus-study fold mismatch census differs")

    query_values: dict[str, list[np.ndarray]] = {head_id: [] for head_id in HEADS}
    audit_rows: list[dict[str, Any]] = []
    max_reproduction_error = 0.0
    for head_id in HEADS:
      for seed in SEEDS:
        for fold in FOLDS:
            shard = heads / "shards" / f"{head_id}__seed{seed}__fold{fold}"
            receipt = json.loads((shard / "prediction_receipt.json").read_text())
            if (
                receipt.get("model_id") != "transcriptformer_tf_sapiens"
                or receipt.get("head_id") != head_id
                or receipt.get("screen_seed") != seed
                or receipt.get("outer_fold") != fold
                or receipt.get("outer_prediction_uses_final_refit") is not True
                or receipt.get("input_artifacts_sha256", {}).get("split") != SPLIT_SHA256
            ):
                raise TranscriptFormerTransferError("source-fit receipt differs")
            state = shard / "states/final_outer_training_refit.npz"
            temperature = float(receipt["temperature"])
            held = np.flatnonzero(source_outer == fold)
            reproduced = _softmax(_head_logits(state, source[held], head_id), temperature)
            fields, observed_rows = _read_tsv(shard / "predictions.tsv")
            expected_fields = [
                "row_id", "donor_id", "dataset", "outer_fold", "predicted_class",
                *(f"probability::{label}" for label in ROSTER),
            ]
            if fields != expected_fields or [row["row_id"] for row in observed_rows] != source_rows[held].tolist():
                raise TranscriptFormerTransferError("source-fit reproduction row contract differs")
            observed = np.asarray(
                [[float(row[f"probability::{label}"]) for label in ROSTER] for row in observed_rows]
            )
            error = float(np.max(np.abs(reproduced - observed)))
            if error > 2.0e-6:
                raise TranscriptFormerTransferError(f"source prediction reproduction failed: {error}")
            max_reproduction_error = max(max_reproduction_error, error)
            query_probabilities = _softmax(_head_logits(state, query, head_id), temperature)
            query_values[head_id].append(query_probabilities)
            audit_rows.append(
                {
                    "head_id": head_id,
                    "seed": seed,
                    "outer_fold": fold,
                    "training_rows": receipt["training_rows"],
                    "training_donors": receipt["training_donors"],
                    "training_studies": ",".join(receipt["training_studies"]),
                    "temperature": format(temperature, ".17g"),
                    "state_sha256": sha256_file(state),
                    "source_prediction_max_abs_error": format(error, ".17g"),
                }
            )
    args.output.mkdir(mode=0o750)
    (args.output / "predictions").mkdir()
    model_ids = []
    for head_id in HEADS:
        if len(query_values[head_id]) != 15:
            raise TranscriptFormerTransferError("head source-fit census differs")
        ensemble = np.mean(np.stack(query_values[head_id]), axis=0)
        ensemble /= ensemble.sum(axis=1, keepdims=True)
        if not np.all(np.isfinite(ensemble)) or not np.allclose(ensemble.sum(axis=1), 1.0, atol=1e-12):
            raise TranscriptFormerTransferError("query ensemble probabilities differ")
        model_id = f"transcriptformer_tf_sapiens__{head_id}"
        model_ids.append(model_id)
        predicted = np.argmax(ensemble, axis=1)
        _write_tsv(
            args.output / "predictions" / f"{model_id}.tsv",
            (
                {
                    "row_id": query_rows[index],
                    "predicted_class": ROSTER[int(predicted[index])],
                    **{f"probability::{label}": format(float(ensemble[index, offset]), ".17g") for offset, label in enumerate(ROSTER)},
                }
                for index in range(len(query_rows))
            ),
        )
    with (args.output / "source_fit_audit.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(audit_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(audit_rows)
    receipt = {
        "schema_version": "masld-bench-gse296875-transcriptformer-outer-ensemble-v1",
        "status": "pass_outcome_blind_predictions",
        "view_id": VIEW_ID,
        "model_id": "transcriptformer_tf_sapiens",
        "head_ids": list(HEADS),
        "rows": 7_500,
        "source_rows": 50_000,
        "source_fit_census_per_head": 15,
        "source_fit_census_total": 30,
        "source_donors": 102,
        "source_studies": 7,
        "source_split_artifacts_sha256": SPLIT_SHA256,
        "source_fold_authority": "resource_atlas_study_outer_5fold_v1",
        "embedding_activation_fold_annotations_used_for_head_reproduction": False,
        "activation_vs_study_fold_mismatch_rows": activation_split_mismatch_rows,
        "ensemble": "unweighted_probability_mean_over_three_seeds_by_five_study_outer_fits",
        "maximum_source_prediction_reproduction_error": max_reproduction_error,
        "query_row_ids_sha256": canonical_sha256(query_rows.tolist()),
        "query_labels_read": False,
        "query_donor_ids_read": False,
        "query_atac_read": False,
        "atac_state": "structurally_missing",
        "target_adaptation_performed": False,
        "all_source_refit_performed": False,
        "metrics_calculated": False,
        "project_exposed_development_only": True,
        "external_or_champion_claim_allowed": False,
        "sealed_outcomes_read": False,
    }
    write_json_exclusive(args.output / "prediction_receipt.json", receipt)
    freeze_tree(
        args.output,
        {
            "artifact_class": "gse296875_cell_state_outcome_blind_prediction_branch",
            "branch_id": "transcriptformer_outer_ensemble",
            "view_id": VIEW_ID,
            "models": model_ids,
            "rows": 7_500,
            "source_fit_census_per_model": 15,
            "source_fit_census_total": 30,
            "query_labels_read": False,
            "query_donor_ids_read": False,
            "query_atac_read": False,
            "target_adaptation_performed": False,
            "metrics_calculated": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--query-embeddings", type=Path, required=True)
    parser.add_argument("--expected-query-embeddings-sha256", required=True)
    parser.add_argument("--source-embeddings", type=Path, required=True)
    parser.add_argument("--expected-source-embeddings-sha256", required=True)
    parser.add_argument("--source-heads", type=Path, required=True)
    parser.add_argument("--expected-source-heads-sha256", required=True)
    parser.add_argument(
        "--split",
        type=Path,
        default=Path(__file__).parents[1] / "executions/model-data-052-21077595/split",
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(run(parse_args()), sort_keys=True))
