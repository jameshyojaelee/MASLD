#!/usr/bin/env python
"""Geneformer common-lane adapter for the frozen cell-state task.

The encoder is frozen.  Only the common head is trained, and it is trained on
outer-training donors alone, so the held-out donors of the fold never influence
any fitted state.  Embeddings follow the registered policy exactly:
``last_hidden_state`` averaged over non-special, non-padding tokens.  The
pooler is deliberately unused; upstream ships no trained pooler and transformers
re-initialises one on load, which would otherwise be silent untrained state.

Rank-value encoding reproduces Geneformer's own contract: counts are scaled to
a fixed library size, divided by the frozen gene median, ranked descending, and
truncated to the model's ``max_position_embeddings``.  Genes absent from either
the median or token dictionary are dropped, and the dropped fraction is
reported rather than hidden.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from hvg_pca_logistic import (  # type: ignore[import-not-found]
    ScientificAdapterError,
    _artifact_record,
    _input_by_role,
    _require_sha256,
    _validate_file_artifact,
    _write_json,
    _write_tsv,
    fold_index,
    join_hash,
)

TASK_ID = "cell_state_mapping"
DATASET_ID = "resource_atlas_current"
RUNTIME_ID = "gpu_geneformer"
RECEIPT_SCHEMA = "masld-bench-adapter-receipt-v1"
ADAPTER_ID = "geneformer_common_lane_v1"
MODEL_IDS = ("geneformer_v1_10m", "geneformer_v2_104m", "geneformer_v2_316m")
HEAD_IDS = ("linear", "two_layer_mlp")
EMBEDDING_POLICY = (
    "last_hidden_state_mean_over_non_special_nonpadding_tokens"
)


def _select_device(requested: str) -> str:
    import torch

    if requested == "cpu":
        return "cpu"
    if torch.cuda.is_available():
        return "cuda"
    if requested == "cuda":
        raise ScientificAdapterError("CUDA was required but is unavailable")
    return "cpu"


def _load_dictionaries(bundle: Path, variant: str) -> tuple[dict, dict]:
    stem = "gene_dictionaries_30m/" if variant == "V1" else ""
    suffix = "gc30M" if variant == "V1" else "gc104M"
    root = bundle / "sanitized" / "geneformer" / stem
    median = json.loads(
        (root / f"gene_median_dictionary_{suffix}.pkl.json").read_text("utf-8")
    )
    token = json.loads(
        (root / f"token_dictionary_{suffix}.pkl.json").read_text("utf-8")
    )
    if not median or not token:
        raise ScientificAdapterError("Geneformer dictionaries are empty")
    return median, token


def _rank_value_encode(
    counts: Any,
    gene_ids: Sequence[str],
    median: Mapping[str, float],
    token: Mapping[str, int],
    max_tokens: int,
) -> tuple[list[list[int]], dict[str, Any]]:
    import numpy as np

    usable = [
        index
        for index, gene in enumerate(gene_ids)
        if gene in median and gene in token and float(median[gene]) > 0.0
    ]
    if not usable:
        raise ScientificAdapterError("no gene maps into both frozen dictionaries")
    medians = np.asarray([float(median[gene_ids[i]]) for i in usable])
    tokens = np.asarray([int(token[gene_ids[i]]) for i in usable])
    encoded: list[list[int]] = []
    truncated = 0
    for row in range(counts.shape[0]):
        values = np.asarray(counts[row, usable], dtype=float).ravel()
        total = values.sum()
        if total <= 0.0:
            raise ScientificAdapterError("a cell has no usable counts")
        normalized = (values / total * 1e4) / medians
        nonzero = np.nonzero(normalized)[0]
        order = nonzero[np.argsort(-normalized[nonzero], kind="stable")]
        if order.size > max_tokens:
            truncated += 1
            order = order[:max_tokens]
        encoded.append([int(value) for value in tokens[order]])
    summary = {
        "schema_version": "masld-bench-geneformer-tokenization-v1",
        "embedding_policy": EMBEDDING_POLICY,
        "genes_total": len(gene_ids),
        "genes_usable": len(usable),
        "genes_dropped": len(gene_ids) - len(usable),
        "max_tokens": max_tokens,
        "cells_truncated": truncated,
        "cells": len(encoded),
    }
    return encoded, summary


def _embed(bundle: Path, subdir: str, encoded, device: str, batch_size: int):
    import torch
    from transformers import AutoModel

    model = AutoModel.from_pretrained((bundle / "upstream" / subdir).as_posix())
    model.eval()
    model.to(device)
    outputs = []
    with torch.no_grad():
        for start in range(0, len(encoded), batch_size):
            chunk = encoded[start : start + batch_size]
            width = max(len(item) for item in chunk)
            ids = torch.zeros(len(chunk), width, dtype=torch.long)
            attention = torch.zeros(len(chunk), width, dtype=torch.long)
            for index, item in enumerate(chunk):
                ids[index, : len(item)] = torch.tensor(item, dtype=torch.long)
                attention[index, : len(item)] = 1
            hidden = model(
                input_ids=ids.to(device), attention_mask=attention.to(device)
            ).last_hidden_state
            mask = attention.to(device).unsqueeze(-1).to(hidden.dtype)
            pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1)
            outputs.append(pooled.detach().cpu())
    return torch.cat(outputs, dim=0)


def _build_head(head_id: str, dimensions: int, classes: int, seed: int):
    import torch
    from torch import nn

    torch.manual_seed(seed)
    if head_id == "linear":
        return nn.Linear(dimensions, classes)
    if head_id == "two_layer_mlp":
        return nn.Sequential(
            nn.Linear(dimensions, 256), nn.ReLU(), nn.Linear(256, classes)
        )
    raise ScientificAdapterError(f"unsupported common head: {head_id}")


def _validate_request(request: Mapping[str, Any], action: str):
    if request.get("schema_version") != "masld-bench-adapter-request-v1":
        raise ScientificAdapterError("unsupported adapter request schema")
    if request.get("action") != action:
        raise ScientificAdapterError("request action differs")
    _require_sha256(request.get("run_id"), "run_id")
    run_spec = request.get("run_spec")
    if not isinstance(run_spec, Mapping):
        raise ScientificAdapterError("request has no RunSpec")
    model_id = str(run_spec.get("model_id", ""))
    if model_id not in MODEL_IDS:
        raise ScientificAdapterError("RunSpec names an unsupported Geneformer model")
    for field, expected in (
        ("task_id", TASK_ID),
        ("split_id", "donor_outer"),
        ("adaptation_regime", "common_lane"),
        ("runtime_id", RUNTIME_ID),
    ):
        if run_spec.get(field) != expected:
            raise ScientificAdapterError(f"RunSpec {field} differs from {expected}")
    parameters = run_spec.get("hyperparameters")
    if not isinstance(parameters, Mapping):
        raise ScientificAdapterError("RunSpec carries no frozen hyperparameters")
    if parameters.get("embedding_policy") != EMBEDDING_POLICY:
        raise ScientificAdapterError("embedding policy differs from the frozen policy")
    head_id = str(parameters.get("common_head_id", ""))
    if head_id not in HEAD_IDS:
        raise ScientificAdapterError(f"unsupported common head: {head_id!r}")
    roster = parameters.get("cell_state_roster")
    if not isinstance(roster, (list, tuple)) or not roster:
        raise ScientificAdapterError("frozen cell-state roster is missing")
    roster = tuple(str(item) for item in roster)
    if list(roster) != sorted(set(roster)):
        raise ScientificAdapterError("cell-state roster must be sorted and unique")
    return parameters, roster, model_id, head_id


def _bundle_root(request: Mapping[str, Any], model_id: str) -> Path:
    role = f"model_checkpoint_bundle:{model_id}"
    manifest = _validate_file_artifact(_input_by_role(request, role), "bundle manifest")
    return manifest.parent


def _variant(model_id: str) -> tuple[str, str]:
    if model_id == "geneformer_v1_10m":
        return "V1", "Geneformer-V1-10M"
    if model_id == "geneformer_v2_104m":
        return "V2", "Geneformer-V2-104M"
    return "V2", "Geneformer-V2-316M"


def _receipt(*, action, request, output, artifacts, extra_metadata) -> None:
    receipt = {
        "schema_version": RECEIPT_SCHEMA,
        "action": action,
        "run_id": request["run_id"],
        "status": "complete",
        "artifacts": [_artifact_record(p, relative_to=output) for p in artifacts],
        "metadata": {
            "adapter": ADAPTER_ID,
            "runtime_id": RUNTIME_ID,
            "fit_dataset_ids": list(request["fit_dataset_ids"]),
            **dict(extra_metadata),
        },
    }
    _write_json(output / "adapter_receipt.json", receipt)


def _atlas_path(request: Mapping[str, Any]) -> Path:
    role = "dataset_view_data:resource_atlas_geneformer_smoke_1000_v1"
    return _validate_file_artifact(_input_by_role(request, role), "Atlas smoke H5AD")


def _split_rows(adata: Any, parameters: Mapping[str, Any], run_spec: Mapping[str, Any]):
    outer_folds = int(parameters["outer_folds"])
    fold = int(run_spec["fold"])
    if not 0 <= fold < outer_folds:
        raise ScientificAdapterError("outer-fold contract differs")
    namespace = str(parameters["join_namespace"])
    rows = []
    for row_id, donor in zip(
        map(str, adata.obs_names), map(str, adata.obs["donor_id"])
    ):
        assigned = fold_index(donor, seed=20260821, outer_folds=outer_folds)
        rows.append(
            {
                "row_hash": join_hash(namespace, "row", row_id),
                "unit_hash": join_hash(namespace, "unit", donor),
                "fold": assigned,
                "held_out": assigned == fold,
            }
        )
    rows.sort(key=lambda item: item["row_hash"])
    if len({r["row_hash"] for r in rows}) != len(rows):
        raise ScientificAdapterError("row identifiers collide")
    return rows


def _load_atlas(path: Path):
    import anndata

    adata = anndata.read_h5ad(path)
    if not {"donor_id", "broad_label"}.issubset(adata.obs.columns):
        raise ScientificAdapterError("Atlas H5AD lacks required observation fields")
    if len(set(map(str, adata.obs_names))) != adata.n_obs:
        raise ScientificAdapterError("Atlas row identifiers are not unique")
    return adata


def _counts(adata: Any):
    import numpy as np

    matrix = adata.X
    return matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)


def prepare(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    parameters, roster, model_id, _ = _validate_request(request, "prepare")
    adata = _load_atlas(_atlas_path(request))
    observed = set(map(str, adata.obs["broad_label"]))
    if observed != set(roster):
        raise ScientificAdapterError("Atlas label set differs from the frozen roster")
    rows = _split_rows(adata, parameters, request["run_spec"])
    split_path = output / "split_rows.tsv"
    _write_tsv(
        split_path,
        ("row_hash", "unit_hash", "fold", "held_out"),
        (
            {**row, "held_out": "true" if row["held_out"] else "false"}
            for row in rows
        ),
    )

    variant, subdir = _variant(model_id)
    bundle = _bundle_root(request, model_id)
    median, token = _load_dictionaries(bundle, variant)
    config = json.loads(
        (bundle / "upstream" / subdir / "config.json").read_text("utf-8")
    )
    max_tokens = int(config["max_position_embeddings"])
    if int(config["vocab_size"]) != len(token):
        raise ScientificAdapterError(
            "token dictionary size differs from the checkpoint vocabulary"
        )
    gene_ids = [str(g) for g in adata.var["ensembl_id"]]
    encoded, summary = _rank_value_encode(
        _counts(adata), gene_ids, median, token, max_tokens
    )
    order = {str(name): i for i, name in enumerate(adata.obs_names)}
    namespace = str(parameters["join_namespace"])
    by_hash = {
        join_hash(namespace, "row", name): encoded[index]
        for name, index in order.items()
    }
    tokens_path = output / "rank_value_tokens.jsonl"
    with tokens_path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(
                json.dumps(
                    {"row_hash": row["row_hash"], "tokens": by_hash[row["row_hash"]]},
                    sort_keys=True,
                )
                + "\n"
            )
    summary_path = output / "tokenization_summary.json"
    _write_json(summary_path, summary)
    _receipt(
        action="prepare",
        request=request,
        output=output,
        artifacts=[split_path, tokens_path, summary_path],
        extra_metadata={
            "model_id": model_id,
            "genes_usable": summary["genes_usable"],
            "genes_dropped": summary["genes_dropped"],
            "cells_truncated": summary["cells_truncated"],
        },
    )


def _read_tokens(path: Path) -> dict[str, list[int]]:
    result: dict[str, list[int]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            result[str(record["row_hash"])] = [int(v) for v in record["tokens"]]
    return result


def _read_split(path: Path) -> list[dict[str, Any]]:
    import csv

    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [
        {
            "row_hash": r["row_hash"],
            "unit_hash": r["unit_hash"],
            "held_out": r["held_out"] == "true",
        }
        for r in rows
    ]


def _prior_output(request: Mapping[str, Any], action: str) -> Path:
    outputs = request.get("prior_action_outputs")
    if not isinstance(outputs, Mapping) or action not in outputs:
        raise ScientificAdapterError(f"request lacks the frozen {action} output")
    return Path(str(outputs[action])).resolve(strict=True)


def _labels_by_hash(request, parameters, roster):
    adata = _load_atlas(_atlas_path(request))
    namespace = str(parameters["join_namespace"])
    index = {label: position for position, label in enumerate(roster)}
    return {
        join_hash(namespace, "row", str(name)): index[str(label)]
        for name, label in zip(adata.obs_names, adata.obs["broad_label"])
    }


def fit(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import torch

    parameters, roster, model_id, head_id = _validate_request(request, "fit")
    prepare_dir = _prior_output(request, "prepare")
    split = _read_split(prepare_dir / "split_rows.tsv")
    tokens = _read_tokens(prepare_dir / "rank_value_tokens.jsonl")
    labels = _labels_by_hash(request, parameters, roster)

    # Only outer-training donors may influence any fitted state.
    training = [row for row in split if not row["held_out"]]
    if not training:
        raise ScientificAdapterError("fold has no training donors")
    held_units = {row["unit_hash"] for row in split if row["held_out"]}
    if held_units & {row["unit_hash"] for row in training}:
        raise ScientificAdapterError("a donor appears in both training and held-out")

    variant, subdir = _variant(model_id)
    bundle = _bundle_root(request, model_id)
    device = _select_device(str(parameters.get("device", "auto")))
    seed = int(request["run_spec"]["seed"])
    torch.manual_seed(seed)

    embeddings = _embed(
        bundle,
        subdir,
        [tokens[row["row_hash"]] for row in training],
        device,
        int(parameters.get("embedding_batch_size", 8)),
    )
    targets = torch.tensor(
        [labels[row["row_hash"]] for row in training], dtype=torch.long
    )

    # Donor- and class-balanced weighting, matching the classical lane.
    from collections import Counter

    donor_counts = Counter(row["unit_hash"] for row in training)
    class_counts = Counter(int(t) for t in targets)
    weights = torch.tensor(
        [
            1.0
            / (donor_counts[row["unit_hash"]] * class_counts[int(target)])
            for row, target in zip(training, targets)
        ],
        dtype=torch.float32,
    )
    weights = weights / weights.sum() * len(weights)

    head = _build_head(head_id, embeddings.shape[1], len(roster), seed).to(device)
    optimizer = torch.optim.Adam(
        head.parameters(), lr=float(parameters.get("head_learning_rate", 1e-3))
    )
    epochs = int(parameters.get("head_epochs", 60))
    features = embeddings.to(device)
    target_device = targets.to(device)
    weight_device = weights.to(device)
    losses = []
    for _ in range(epochs):
        optimizer.zero_grad()
        logits = head(features)
        per_row = torch.nn.functional.cross_entropy(
            logits, target_device, reduction="none"
        )
        loss = (per_row * weight_device).mean()
        loss.backward()
        optimizer.step()
        losses.append(float(loss.detach().cpu()))

    state_path = output / "head_state.pt"
    torch.save(
        {
            "head_id": head_id,
            "state_dict": {k: v.cpu() for k, v in head.state_dict().items()},
            "embedding_dimensions": int(embeddings.shape[1]),
            "roster": list(roster),
            "seed": seed,
        },
        state_path,
    )
    summary_path = output / "fit_summary.json"
    _write_json(
        summary_path,
        {
            "schema_version": "masld-bench-geneformer-fit-v1",
            "model_id": model_id,
            "common_head_id": head_id,
            "embedding_policy": EMBEDDING_POLICY,
            "encoder_frozen": True,
            "device": device,
            "training_rows": len(training),
            "training_donors": len(donor_counts),
            "held_out_donors": len(held_units),
            "epochs": epochs,
            "first_loss": losses[0],
            "final_loss": losses[-1],
        },
    )
    _receipt(
        action="fit",
        request=request,
        output=output,
        artifacts=[state_path, summary_path],
        extra_metadata={
            "model_id": model_id,
            "common_head_id": head_id,
            "encoder_frozen": True,
            "device": device,
            "training_donors": len(donor_counts),
        },
    )


def predict(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import torch

    parameters, roster, model_id, head_id = _validate_request(request, "predict")
    prepare_dir = _prior_output(request, "prepare")
    fit_dir = _prior_output(request, "fit")
    split = _read_split(prepare_dir / "split_rows.tsv")
    tokens = _read_tokens(prepare_dir / "rank_value_tokens.jsonl")

    held = [row for row in split if row["held_out"]]
    if not held:
        raise ScientificAdapterError("fold has no held-out rows to predict")

    state = torch.load(fit_dir / "head_state.pt", map_location="cpu", weights_only=True)
    if state["head_id"] != head_id or list(state["roster"]) != list(roster):
        raise ScientificAdapterError("fitted head does not match this run's contract")

    variant, subdir = _variant(model_id)
    bundle = _bundle_root(request, model_id)
    device = _select_device(str(parameters.get("device", "auto")))
    embeddings = _embed(
        bundle,
        subdir,
        [tokens[row["row_hash"]] for row in held],
        device,
        int(parameters.get("embedding_batch_size", 8)),
    )
    if embeddings.shape[1] != int(state["embedding_dimensions"]):
        raise ScientificAdapterError("embedding width differs from the fitted head")

    head = _build_head(head_id, embeddings.shape[1], len(roster), int(state["seed"]))
    head.load_state_dict(state["state_dict"])
    head.eval()
    with torch.no_grad():
        probabilities = torch.softmax(head(embeddings), dim=1).cpu()

    prediction_rows = []
    probability_rows = []
    for position, row in enumerate(held):
        values = [float(v) for v in probabilities[position]]
        # Deterministic tie-break identical to the sealed evaluator.
        winner = min(range(len(roster)), key=lambda i: (-values[i], roster[i]))
        prediction_rows.append(
            {
                "row_hash": row["row_hash"],
                "unit_hash": row["unit_hash"],
                "predicted_class": roster[winner],
            }
        )
        probability_rows.append(
            {
                "row_hash": row["row_hash"],
                "unit_hash": row["unit_hash"],
                **{
                    f"probability::{label}": repr(values[index])
                    for index, label in enumerate(roster)
                },
            }
        )
    prediction_rows.sort(key=lambda item: item["row_hash"])
    probability_rows.sort(key=lambda item: item["row_hash"])

    predictions_path = output / "predictions.tsv"
    _write_tsv(
        predictions_path, ("row_hash", "unit_hash", "predicted_class"), prediction_rows
    )
    row_ids_path = output / "row_ids.tsv"
    _write_tsv(
        row_ids_path,
        ("row_hash", "unit_hash"),
        ({"row_hash": r["row_hash"], "unit_hash": r["unit_hash"]} for r in prediction_rows),
    )
    probabilities_path = output / "class_probabilities.tsv"
    _write_tsv(
        probabilities_path,
        ("row_hash", "unit_hash", *(f"probability::{label}" for label in roster)),
        probability_rows,
    )
    _receipt(
        action="predict",
        request=request,
        output=output,
        artifacts=[predictions_path, row_ids_path, probabilities_path],
        extra_metadata={
            "model_id": model_id,
            "common_head_id": head_id,
            "device": device,
            "predicted_rows": len(prediction_rows),
            "class_roster": list(roster),
        },
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--action", required=True, choices=("prepare", "fit", "predict"))
    parser.add_argument("--request", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args(argv)
    try:
        request = json.loads(arguments.request.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScientificAdapterError("adapter request is invalid JSON") from error
    if not isinstance(request, Mapping):
        raise ScientificAdapterError("adapter request must be an object")
    if arguments.output.exists():
        raise ScientificAdapterError(f"adapter output already exists: {arguments.output}")
    arguments.output.mkdir(parents=True, exist_ok=False)
    {"prepare": prepare, "fit": fit, "predict": predict}[arguments.action](
        arguments.request.resolve(strict=True), request, arguments.output
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
