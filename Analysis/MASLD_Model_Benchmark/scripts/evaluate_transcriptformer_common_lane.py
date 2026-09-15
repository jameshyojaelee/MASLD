#!/usr/bin/env python
"""Evaluate frozen cell-foundation embeddings with the common cell-state heads.

This is a development-only evaluator. It joins labels only after embedding
extraction, keeps every donor in one outer and inner role, fits every transform
inside its inner training partition, and never reads held-back outcomes.
"""

from __future__ import annotations

import argparse
import copy
import csv
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import random
from typing import Any, Iterable, Mapping, Sequence


ROSTER = (
    "cholangiocyte",
    "endothelial",
    "hepatocyte",
    "immune",
    "mesenchymal_stromal",
)
HEADS = ("linear", "two_layer_mlp")
BASELINES = (
    "hvg_pca_nearest_centroid",
    "hvg_pca_knn",
    "hvg_pca_logistic",
    "hvg_pca_elastic_net",
    "hvg_pca_linear_svm",
)
EXISTING_CELL_CANDIDATES = (
    "geneformer_v1_10m",
    "geneformer_v2_104m",
    "geneformer_v2_316m",
    "scimilarity_v1_1",
    "scgpt_continual",
    "uce_4l",
)
NAMESPACE = "resource_atlas_current:cell_state_mapping:donor_outer:v1"
OUTER_SEED = 20260821


class EvaluationError(RuntimeError):
    """Raised when the development evaluation requirements are not met."""


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _join_hash(kind: str, identifier: str) -> str:
    return sha256(f"{NAMESPACE}\0{kind}\0{identifier}".encode()).hexdigest()


def _fold_index(unit_id: str, folds: int = 5) -> int:
    digest = sha256(f"{OUTER_SEED}\0{unit_id}".encode()).digest()
    return int.from_bytes(digest[:8], "big") % folds


def _inner_assignment(donors: Sequence[str], outer_fold: int) -> dict[str, int]:
    unique = sorted(
        set(donors),
        key=lambda donor: sha256(
            f"1103\0{outer_fold}\0inner\0{donor}".encode()
        ).hexdigest(),
    )
    if len(unique) < 5:
        raise EvaluationError("outer training partition has fewer than five donors")
    return {donor: index % 5 for index, donor in enumerate(unique)}


def _donor_class_weights(donors, labels, *, require_full_roster: bool = True):
    import numpy as np

    counts: dict[tuple[str, str], int] = {}
    donors_by_class: dict[str, set[str]] = {}
    for donor, label in zip(donors, labels, strict=True):
        key = (str(donor), str(label))
        counts[key] = counts.get(key, 0) + 1
        donors_by_class.setdefault(str(label), set()).add(str(donor))
    if require_full_roster and len(donors_by_class) != len(ROSTER):
        raise EvaluationError("training partition lacks a frozen class")
    class_count = len(donors_by_class)
    weights = [
        1.0
        / (
            class_count
            * len(donors_by_class[str(label)])
            * counts[(str(donor), str(label))]
        )
        for donor, label in zip(donors, labels, strict=True)
    ]
    values = np.asarray(weights, dtype=np.float32)
    return values / values.sum() * len(values)


def _build_head(head_id: str, width: int, classes: int):
    from torch import nn

    if head_id == "linear":
        return nn.Linear(width, classes)
    if head_id == "two_layer_mlp":
        return nn.Sequential(
            nn.Linear(width, min(256, width)),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(min(256, width), classes),
        )
    raise EvaluationError(f"unknown common head: {head_id}")


def _weighted_cross_entropy(logits, targets, weights):
    import torch

    losses = torch.nn.functional.cross_entropy(logits, targets, reduction="none")
    return (losses * weights).sum() / weights.sum()


def _fit_inner_head(
    *,
    head_id: str,
    features,
    targets,
    donors,
    fit_indices,
    validation_indices,
    seed: int,
):
    import numpy as np
    import torch

    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.use_deterministic_algorithms(True)
    torch.set_num_threads(max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))))

    fit_x = features[fit_indices].astype(np.float32, copy=False)
    validation_x = features[validation_indices].astype(np.float32, copy=False)
    mean = fit_x.mean(axis=0, dtype=np.float64).astype(np.float32)
    scale = fit_x.std(axis=0, dtype=np.float64).astype(np.float32)
    scale[scale == 0.0] = 1.0
    fit_x = (fit_x - mean) / scale
    validation_x = (validation_x - mean) / scale

    fit_y = targets[fit_indices]
    validation_y = targets[validation_indices]
    fit_donors = donors[fit_indices]
    validation_donors = donors[validation_indices]
    fit_w = _donor_class_weights(fit_donors, fit_y, require_full_roster=True)
    validation_w = _donor_class_weights(
        validation_donors, validation_y, require_full_roster=False
    )

    model = _build_head(head_id, features.shape[1], len(ROSTER))
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=0.001, weight_decay=0.01
    )
    x = torch.from_numpy(fit_x)
    y = torch.from_numpy(fit_y.astype(np.int64, copy=False))
    w = torch.from_numpy(fit_w)
    val_x = torch.from_numpy(validation_x)
    val_y = torch.from_numpy(validation_y.astype(np.int64, copy=False))
    val_w = torch.from_numpy(validation_w)
    generator = torch.Generator().manual_seed(seed)
    best_loss = math.inf
    best_epoch = -1
    best_state = None
    stale = 0
    for epoch in range(200):
        model.train()
        order = torch.randperm(len(x), generator=generator)
        for start in range(0, len(x), 128):
            batch = order[start : start + 128]
            optimizer.zero_grad(set_to_none=True)
            loss = _weighted_cross_entropy(model(x[batch]), y[batch], w[batch])
            loss.backward()
            optimizer.step()
        model.eval()
        with torch.no_grad():
            validation_loss = float(
                _weighted_cross_entropy(model(val_x), val_y, val_w)
            )
        if validation_loss < best_loss - 0.0001:
            best_loss = validation_loss
            best_epoch = epoch + 1
            best_state = copy.deepcopy(model.state_dict())
            stale = 0
        else:
            stale += 1
        if stale >= 15:
            break
    if best_state is None:
        raise EvaluationError("inner head produced no accepted epoch")
    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        validation_logits = model(val_x).numpy()
    return model, mean, scale, validation_logits, {
        "best_epoch": best_epoch,
        "best_validation_weighted_cross_entropy": best_loss,
        "fit_rows": len(fit_indices),
        "validation_rows": len(validation_indices),
        "fit_donors": len(set(map(str, fit_donors))),
        "validation_donors": len(set(map(str, validation_donors))),
    }


def _fit_temperature(logits, targets, weights) -> float:
    import torch

    values = torch.tensor(logits, dtype=torch.float64)
    outcomes = torch.tensor(targets, dtype=torch.long)
    sample_weights = torch.tensor(weights, dtype=torch.float64)
    log_temperature = torch.tensor(0.0, dtype=torch.float64, requires_grad=True)
    optimizer = torch.optim.LBFGS(
        [log_temperature], lr=0.1, max_iter=100, line_search_fn="strong_wolfe"
    )

    def closure():
        optimizer.zero_grad()
        temperature = log_temperature.exp().clamp(0.05, 20.0)
        loss = _weighted_cross_entropy(
            values / temperature, outcomes, sample_weights
        )
        loss.backward()
        return loss

    optimizer.step(closure)
    temperature = float(log_temperature.exp().clamp(0.05, 20.0).detach())
    if not math.isfinite(temperature) or temperature <= 0:
        raise EvaluationError("temperature calibration failed")
    return temperature


def _train_predict_head(head_id, features, targets, donors, outer_folds):
    import numpy as np
    import torch

    probabilities = np.zeros((len(features), len(ROSTER)), dtype=np.float64)
    fold_records = []
    for outer_fold in range(5):
        outer_training = np.flatnonzero(outer_folds != outer_fold)
        outer_test = np.flatnonzero(outer_folds == outer_fold)
        if not len(outer_test):
            raise EvaluationError(f"outer fold {outer_fold} is empty")
        assignment = _inner_assignment(donors[outer_training], outer_fold)
        inner_models = []
        oof_logits = []
        oof_targets = []
        oof_donors = []
        inner_records = []
        for inner_fold in range(5):
            validation_mask = np.asarray(
                [assignment[str(donors[index])] == inner_fold for index in outer_training]
            )
            fit_indices = outer_training[~validation_mask]
            validation_indices = outer_training[validation_mask]
            seed = 1103 + outer_fold * 101 + inner_fold * 17 + HEADS.index(head_id)
            model, mean, scale, logits, record = _fit_inner_head(
                head_id=head_id,
                features=features,
                targets=targets,
                donors=donors,
                fit_indices=fit_indices,
                validation_indices=validation_indices,
                seed=seed,
            )
            inner_models.append((model, mean, scale))
            oof_logits.append(logits)
            oof_targets.append(targets[validation_indices])
            oof_donors.append(donors[validation_indices])
            inner_records.append({"inner_fold": inner_fold, "seed": seed, **record})
        all_logits = np.concatenate(oof_logits)
        all_targets = np.concatenate(oof_targets)
        all_donors = np.concatenate(oof_donors)
        calibration_weights = _donor_class_weights(all_donors, all_targets)
        temperature = _fit_temperature(
            all_logits, all_targets, calibration_weights
        )
        fold_probabilities = []
        for model, mean, scale in inner_models:
            held = torch.from_numpy(
                ((features[outer_test] - mean) / scale).astype(np.float32)
            )
            model.eval()
            with torch.no_grad():
                fold_probabilities.append(
                    torch.softmax(model(held) / temperature, dim=1).numpy()
                )
        probabilities[outer_test] = np.mean(fold_probabilities, axis=0)
        fold_records.append(
            {
                "outer_fold": outer_fold,
                "outer_training_rows": len(outer_training),
                "outer_training_donors": len(set(map(str, donors[outer_training]))),
                "outer_test_rows": len(outer_test),
                "outer_test_donors": len(set(map(str, donors[outer_test]))),
                "temperature": temperature,
                "inner_models": inner_records,
            }
        )
    if not np.isfinite(probabilities).all():
        raise EvaluationError("common-head probabilities are non-finite")
    if not np.allclose(probabilities.sum(axis=1), 1.0, atol=1e-6):
        raise EvaluationError("common-head probabilities do not sum to one")
    return probabilities, fold_records


def _read_table(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _load_baseline(root: Path, model_id: str, expected_rows: set[str]):
    prediction_paths = sorted(
        (root / model_id).glob("fold*/adapter_actions/003-predict/predictions.tsv")
    )
    probability_paths = sorted(
        (root / model_id).glob(
            "fold*/adapter_actions/003-predict/class_probabilities.tsv"
        )
    )
    if len(prediction_paths) != 5 or len(probability_paths) != 5:
        raise EvaluationError(f"baseline {model_id} lacks five frozen folds")
    predictions: dict[str, str] = {}
    probabilities: dict[str, dict[str, float]] = {}
    for path in prediction_paths:
        for row in _read_table(path):
            row_hash = row["row_hash"]
            if row_hash in predictions:
                raise EvaluationError(f"baseline {model_id} repeats a row")
            predictions[row_hash] = row["predicted_class"]
    for path in probability_paths:
        for row in _read_table(path):
            row_hash = row["row_hash"]
            if row_hash in probabilities:
                raise EvaluationError(f"baseline {model_id} repeats probabilities")
            probabilities[row_hash] = {
                label: float(row[f"probability::{label}"]) for label in ROSTER
            }
    if set(predictions) != expected_rows or set(probabilities) != expected_rows:
        raise EvaluationError(f"baseline {model_id} row universe differs")
    records = [
        {
            "path": str(path.resolve()),
            "sha256": _sha256_file(path),
            "size_bytes": path.stat().st_size,
        }
        for path in (*prediction_paths, *probability_paths)
    ]
    return predictions, probabilities, records


def _load_existing_candidate(
    root: Path, model_id: str, head_id: str, expected_rows: set[str]
):
    prediction_paths = sorted(
        root.glob(f"{model_id}*/{head_id}/fold*/predict/predictions.tsv")
    )
    probability_paths = sorted(
        root.glob(f"{model_id}*/{head_id}/fold*/predict/class_probabilities.tsv")
    )
    if len(prediction_paths) != 5 or len(probability_paths) != 5:
        raise EvaluationError(
            f"existing candidate {model_id}/{head_id} lacks five frozen folds"
        )
    predictions: dict[str, str] = {}
    probabilities: dict[str, dict[str, float]] = {}
    for path in prediction_paths:
        for row in _read_table(path):
            row_hash = row["row_hash"]
            if row_hash in predictions:
                raise EvaluationError(f"existing candidate {model_id} repeats a row")
            predictions[row_hash] = row["predicted_class"]
    for path in probability_paths:
        for row in _read_table(path):
            row_hash = row["row_hash"]
            if row_hash in probabilities:
                raise EvaluationError(
                    f"existing candidate {model_id} repeats probabilities"
                )
            probabilities[row_hash] = {
                label: float(row[f"probability::{label}"]) for label in ROSTER
            }
    if set(predictions) != expected_rows or set(probabilities) != expected_rows:
        raise EvaluationError(
            f"existing candidate {model_id}/{head_id} row universe differs"
        )
    records = [
        {
            "path": str(path.resolve()),
            "sha256": _sha256_file(path),
            "size_bytes": path.stat().st_size,
        }
        for path in (*prediction_paths, *probability_paths)
    ]
    return predictions, probabilities, records


def _score(labels, predictions, donors, probabilities):
    from masld_bench.evaluators.cell_state_development import score_rows

    return score_rows(
        observed=list(labels),
        predicted=list(predictions),
        donors=list(donors),
        probabilities=list(probabilities),
        class_roster=ROSTER,
    )


def _bootstrap_delta(
    *, labels, donors, candidate_predictions, candidate_probabilities,
    baseline_predictions, baseline_probabilities, seed: int, resamples: int,
):
    import numpy as np

    unit_roster = sorted(set(map(str, donors)))
    by_unit = {
        unit: np.flatnonzero(donors == unit).tolist() for unit in unit_roster
    }
    rng = np.random.default_rng(seed)
    f1_deltas = []
    brier_deltas = []
    for _ in range(resamples):
        indices = []
        sampled_donors = []
        for draw, unit in enumerate(
            rng.choice(unit_roster, size=len(unit_roster), replace=True)
        ):
            indices.extend(by_unit[str(unit)])
            sampled_donors.extend([f"{draw}:{unit}"] * len(by_unit[str(unit)]))
        candidate = _score(
            labels[indices], candidate_predictions[indices], sampled_donors,
            [candidate_probabilities[index] for index in indices],
        )
        baseline = _score(
            labels[indices], baseline_predictions[indices], sampled_donors,
            [baseline_probabilities[index] for index in indices],
        )
        f1_deltas.append(
            candidate["donor_class_balanced_macro_f1"]
            - baseline["donor_class_balanced_macro_f1"]
        )
        brier_deltas.append(
            candidate["multiclass_brier_score"]
            - baseline["multiclass_brier_score"]
        )
    return {
        "resamples": resamples,
        "seed": seed,
        "f1_delta_ci95": [
            float(np.quantile(f1_deltas, 0.025)),
            float(np.quantile(f1_deltas, 0.975)),
        ],
        "f1_delta_positive_fraction": float(np.mean(np.asarray(f1_deltas) > 0)),
        "brier_delta_ci95": [
            float(np.quantile(brier_deltas, 0.025)),
            float(np.quantile(brier_deltas, 0.975)),
        ],
        "f1_delta_sha256": sha256(
            json.dumps(f1_deltas, separators=(",", ":")).encode()
        ).hexdigest(),
        "brier_delta_sha256": sha256(
            json.dumps(brier_deltas, separators=(",", ":")).encode()
        ).hexdigest(),
    }


def _write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]):
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> int:
    import anndata
    import numpy as np

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings", required=True, type=Path)
    parser.add_argument("--atlas", required=True, type=Path)
    parser.add_argument("--baseline-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--bootstrap-resamples", type=int, default=10000)
    parser.add_argument("--model-id", default="transcriptformer_tf_sapiens")
    parser.add_argument("--representation-id", default="common")
    parser.add_argument("--expected-width", type=int, default=2048)
    parser.add_argument(
        "--row-id-source",
        type=Path,
        help="Optional frozen h5ad supplying row_id and outer_fold without loading an object array",
    )
    parser.add_argument(
        "--schema-version",
        default="masld-bench-transcriptformer-common-lane-evaluation-v1",
    )
    parser.add_argument(
        "--activation-fold-namespace",
        default="transcriptformer-resource-atlas-smoke-v1",
    )
    args = parser.parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=False)

    bundle = np.load(args.embeddings, allow_pickle=False)
    features = np.asarray(bundle["embeddings"], dtype=np.float32)
    activation_folds = np.asarray(bundle["outer_folds"], dtype=np.int64)
    if args.row_id_source is None:
        row_ids = np.asarray(bundle["row_ids"]).astype(str)
    else:
        row_source = anndata.read_h5ad(args.row_id_source, backed="r")
        if "row_id" not in row_source.obs or "outer_fold" not in row_source.obs:
            raise EvaluationError("row-ID source columns differ")
        row_ids = row_source.obs["row_id"].astype(str).to_numpy()
        source_folds = row_source.obs["outer_fold"].astype(int).to_numpy()
        row_source.file.close()
        if not np.array_equal(activation_folds, source_folds):
            raise EvaluationError("embedding and row-ID source orders differ")
    if features.shape != (1000, args.expected_width) or len(set(row_ids)) != 1000:
        raise EvaluationError("cell-foundation embedding shape or row IDs differ")
    if not np.isfinite(features).all() or set(activation_folds) != set(range(5)):
        raise EvaluationError("cell-foundation embeddings or folds differ")

    atlas = anndata.read_h5ad(args.atlas, backed="r")
    atlas_index = {str(row_id): index for index, row_id in enumerate(atlas.obs_names)}
    if set(row_ids) != set(atlas_index):
        raise EvaluationError("embedding and Atlas row universes differ")
    donors = np.asarray(
        [str(atlas.obs.iloc[atlas_index[row_id]]["donor_id"]) for row_id in row_ids]
    )
    labels = np.asarray(
        [str(atlas.obs.iloc[atlas_index[row_id]]["broad_label"]) for row_id in row_ids]
    )
    atlas.file.close()
    if set(labels) != set(ROSTER):
        raise EvaluationError("development labels differ from the frozen roster")
    outer_folds = np.asarray([_fold_index(donor) for donor in donors])
    if any(
        len(set(activation_folds[donors == donor])) != 1
        or len(set(outer_folds[donors == donor])) != 1
        for donor in set(donors)
    ):
        raise EvaluationError("a donor spans activation or evaluation folds")

    class_index = {label: index for index, label in enumerate(ROSTER)}
    targets = np.asarray([class_index[label] for label in labels], dtype=np.int64)
    head_probabilities = {}
    head_metrics = {}
    head_fold_records = {}
    for head_id in HEADS:
        probabilities, records = _train_predict_head(
            head_id, features, targets, donors, outer_folds
        )
        predicted = np.asarray(
            [ROSTER[int(np.argmax(row))] for row in probabilities]
        )
        probability_rows = [
            {label: float(row[index]) for index, label in enumerate(ROSTER)}
            for row in probabilities
        ]
        head_probabilities[head_id] = probabilities
        head_metrics[head_id] = _score(
            labels, predicted, donors, probability_rows
        )
        head_fold_records[head_id] = records
        _write_tsv(
            args.output / f"{head_id}_predictions.tsv",
            ("row_hash", "unit_hash", "predicted_class"),
            (
                {
                    "row_hash": _join_hash("row", row_id),
                    "unit_hash": _join_hash("unit", donor),
                    "predicted_class": prediction,
                }
                for row_id, donor, prediction in zip(
                    row_ids, donors, predicted, strict=True
                )
            ),
        )
        _write_tsv(
            args.output / f"{head_id}_class_probabilities.tsv",
            (
                "row_hash", "unit_hash",
                *(f"probability::{label}" for label in ROSTER),
            ),
            (
                {
                    "row_hash": _join_hash("row", row_id),
                    "unit_hash": _join_hash("unit", donor),
                    **{
                        f"probability::{label}": repr(float(row[index]))
                        for index, label in enumerate(ROSTER)
                    },
                }
                for row_id, donor, row in zip(
                    row_ids, donors, probabilities, strict=True
                )
            ),
        )

    row_hashes = [_join_hash("row", row_id) for row_id in row_ids]
    expected_rows = set(row_hashes)
    baseline_metrics = {}
    baseline_data = {}
    baseline_inputs = []
    for model_id in BASELINES:
        predictions, probabilities, records = _load_baseline(
            args.baseline_root, model_id, expected_rows
        )
        prediction_vector = np.asarray([predictions[row] for row in row_hashes])
        probability_vector = [probabilities[row] for row in row_hashes]
        baseline_metrics[model_id] = _score(
            labels, prediction_vector, donors, probability_vector
        )
        baseline_data[model_id] = (prediction_vector, probability_vector)
        baseline_inputs.extend({"model_id": model_id, **record} for record in records)

    existing_candidate_metrics = {}
    existing_candidate_inputs = []
    for model_id in EXISTING_CELL_CANDIDATES:
        for head_id in HEADS:
            comparison_id = f"{model_id}::{head_id}"
            predictions, probabilities, records = _load_existing_candidate(
                args.baseline_root, model_id, head_id, expected_rows
            )
            prediction_vector = np.asarray([predictions[row] for row in row_hashes])
            probability_vector = [probabilities[row] for row in row_hashes]
            existing_candidate_metrics[comparison_id] = _score(
                labels, prediction_vector, donors, probability_vector
            )
            existing_candidate_inputs.extend(
                {"comparison_id": comparison_id, **record} for record in records
            )
    strongest = max(
        BASELINES,
        key=lambda model_id: (
            baseline_metrics[model_id]["donor_class_balanced_macro_f1"],
            -baseline_metrics[model_id]["multiclass_brier_score"],
            model_id,
        ),
    )
    comparisons = {}
    baseline_predictions, baseline_probabilities = baseline_data[strongest]
    for head_id in HEADS:
        probabilities = head_probabilities[head_id]
        predictions = np.asarray(
            [ROSTER[int(np.argmax(row))] for row in probabilities]
        )
        probability_rows = [
            {label: float(row[index]) for index, label in enumerate(ROSTER)}
            for row in probabilities
        ]
        comparisons[head_id] = {
            "strongest_registered_baseline": strongest,
            "delta_donor_class_balanced_macro_f1": (
                head_metrics[head_id]["donor_class_balanced_macro_f1"]
                - baseline_metrics[strongest]["donor_class_balanced_macro_f1"]
            ),
            "delta_multiclass_brier_score": (
                head_metrics[head_id]["multiclass_brier_score"]
                - baseline_metrics[strongest]["multiclass_brier_score"]
            ),
            "paired_donor_bootstrap": _bootstrap_delta(
                labels=labels,
                donors=donors,
                candidate_predictions=predictions,
                candidate_probabilities=probability_rows,
                baseline_predictions=baseline_predictions,
                baseline_probabilities=baseline_probabilities,
                seed=20260824 + HEADS.index(head_id),
                resamples=args.bootstrap_resamples,
            ),
            "delta_vs_existing_cell_candidates": {
                comparison_id: {
                    "delta_donor_class_balanced_macro_f1": (
                        head_metrics[head_id]["donor_class_balanced_macro_f1"]
                        - metrics["donor_class_balanced_macro_f1"]
                    ),
                    "delta_multiclass_brier_score": (
                        head_metrics[head_id]["multiclass_brier_score"]
                        - metrics["multiclass_brier_score"]
                    ),
                }
                for comparison_id, metrics in existing_candidate_metrics.items()
            },
        }

    result = {
        "schema_version": args.schema_version,
        "status": "pass_development_smoke",
        "task_id": "cell_state_mapping",
        "model_id": args.model_id,
        "representation_id": args.representation_id,
        "dataset_view_id": "resource_atlas_geneformer_smoke_1000_v1",
        "development_rows": len(row_ids),
        "development_donors": len(set(donors)),
        "outer_folds": 5,
        "inner_folds": 5,
        "activation_fold_namespace": args.activation_fold_namespace,
        "evaluation_fold_namespace": "resource_atlas_current:cell_state_mapping:donor_outer:v1",
        "activation_fold_annotations_used_for_head_fitting": False,
        "fold_annotation_mismatch_rows": int(
            np.sum(activation_folds != outer_folds)
        ),
        "fold_annotation_note": "The outcome-blind activation hash is donor-safe but differs from the shared benchmark hash. Frozen per-cell embeddings contain no fold-fitted state, so all common heads use the shared registered donor_outer folds for exact comparison.",
        "seed": 1103,
        "heads": head_metrics,
        "head_fold_records": head_fold_records,
        "registered_classical_baselines": baseline_metrics,
        "existing_cell_foundation_common_lane": existing_candidate_metrics,
        "existing_candidate_comparison_scope": "retrospective exact frozen prediction files; no sealed claim and no assertion that earlier head implementations satisfy the finalized nested common-head contract",
        "comparisons": comparisons,
        "development_only": True,
        "smoke_only": True,
        "exposure_status": "unknown",
        "sealed_champion_eligible": False,
        "evaluation_labels_read": ["broad_label"],
        "sealed_outcomes_read": False,
        "label_join_stage": "after_outcome_blind_embedding_extraction",
        "baseline_inputs": baseline_inputs,
        "existing_candidate_inputs": existing_candidate_inputs,
        "input_artifacts": {
            "embeddings_sha256": _sha256_file(args.embeddings),
            "atlas_sha256": _sha256_file(args.atlas),
        },
    }
    (args.output / "evaluation.json").write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": result["status"],
        "head_metrics": head_metrics,
        "comparisons": comparisons,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
