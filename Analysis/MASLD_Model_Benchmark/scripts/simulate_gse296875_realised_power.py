#!/usr/bin/env python3
"""Realised power of the frozen GSE296875 pipeline against a planted effect.

The information ceiling says what the metric could detect if a perfect predictor
were handed to the evaluator. This says what is left after the predictor has to
be *learned* from 39 donors by the frozen pipeline. The gap between the two is
the price of estimation.

Non-circularity. The effect is planted at a stated size and the detection rate
is measured. No observed model performance is read: this script never opens the
prediction table, the score output file, or any reported metric. The frozen inputs
it does use are design, not result — the donor folds, the endpoint vectors and
their denominators, and the real donor-by-lineage expression matrix, which is
used precisely so the simulation inherits the true feature covariance and the
true n-to-p ratio rather than a flattering synthetic substitute.

The one modelling choice not fixed by the frozen design is how much donor-level
expression variance the planted biological program occupies. It is swept at two
prespecified values so its influence is visible rather than buried.
"""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression, Ridge


class PowerSimulationError(RuntimeError):
    """Raised when a frozen input does not verify."""


ALPHA_GRID = (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)
C_GRID = (0.001, 0.01, 0.1, 1.0, 10.0)
DETECTION_FRACTION = 0.5
N_HIGHLY_VARIABLE = 2_000
MAX_COMPONENTS = 20

# Hepatocyte is the deepest lineage: all 39 donor units clear the frozen
# minimum-cell threshold, so the realised curve is not confounded by imputed
# blocks. It is the most favourable scope, which makes this an upper bound on
# realised power across scopes.
REPRESENTATIVE_LINEAGE = "hepatocyte"

BOUND = {
    "phenotype_endpoints": "6518494a7249cb1eac81fbad75989e7dfe0269b611254d5097736b9f5b1da140",
    "donor_folds": "e6bc9161acaab508c298c49036af1ad10a720949178a14cba42282933195d416",
}


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def average_precision(labels: np.ndarray, scores: np.ndarray) -> float:
    order = np.argsort(-scores, kind="stable")
    ranked = labels[order]
    values = scores[order]
    positives = int(ranked.sum())
    if positives == 0:
        raise PowerSimulationError("average precision needs a positive")
    true_positive = 0
    false_positive = 0
    previous_recall = 0.0
    area = 0.0
    index = 0
    while index < len(ranked):
        stop = index
        while stop < len(ranked) and values[stop] == values[index]:
            stop += 1
        group = ranked[index:stop]
        true_positive += int(group.sum())
        false_positive += len(group) - int(group.sum())
        recall = true_positive / positives
        precision = true_positive / (true_positive + false_positive)
        area += (recall - previous_recall) * precision
        previous_recall = recall
        index = stop
    return float(area)


def spearman(left: np.ndarray, right: np.ndarray) -> float:
    def ranks(values: np.ndarray) -> np.ndarray:
        order = np.argsort(values, kind="stable")
        output = np.empty(len(values), dtype=np.float64)
        index = 0
        while index < len(values):
            stop = index
            while stop < len(values) and values[order[stop]] == values[order[index]]:
                stop += 1
            output[order[index:stop]] = (index + stop + 1) / 2.0
            index = stop
        return output

    a = ranks(left)
    b = ranks(right)
    a = a - a.mean()
    b = b - b.mean()
    denominator = np.sqrt((a * a).sum() * (b * b).sum())
    if denominator == 0:
        return 0.0
    return float((a * b).sum() / denominator)


def transform_block(
    expression: np.ndarray, train_rows: np.ndarray, all_rows: np.ndarray
) -> np.ndarray:
    """The frozen transform chain, fitted on training rows only."""

    fit = expression[train_rows]
    detected = (fit > 0).mean(axis=0) >= DETECTION_FRACTION
    if detected.sum() < 10:
        detected = (fit > 0).mean(axis=0) > 0
    fit = fit[:, detected]
    variance = fit.var(axis=0)
    keep = np.argsort(-variance)[: min(N_HIGHLY_VARIABLE, fit.shape[1])]
    keep.sort()
    fit = fit[:, keep]
    mean = fit.mean(axis=0)
    deviation = fit.std(axis=0)
    deviation[deviation == 0.0] = 1.0
    components = min(MAX_COMPONENTS, train_rows.size - 1, fit.shape[1])
    pca = PCA(n_components=components, svd_solver="full", random_state=0)
    pca.fit((fit - mean) / deviation)
    applied = expression[all_rows][:, detected][:, keep]
    return pca.transform((applied - mean) / deviation)


def fit_predict(
    endpoint: str,
    features_train: np.ndarray,
    labels_train: np.ndarray,
    features_test: np.ndarray,
    penalty: float,
) -> np.ndarray:
    if endpoint == "steatosis":
        model = Ridge(alpha=penalty, random_state=0)
        model.fit(features_train, labels_train)
        return model.predict(features_test)
    model = LogisticRegression(
        C=penalty, penalty="l2", solver="lbfgs", max_iter=2_000, random_state=0
    )
    model.fit(features_train, labels_train)
    return model.predict_proba(features_test)[:, 1]


def cross_fitted_predictions(
    endpoint: str,
    expression: np.ndarray,
    labels: np.ndarray,
    folds: np.ndarray,
) -> np.ndarray:
    """Pooled out-of-fold predictions under the frozen donor-grouped folds."""

    grid = ALPHA_GRID if endpoint == "steatosis" else C_GRID
    out = np.empty(len(labels), dtype=np.float64)
    for fold in np.unique(folds):
        test = np.flatnonzero(folds == fold)
        train = np.flatnonzero(folds != fold)
        if endpoint == "fibrosis" and len(np.unique(labels[train])) < 2:
            raise PowerSimulationError("a training fold has one class only")
        # Nested leave-one-donor-out penalty selection inside the training fold.
        totals = np.zeros(len(grid))
        for position in range(train.size):
            inner_train = np.delete(train, position)
            held = train[position : position + 1]
            inner_labels = labels[inner_train]
            if endpoint == "fibrosis" and len(np.unique(inner_labels)) < 2:
                continue
            combined = np.concatenate([inner_train, held])
            block = transform_block(expression, inner_train, combined)
            for index, penalty in enumerate(grid):
                predicted = fit_predict(
                    endpoint,
                    block[: inner_train.size],
                    inner_labels,
                    block[inner_train.size :],
                    penalty,
                )
                truth = labels[train[position]]
                if endpoint == "steatosis":
                    totals[index] += (truth - predicted[0]) ** 2
                else:
                    clipped = min(max(float(predicted[0]), 1e-9), 1 - 1e-9)
                    totals[index] += -(
                        truth * np.log(clipped) + (1 - truth) * np.log(1 - clipped)
                    )
        penalty = grid[int(np.argmin(totals))]
        combined = np.concatenate([train, test])
        block = transform_block(expression, train, combined)
        out[test] = fit_predict(
            endpoint, block[: train.size], labels[train], block[train.size :], penalty
        )
    return out


def plant(
    base: np.ndarray,
    latent: np.ndarray,
    variance_fraction: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Add a donor-level program of stated variance share along a random axis."""

    loading = rng.normal(size=base.shape[1])
    loading /= np.linalg.norm(loading)
    centred = latent - latent.mean()
    if centred.std() > 0:
        centred = centred / centred.std()
    signal = np.outer(centred, loading)
    base_scale = base.var(axis=0).sum()
    signal_scale = signal.var(axis=0).sum()
    if signal_scale <= 0:
        raise PowerSimulationError("planted signal has no variance")
    factor = np.sqrt(variance_fraction * base_scale / signal_scale)
    return base + factor * signal


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--pseudobulk", type=Path, required=True)
    parser.add_argument("--pseudobulk-sha256", required=True)
    parser.add_argument("--endpoint", choices=("steatosis", "fibrosis"), required=True)
    parser.add_argument("--variance-fraction", type=float, required=True)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260825)
    parser.add_argument("--critical-nominal", type=float, required=True)
    parser.add_argument("--critical-bh", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    root = arguments.root
    for name, expected in BOUND.items():
        folder = {
            "phenotype_endpoints": "executions/gse296875-phenotype-endpoints-20260825",
            "donor_folds": "executions/gse296875-donor-folds-20260825",
        }[name]
        if digest(root / folder / "ARTIFACTS.json") != expected:
            raise PowerSimulationError(f"bound fixture changed: {name}")
    if digest(arguments.pseudobulk / "ARTIFACTS.json") != arguments.pseudobulk_sha256:
        raise PowerSimulationError("pseudobulk fixture changed")

    units = read_tsv(arguments.pseudobulk / "pseudobulk" / "units.tsv")
    with np.load(arguments.pseudobulk / "pseudobulk" / "donor_lineage_counts.npz") as p:
        counts = p["counts"]
    endpoint_rows = read_tsv(
        root
        / "executions/gse296875-phenotype-endpoints-20260825/endpoint_lock/donor_endpoints.tsv"
    )
    fold_rows = read_tsv(
        root / "executions/gse296875-donor-folds-20260825/split_lock/donor_folds.tsv"
    )

    if arguments.endpoint == "steatosis":
        observed = {
            row["donor_id"]: float(row["steatosis_numeric"])
            for row in endpoint_rows
            if row["steatosis_observed"] == "true"
        }
        grid = (0.0, 0.2, 0.35, 0.5, 0.65, 0.8)
    else:
        observed = {
            row["donor_id"]: 1.0 if row["fibrosis_any"] == "true" else 0.0
            for row in endpoint_rows
            if row["fibrosis_observed"] == "true"
        }
        grid = (0.0, 0.4, 0.8, 1.2, 1.6, 2.0)

    donors = [row["donor_id"] for row in fold_rows if row["donor_id"] in observed]
    fold_of = {row["donor_id"]: int(row["outer_fold"]) for row in fold_rows}
    folds = np.asarray([fold_of[d] for d in donors], dtype=np.int64)
    labels = np.asarray([observed[d] for d in donors], dtype=np.float64)

    row_of = {
        (row["donor_id"], row["lineage_id"]): int(row["unit_index"]) for row in units
    }
    indices = np.asarray(
        [row_of[(d, REPRESENTATIVE_LINEAGE)] for d in donors], dtype=np.int64
    )
    raw = counts[indices].astype(np.float64)
    library = raw.sum(axis=1, keepdims=True)
    base = np.log1p(raw / library * 1e6)

    ranks = np.argsort(np.argsort(labels)).astype(np.float64)
    anchor = (ranks - ranks.mean()) / ranks.std()

    rng = np.random.default_rng(arguments.seed)
    rows: list[dict[str, object]] = []
    for effect in grid:
        hits_nominal = 0
        hits_bh = 0
        samples: list[float] = []
        for _ in range(arguments.replicates):
            if arguments.endpoint == "steatosis":
                residual = np.sqrt(max(0.0, 1.0 - effect * effect))
                latent = effect * anchor + residual * rng.normal(size=len(donors))
            else:
                latent = rng.normal(size=len(donors)) + effect * labels
            planted = plant(base, latent, arguments.variance_fraction, rng)
            predicted = cross_fitted_predictions(
                arguments.endpoint, planted, labels, folds
            )
            if arguments.endpoint == "steatosis":
                sample = abs(spearman(predicted, labels))
            else:
                sample = average_precision(labels.astype(np.int64), predicted)
            samples.append(sample)
            if sample >= arguments.critical_nominal:
                hits_nominal += 1
            if sample >= arguments.critical_bh:
                hits_bh += 1
        rows.append(
            {
                "endpoint": arguments.endpoint,
                "lineage_scope": REPRESENTATIVE_LINEAGE,
                "planted_effect": effect,
                "planted_variance_fraction": arguments.variance_fraction,
                "n": len(donors),
                "replicates": arguments.replicates,
                "mean_realised_metric": float(np.mean(samples)),
                "realised_power_nominal": hits_nominal / arguments.replicates,
                "realised_power_bh_single_true": hits_bh / arguments.replicates,
            }
        )
        print(json.dumps(rows[-1], sort_keys=True), flush=True)

    arguments.output.mkdir(parents=True, exist_ok=False)
    write_tsv(arguments.output / "realised_power.tsv", rows)
    audit = {
        "schema_version": "masld-bench-gse296875-realised-power-v1",
        "endpoint": arguments.endpoint,
        "lineage_scope": REPRESENTATIVE_LINEAGE,
        "planted_variance_fraction": arguments.variance_fraction,
        "replicates_per_point": arguments.replicates,
        "grid": list(grid),
        "n": len(donors),
        "seed": arguments.seed,
        "critical_nominal": arguments.critical_nominal,
        "critical_bh_single_true": arguments.critical_bh,
        "effect_is_planted_not_observed": True,
        "observed_model_performance_read": False,
        "prediction_artifact_opened": False,
        "score_artifact_opened": False,
        "status": "pass_realised_power",
    }
    (arguments.output / "audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
