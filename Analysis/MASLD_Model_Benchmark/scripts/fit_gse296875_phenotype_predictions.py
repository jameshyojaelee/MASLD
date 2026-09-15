#!/usr/bin/env python3
"""Fit the three arms on frozen donor folds and emit out-of-fold predictions.

This process writes predictions and nothing else.  It computes no metric, and
it runs on a Python that cannot import the benchmark package, so the scoring
module is unreachable from here.  Scoring happens in a separate job against the
hash-frozen prediction table.

Label handling.  A cross-fitted supervised model must see training-fold labels;
that is what fitting is.  What the separation guarantees is narrower and is the
part that matters: no test-fold label enters any fit, any hyperparameter
choice, or any transform.  Every normalisation, filter, variance ranking,
scaling, principal-component basis, and imputation mean is estimated on
training donors alone, and the inner leave-one-donor-out selection is nested
inside the training fold.

Donors whose lineage unit falls below the frozen minimum-cell threshold keep a
prediction.  Their block is filled with the training-fold lineage mean and
flagged by a binary observed indicator, so the model is told the block is
imputed rather than being handed a fabricated zero.  The imputed count is
reported per scope because a scope with many imputed donors is predicting
mostly from the indicator.
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


class FitError(RuntimeError):
    """Raised when a frozen input or a structural expectation does not hold."""


PRIMARY_LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")
SECONDARY_LINEAGES = ("endothelial_cell", "b_cell")
SCOPES = ("all_lineage",) + PRIMARY_LINEAGES + SECONDARY_LINEAGES
ARMS = ("molecular", "metadata", "molecular_metadata")

ENDPOINTS = ("steatosis", "fibrosis")
ALPHA_GRID = (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)
C_GRID = (0.001, 0.01, 0.1, 1.0, 10.0)

DETECTION_FRACTION = 0.5
N_HIGHLY_VARIABLE = 2_000
MAX_COMPONENTS = 20
SEED = 20260825

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


def log_cpm(counts: np.ndarray) -> np.ndarray:
    """Per-unit counts per million, then log1p.  No parameter is fitted here.

    A unit with no nuclei has no library to normalise against.  Donor 733 has
    zero B cells, so its b_cell unit is structurally missing rather than thin.
    Such a row becomes NaN, a typed missing value, and never zero: a zero row
    would be a claim that the lineage was sampled and found silent.  NaN also
    fails loudly if it ever reaches a fit, which is the point.
    """

    library = counts.sum(axis=1, keepdims=True).astype(np.float64)
    empty = (library <= 0).ravel()
    safe = np.where(library <= 0, 1.0, library)
    transformed = np.log1p(counts.astype(np.float64) / safe * 1e6)
    transformed[empty] = np.nan
    return transformed


def molecular_block(
    expression: np.ndarray,
    observed: np.ndarray,
    train_rows: np.ndarray,
    all_rows: np.ndarray,
) -> np.ndarray:
    """Fit every transform on observed training units, then apply to all rows.

    Returns the principal-component block for ``all_rows`` with one extra
    column carrying the observed indicator.
    """

    fit_rows = train_rows[observed[train_rows]]
    if fit_rows.size < 3:
        raise FitError("too few observed training units to fit a transform")
    fit_matrix = expression[fit_rows]
    if not np.isfinite(fit_matrix).all():
        raise FitError("a structurally missing unit reached the fitted transform")

    detected = (fit_matrix > 0).mean(axis=0) >= DETECTION_FRACTION
    if detected.sum() < 10:
        detected = (fit_matrix > 0).mean(axis=0) > 0
    fit_matrix = fit_matrix[:, detected]

    variance = fit_matrix.var(axis=0)
    keep = np.argsort(-variance)[: min(N_HIGHLY_VARIABLE, fit_matrix.shape[1])]
    keep.sort()
    fit_matrix = fit_matrix[:, keep]

    mean = fit_matrix.mean(axis=0)
    deviation = fit_matrix.std(axis=0)
    deviation[deviation == 0.0] = 1.0

    components = min(MAX_COMPONENTS, fit_rows.size - 1, fit_matrix.shape[1])
    pca = PCA(n_components=components, svd_solver="full", random_state=SEED)
    pca.fit((fit_matrix - mean) / deviation)

    # Masked units are filled with the training-fold mean, never with zero.
    applied = expression[all_rows][:, detected][:, keep]
    masked = ~observed[all_rows]
    if masked.any():
        applied[masked] = mean
    projected = pca.transform((applied - mean) / deviation)
    if not np.isfinite(projected).all():
        raise FitError("a non-finite value survived into the feature block")
    indicator = observed[all_rows].astype(np.float64).reshape(-1, 1)
    return np.hstack([projected, indicator])


def metadata_block(
    metadata: np.ndarray, train_rows: np.ndarray, all_rows: np.ndarray
) -> np.ndarray:
    fit_matrix = metadata[train_rows]
    mean = fit_matrix.mean(axis=0)
    deviation = fit_matrix.std(axis=0)
    deviation[deviation == 0.0] = 1.0
    return (metadata[all_rows] - mean) / deviation


def build_features(
    arm: str,
    scope: str,
    expression: dict[str, np.ndarray],
    observed: dict[str, np.ndarray],
    metadata: np.ndarray,
    train_index: np.ndarray,
    all_index: np.ndarray,
) -> np.ndarray:
    blocks: list[np.ndarray] = []
    if arm in {"molecular", "molecular_metadata"}:
        lineages = PRIMARY_LINEAGES if scope == "all_lineage" else (scope,)
        for lineage in lineages:
            blocks.append(
                molecular_block(
                    expression[lineage], observed[lineage], train_index, all_index
                )
            )
    if arm in {"metadata", "molecular_metadata"}:
        blocks.append(metadata_block(metadata, train_index, all_index))
    if not blocks:
        raise FitError(f"arm produced no features: {arm}")
    return np.hstack(blocks)


def fit_predict(
    endpoint: str,
    features_train: np.ndarray,
    labels_train: np.ndarray,
    features_test: np.ndarray,
    penalty: float,
) -> np.ndarray:
    if endpoint == "steatosis":
        model = Ridge(alpha=penalty, random_state=SEED)
        model.fit(features_train, labels_train)
        return model.predict(features_test)
    model = LogisticRegression(
        C=penalty, penalty="l2", solver="lbfgs", max_iter=5_000, random_state=SEED
    )
    model.fit(features_train, labels_train)
    return model.predict_proba(features_test)[:, 1]


def inner_loss(endpoint: str, truth: float, prediction: float) -> float:
    if endpoint == "steatosis":
        return (truth - prediction) ** 2
    clipped = min(max(prediction, 1e-9), 1.0 - 1e-9)
    return -(truth * np.log(clipped) + (1.0 - truth) * np.log(1.0 - clipped))


def select_penalty(
    endpoint: str,
    arm: str,
    scope: str,
    expression: dict[str, np.ndarray],
    observed: dict[str, np.ndarray],
    metadata: np.ndarray,
    train_index: np.ndarray,
    labels: np.ndarray,
) -> float:
    """Leave one training donor out, refitting every transform each time."""

    grid = ALPHA_GRID if endpoint == "steatosis" else C_GRID
    totals = np.zeros(len(grid), dtype=np.float64)
    usable = 0
    for position in range(train_index.size):
        inner_train = np.delete(train_index, position)
        held = train_index[position : position + 1]
        inner_labels = np.delete(labels, position)
        if endpoint == "fibrosis" and len(set(inner_labels.tolist())) < 2:
            continue
        combined = np.concatenate([inner_train, held])
        try:
            features = build_features(
                arm, scope, expression, observed, metadata, inner_train, combined
            )
        except FitError:
            continue
        train_features = features[: inner_train.size]
        test_features = features[inner_train.size :]
        for index, penalty in enumerate(grid):
            prediction = fit_predict(
                endpoint, train_features, inner_labels, test_features, penalty
            )
            totals[index] += inner_loss(endpoint, float(labels[position]), float(prediction[0]))
        usable += 1
    if usable == 0:
        return float(grid[len(grid) // 2])
    return float(grid[int(np.argmin(totals))])


def run_partition(
    partition: str,
    fold_of: dict[str, int],
    donors_in_scope: list[str],
    endpoint_values: dict[str, dict[str, float]],
    row_of: dict[tuple[str, str], int],
    expression: dict[str, np.ndarray],
    observed: dict[str, np.ndarray],
    metadata_of: dict[str, np.ndarray],
    *,
    arms: tuple[str, ...] = ARMS,
    scopes: tuple[str, ...] = SCOPES,
    only_fold: int | None = None,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Fit one partition.  Nothing outside the reported grid is computed."""

    predictions: list[dict[str, object]] = []
    receipts: list[dict[str, object]] = []

    for endpoint in ENDPOINTS:
        values = endpoint_values[endpoint]
        donors = [donor for donor in donors_in_scope if donor in values]
        metadata = np.vstack([metadata_of[donor] for donor in donors])
        index_of = {donor: position for position, donor in enumerate(donors)}
        rows = {
            lineage: np.asarray(
                [row_of[(donor, lineage)] for donor in donors], dtype=np.int64
            )
            for lineage in PRIMARY_LINEAGES + SECONDARY_LINEAGES
        }
        expression_by_lineage = {
            lineage: expression["counts"][rows[lineage]] for lineage in rows
        }
        observed_by_lineage = {
            lineage: observed["unit"][rows[lineage]] for lineage in rows
        }
        labels = np.asarray([values[donor] for donor in donors], dtype=np.float64)
        folds = sorted({fold_of[donor] for donor in donors})

        for arm in arms:
            for scope in scopes:
                if arm == "metadata" and scope != "all_lineage":
                    continue  # metadata does not vary by lineage
                for fold in folds:
                    if only_fold is not None and fold != only_fold:
                        continue
                    test_positions = np.asarray(
                        [index_of[d] for d in donors if fold_of[d] == fold],
                        dtype=np.int64,
                    )
                    train_positions = np.asarray(
                        [index_of[d] for d in donors if fold_of[d] != fold],
                        dtype=np.int64,
                    )
                    if test_positions.size == 0 or train_positions.size < 4:
                        continue
                    if endpoint == "fibrosis" and len(
                        set(labels[train_positions].tolist())
                    ) < 2:
                        raise FitError(
                            f"training fold {fold} has one fibrosis class only"
                        )
                    penalty = select_penalty(
                        endpoint,
                        arm,
                        scope,
                        expression_by_lineage,
                        observed_by_lineage,
                        metadata,
                        train_positions,
                        labels[train_positions],
                    )
                    combined = np.concatenate([train_positions, test_positions])
                    features = build_features(
                        arm,
                        scope,
                        expression_by_lineage,
                        observed_by_lineage,
                        metadata,
                        train_positions,
                        combined,
                    )
                    predicted = fit_predict(
                        endpoint,
                        features[: train_positions.size],
                        labels[train_positions],
                        features[train_positions.size :],
                        penalty,
                    )
                    for offset, position in enumerate(test_positions):
                        predictions.append(
                            {
                                "donor_id": donors[position],
                                "arm": arm,
                                "lineage_scope": scope,
                                "endpoint": endpoint,
                                "prediction": f"{float(predicted[offset]):.10f}",
                                "outer_fold": fold,
                            }
                        )
                    imputed = 0
                    if arm != "metadata":
                        lineages = (
                            PRIMARY_LINEAGES if scope == "all_lineage" else (scope,)
                        )
                        imputed = int(
                            sum(
                                (~observed_by_lineage[lineage][test_positions]).sum()
                                for lineage in lineages
                            )
                        )
                    receipts.append(
                        {
                            "partition": partition,
                            "endpoint": endpoint,
                            "arm": arm,
                            "lineage_scope": scope,
                            "outer_fold": fold,
                            "train_donors": int(train_positions.size),
                            "test_donors": int(test_positions.size),
                            "selected_penalty": penalty,
                            "n_features": int(features.shape[1]),
                            "test_units_imputed": imputed,
                        }
                    )
    return predictions, receipts


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--pseudobulk", type=Path, required=True)
    parser.add_argument("--pseudobulk-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    root = arguments.root
    endpoints_root = root / "executions/gse296875-phenotype-endpoints-20260825"
    folds_root = root / "executions/gse296875-donor-folds-20260825"
    if digest(endpoints_root / "ARTIFACTS.json") != BOUND["phenotype_endpoints"]:
        raise FitError("endpoint lock changed")
    if digest(folds_root / "ARTIFACTS.json") != BOUND["donor_folds"]:
        raise FitError("donor fold lock changed")
    if digest(arguments.pseudobulk / "ARTIFACTS.json") != arguments.pseudobulk_sha256:
        raise FitError("pseudobulk fixture changed")

    units = read_tsv(arguments.pseudobulk / "pseudobulk" / "units.tsv")
    with np.load(arguments.pseudobulk / "pseudobulk" / "donor_lineage_counts.npz") as p:
        counts = p["counts"]
        unit_observed = p["unit_observed"]
    row_of = {
        (row["donor_id"], row["lineage_id"]): int(row["unit_index"]) for row in units
    }

    endpoint_rows = read_tsv(
        endpoints_root / "endpoint_lock" / "donor_endpoints.tsv"
    )
    fold_rows = read_tsv(folds_root / "split_lock" / "donor_folds.tsv")

    values = {
        "steatosis": {
            row["donor_id"]: float(row["steatosis_numeric"])
            for row in endpoint_rows
            if row["steatosis_observed"] == "true"
        },
        "fibrosis": {
            row["donor_id"]: 1.0 if row["fibrosis_any"] == "true" else 0.0
            for row in endpoint_rows
            if row["fibrosis_observed"] == "true"
        },
    }
    metadata_of = {
        row["donor_id"]: np.asarray(
            [
                float(row["age_in_yr"]),
                1.0 if row["reported_sex"] == "M" else 0.0,
                float(row["BMI"]),
            ],
            dtype=np.float64,
        )
        for row in endpoint_rows
    }
    fold_of = {row["donor_id"]: int(row["outer_fold"]) for row in fold_rows}
    well_of = {row["donor_id"]: row["well_id"] for row in fold_rows}
    adults = [row["donor_id"] for row in fold_rows if row["is_adult"] == "true"]
    all_donors = [row["donor_id"] for row in fold_rows]

    expression = {"counts": log_cpm(counts)}
    observed = {"unit": unit_observed}

    arguments.output.mkdir(parents=True, exist_ok=False)
    receipts: list[dict[str, object]] = []

    predictions, fold_receipts = run_partition(
        "all_age", fold_of, all_donors, values, row_of, expression, observed, metadata_of
    )
    receipts.extend(fold_receipts)
    write_tsv(arguments.output / "predictions.tsv", predictions)

    adult_folds = {donor: fold_of[donor] for donor in adults}
    adult_predictions, adult_receipts = run_partition(
        "adult_only_refit",
        adult_folds,
        adults,
        values,
        row_of,
        expression,
        observed,
        metadata_of,
        arms=("molecular",),
        scopes=("all_lineage", *PRIMARY_LINEAGES),
    )
    receipts.extend(adult_receipts)
    write_tsv(arguments.output / "predictions_adult_refit.tsv", adult_predictions)

    well_predictions: list[dict[str, object]] = []
    for held_out in sorted({well_of[donor] for donor in all_donors}):
        pseudo_fold = {
            donor: (0 if well_of[donor] == held_out else 1) for donor in all_donors
        }
        rows, well_receipts = run_partition(
            f"leave_one_well_out_{held_out}",
            pseudo_fold,
            all_donors,
            values,
            row_of,
            expression,
            observed,
            metadata_of,
            arms=("molecular",),
            scopes=("all_lineage",),
            only_fold=0,
        )
        for row in rows:
            if row["outer_fold"] != 0:
                continue  # keep only the held-out well's predictions
            row = dict(row)
            row["outer_fold"] = held_out
            well_predictions.append(row)
        receipts.extend(
            r for r in well_receipts if r["outer_fold"] == 0
        )
    write_tsv(arguments.output / "predictions_leave_one_well_out.tsv", well_predictions)
    write_tsv(arguments.output / "fit_receipts.tsv", receipts)

    audit = {
        "schema_version": "masld-bench-gse296875-phenotype-predictions-v1",
        "dataset_id": "gse296875",
        "unit_of_inference": "donor",
        "arms": list(ARMS),
        "scopes": list(SCOPES),
        "confirmatory_scopes": ["all_lineage", *PRIMARY_LINEAGES],
        "secondary_scopes": list(SECONDARY_LINEAGES),
        "endpoints": list(ENDPOINTS),
        "all_age_predictions": len(predictions),
        "adult_refit_predictions": len(adult_predictions),
        "leave_one_well_out_predictions": len(well_predictions),
        "test_fold_labels_never_entered_a_fit_or_a_transform": True,
        "transforms_fit_inside_training_folds_only": True,
        "hyperparameters_selected_by_nested_leave_one_donor_out": True,
        "masked_units_imputed_at_training_mean_with_indicator": True,
        "masked_units_never_zero": True,
        "metrics_computed": 0,
        "metrics_are_computed_in_a_separate_process": True,
        "seed": SEED,
        "status": "pass_predictions_frozen",
    }
    (arguments.output / "audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
