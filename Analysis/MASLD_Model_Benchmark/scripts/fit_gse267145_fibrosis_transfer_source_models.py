#!/usr/bin/env python3
"""Fit GSE267145 mild-versus-advanced fibrosis models in per-array representations.

GSE49541 is log2 GPL570 intensity; GSE267145 is RNA-seq abundance.  The frozen
GSE267145 stage3 model is HVG/PCA on log1p CPM of sequencing counts and is not
transferable to an array, so this lane fits new source models in the two
representations the GSE49541 TaskSpec names by its baseline identifiers:

  per_array_rank  normalised within-array rank over the shared gene axis
  gene_median     log2 value minus that same array's own median gene

Both are computed strictly inside one participant column, so no statistic ever
pools across arrays and none can cross the platform boundary.  Everything
learned downstream - eligibility, variance ranking, standardisation, PCA,
hyperparameters, coefficients, calibration - is fit on GSE267145 participants
only and is frozen here before GSE49541 expression is ever opened.

The source fibrosis endpoint is the recorded NASH-CRN 0-to-3 scale.  GSE49541
deposits mild F0-F1 against advanced F3-F4 and excludes F2 from both arms, so
the mapping applies that same definition to the source rather than inventing a
new one; the F2 participants are dropped, not folded into either group.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence
import warnings

import numpy as np


CLASSES = ("mild_f0_f1", "advanced_f3_f4")
POSITIVE_CLASS = "advanced_f3_f4"
REPRESENTATIONS = ("per_array_rank", "gene_median")
LEARNED_MODEL_IDS = (
    "gene_median_elastic_net",
    "gene_median_linear_svm",
    "gene_median_pca_elastic_net",
    "per_array_rank_elastic_net",
)
PRIOR_MODEL_ID = "training_prevalence"
MODEL_IDS = (PRIOR_MODEL_ID,) + LEARNED_MODEL_IDS


class SourceFibrosisFitError(RuntimeError):
    """Raised when the source fit would violate its frozen recipe."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SourceFibrosisFitError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


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


# ---------------------------------------------------------------------------
# Per-array representations.  The predictor imports these so the source and the
# external array can never drift onto two different definitions.
# ---------------------------------------------------------------------------


def log2_cpm_complete_axis(matrix: np.ndarray) -> np.ndarray:
    """log2(1 + CPM) with the library total taken over the complete source axis.

    log2 rather than natural log because GSE49541 arrives as log2 array
    intensity; a mismatched log base would silently rescale every external
    feature by 1/ln(2) relative to the source-fitted standardiser.
    """

    values = np.asarray(matrix, dtype=np.float64)
    if values.ndim != 2 or not np.all(np.isfinite(values)) or np.any(values < 0):
        raise SourceFibrosisFitError("source RNA matrix is not finite and nonnegative")
    totals = values.sum(axis=1)
    if np.any(totals <= 0):
        raise SourceFibrosisFitError("source participant RNA library is empty")
    return np.log2(1.0 + values * (1_000_000.0 / totals[:, None]))


def per_array_rank(values: np.ndarray) -> np.ndarray:
    """Normalised within-array rank of every gene on the shared axis."""

    from scipy.stats import rankdata

    matrix = np.asarray(values, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[1] < 2 or not np.all(np.isfinite(matrix)):
        raise SourceFibrosisFitError("rank representation input is invalid")
    ranks = rankdata(matrix, method="average", axis=1)
    return (ranks - 0.5) / matrix.shape[1]


def gene_median_center(values: np.ndarray) -> np.ndarray:
    """Subtract each array's own median gene on the shared log2 axis."""

    matrix = np.asarray(values, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[1] < 2 or not np.all(np.isfinite(matrix)):
        raise SourceFibrosisFitError("gene-median representation input is invalid")
    return matrix - np.median(matrix, axis=1, keepdims=True)


REPRESENTATION_FUNCTIONS = {
    "per_array_rank": per_array_rank,
    "gene_median": gene_median_center,
}


def assert_row_independent(name: str, values: np.ndarray) -> None:
    """Prove the representation is a per-array function, not a pooled one.

    A transform that leaked information across arrays would change when its
    neighbours change.  Recomputing one row on its own and recomputing a
    reordered subset must both reproduce the full-matrix answer exactly.
    """

    function = REPRESENTATION_FUNCTIONS[name]
    matrix = np.asarray(values, dtype=np.float64)
    if matrix.shape[0] < 3:
        return
    full = function(matrix)
    alone = function(matrix[:1])
    if not np.array_equal(alone[0], full[0]):
        raise SourceFibrosisFitError(f"{name} is not a per-array transform")
    order = np.asarray([2, 0, 1], dtype=np.int64)
    subset = function(matrix[order])
    if not np.array_equal(subset, full[order]):
        raise SourceFibrosisFitError(f"{name} depends on array ordering")


def build_representation(
    *, name: str, log2_shared: np.ndarray
) -> np.ndarray:
    if name not in REPRESENTATION_FUNCTIONS:
        raise SourceFibrosisFitError(f"unregistered representation: {name}")
    assert_row_independent(name, log2_shared)
    return REPRESENTATION_FUNCTIONS[name](log2_shared)


# ---------------------------------------------------------------------------
# Source-only learned steps.
# ---------------------------------------------------------------------------


def class_weights(labels: np.ndarray) -> dict[int, float]:
    counts = Counter(int(value) for value in labels)
    if set(counts) != {0, 1}:
        raise SourceFibrosisFitError("source fitting partition lacks a fibrosis class")
    return {key: len(labels) / (2.0 * count) for key, count in counts.items()}


def select_feature_indices(
    *,
    eligibility_matrix: np.ndarray,
    representation: np.ndarray,
    gene_ids: Sequence[str],
    feature_count: int,
) -> np.ndarray:
    """Rank eligible genes by fitting-partition variance, tie-broken by gene ID."""

    ids = np.asarray(gene_ids, dtype=str)
    if (
        eligibility_matrix.shape != representation.shape
        or representation.shape[1] != len(ids)
    ):
        raise SourceFibrosisFitError("feature-selection axes differ")
    minimum_positive = max(3, math.ceil(0.10 * representation.shape[0]))
    eligible = np.flatnonzero(
        np.sum(eligibility_matrix > 0, axis=0) >= minimum_positive
    )
    if len(eligible) < feature_count:
        raise SourceFibrosisFitError("shared axis has too few eligible genes")
    variances = np.var(representation[:, eligible], axis=0, ddof=0)
    order = np.lexsort((ids[eligible], -variances))
    return np.asarray(np.sort(eligible[order[:feature_count]]), dtype=np.int64)


ELASTIC_NET_MAX_ITER = 20_000


def fit_elastic_net(
    x: np.ndarray,
    labels: np.ndarray,
    *,
    c_value: float,
    l1_ratio: float,
    seed: int,
    max_iter: int = ELASTIC_NET_MAX_ITER,
):
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        model = LogisticRegression(
            C=c_value,
            penalty="elasticnet",
            solver="saga",
            l1_ratio=l1_ratio,
            max_iter=max_iter,
            tol=1e-4,
            class_weight=class_weights(labels),
            random_state=seed,
        )
        model.fit(x, labels)
    converged = not any(
        isinstance(item.message, ConvergenceWarning) for item in caught
    ) and int(np.max(model.n_iter_)) < max_iter
    return model, converged


def classes_are_linearly_separable(design: np.ndarray, labels: np.ndarray) -> bool:
    """Exact feasibility test for a separating hyperplane.

    When the two classes can be separated, the logistic likelihood has no finite
    maximum: the coefficient norm grows without bound and the fit either never
    converges or is held together only by a penalty strong enough to erase it.
    That is the regime four positives put this source in, and it is a property
    of the data, not of the solver, so it is recorded rather than worked around.
    """

    from scipy.optimize import linprog

    x = np.asarray(design, dtype=np.float64)
    signs = np.where(np.asarray(labels, dtype=np.int64) == 1, 1.0, -1.0)
    augmented = np.hstack([x, np.ones((x.shape[0], 1))])
    result = linprog(
        c=np.zeros(augmented.shape[1]),
        A_ub=-signs[:, None] * augmented,
        b_ub=-np.ones(x.shape[0]),
        bounds=[(None, None)] * augmented.shape[1],
        method="highs",
    )
    return bool(result.status == 0 and result.success)


def fit_linear_svm(x: np.ndarray, labels: np.ndarray, *, c_value: float, seed: int):
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.svm import LinearSVC

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        model = LinearSVC(
            C=c_value,
            dual="auto",
            class_weight=class_weights(labels),
            max_iter=50_000,
            random_state=seed,
        )
        model.fit(x, labels)
    converged = not any(
        isinstance(item.message, ConvergenceWarning) for item in caught
    ) and int(np.max(np.atleast_1d(model.n_iter_))) < 50_000
    return model, converged


def fit_score_calibrator(raw: np.ndarray, labels: np.ndarray, *, seed: int):
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        model = LogisticRegression(
            C=1.0,
            max_iter=50_000,
            class_weight=class_weights(labels),
            random_state=seed,
        )
        model.fit(np.asarray(raw, dtype=np.float64).reshape(-1, 1), labels)
    converged = not any(
        isinstance(item.message, ConvergenceWarning) for item in caught
    )
    return model, converged


def average_precision(labels: np.ndarray, scores: np.ndarray) -> float:
    from sklearn.metrics import average_precision_score

    y = np.asarray(labels, dtype=np.int64)
    if len(np.unique(y)) != 2:
        raise SourceFibrosisFitError("average precision needs both classes")
    return float(average_precision_score(y, np.asarray(scores, dtype=np.float64)))


def average_precision_null(
    labels: np.ndarray, scores: np.ndarray, *, replicates: int, seed: int
) -> dict[str, float]:
    """Null for this exact score vector, built by permuting the source outcome.

    Average precision is not comparable across scorers with different tie
    structures.  A model whose coefficients have all been shrunk to zero emits a
    handful of distinct values, and its average precision sits near prevalence
    for reasons that have nothing to do with signal - which is precisely how a
    collapsed model outranks a working one when raw average precision is the
    selection metric.  Referencing every configuration to its own null removes
    that output file.
    """

    rng = np.random.default_rng(seed)
    permuted = np.asarray(labels, dtype=np.int64).copy()
    draws = np.empty(replicates, dtype=np.float64)
    for index in range(replicates):
        rng.shuffle(permuted)
        draws[index] = average_precision(permuted, scores)
    point = average_precision(labels, scores)
    return {
        "average_precision": point,
        "null_mean": float(np.mean(draws)),
        "null_p95": float(np.quantile(draws, 0.95)),
        "excess_over_null_mean": float(point - np.mean(draws)),
        "permutation_p_value": float(
            (1.0 + float(np.sum(draws >= point))) / (replicates + 1.0)
        ),
    }


def _sigmoid(values: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.asarray(values, dtype=np.float64)))


def _reduce(
    *,
    selected: np.ndarray,
    center: np.ndarray,
    scale: np.ndarray,
    pca_mean: np.ndarray | None,
    pca_components: np.ndarray | None,
    values: np.ndarray,
) -> np.ndarray:
    scaled = (np.asarray(values, dtype=np.float64)[:, selected] - center) / scale
    if pca_components is None:
        return scaled
    assert pca_mean is not None
    return (scaled - pca_mean) @ pca_components.T


def fit_pipeline(
    *,
    model_id: str,
    config: Mapping[str, Any],
    representation: np.ndarray,
    eligibility: np.ndarray,
    gene_ids: Sequence[str],
    fitting: np.ndarray,
    labels: np.ndarray,
    feature_count: int,
    pca_components: int,
    hyperparameters: Mapping[str, Any],
    seed: int,
) -> dict[str, Any]:
    """Fit one model on one partition and return its applicable linear state."""

    from sklearn.decomposition import PCA

    selected = select_feature_indices(
        eligibility_matrix=eligibility[fitting],
        representation=representation[fitting],
        gene_ids=gene_ids,
        feature_count=feature_count,
    )
    fit_dense = representation[fitting][:, selected]
    center = fit_dense.mean(axis=0)
    scale = fit_dense.std(axis=0)
    scale[scale <= np.finfo(np.float64).eps] = 1.0
    fit_scaled = (fit_dense - center) / scale
    pca_mean = None
    components = None
    if config["reducer"] == "pca":
        if pca_components >= min(fit_scaled.shape):
            raise SourceFibrosisFitError("partition cannot support the locked PCA size")
        pca = PCA(
            n_components=pca_components,
            svd_solver="randomized",
            whiten=False,
            random_state=seed + 90_000,
        )
        design = pca.fit_transform(fit_scaled)
        pca_mean = np.asarray(pca.mean_, dtype=np.float64)
        components = np.asarray(pca.components_, dtype=np.float64)
    elif config["reducer"] == "none":
        design = fit_scaled
    else:
        raise SourceFibrosisFitError(f"unregistered reducer: {config['reducer']}")
    state: dict[str, Any] = {
        "selected_indices": selected,
        "center": center,
        "scale": scale,
    }
    if components is not None:
        state["pca_mean"] = pca_mean
        state["pca_components"] = components
    converged = True
    if config["classifier"] == "elastic_net":
        model, ok = fit_elastic_net(
            design,
            labels[fitting],
            c_value=float(hyperparameters["c"]),
            l1_ratio=float(hyperparameters["l1_ratio"]),
            seed=seed + 91_100,
        )
        converged = converged and ok
        state["coef"] = np.asarray(model.coef_, dtype=np.float64)
        state["intercept"] = np.asarray(model.intercept_, dtype=np.float64)
        state["nonzero_coefficients"] = np.asarray(
            [int(np.sum(np.abs(model.coef_) > 0))], dtype=np.int64
        )
    elif config["classifier"] == "linear_svm":
        model, ok = fit_linear_svm(
            design,
            labels[fitting],
            c_value=float(hyperparameters["c"]),
            seed=seed + 91_000,
        )
        converged = converged and ok
        state["svm_coef"] = np.asarray(model.coef_, dtype=np.float64)
        state["svm_intercept"] = np.asarray(model.intercept_, dtype=np.float64)
        state["nonzero_coefficients"] = np.asarray(
            [int(np.sum(np.abs(model.coef_) > 0))], dtype=np.int64
        )
    else:
        raise SourceFibrosisFitError(f"unregistered classifier: {config['classifier']}")
    state["_converged"] = converged
    state["_design_columns"] = design.shape[1]
    state["_design"] = design
    return state


def apply_pipeline(
    *, config: Mapping[str, Any], state: Mapping[str, Any], representation: np.ndarray
) -> np.ndarray:
    """Apply a frozen pipeline state and return P(advanced_f3_f4)."""

    design = _reduce(
        selected=np.asarray(state["selected_indices"], dtype=np.int64),
        center=np.asarray(state["center"], dtype=np.float64),
        scale=np.asarray(state["scale"], dtype=np.float64),
        pca_mean=(
            np.asarray(state["pca_mean"], dtype=np.float64)
            if "pca_mean" in state
            else None
        ),
        pca_components=(
            np.asarray(state["pca_components"], dtype=np.float64)
            if "pca_components" in state
            else None
        ),
        values=representation,
    )
    if config["classifier"] == "elastic_net":
        logits = (
            design @ np.asarray(state["coef"], dtype=np.float64).T
            + np.asarray(state["intercept"], dtype=np.float64)
        )
        return _sigmoid(logits).ravel()
    if config["classifier"] == "linear_svm":
        raw = (
            design @ np.asarray(state["svm_coef"], dtype=np.float64).T
            + np.asarray(state["svm_intercept"], dtype=np.float64)
        )
        logits = (
            raw @ np.asarray(state["calibrator_coef"], dtype=np.float64).T
            + np.asarray(state["calibrator_intercept"], dtype=np.float64)
        )
        return _sigmoid(logits).ravel()
    raise SourceFibrosisFitError(f"unregistered classifier: {config['classifier']}")


def raw_decision(
    *, config: Mapping[str, Any], state: Mapping[str, Any], representation: np.ndarray
) -> np.ndarray:
    design = _reduce(
        selected=np.asarray(state["selected_indices"], dtype=np.int64),
        center=np.asarray(state["center"], dtype=np.float64),
        scale=np.asarray(state["scale"], dtype=np.float64),
        pca_mean=(
            np.asarray(state["pca_mean"], dtype=np.float64)
            if "pca_mean" in state
            else None
        ),
        pca_components=(
            np.asarray(state["pca_components"], dtype=np.float64)
            if "pca_components" in state
            else None
        ),
        values=representation,
    )
    if config["classifier"] == "linear_svm":
        return (
            design @ np.asarray(state["svm_coef"], dtype=np.float64).T
            + np.asarray(state["svm_intercept"], dtype=np.float64)
        ).ravel()
    logits = (
        design @ np.asarray(state["coef"], dtype=np.float64).T
        + np.asarray(state["intercept"], dtype=np.float64)
    )
    return logits.ravel()


def out_of_fold_scores(
    *,
    model_id: str,
    config: Mapping[str, Any],
    representation: np.ndarray,
    eligibility: np.ndarray,
    gene_ids: Sequence[str],
    outer_folds: np.ndarray,
    labels: np.ndarray,
    feature_count: int,
    pca_components: int,
    hyperparameters: Mapping[str, Any],
    seed: int,
) -> tuple[np.ndarray, bool, bool]:
    """Pooled out-of-fold decision scores over the frozen source outer folds.

    Returns the scores, whether every fold converged, and whether every fold's
    fitted model kept at least one non-zero coefficient.  A fold whose
    coefficients were all shrunk to zero contributes only its intercept, so the
    pooled vector degenerates into one value per fold and stops being a score.
    """

    scores = np.full(len(labels), np.nan, dtype=np.float64)
    converged = True
    every_fold_uses_the_data = True
    for fold in sorted(set(int(value) for value in outer_folds)):
        fitting = np.flatnonzero(outer_folds != fold)
        held = np.flatnonzero(outer_folds == fold)
        state = fit_pipeline(
            model_id=model_id,
            config=config,
            representation=representation,
            eligibility=eligibility,
            gene_ids=gene_ids,
            fitting=fitting,
            labels=labels,
            feature_count=feature_count,
            pca_components=pca_components,
            hyperparameters=hyperparameters,
            seed=seed + 10_000 * fold,
        )
        converged = converged and bool(state.pop("_converged"))
        state.pop("_design_columns", None)
        state.pop("_design", None)
        every_fold_uses_the_data = every_fold_uses_the_data and int(
            state["nonzero_coefficients"][0]
        ) > 0
        scores[held] = raw_decision(
            config=config, state=state, representation=representation[held]
        )
    if not np.all(np.isfinite(scores)):
        raise SourceFibrosisFitError("source out-of-fold scores are incomplete")
    return scores, converged, every_fold_uses_the_data


def select_hyperparameters(
    *,
    model_id: str,
    config: Mapping[str, Any],
    grid: Sequence[Mapping[str, Any]],
    representation: np.ndarray,
    eligibility: np.ndarray,
    gene_ids: Sequence[str],
    outer_folds: np.ndarray,
    labels: np.ndarray,
    feature_count: int,
    pca_components: int,
    seed: int,
    null_replicates: int,
    full_indices: np.ndarray,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Pick hyperparameters on GSE267145 alone, referenced to each config's null.

    GSE49541 plays no part: the grid, the folds, the metric, and the tie-break
    are all source-internal.  Two guards matter more than the grid itself.

    A configuration is eligible only if every outer-fold model and the full-fit
    model kept at least one non-zero coefficient and the full-fit model varies
    across participants.  An L1 fit that shrinks every coefficient to zero still
    returns a model object and still scores, but it is an intercept, not a
    classifier, and it must not win.

    Eligible configurations are then ranked by average precision *in excess of
    their own permutation null*, never by raw average precision, because scorers
    with different tie structures have different nulls and the raw number is not
    comparable across them.

    With four advanced-fibrosis source participants this selection is
    low-confidence whatever it returns, and the receipt records that.
    """

    rows: list[dict[str, Any]] = []
    best: tuple[float, float, float] | None = None
    chosen: dict[str, Any] | None = None
    for candidate in grid:
        scores, converged, folds_use_data = out_of_fold_scores(
            model_id=model_id,
            config=config,
            representation=representation,
            eligibility=eligibility,
            gene_ids=gene_ids,
            outer_folds=outer_folds,
            labels=labels,
            feature_count=feature_count,
            pca_components=pca_components,
            hyperparameters=candidate,
            seed=seed,
        )
        full_state = fit_pipeline(
            model_id=model_id,
            config=config,
            representation=representation,
            eligibility=eligibility,
            gene_ids=gene_ids,
            fitting=full_indices,
            labels=labels,
            feature_count=feature_count,
            pca_components=pca_components,
            hyperparameters=candidate,
            seed=seed,
        )
        full_converged = bool(full_state.pop("_converged"))
        full_state.pop("_design_columns", None)
        full_state.pop("_design", None)
        full_nonzero = int(full_state["nonzero_coefficients"][0])
        full_raw = raw_decision(
            config=config, state=full_state, representation=representation
        )
        full_distinct = int(len(np.unique(np.round(full_raw, 12))))
        null = average_precision_null(
            labels, scores, replicates=null_replicates, seed=seed + 55_000
        )
        distinct = int(len(np.unique(np.round(scores, 12))))
        eligible = bool(
            converged
            and full_converged
            and folds_use_data
            and full_nonzero > 0
            and full_distinct > 1
            and distinct > 1
        )
        rows.append(
            {
                "model_id": model_id,
                "c": float(candidate["c"]),
                "l1_ratio": float(candidate.get("l1_ratio", float("nan"))),
                "source_pooled_oof_average_precision": null["average_precision"],
                "source_pooled_oof_null_mean": null["null_mean"],
                "source_pooled_oof_null_p95": null["null_p95"],
                "source_pooled_oof_excess_over_null_mean": null["excess_over_null_mean"],
                "source_pooled_oof_permutation_p_value": null["permutation_p_value"],
                "distinct_oof_scores": distinct,
                "distinct_full_fit_scores": full_distinct,
                "full_fit_nonzero_coefficients": full_nonzero,
                "every_outer_fold_uses_the_data": folds_use_data,
                "converged": converged and full_converged,
                "eligible": eligible,
            }
        )
        if not eligible:
            continue
        key = (
            -null["excess_over_null_mean"],
            float(candidate["c"]),
            float(candidate.get("l1_ratio", 0.0)),
        )
        if best is None or key < best:
            best = key
            chosen = dict(candidate)
    eligible_count = int(sum(1 for row in rows if row["eligible"]))
    for row in rows:
        row["eligible_configurations_for_this_model"] = eligible_count
        row["selected"] = bool(
            chosen is not None
            and row["c"] == float(chosen["c"])
            and (
                np.isnan(row["l1_ratio"])
                or row["l1_ratio"] == float(chosen.get("l1_ratio", 0.0))
            )
        )
    # `chosen is None` is a result, not a crash: it means no configuration in the
    # registered grid is well posed on this source.  The caller records the model
    # as unfittable with its evidence rather than widening the grid until
    # something survives.
    return chosen, rows


def build_labels(
    *, endpoint_rows: Sequence[Mapping[str, str]], mapping: Mapping[str, Any]
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Map the recorded source fibrosis scale onto GSE49541's deposited groups."""

    field = str(mapping["source_field"])
    scale = [int(value) for value in mapping["source_scale"]]
    mild = {int(value) for value in mapping["mild_f0_f1"]}
    advanced = {int(value) for value in mapping["advanced_f3_f4"]}
    excluded = {int(value) for value in mapping["excluded"]}
    if (
        mild & advanced
        or mild & excluded
        or advanced & excluded
        or (mild | advanced | excluded) != set(scale)
    ):
        raise SourceFibrosisFitError("fibrosis mapping does not partition the scale")
    labels: list[int] = []
    keep: list[int] = []
    observed: Counter[int] = Counter()
    for index, row in enumerate(endpoint_rows):
        try:
            value = int(row[field])
        except (KeyError, ValueError) as error:
            raise SourceFibrosisFitError("source fibrosis value is not an integer") from error
        if value not in set(scale):
            raise SourceFibrosisFitError("source fibrosis value is off its recorded scale")
        observed[value] += 1
        if value in excluded:
            continue
        keep.append(index)
        labels.append(1 if value in advanced else 0)
    label_array = np.asarray(labels, dtype=np.int64)
    keep_array = np.asarray(keep, dtype=np.int64)
    if set(int(value) for value in label_array) != {0, 1}:
        raise SourceFibrosisFitError("source fibrosis mapping lacks a class")
    audit = {
        "source_field": field,
        "source_scale": scale,
        "source_scale_census": {str(key): observed[key] for key in sorted(observed)},
        "mild_f0_f1_source_values": sorted(mild),
        "advanced_f3_f4_source_values": sorted(advanced),
        "excluded_source_values": sorted(excluded),
        "excluded_participants": int(len(endpoint_rows) - len(keep_array)),
        "fitting_participants": int(len(keep_array)),
        "class_counts": {
            CLASSES[0]: int(np.sum(label_array == 0)),
            CLASSES[1]: int(np.sum(label_array == 1)),
        },
        "training_prevalence_advanced_f3_f4": float(np.mean(label_array)),
        "source_scale_has_stage_4": bool(4 in set(scale) and observed.get(4, 0) > 0),
        "mapping_rule": (
            "GSE49541's own deposited definition applied to the recorded source "
            "scale: F0-F1 mild, F3-F4 advanced, F2 excluded from both arms."
        ),
    }
    return label_array, keep_array, audit


def _check_hash(path: Path, expected: str, label: str) -> None:
    if sha256_file(path) != expected:
        raise SourceFibrosisFitError(f"{label} SHA-256 differs")


def run_fit(
    *, benchmark_root: Path, contract_path: Path, contract_sha256: str, output: Path
) -> dict[str, Any]:
    from masld_bench.artifacts import verify_frozen_tree

    if output.exists():
        raise SourceFibrosisFitError(f"refusing to overwrite source fit: {output}")
    _check_hash(contract_path, contract_sha256, "source fibrosis fit contract")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if (
        contract.get("status") != "registered_source_fit_pending"
        or contract.get("firewall", {}).get("external_labels_read") is not False
        or contract.get("firewall", {}).get("external_expression_values_read") is not False
    ):
        raise SourceFibrosisFitError("source fibrosis fit contract differs")

    activation_spec = contract["activation"]
    activation_root = benchmark_root / activation_spec["path"]
    _check_hash(
        activation_root / "ARTIFACTS.json",
        activation_spec["artifacts_sha256"],
        "transfer activation",
    )
    verify_frozen_tree(activation_root)
    axis_receipt = json.loads(
        (activation_root / "gene_axis_receipt.json").read_text(encoding="utf-8")
    )
    if (
        axis_receipt.get("status") != "pass_label_blind_shared_gene_axis"
        or axis_receipt.get("external_expression_values_read") is not False
        or axis_receipt.get("shared_genes") != activation_spec["shared_genes"]
    ):
        raise SourceFibrosisFitError("transfer activation receipt differs")

    inputs = contract["source_inputs"]
    molecular = benchmark_root / inputs["molecular_path"]
    outcomes = benchmark_root / inputs["outcomes_path"]
    folds = benchmark_root / inputs["folds_path"]
    for root, expected, label in (
        (molecular, inputs["molecular_artifacts_sha256"], "source molecular"),
        (outcomes, inputs["outcomes_artifacts_sha256"], "source outcomes"),
        (folds, inputs["folds_artifacts_sha256"], "source folds"),
    ):
        _check_hash(root / "ARTIFACTS.json", expected, label)
        verify_frozen_tree(root)

    _, participants = read_tsv(molecular / "participant_axis.tsv")
    _, source_axis = read_tsv(molecular / "rna_feature_axis.tsv")
    _, endpoint_rows = read_tsv(outcomes / "participant_endpoints.tsv")
    _, fold_rows = read_tsv(folds / "participant_outer_folds.tsv")
    participant_ids = [row["participant_id"] for row in participants]
    if (
        len(participant_ids) != inputs["participants"]
        or participant_ids != [row["participant_id"] for row in endpoint_rows]
        or participant_ids != [row["participant_id"] for row in fold_rows]
    ):
        raise SourceFibrosisFitError("source participant axes differ")

    labels, keep, endpoint_audit = build_labels(
        endpoint_rows=endpoint_rows, mapping=contract["endpoint_mapping"]
    )
    if endpoint_audit["class_counts"] != contract["endpoint_mapping"]["expected_class_counts"]:
        raise SourceFibrosisFitError("source fibrosis class census differs")
    outer_folds = np.asarray(
        [int(fold_rows[index]["outer_fold"]) for index in keep], dtype=np.int64
    )
    for fold in sorted(set(int(value) for value in outer_folds)):
        fitting = labels[outer_folds != fold]
        if set(int(value) for value in fitting) != {0, 1}:
            raise SourceFibrosisFitError(
                "a source outer-fold training partition lacks a fibrosis class"
            )

    _, axis_rows = read_tsv(activation_root / "common_stable_gene_axis.tsv")
    shared_ids = [row["stable_gene_id"] for row in axis_rows]
    source_ids = [row["stable_gene_id"] for row in source_axis]
    lookup = {value: index for index, value in enumerate(source_ids)}
    if len(lookup) != len(source_ids) or any(value not in lookup for value in shared_ids):
        raise SourceFibrosisFitError("shared gene is absent from the source RNA axis")
    raw = np.load(molecular / "rna_values.npy", mmap_mode="r", allow_pickle=False)
    if raw.shape != (inputs["participants"], len(source_ids)):
        raise SourceFibrosisFitError("source RNA matrix shape differs")
    columns = [lookup[value] for value in shared_ids]
    log2_full = log2_cpm_complete_axis(raw)
    log2_shared = log2_full[np.ix_(keep, columns)]
    eligibility = np.asarray(raw, dtype=np.float64)[np.ix_(keep, columns)]

    recipe = contract["recipe"]
    feature_count = int(recipe["feature_count"])
    pca_components = int(recipe["pca_components"])
    seeds = [int(value) for value in recipe["model_seeds"]]
    selection_seed = int(recipe["hyperparameter_selection_seed"])
    null_replicates = int(recipe["source_null_replicates"])
    representations = {
        name: build_representation(name=name, log2_shared=log2_shared)
        for name in REPRESENTATIONS
    }
    spread = {
        name: {
            "source_mean_interquartile_range": float(
                np.mean(
                    np.quantile(values, 0.75, axis=1) - np.quantile(values, 0.25, axis=1)
                )
            ),
            "source_mean_standard_deviation": float(np.mean(np.std(values, axis=1))),
        }
        for name, values in representations.items()
    }

    output.mkdir(parents=True)
    model_root = output / "models"
    model_root.mkdir()
    grid_rows: list[dict[str, Any]] = []
    oof_rows: list[dict[str, Any]] = []
    seed_rows: list[dict[str, Any]] = []
    selected_hyperparameters: dict[str, dict[str, Any]] = {}
    eligible_counts: dict[str, int] = {}
    separability: dict[str, bool] = {}
    unfittable: list[dict[str, Any]] = []
    source_model_ranking: list[dict[str, Any]] = []
    selected_axis_by_representation: dict[str, list[str]] = {}
    full_indices = np.arange(len(labels), dtype=np.int64)

    for model_id in LEARNED_MODEL_IDS:
        config = contract["models"][model_id]
        representation = representations[config["representation"]]
        grid = contract["hyperparameter_grid"][config["classifier"]]
        candidates = [
            {"c": c_value, "l1_ratio": l1_ratio}
            for c_value in grid["c"]
            for l1_ratio in grid.get("l1_ratio", [0.0])
        ]
        chosen, rows = select_hyperparameters(
            model_id=model_id,
            config=config,
            grid=candidates,
            representation=representation,
            eligibility=eligibility,
            gene_ids=shared_ids,
            outer_folds=outer_folds,
            labels=labels,
            feature_count=feature_count,
            pca_components=pca_components,
            seed=selection_seed,
            null_replicates=null_replicates,
            full_indices=full_indices,
        )
        grid_rows.extend(rows)
        eligible_counts[model_id] = int(
            rows[0]["eligible_configurations_for_this_model"]
        )
        if chosen is None:
            unfittable.append(
                {
                    "model_id": model_id,
                    "representation": config["representation"],
                    "reducer": config["reducer"],
                    "classifier": config["classifier"],
                    "grid_configurations": len(rows),
                    "eligible_grid_configurations": 0,
                    "configurations_that_converged": int(
                        sum(1 for row in rows if row["converged"])
                    ),
                    "configurations_with_a_nonzero_coefficient": int(
                        sum(1 for row in rows if row["full_fit_nonzero_coefficients"] > 0)
                    ),
                    "reason": (
                        "no registered configuration both converged in every outer "
                        "fold and in the full fit and kept a non-zero coefficient; "
                        "the classes are separable at four positives, so weakening "
                        "the penalty diverges and strengthening it erases the model"
                    ),
                }
            )
            continue
        selected_hyperparameters[model_id] = chosen

        for seed in seeds:
            oof_raw, oof_converged, oof_folds_use_data = out_of_fold_scores(
                model_id=model_id,
                config=config,
                representation=representation,
                eligibility=eligibility,
                gene_ids=shared_ids,
                outer_folds=outer_folds,
                labels=labels,
                feature_count=feature_count,
                pca_components=pca_components,
                hyperparameters=chosen,
                seed=seed,
            )
            state = fit_pipeline(
                model_id=model_id,
                config=config,
                representation=representation,
                eligibility=eligibility,
                gene_ids=shared_ids,
                fitting=full_indices,
                labels=labels,
                feature_count=feature_count,
                pca_components=pca_components,
                hyperparameters=chosen,
                seed=seed,
            )
            converged = bool(state.pop("_converged")) and oof_converged
            design_columns = int(state.pop("_design_columns"))
            design = state.pop("_design")
            separable = classes_are_linearly_separable(design, labels)
            if config["classifier"] == "linear_svm":
                calibrator, calibrator_ok = fit_score_calibrator(
                    oof_raw, labels, seed=seed + 92_000
                )
                converged = converged and calibrator_ok
                state["calibrator_coef"] = np.asarray(
                    calibrator.coef_, dtype=np.float64
                )
                state["calibrator_intercept"] = np.asarray(
                    calibrator.intercept_, dtype=np.float64
                )
            if not converged:
                raise SourceFibrosisFitError(
                    f"selected configuration for {model_id} did not converge at seed {seed}"
                )
            selected_ids = [shared_ids[index] for index in state["selected_indices"]]
            prior = selected_axis_by_representation.setdefault(
                config["representation"], selected_ids
            )
            if prior != selected_ids:
                raise SourceFibrosisFitError(
                    "feature selection is not a function of the fitting partition alone"
                )
            fitted = apply_pipeline(
                config=config, state=state, representation=representation
            )
            if not np.all(np.isfinite(fitted)) or np.any(fitted < 0) or np.any(fitted > 1):
                raise SourceFibrosisFitError("source fitted probabilities are invalid")
            oof_probability = (
                _sigmoid(
                    np.asarray(state["calibrator_coef"], dtype=np.float64).ravel()[0]
                    * oof_raw
                    + np.asarray(state["calibrator_intercept"], dtype=np.float64).ravel()[0]
                )
                if config["classifier"] == "linear_svm"
                else _sigmoid(oof_raw)
            )
            oof_null = average_precision_null(
                labels, oof_raw, replicates=null_replicates, seed=seed + 55_000
            )
            oof_ap = oof_null["average_precision"]
            full_distinct = int(len(np.unique(np.round(fitted, 12))))
            if (
                int(state["nonzero_coefficients"][0]) == 0
                or full_distinct < 2
                or not oof_folds_use_data
            ):
                raise SourceFibrosisFitError(
                    f"selected configuration for {model_id} is degenerate at seed {seed}"
                )
            for index, participant_index in enumerate(keep):
                oof_rows.append(
                    {
                        "model_id": model_id,
                        "model_seed": seed,
                        "participant_id": participant_ids[participant_index],
                        "outer_fold": int(outer_folds[index]),
                        "fibrosis_group": CLASSES[int(labels[index])],
                        "oof_raw_score": format(float(oof_raw[index]), ".17g"),
                        "oof_probability_advanced_f3_f4": format(
                            float(oof_probability[index]), ".17g"
                        ),
                    }
                )
            current = model_root / model_id / f"seed_{seed}"
            current.mkdir(parents=True)
            np.savez_compressed(
                current / "model_state.npz",
                **{
                    key: np.asarray(value)
                    for key, value in state.items()
                },
            )
            state_sha256 = sha256_file(current / "model_state.npz")
            seed_receipt = {
                "schema_version": "masld-bench-gse267145-fibrosis-transfer-seed-v1",
                "model_id": model_id,
                "model_seed": int(seed),
                "representation": config["representation"],
                "reducer": config["reducer"],
                "classifier": config["classifier"],
                "hyperparameters": chosen,
                "fitting_participants": int(len(labels)),
                "selected_features": len(selected_ids),
                "design_columns": design_columns,
                "training_partition_linearly_separable": separable,
                "eligible_grid_configurations": eligible_counts[model_id],
                "selected_configuration_is_isolated": eligible_counts[model_id] <= 1,
                "nonzero_coefficients": int(state["nonzero_coefficients"][0]),
                "source_pooled_oof_average_precision": oof_ap,
                "source_pooled_oof_null_mean": oof_null["null_mean"],
                "source_pooled_oof_null_p95": oof_null["null_p95"],
                "source_pooled_oof_excess_over_null_mean": oof_null[
                    "excess_over_null_mean"
                ],
                "source_pooled_oof_permutation_p_value": oof_null["permutation_p_value"],
                "distinct_full_fit_scores": full_distinct,
                "model_state_sha256": state_sha256,
                "external_expression_values_read": False,
                "external_labels_read": False,
            }
            with (current / "receipt.json").open("x", encoding="utf-8") as handle:
                handle.write(json.dumps(seed_receipt, indent=2, sort_keys=True) + "\n")
            seed_rows.append(
                {
                    "model_id": model_id,
                    "model_seed": seed,
                    "representation": config["representation"],
                    "model_state_path": f"models/{model_id}/seed_{seed}/model_state.npz",
                    "model_state_sha256": state_sha256,
                    "source_pooled_oof_average_precision": oof_ap,
                    "source_pooled_oof_null_mean": oof_null["null_mean"],
                    "source_pooled_oof_excess_over_null_mean": oof_null[
                        "excess_over_null_mean"
                    ],
                    "source_pooled_oof_permutation_p_value": oof_null[
                        "permutation_p_value"
                    ],
                    "nonzero_coefficients": int(state["nonzero_coefficients"][0]),
                    "distinct_full_fit_scores": full_distinct,
                    "training_partition_linearly_separable": separable,
                    "eligible_grid_configurations": eligible_counts[model_id],
                    "selected_configuration_is_isolated": eligible_counts[model_id] <= 1,
                }
            )
        model_rows = [row for row in seed_rows if row["model_id"] == model_id]
        model_seed_ap = [
            row["source_pooled_oof_average_precision"] for row in model_rows
        ]
        model_seed_excess = [
            row["source_pooled_oof_excess_over_null_mean"] for row in model_rows
        ]
        model_seed_p = [
            row["source_pooled_oof_permutation_p_value"] for row in model_rows
        ]
        source_model_ranking.append(
            {
                "model_id": model_id,
                "mean_source_pooled_oof_average_precision": float(
                    np.mean(model_seed_ap)
                ),
                "mean_source_pooled_oof_null_mean": float(
                    np.mean([row["source_pooled_oof_null_mean"] for row in model_rows])
                ),
                "mean_source_pooled_oof_excess_over_null_mean": float(
                    np.mean(model_seed_excess)
                ),
                "min_source_pooled_oof_permutation_p_value": float(np.min(model_seed_p)),
                "exceeds_its_own_null_on_every_seed": bool(
                    all(value > 0.0 for value in model_seed_excess)
                ),
                "selected_c": float(chosen["c"]),
                "selected_l1_ratio": float(chosen.get("l1_ratio", float("nan"))),
                "eligible_grid_configurations": eligible_counts[model_id],
                "selected_configuration_is_isolated": eligible_counts[model_id] <= 1,
                "training_partition_linearly_separable": bool(
                    model_rows[0]["training_partition_linearly_separable"]
                ),
            }
        )
        separability[model_id] = bool(
            model_rows[0]["training_partition_linearly_separable"]
        )

    if not source_model_ranking:
        raise SourceFibrosisFitError(
            "every learned model is unfittable on this source; there is nothing to transfer"
        )
    source_model_ranking.sort(
        key=lambda row: (
            -row["mean_source_pooled_oof_excess_over_null_mean"],
            row["model_id"],
        )
    )
    selected_source_model_id = source_model_ranking[0]["model_id"]
    no_source_model_exceeds_its_own_null = not any(
        row["mean_source_pooled_oof_excess_over_null_mean"] > 0.0
        for row in source_model_ranking
    )

    for representation, selected_ids in sorted(selected_axis_by_representation.items()):
        write_tsv(
            output / f"selected_feature_axis__{representation}.tsv",
            ("feature_index", "shared_axis_index", "stable_gene_id"),
            [
                {
                    "feature_index": index,
                    "shared_axis_index": shared_ids.index(stable_id),
                    "stable_gene_id": stable_id,
                }
                for index, stable_id in enumerate(selected_ids)
            ],
        )
    write_tsv(
        output / "hyperparameter_selection.tsv",
        tuple(grid_rows[0]),
        grid_rows,
    )
    write_tsv(
        output / "source_out_of_fold_scores.tsv",
        (
            "model_id",
            "model_seed",
            "participant_id",
            "outer_fold",
            "fibrosis_group",
            "oof_raw_score",
            "oof_probability_advanced_f3_f4",
        ),
        oof_rows,
    )
    write_tsv(output / "seed_model_index.tsv", tuple(seed_rows[0]), seed_rows)
    write_tsv(
        output / "source_model_ranking.tsv",
        tuple(source_model_ranking[0]),
        source_model_ranking,
    )

    receipt = {
        "schema_version": "masld-bench-gse267145-fibrosis-transfer-fit-v1",
        "status": "passed_source_fibrosis_transfer_fit_preprocessing_locked",
        "source_series": "GSE267145",
        "external_series": "GSE49541",
        "arm_id": contract["arm_id"],
        "gate_eligible_arm": bool(contract["gate_eligible_arm"]),
        "source_participants_deposited": int(inputs["participants"]),
        "fitting_participants": int(len(labels)),
        "endpoint_mapping_audit": endpoint_audit,
        "training_prevalence_advanced_f3_f4": endpoint_audit[
            "training_prevalence_advanced_f3_f4"
        ],
        "shared_genes": len(shared_ids),
        "selected_features": feature_count,
        "pca_components": pca_components,
        "model_seeds": seeds,
        "seeds_are_biological_replicates": False,
        "hyperparameter_selection": "source_pooled_out_of_fold_average_precision",
        "hyperparameter_selection_seed": selection_seed,
        "selected_hyperparameters": selected_hyperparameters,
        "source_model_ranking": source_model_ranking,
        "selected_source_model_id": selected_source_model_id,
        "fitted_model_ids": sorted(selected_hyperparameters),
        "unfittable_model_ids": sorted(row["model_id"] for row in unfittable),
        "unfittable_models": unfittable,
        "selected_source_model_selection_basis": "source_only_pooled_oof_average_precision_in_excess_of_its_own_permutation_null",
        "source_null_replicates": null_replicates,
        "elastic_net_max_iter": ELASTIC_NET_MAX_ITER,
        "training_partition_linearly_separable": dict(sorted(separability.items())),
        "eligible_grid_configurations": dict(sorted(eligible_counts.items())),
        "models_with_an_isolated_surviving_configuration": sorted(
            model_id for model_id, count in eligible_counts.items() if count <= 1
        ),
        "complete_separation_diagnosis": (
            "With four advanced-fibrosis participants the classes are linearly "
            "separable in every fitted design, so the logistic likelihood has no "
            "finite maximum. Configurations either shrink every coefficient to "
            "zero or fail to converge; a model that survives at a single grid "
            "point sits on the boundary between those two failures and is not a "
            "well-posed fit."
        ),
        "no_source_model_exceeds_its_own_null": no_source_model_exceeds_its_own_null,
        "source_signal_present": not no_source_model_exceeds_its_own_null,
        "representation_spread_diagnostics": spread,
        "per_array_transform_is_row_independent": True,
        "cross_array_pooling_performed": False,
        "external_expression_values_read": False,
        "external_labels_read": False,
        "external_fit_or_calibration_performed": False,
        "model_selected_or_repaired_from_external_outcomes": False,
        "low_positive_count_warning": bool(
            endpoint_audit["class_counts"][POSITIVE_CLASS] < 10
        ),
        "project_sealed": False,
        "champion_claim_allowed": False,
        "diagnostic_or_prognostic_claim_allowed": False,
        "clinical_claim_allowed": False,
    }
    with (output / "fit_receipt.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-root", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--contract-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    receipt = run_fit(
        benchmark_root=arguments.benchmark_root,
        contract_path=arguments.contract,
        contract_sha256=arguments.contract_sha256,
        output=arguments.output,
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
