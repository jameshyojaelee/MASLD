"""Fold-local modality controls for the observed-multiome task.

These helpers fit profile predictors but never calculate benchmark metrics.
The caller must remove scored ATAC bins and their buffer before passing the
``observed_atac`` matrix. Missing modalities are represented by ``None`` and
are never converted to all-zero observations.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np


MODEL_IDS = (
    "masked_modality",
    "observed_atac_glm",
    "observed_atac_only",
    "rna_only",
    "shuffled_modality",
)


class ObservedMultiomeBaselineError(ValueError):
    """Raised when a modality, shape, missingness, or fit requirement differs."""


@dataclass(frozen=True)
class BaselineState:
    model_id: str
    seed: int
    feature_mean: np.ndarray
    feature_scale: np.ndarray
    coefficients: np.ndarray
    intercept: np.ndarray
    link: str
    training_pool: np.ndarray | None = None
    training_pool_strata: tuple[str, ...] = ()


def _matrix(value: np.ndarray | None, *, rows: int, label: str) -> np.ndarray:
    if value is None:
        raise ObservedMultiomeBaselineError(f"{label} is structurally missing")
    array = np.asarray(value, dtype=np.float64)
    if array.ndim != 2 or array.shape[0] != rows or array.shape[1] < 1:
        raise ObservedMultiomeBaselineError(f"{label} shape differs")
    if not np.isfinite(array).all():
        raise ObservedMultiomeBaselineError(f"{label} contains non-finite values")
    return array


def _validate_rows(row_ids: Sequence[str], strata: Sequence[str]) -> tuple[str, ...]:
    rows = tuple(str(item) for item in row_ids)
    groups = tuple(str(item) for item in strata)
    if not rows or len(rows) != len(groups) or len(set(rows)) != len(rows):
        raise ObservedMultiomeBaselineError("row IDs or strata differ")
    if any(not item for item in rows) or any(not item for item in groups):
        raise ObservedMultiomeBaselineError("row IDs and strata must be non-empty")
    return groups


def _feature_matrix(
    model_id: str,
    *,
    row_ids: Sequence[str],
    strata: Sequence[str],
    rna: np.ndarray | None,
    sequence: np.ndarray | None,
    observed_atac: np.ndarray | None,
    seed: int,
    training_pool: np.ndarray | None = None,
    training_pool_strata: Sequence[str] = (),
) -> np.ndarray:
    if model_id not in MODEL_IDS:
        raise ObservedMultiomeBaselineError(f"unsupported model_id: {model_id}")
    groups = _validate_rows(row_ids, strata)
    rows = len(row_ids)
    if model_id == "masked_modality":
        return np.concatenate(
            (
                _matrix(rna, rows=rows, label="RNA"),
                _matrix(sequence, rows=rows, label="sequence"),
            ),
            axis=1,
        )
    if model_id == "rna_only":
        return _matrix(rna, rows=rows, label="RNA")
    if model_id in {"observed_atac_glm", "observed_atac_only"}:
        atac = _matrix(observed_atac, rows=rows, label="observed ATAC outside mask")
        if np.any(atac < 0):
            raise ObservedMultiomeBaselineError("observed ATAC counts must be nonnegative")
        return atac

    atac = _matrix(observed_atac, rows=rows, label="observed ATAC outside mask")
    if training_pool is None:
        output = np.empty_like(atac)
        for group in sorted(set(groups)):
            indices = [index for index, value in enumerate(groups) if value == group]
            if len(indices) < 2:
                raise ObservedMultiomeBaselineError(
                    "shuffled-modality fit requires at least two rows per stratum"
                )
            ranked = sorted(
                indices,
                key=lambda index: sha256(
                    f"{seed}\0fit\0{row_ids[index]}\0{group}".encode("utf-8")
                ).digest(),
            )
            for left, right in zip(ranked, ranked[1:] + ranked[:1]):
                output[left] = atac[right]
        return output

    pool = np.asarray(training_pool, dtype=np.float64)
    pool_groups = tuple(str(item) for item in training_pool_strata)
    if pool.ndim != 2 or pool.shape[0] != len(pool_groups) or pool.shape[1] != atac.shape[1]:
        raise ObservedMultiomeBaselineError("training-frozen shuffle pool differs")
    output = np.empty_like(atac)
    by_group = {
        group: np.asarray([index for index, value in enumerate(pool_groups) if value == group])
        for group in sorted(set(pool_groups))
    }
    for index, (row_id, group) in enumerate(zip(row_ids, groups)):
        eligible = by_group.get(group)
        if eligible is None or eligible.size == 0:
            raise ObservedMultiomeBaselineError("query stratum is absent from the training pool")
        digest = sha256(f"{seed}\0predict\0{row_id}\0{group}".encode("utf-8")).digest()
        output[index] = pool[int(eligible[int.from_bytes(digest[:8], "big") % eligible.size])]
    return output


def _standardize(features: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = features.mean(axis=0)
    scale = features.std(axis=0)
    scale[scale < 1e-12] = 1.0
    return (features - mean) / scale, mean, scale


def _ridge_fit(features: np.ndarray, target: np.ndarray, alpha: float) -> tuple[np.ndarray, np.ndarray]:
    design = np.column_stack((np.ones(features.shape[0]), features))
    penalty = np.eye(design.shape[1]) * alpha
    penalty[0, 0] = 0.0
    coefficients = np.linalg.solve(
        design.T @ design + penalty,
        design.T @ np.log1p(target),
    )
    return coefficients[1:], coefficients[0]


def _poisson_fit(
    features: np.ndarray,
    target: np.ndarray,
    *,
    alpha: float,
    max_iter: int,
    tolerance: float,
) -> tuple[np.ndarray, np.ndarray]:
    design = np.column_stack((np.ones(features.shape[0]), features))
    penalty = np.eye(design.shape[1]) * alpha
    penalty[0, 0] = 0.0
    fitted = np.zeros((design.shape[1], target.shape[1]), dtype=np.float64)
    for column in range(target.shape[1]):
        response = target[:, column]
        beta = np.zeros(design.shape[1], dtype=np.float64)
        beta[0] = np.log(max(float(response.mean()), 1e-6))
        for _ in range(max_iter):
            eta = np.clip(design @ beta, -12.0, 12.0)
            mean = np.exp(eta)
            gradient = design.T @ (mean - response) + penalty @ beta
            hessian = design.T @ (mean[:, None] * design) + penalty
            step = np.linalg.solve(hessian, gradient)
            beta -= step
            if float(np.max(np.abs(step))) < tolerance:
                break
        else:
            raise ObservedMultiomeBaselineError("Poisson GLM did not converge")
        fitted[:, column] = beta
    return fitted[1:], fitted[0]


def fit(
    model_id: str,
    *,
    row_ids: Sequence[str],
    strata: Sequence[str],
    rna: np.ndarray | None,
    sequence: np.ndarray | None,
    observed_atac: np.ndarray | None,
    target_counts: np.ndarray,
    seed: int,
    ridge_alpha: float = 1.0,
    poisson_alpha: float = 0.5,
    poisson_max_iter: int = 100,
    poisson_tolerance: float = 1e-8,
) -> BaselineState:
    """Fit one baseline using outer-training rows and targets only."""

    target = _matrix(target_counts, rows=len(row_ids), label="training target")
    if np.any(target < 0):
        raise ObservedMultiomeBaselineError("training target counts must be nonnegative")
    features = _feature_matrix(
        model_id,
        row_ids=row_ids,
        strata=strata,
        rna=rna,
        sequence=sequence,
        observed_atac=observed_atac,
        seed=seed,
    )
    standardized, mean, scale = _standardize(features)
    if model_id == "observed_atac_glm":
        coefficients, intercept = _poisson_fit(
            standardized,
            target,
            alpha=poisson_alpha,
            max_iter=poisson_max_iter,
            tolerance=poisson_tolerance,
        )
        link = "log"
    else:
        coefficients, intercept = _ridge_fit(standardized, target, ridge_alpha)
        link = "log1p_identity"
    pool = None
    pool_strata: tuple[str, ...] = ()
    if model_id == "shuffled_modality":
        pool = _matrix(
            observed_atac,
            rows=len(row_ids),
            label="observed ATAC outside mask",
        ).copy()
        pool_strata = tuple(str(item) for item in strata)
    return BaselineState(
        model_id=model_id,
        seed=int(seed),
        feature_mean=mean,
        feature_scale=scale,
        coefficients=coefficients,
        intercept=intercept,
        link=link,
        training_pool=pool,
        training_pool_strata=pool_strata,
    )


def predict(
    state: BaselineState,
    *,
    row_ids: Sequence[str],
    strata: Sequence[str],
    rna: np.ndarray | None,
    sequence: np.ndarray | None,
    observed_atac: np.ndarray | None,
) -> np.ndarray:
    """Predict masked target counts without accepting a target argument."""

    features = _feature_matrix(
        state.model_id,
        row_ids=row_ids,
        strata=strata,
        rna=rna,
        sequence=sequence,
        observed_atac=observed_atac,
        seed=state.seed,
        training_pool=state.training_pool,
        training_pool_strata=state.training_pool_strata,
    )
    if features.shape[1] != state.feature_mean.shape[0]:
        raise ObservedMultiomeBaselineError("query feature axis differs from training")
    linear = ((features - state.feature_mean) / state.feature_scale) @ state.coefficients + state.intercept
    if state.link == "log":
        output = np.exp(np.clip(linear, -12.0, 12.0))
    elif state.link == "log1p_identity":
        output = np.expm1(np.clip(linear, 0.0, 12.0))
    else:
        raise ObservedMultiomeBaselineError("unknown fitted link")
    if output.ndim != 2 or not np.isfinite(output).all() or np.any(output < 0):
        raise ObservedMultiomeBaselineError("prediction output is invalid")
    return output


def export(state: BaselineState, path: Path) -> None:
    """Export a project-owned state as non-pickle arrays plus JSON metadata."""

    path.mkdir(parents=True, exist_ok=False)
    arrays = {
        "feature_mean": state.feature_mean,
        "feature_scale": state.feature_scale,
        "coefficients": state.coefficients,
        "intercept": state.intercept,
    }
    if state.training_pool is not None:
        arrays["training_pool"] = state.training_pool
    with (path / "state.npz").open("xb") as handle:
        np.savez_compressed(handle, **arrays)
    metadata = {
        "schema_version": "masld-bench-observed-multiome-baseline-state-v1",
        "model_id": state.model_id,
        "seed": state.seed,
        "link": state.link,
        "training_pool_strata": list(state.training_pool_strata),
        "pickle_used": False,
        "benchmark_metrics_calculated": False,
    }
    with (path / "metadata.json").open("x", encoding="utf-8") as handle:
        json.dump(metadata, handle, sort_keys=True, separators=(",", ":"))
        handle.write("\n")
