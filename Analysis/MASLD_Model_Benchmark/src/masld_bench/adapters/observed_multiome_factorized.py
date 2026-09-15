"""Target-generalizing factorized baselines for observed multiome.

The models learn one shared coefficient vector over row features, target
features, and their interactions. They cannot learn a per-target parameter.
Prediction rejects every target ID seen during fit and never accepts outcomes.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Sequence

import numpy as np


MODEL_IDS = (
    "masked_modality",
    "observed_atac_glm",
    "observed_atac_only",
    "rna_only",
    "shuffled_modality",
)


class FactorizedBaselineError(ValueError):
    """Raised when target, row, missingness, or crossed-fold evidence differs."""


@dataclass(frozen=True)
class FactorizedState:
    model_id: str
    seed: int
    row_mean: np.ndarray
    row_scale: np.ndarray
    target_mean: np.ndarray
    target_scale: np.ndarray
    coefficients: np.ndarray
    link: str
    target_mode: str
    training_target_ids: tuple[str, ...]
    training_pool: np.ndarray | None = None
    training_pool_strata: tuple[str, ...] = ()


def _matrix(value: np.ndarray | None, rows: int, label: str) -> np.ndarray:
    if value is None:
        raise FactorizedBaselineError(f"{label} is structurally missing")
    array = np.asarray(value, dtype=np.float64)
    if array.ndim != 2 or array.shape[0] != rows or array.shape[1] < 1:
        raise FactorizedBaselineError(f"{label} shape differs")
    if not np.isfinite(array).all():
        raise FactorizedBaselineError(f"{label} contains non-finite values")
    return array


def _identifiers(values: Sequence[str], label: str) -> tuple[str, ...]:
    identifiers = tuple(str(item) for item in values)
    if not identifiers or any(not item for item in identifiers) or len(set(identifiers)) != len(identifiers):
        raise FactorizedBaselineError(f"{label} identifiers differ")
    return identifiers


def _strata(row_ids: Sequence[str], strata: Sequence[str]) -> tuple[str, ...]:
    rows = _identifiers(row_ids, "row")
    groups = tuple(str(item) for item in strata)
    if len(groups) != len(rows) or any(not item for item in groups):
        raise FactorizedBaselineError("row strata differ")
    return groups


def _shuffled_fit(
    observed_atac: np.ndarray,
    *,
    row_ids: Sequence[str],
    strata: Sequence[str],
    seed: int,
) -> np.ndarray:
    groups = _strata(row_ids, strata)
    output = np.empty_like(observed_atac)
    for group in sorted(set(groups)):
        indices = [index for index, value in enumerate(groups) if value == group]
        if len(indices) < 2:
            raise FactorizedBaselineError("shuffle fit requires two rows per stratum")
        ranked = sorted(
            indices,
            key=lambda index: sha256(
                f"{seed}\0fit\0{row_ids[index]}\0{group}".encode("utf-8")
            ).digest(),
        )
        for left, right in zip(ranked, ranked[1:] + ranked[:1]):
            output[left] = observed_atac[right]
    return output


def _shuffled_query(
    training_pool: np.ndarray,
    training_pool_strata: Sequence[str],
    *,
    row_ids: Sequence[str],
    strata: Sequence[str],
    seed: int,
) -> np.ndarray:
    groups = _strata(row_ids, strata)
    pool_groups = tuple(str(item) for item in training_pool_strata)
    if training_pool.ndim != 2 or training_pool.shape[0] != len(pool_groups):
        raise FactorizedBaselineError("training shuffle pool differs")
    by_group = {
        group: np.asarray([index for index, value in enumerate(pool_groups) if value == group])
        for group in sorted(set(pool_groups))
    }
    output = np.empty((len(row_ids), training_pool.shape[1]), dtype=np.float64)
    for index, (row_id, group) in enumerate(zip(row_ids, groups)):
        eligible = by_group.get(group)
        if eligible is None or eligible.size == 0:
            raise FactorizedBaselineError("query stratum is absent from training pool")
        digest = sha256(f"{seed}\0predict\0{row_id}\0{group}".encode("utf-8")).digest()
        output[index] = training_pool[int(eligible[int.from_bytes(digest[:8], "big") % eligible.size])]
    return output


def _row_features(
    model_id: str,
    *,
    row_ids: Sequence[str],
    strata: Sequence[str],
    rna: np.ndarray | None,
    observed_atac: np.ndarray | None,
    seed: int,
    training_pool: np.ndarray | None = None,
    training_pool_strata: Sequence[str] = (),
) -> np.ndarray:
    rows = len(_identifiers(row_ids, "row"))
    _strata(row_ids, strata)
    if model_id in {"masked_modality", "rna_only"}:
        return _matrix(rna, rows, "RNA context")
    atac = _matrix(observed_atac, rows, "observed ATAC outside mask and buffer")
    if np.any(atac < 0):
        raise FactorizedBaselineError("observed ATAC must be nonnegative")
    if model_id != "shuffled_modality":
        return atac
    if training_pool is None:
        return _shuffled_fit(atac, row_ids=row_ids, strata=strata, seed=seed)
    return _shuffled_query(
        training_pool,
        training_pool_strata,
        row_ids=row_ids,
        strata=strata,
        seed=seed,
    )


def _target_features(
    model_id: str,
    *,
    target_ids: Sequence[str],
    sequence_features: np.ndarray | None,
    target_covariates: np.ndarray | None,
) -> tuple[np.ndarray, str]:
    targets = _identifiers(target_ids, "target")
    if model_id in {"masked_modality", "observed_atac_glm"}:
        return _matrix(sequence_features, len(targets), "target sequence features"), "sequence"
    return _matrix(target_covariates, len(targets), "target covariates"), "covariate"


def _standardize(features: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = features.mean(axis=0)
    scale = features.std(axis=0)
    scale[scale < 1e-12] = 1.0
    return (features - mean) / scale, mean, scale


def _design(rows: np.ndarray, targets: np.ndarray, maximum_pairs: int) -> np.ndarray:
    pair_count = rows.shape[0] * targets.shape[0]
    if pair_count < 1 or pair_count > maximum_pairs:
        raise FactorizedBaselineError("expanded pair count exceeds its fixture ceiling")
    row_expanded = np.repeat(rows, targets.shape[0], axis=0)
    target_expanded = np.tile(targets, (rows.shape[0], 1))
    interactions = np.einsum("ni,nj->nij", row_expanded, target_expanded).reshape(pair_count, -1)
    return np.column_stack(
        (np.ones(pair_count), row_expanded, target_expanded, interactions)
    )


def _ridge(design: np.ndarray, target: np.ndarray, alpha: float) -> np.ndarray:
    penalty = np.eye(design.shape[1]) * alpha
    penalty[0, 0] = 0.0
    return np.linalg.solve(design.T @ design + penalty, design.T @ np.log1p(target))


def _poisson(
    design: np.ndarray,
    target: np.ndarray,
    *,
    alpha: float,
    max_iter: int,
    tolerance: float,
) -> np.ndarray:
    penalty = np.eye(design.shape[1]) * alpha
    penalty[0, 0] = 0.0
    beta = np.zeros(design.shape[1], dtype=np.float64)
    beta[0] = np.log(max(float(target.mean()), 1e-6))
    for _ in range(max_iter):
        eta = np.clip(design @ beta, -12.0, 12.0)
        mean = np.exp(eta)
        gradient = design.T @ (mean - target) + penalty @ beta
        hessian = design.T @ (mean[:, None] * design) + penalty
        step = np.linalg.solve(hessian, gradient)
        beta -= step
        if float(np.max(np.abs(step))) < tolerance:
            return beta
    raise FactorizedBaselineError("factorized Poisson GLM did not converge")


def fit(
    model_id: str,
    *,
    row_ids: Sequence[str],
    strata: Sequence[str],
    target_ids: Sequence[str],
    rna: np.ndarray | None,
    observed_atac: np.ndarray | None,
    sequence_features: np.ndarray | None,
    target_covariates: np.ndarray | None,
    target_counts: np.ndarray,
    seed: int,
    ridge_alpha: float = 1.0,
    poisson_alpha: float = 0.5,
    poisson_max_iter: int = 100,
    poisson_tolerance: float = 1e-8,
    maximum_expanded_pairs: int = 100_000,
) -> FactorizedState:
    if model_id not in MODEL_IDS:
        raise FactorizedBaselineError(f"unsupported model_id: {model_id}")
    rows = _identifiers(row_ids, "row")
    targets = _identifiers(target_ids, "target")
    response = _matrix(target_counts, len(rows), "training target count matrix")
    if response.shape[1] != len(targets) or np.any(response < 0):
        raise FactorizedBaselineError("training target matrix differs")
    row_features = _row_features(
        model_id,
        row_ids=rows,
        strata=strata,
        rna=rna,
        observed_atac=observed_atac,
        seed=seed,
    )
    target_features, target_mode = _target_features(
        model_id,
        target_ids=targets,
        sequence_features=sequence_features,
        target_covariates=target_covariates,
    )
    row_scaled, row_mean, row_scale = _standardize(row_features)
    target_scaled, target_mean, target_scale = _standardize(target_features)
    design = _design(row_scaled, target_scaled, maximum_expanded_pairs)
    flattened = response.reshape(-1)
    if model_id == "observed_atac_glm":
        coefficients = _poisson(
            design,
            flattened,
            alpha=poisson_alpha,
            max_iter=poisson_max_iter,
            tolerance=poisson_tolerance,
        )
        link = "log"
    else:
        coefficients = _ridge(design, flattened, ridge_alpha)
        link = "log1p_identity"
    pool = None
    pool_strata: tuple[str, ...] = ()
    if model_id == "shuffled_modality":
        pool = _matrix(observed_atac, len(rows), "observed ATAC outside mask and buffer").copy()
        pool_strata = tuple(str(item) for item in strata)
    return FactorizedState(
        model_id=model_id,
        seed=int(seed),
        row_mean=row_mean,
        row_scale=row_scale,
        target_mean=target_mean,
        target_scale=target_scale,
        coefficients=coefficients,
        link=link,
        target_mode=target_mode,
        training_target_ids=targets,
        training_pool=pool,
        training_pool_strata=pool_strata,
    )


def predict(
    state: FactorizedState,
    *,
    row_ids: Sequence[str],
    strata: Sequence[str],
    target_ids: Sequence[str],
    rna: np.ndarray | None,
    observed_atac: np.ndarray | None,
    sequence_features: np.ndarray | None,
    target_covariates: np.ndarray | None,
    maximum_expanded_pairs: int = 100_000,
) -> np.ndarray:
    rows = _identifiers(row_ids, "row")
    targets = _identifiers(target_ids, "target")
    overlap = sorted(set(targets).intersection(state.training_target_ids))
    if overlap:
        raise FactorizedBaselineError("prediction target overlaps a fitted genomic target")
    row_features = _row_features(
        state.model_id,
        row_ids=rows,
        strata=strata,
        rna=rna,
        observed_atac=observed_atac,
        seed=state.seed,
        training_pool=state.training_pool,
        training_pool_strata=state.training_pool_strata,
    )
    target_features, target_mode = _target_features(
        state.model_id,
        target_ids=targets,
        sequence_features=sequence_features,
        target_covariates=target_covariates,
    )
    if target_mode != state.target_mode:
        raise FactorizedBaselineError("prediction target feature mode differs")
    if row_features.shape[1] != state.row_mean.shape[0] or target_features.shape[1] != state.target_mean.shape[0]:
        raise FactorizedBaselineError("prediction feature axis differs")
    design = _design(
        (row_features - state.row_mean) / state.row_scale,
        (target_features - state.target_mean) / state.target_scale,
        maximum_expanded_pairs,
    )
    linear = design @ state.coefficients
    if state.link == "log":
        output = np.exp(np.clip(linear, -12.0, 12.0))
    elif state.link == "log1p_identity":
        output = np.expm1(np.clip(linear, 0.0, 12.0))
    else:
        raise FactorizedBaselineError("fitted link differs")
    result = output.reshape(len(rows), len(targets))
    if not np.isfinite(result).all() or np.any(result < 0):
        raise FactorizedBaselineError("factorized prediction is invalid")
    return result


def export(state: FactorizedState, path: Path) -> None:
    path.mkdir(parents=True, exist_ok=False)
    arrays = {
        "row_mean": state.row_mean,
        "row_scale": state.row_scale,
        "target_mean": state.target_mean,
        "target_scale": state.target_scale,
        "coefficients": state.coefficients,
        "training_target_ids": np.asarray(state.training_target_ids, dtype=str),
    }
    if state.training_pool is not None:
        arrays["training_pool"] = state.training_pool
    with (path / "state.npz").open("xb") as handle:
        np.savez_compressed(handle, **arrays)
    target_hash = sha256(
        "\n".join(sorted(state.training_target_ids)).encode("utf-8")
    ).hexdigest()
    metadata = {
        "schema_version": "masld-bench-observed-multiome-factorized-state-v1",
        "model_id": state.model_id,
        "seed": state.seed,
        "link": state.link,
        "target_mode": state.target_mode,
        "training_target_count": len(state.training_target_ids),
        "training_target_set_sha256": target_hash,
        "training_pool_strata": list(state.training_pool_strata),
        "target_specific_coefficients": False,
        "pickle_used": False,
        "benchmark_metrics_calculated": False,
    }
    with (path / "metadata.json").open("x", encoding="utf-8") as handle:
        json.dump(metadata, handle, sort_keys=True, separators=(",", ":"))
        handle.write("\n")


def load(path: Path) -> FactorizedState:
    """Load a non-pickle state while preserving the fitted-target separation."""
    metadata_path = path / "metadata.json"
    state_path = path / "state.npz"
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise FactorizedBaselineError("factorized state metadata is invalid") from error
    if not isinstance(metadata, dict) or metadata.get("schema_version") != "masld-bench-observed-multiome-factorized-state-v1":
        raise FactorizedBaselineError("factorized state schema differs")
    model_id = metadata.get("model_id")
    if model_id not in MODEL_IDS or metadata.get("pickle_used") is not False:
        raise FactorizedBaselineError("factorized state identity differs")
    required = {
        "row_mean",
        "row_scale",
        "target_mean",
        "target_scale",
        "coefficients",
        "training_target_ids",
    }
    try:
        with np.load(state_path, allow_pickle=False) as archive:
            names = set(archive.files)
            if names not in (required, required | {"training_pool"}):
                raise FactorizedBaselineError("factorized state array roster differs")
            arrays = {name: np.asarray(archive[name]).copy() for name in names}
    except (OSError, ValueError) as error:
        raise FactorizedBaselineError("factorized state arrays are invalid") from error
    row_mean = np.asarray(arrays["row_mean"], dtype=np.float64)
    row_scale = np.asarray(arrays["row_scale"], dtype=np.float64)
    target_mean = np.asarray(arrays["target_mean"], dtype=np.float64)
    target_scale = np.asarray(arrays["target_scale"], dtype=np.float64)
    coefficients = np.asarray(arrays["coefficients"], dtype=np.float64)
    if (
        row_mean.ndim != 1
        or row_mean.shape != row_scale.shape
        or target_mean.ndim != 1
        or target_mean.shape != target_scale.shape
        or coefficients.ndim != 1
        or not all(np.isfinite(value).all() for value in (row_mean, row_scale, target_mean, target_scale, coefficients))
        or np.any(row_scale <= 0)
        or np.any(target_scale <= 0)
        or coefficients.size != 1 + row_mean.size + target_mean.size + row_mean.size * target_mean.size
    ):
        raise FactorizedBaselineError("factorized state dimensions differ")
    training_target_ids = _identifiers(
        [str(value) for value in arrays["training_target_ids"].tolist()],
        "training target",
    )
    if (
        len(training_target_ids) != metadata.get("training_target_count")
        or sha256("\n".join(sorted(training_target_ids)).encode("utf-8")).hexdigest()
        != metadata.get("training_target_set_sha256")
    ):
        raise FactorizedBaselineError("factorized training-target identity differs")
    strata_value = metadata.get("training_pool_strata")
    if not isinstance(strata_value, list) or any(not isinstance(value, str) or not value for value in strata_value):
        raise FactorizedBaselineError("factorized training-pool strata differ")
    pool = arrays.get("training_pool")
    if model_id == "shuffled_modality":
        if pool is None or pool.ndim != 2 or pool.shape != (len(strata_value), row_mean.size):
            raise FactorizedBaselineError("shuffled training pool differs")
        pool = np.asarray(pool, dtype=np.float64)
    elif pool is not None or strata_value:
        raise FactorizedBaselineError("unexpected training pool in factorized state")
    link = metadata.get("link")
    target_mode = metadata.get("target_mode")
    if link not in {"log", "log1p_identity"} or target_mode not in {"sequence", "covariate"}:
        raise FactorizedBaselineError("factorized state link or target mode differs")
    seed = metadata.get("seed")
    if not isinstance(seed, int) or isinstance(seed, bool):
        raise FactorizedBaselineError("factorized state seed differs")
    return FactorizedState(
        model_id=model_id,
        seed=seed,
        row_mean=row_mean,
        row_scale=row_scale,
        target_mean=target_mean,
        target_scale=target_scale,
        coefficients=coefficients,
        link=link,
        target_mode=target_mode,
        training_target_ids=training_target_ids,
        training_pool=pool,
        training_pool_strata=tuple(strata_value),
    )
