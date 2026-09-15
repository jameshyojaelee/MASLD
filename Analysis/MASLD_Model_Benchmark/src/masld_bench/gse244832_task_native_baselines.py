"""Source-only task-native baselines for blind GSE244832 ATAC transport.

The estimators in this module never accept target-role GSE244832 values.  They
fit on donor-lineage pseudobulk counts from GSE296875 and transform only the
observed rotation-specific context ATAC at prediction time.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np


MODEL_IDS = ("lsi", "observed_atac_glm", "observed_atac_only")
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage")
STATE_SCHEMA = "masld-bench-gse244832-task-native-baseline-state-v1"


class TaskNativeBaselineError(ValueError):
    """Raised when a source fit or blind query does not meet the frozen requirements."""


@dataclass(frozen=True)
class BaselineState:
    model_id: str
    arrays: Mapping[str, np.ndarray]
    metadata: Mapping[str, Any]


def _counts(value: np.ndarray, *, rows: int | None = None, label: str) -> np.ndarray:
    array = np.asarray(value, dtype=np.float64)
    if (
        array.ndim != 2
        or array.shape[0] < 1
        or array.shape[1] < 2
        or rows is not None
        and array.shape[0] != rows
        or not np.isfinite(array).all()
        or np.any(array < 0)
    ):
        raise TaskNativeBaselineError(f"{label} count matrix differs")
    totals = array.sum(axis=1)
    if np.any(totals <= 0):
        raise TaskNativeBaselineError(f"{label} contains an observed zero-library row")
    return array


def _lineages(
    values: Sequence[int], rows: int, *, require_all: bool
) -> np.ndarray:
    array = np.asarray(values, dtype=np.int64)
    observed = set(np.unique(array))
    allowed = set(range(len(LINEAGES)))
    if (
        array.shape != (rows,)
        or not observed
        or not observed.issubset(allowed)
        or require_all
        and observed != allowed
    ):
        raise TaskNativeBaselineError("lineage axis differs")
    return array


def _one_hot(lineages: np.ndarray) -> np.ndarray:
    return np.eye(len(LINEAGES), dtype=np.float64)[lineages]


def _log_cpm(counts: np.ndarray, scale: float = 10_000.0) -> np.ndarray:
    totals = counts.sum(axis=1)
    return np.log1p(scale * counts / totals[:, None])


def _fit_depth(context: np.ndarray, target: np.ndarray, lineages: np.ndarray) -> np.ndarray:
    design = np.column_stack(
        (
            np.ones(context.shape[0], dtype=np.float64),
            _one_hot(lineages)[:, 1:],
            np.log1p(context.sum(axis=1)),
        )
    )
    response = np.log1p(target.sum(axis=1))
    penalty = np.eye(design.shape[1], dtype=np.float64) * 1.0e-3
    penalty[0, 0] = 0.0
    return np.linalg.solve(design.T @ design + penalty, design.T @ response)


def _predict_depth(
    coefficients: np.ndarray, context: np.ndarray, lineages: np.ndarray
) -> np.ndarray:
    design = np.column_stack(
        (
            np.ones(context.shape[0], dtype=np.float64),
            _one_hot(lineages)[:, 1:],
            np.log1p(context.sum(axis=1)),
        )
    )
    return np.expm1(np.clip(design @ coefficients, 0.0, np.log1p(1.0e12)))


def _lineage_profiles(target: np.ndarray, lineages: np.ndarray) -> np.ndarray:
    output = np.empty((len(LINEAGES), target.shape[1]), dtype=np.float64)
    global_profile = target.sum(axis=0)
    if global_profile.sum() <= 0:
        raise TaskNativeBaselineError("source target is empty")
    global_profile /= global_profile.sum()
    for lineage in range(len(LINEAGES)):
        profile = target[lineages == lineage].sum(axis=0)
        output[lineage] = profile / profile.sum() if profile.sum() > 0 else global_profile
    return output


def _profiles_to_counts(
    log_profiles: np.ndarray,
    depth: np.ndarray,
    lineages: np.ndarray,
    fallback: np.ndarray,
) -> np.ndarray:
    raw = np.expm1(np.clip(log_profiles, 0.0, np.log1p(10_000.0)))
    totals = raw.sum(axis=1)
    for row in np.flatnonzero(totals <= 0):
        raw[row] = fallback[lineages[row]]
    totals = raw.sum(axis=1)
    output = raw / totals[:, None] * depth[:, None]
    if not np.isfinite(output).all() or np.any(output < 0):
        raise TaskNativeBaselineError("predicted count matrix is invalid")
    return output


def _fit_lsi(
    context: np.ndarray,
    target: np.ndarray,
    lineages: np.ndarray,
    *,
    components: int,
    ridge_alpha: float,
) -> BaselineState:
    peak_sum = context.sum(axis=0)
    active = peak_sum > 0
    if active.sum() <= components:
        raise TaskNativeBaselineError("source context has too few nonempty windows for LSI")
    idf = context.shape[0] / peak_sum[active]
    tfidf = np.log1p(
        10_000.0
        * context[:, active]
        / context.sum(axis=1, keepdims=True)
        * idf[None, :]
    )
    _left, _singular, right = np.linalg.svd(tfidf, full_matrices=False)
    basis = right[:components]
    embedding = tfidf @ basis.T
    features = np.column_stack((embedding, _one_hot(lineages)))
    feature_mean = features.mean(axis=0)
    feature_scale = features.std(axis=0)
    feature_scale[feature_scale < 1.0e-12] = 1.0
    design = (features - feature_mean) / feature_scale
    response = _log_cpm(target)
    response_mean = response.mean(axis=0)
    penalty = np.eye(design.shape[1], dtype=np.float64) * ridge_alpha
    coefficients = np.linalg.solve(
        design.T @ design + penalty,
        design.T @ (response - response_mean),
    )
    arrays = {
        "active_context": active.astype(np.uint8),
        "idf": idf,
        "lsi_basis": basis,
        "feature_mean": feature_mean,
        "feature_scale": feature_scale,
        "coefficients": coefficients,
        "response_mean": response_mean,
        "depth_coefficients": _fit_depth(context, target, lineages),
        "fallback_profiles": _lineage_profiles(target, lineages),
    }
    return BaselineState(
        model_id="lsi",
        arrays=arrays,
        metadata={
            "algorithm": "source_fitted_signac_method1_tfidf_exact_svd_ridge_decoder",
            "components": components,
            "ridge_alpha": ridge_alpha,
            "idf_formula": "source_units_divided_by_source_window_sum",
            "empty_source_context_windows": int((~active).sum()),
        },
    )


def _fit_glm(
    context: np.ndarray, target: np.ndarray, lineages: np.ndarray
) -> BaselineState:
    rates = np.empty((len(LINEAGES), target.shape[1]), dtype=np.float64)
    for lineage in range(len(LINEAGES)):
        selected = lineages == lineage
        exposure = context[selected].sum()
        if exposure <= 0:
            raise TaskNativeBaselineError("source lineage has zero context exposure")
        rates[lineage] = target[selected].sum(axis=0) / exposure
    return BaselineState(
        model_id="observed_atac_glm",
        arrays={"lineage_window_rates": rates},
        metadata={
            "algorithm": "poisson_log_link_lineage_window_rate_with_context_library_offset",
            "likelihood": "Poisson",
            "fit": "closed_form_maximum_likelihood",
            "offset": "log_observed_context_atac_library",
        },
    )


def _fit_observed_only(
    context: np.ndarray,
    target: np.ndarray,
    lineages: np.ndarray,
    *,
    ridge_alpha: float,
    clip: float,
) -> BaselineState:
    transformed = _log_cpm(context)
    active = context.sum(axis=0) > 0
    feature_mean = transformed[:, active].mean(axis=0)
    feature_scale = transformed[:, active].std(axis=0)
    feature_scale[feature_scale < 1.0e-3] = 1.0
    standardized = np.clip(
        (transformed[:, active] - feature_mean) / feature_scale,
        -clip,
        clip,
    ) / np.sqrt(float(active.sum()))
    features = np.column_stack((standardized, _one_hot(lineages)))
    feature_center = features.mean(axis=0)
    centered = features - feature_center
    response = _log_cpm(target)
    response_mean = response.mean(axis=0)
    kernel = centered @ centered.T
    dual = np.linalg.solve(
        kernel + ridge_alpha * np.eye(kernel.shape[0], dtype=np.float64),
        response - response_mean,
    )
    arrays = {
        "active_context": active.astype(np.uint8),
        "feature_mean": feature_mean,
        "feature_scale": feature_scale,
        "feature_center": feature_center,
        "source_centered_features": centered,
        "dual_coefficients": dual,
        "response_mean": response_mean,
        "depth_coefficients": _fit_depth(context, target, lineages),
        "fallback_profiles": _lineage_profiles(target, lineages),
    }
    return BaselineState(
        model_id="observed_atac_only",
        arrays=arrays,
        metadata={
            "algorithm": "source_only_full_context_dual_ridge",
            "ridge_alpha": ridge_alpha,
            "source_standardized_clip": clip,
            "empty_source_context_windows": int((~active).sum()),
        },
    )


def fit(
    model_id: str,
    *,
    context_counts: np.ndarray,
    target_counts: np.ndarray,
    lineage_indices: Sequence[int],
    lsi_components: int = 16,
    ridge_alpha: float = 1.0,
    standardized_clip: float = 8.0,
) -> BaselineState:
    """Fit one model using GSE296875 source counts only."""

    if model_id not in MODEL_IDS:
        raise TaskNativeBaselineError(f"unsupported model ID: {model_id}")
    context = _counts(context_counts, label="source context")
    target = _counts(target_counts, rows=context.shape[0], label="source target")
    if context.shape[1] != target.shape[1]:
        raise TaskNativeBaselineError("source context and target window counts differ")
    lineages = _lineages(lineage_indices, context.shape[0], require_all=True)
    if model_id == "lsi":
        return _fit_lsi(
            context,
            target,
            lineages,
            components=lsi_components,
            ridge_alpha=ridge_alpha,
        )
    if model_id == "observed_atac_glm":
        return _fit_glm(context, target, lineages)
    return _fit_observed_only(
        context,
        target,
        lineages,
        ridge_alpha=ridge_alpha,
        clip=standardized_clip,
    )


def predict(
    state: BaselineState,
    *,
    context_counts: np.ndarray,
    lineage_indices: Sequence[int],
) -> np.ndarray:
    """Predict target-role counts from observed context ATAC only."""

    context = _counts(context_counts, label="query context")
    lineages = _lineages(lineage_indices, context.shape[0], require_all=False)
    arrays = state.arrays
    if state.model_id == "observed_atac_glm":
        rates = np.asarray(arrays["lineage_window_rates"], dtype=np.float64)
        output = context.sum(axis=1, keepdims=True) * rates[lineages]
    elif state.model_id == "lsi":
        active = np.asarray(arrays["active_context"], dtype=np.uint8).astype(bool)
        if active.shape != (context.shape[1],):
            raise TaskNativeBaselineError("query context feature axis differs")
        idf = np.asarray(arrays["idf"], dtype=np.float64)
        transformed = np.log1p(
            10_000.0
            * context[:, active]
            / context.sum(axis=1, keepdims=True)
            * idf[None, :]
        )
        embedding = transformed @ np.asarray(arrays["lsi_basis"], dtype=np.float64).T
        features = np.column_stack((embedding, _one_hot(lineages)))
        design = (
            features - np.asarray(arrays["feature_mean"], dtype=np.float64)
        ) / np.asarray(arrays["feature_scale"], dtype=np.float64)
        log_profiles = (
            np.asarray(arrays["response_mean"], dtype=np.float64)
            + design @ np.asarray(arrays["coefficients"], dtype=np.float64)
        )
        output = _profiles_to_counts(
            log_profiles,
            _predict_depth(
                np.asarray(arrays["depth_coefficients"], dtype=np.float64),
                context,
                lineages,
            ),
            lineages,
            np.asarray(arrays["fallback_profiles"], dtype=np.float64),
        )
    elif state.model_id == "observed_atac_only":
        active = np.asarray(arrays["active_context"], dtype=np.uint8).astype(bool)
        if active.shape != (context.shape[1],):
            raise TaskNativeBaselineError("query context feature axis differs")
        clip = float(state.metadata["source_standardized_clip"])
        transformed = _log_cpm(context)[:, active]
        standardized = np.clip(
            (
                transformed
                - np.asarray(arrays["feature_mean"], dtype=np.float64)
            )
            / np.asarray(arrays["feature_scale"], dtype=np.float64),
            -clip,
            clip,
        ) / np.sqrt(float(active.sum()))
        features = np.column_stack((standardized, _one_hot(lineages)))
        centered = features - np.asarray(arrays["feature_center"], dtype=np.float64)
        kernel = centered @ np.asarray(
            arrays["source_centered_features"], dtype=np.float64
        ).T
        log_profiles = (
            np.asarray(arrays["response_mean"], dtype=np.float64)
            + kernel @ np.asarray(arrays["dual_coefficients"], dtype=np.float64)
        )
        output = _profiles_to_counts(
            log_profiles,
            _predict_depth(
                np.asarray(arrays["depth_coefficients"], dtype=np.float64),
                context,
                lineages,
            ),
            lineages,
            np.asarray(arrays["fallback_profiles"], dtype=np.float64),
        )
    else:
        raise TaskNativeBaselineError("fitted model ID differs")
    if output.ndim != 2 or not np.isfinite(output).all() or np.any(output < 0):
        raise TaskNativeBaselineError("blind query prediction is invalid")
    return output


def export_state(state: BaselineState, path: Path) -> None:
    if path.exists():
        raise TaskNativeBaselineError("refusing to overwrite fitted state")
    path.mkdir(parents=True, mode=0o750)
    with (path / "state.npz").open("xb") as handle:
        np.savez_compressed(handle, **state.arrays)
    payload = {
        "schema_version": STATE_SCHEMA,
        "model_id": state.model_id,
        "pickle_used": False,
        "metrics_calculated": False,
        **state.metadata,
    }
    (path / "metadata.json").write_text(
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )


def load_state(path: Path) -> BaselineState:
    metadata = json.loads((path / "metadata.json").read_text(encoding="utf-8"))
    if (
        not isinstance(metadata, dict)
        or metadata.get("schema_version") != STATE_SCHEMA
        or metadata.get("model_id") not in MODEL_IDS
        or metadata.get("pickle_used") is not False
        or metadata.get("metrics_calculated") is not False
    ):
        raise TaskNativeBaselineError("fitted-state metadata differs")
    with np.load(path / "state.npz", allow_pickle=False) as archive:
        arrays = {name: np.asarray(archive[name]).copy() for name in archive.files}
    return BaselineState(
        model_id=str(metadata["model_id"]),
        arrays=arrays,
        metadata={
            key: value
            for key, value in metadata.items()
            if key
            not in {"schema_version", "model_id", "pickle_used", "metrics_calculated"}
        },
    )
