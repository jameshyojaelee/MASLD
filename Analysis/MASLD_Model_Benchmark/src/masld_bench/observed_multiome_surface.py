"""Fold-local preprocessing and descriptive metrics for observed multiome."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from scipy import sparse
from scipy.optimize import nnls
from sklearn.decomposition import NMF, TruncatedSVD


class ObservedMultiomeSurfaceError(ValueError):
    """Raised when a fold-local axis, transform, or profile differs."""


@dataclass(frozen=True)
class SurfacePreprocessState:
    rna_components: np.ndarray
    atac_components: np.ndarray
    target_mean: np.ndarray
    target_scale: np.ndarray
    target_components: np.ndarray
    components: int
    seed: int


def read_csr(group: Any) -> sparse.csr_matrix:
    shape = tuple(int(value) for value in group["shape"][:])
    matrix = sparse.csr_matrix(
        (
            np.asarray(group["data"][:]),
            np.asarray(group["indices"][:], dtype=np.int64),
            np.asarray(group["indptr"][:], dtype=np.int64),
        ),
        shape=shape,
        dtype=np.float64,
    )
    if matrix.ndim != 2 or np.any(~np.isfinite(matrix.data)) or np.any(matrix.data < 0):
        raise ObservedMultiomeSurfaceError("count matrix differs")
    return matrix


def decode(values: Any, *, unique: bool) -> list[str]:
    output = [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]
    if not output or any(not value for value in output) or (unique and len(set(output)) != len(output)):
        raise ObservedMultiomeSurfaceError("string axis differs")
    return output


def log_cpm(matrix: sparse.spmatrix, scale: float = 10_000.0) -> tuple[sparse.csr_matrix, int]:
    value = sparse.csr_matrix(matrix, dtype=np.float64, copy=True)
    if value.ndim != 2 or value.shape[0] < 2 or value.shape[1] < 2 or np.any(value.data < 0):
        raise ObservedMultiomeSurfaceError("log-CPM input differs")
    totals = np.asarray(value.sum(axis=1)).ravel()
    zero = int(np.count_nonzero(totals == 0))
    factors = np.divide(scale, totals, out=np.zeros_like(totals), where=totals > 0)
    value = sparse.diags(factors) @ value
    value.data = np.log1p(value.data)
    if np.any(~np.isfinite(value.data)):
        raise ObservedMultiomeSurfaceError("log-CPM output differs")
    return value.tocsr(), zero


def fit_sparse_reduction(matrix: sparse.spmatrix, *, components: int, seed: int, iterations: int) -> tuple[np.ndarray, np.ndarray]:
    value = sparse.csr_matrix(matrix, dtype=np.float64)
    if components < 1 or components >= min(value.shape) or iterations < 1:
        raise ObservedMultiomeSurfaceError("sparse reduction contract differs")
    estimator = TruncatedSVD(n_components=components, n_iter=iterations, random_state=seed, algorithm="randomized")
    transformed = estimator.fit_transform(value)
    basis = np.asarray(estimator.components_, dtype=np.float64)
    if transformed.shape != (value.shape[0], components) or basis.shape != (components, value.shape[1]) or not np.isfinite(transformed).all() or not np.isfinite(basis).all():
        raise ObservedMultiomeSurfaceError("sparse reduction output differs")
    return transformed, basis


def apply_sparse_reduction(matrix: sparse.spmatrix, basis: np.ndarray) -> np.ndarray:
    value = sparse.csr_matrix(matrix, dtype=np.float64)
    components = np.asarray(basis, dtype=np.float64)
    if components.ndim != 2 or components.shape[1] != value.shape[1] or not np.isfinite(components).all():
        raise ObservedMultiomeSurfaceError("sparse reduction basis differs")
    output = np.asarray(value @ components.T, dtype=np.float64)
    if output.shape != (value.shape[0], components.shape[0]) or not np.isfinite(output).all():
        raise ObservedMultiomeSurfaceError("sparse reduction transform differs")
    return output


def apply_nonnegative_reduction(matrix: sparse.spmatrix, basis: np.ndarray) -> np.ndarray:
    value = sparse.csr_matrix(matrix, dtype=np.float64)
    components = np.asarray(basis, dtype=np.float64)
    if components.ndim != 2 or components.shape[1] != value.shape[1] or np.any(components < 0) or not np.isfinite(components).all():
        raise ObservedMultiomeSurfaceError("nonnegative reduction basis differs")
    output = np.vstack([nnls(components.T, value.getrow(index).toarray().ravel())[0] for index in range(value.shape[0])])
    if output.shape != (value.shape[0], components.shape[0]) or np.any(output < 0) or not np.isfinite(output).all():
        raise ObservedMultiomeSurfaceError("nonnegative reduction transform differs")
    return output


def fit_nonnegative_reduction(matrix: sparse.spmatrix, *, components: int, seed: int, max_iter: int) -> tuple[np.ndarray, np.ndarray]:
    value = sparse.csr_matrix(matrix, dtype=np.float64)
    if components < 1 or components >= min(value.shape) or max_iter < 1 or np.any(value.data < 0):
        raise ObservedMultiomeSurfaceError("nonnegative reduction contract differs")
    estimator = NMF(n_components=components, init="nndsvda", solver="cd", beta_loss="frobenius", max_iter=max_iter, tol=1e-4, random_state=seed)
    estimator.fit(value)
    basis = np.asarray(estimator.components_, dtype=np.float64)
    transformed = apply_nonnegative_reduction(value, basis)
    if np.any(basis < 0):
        raise ObservedMultiomeSurfaceError("nonnegative reduction fit differs")
    return transformed, basis


def fit_target_reduction(features: np.ndarray, *, components: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    value = np.asarray(features, dtype=np.float64)
    if value.ndim != 2 or components < 1 or components >= min(value.shape) or not np.isfinite(value).all():
        raise ObservedMultiomeSurfaceError("target reduction input differs")
    mean = value.mean(axis=0)
    scale = value.std(axis=0)
    scale[scale < 1e-12] = 1.0
    standardized = (value - mean) / scale
    _, _, right = np.linalg.svd(standardized, full_matrices=False)
    basis = right[:components]
    transformed = standardized @ basis.T
    if transformed.shape != (value.shape[0], components) or not np.isfinite(transformed).all():
        raise ObservedMultiomeSurfaceError("target reduction output differs")
    return transformed, mean, scale, basis


def apply_target_reduction(features: np.ndarray, mean: np.ndarray, scale: np.ndarray, basis: np.ndarray) -> np.ndarray:
    value = np.asarray(features, dtype=np.float64)
    center = np.asarray(mean, dtype=np.float64)
    spread = np.asarray(scale, dtype=np.float64)
    components = np.asarray(basis, dtype=np.float64)
    if value.ndim != 2 or center.shape != spread.shape or center.shape != (value.shape[1],) or components.ndim != 2 or components.shape[1] != value.shape[1] or np.any(spread <= 0):
        raise ObservedMultiomeSurfaceError("target transform contract differs")
    output = ((value - center) / spread) @ components.T
    if not np.isfinite(output).all():
        raise ObservedMultiomeSurfaceError("target transform output differs")
    return output


def export_preprocess(state: SurfacePreprocessState, path: Path, metadata: dict[str, Any]) -> None:
    path.mkdir(parents=True, exist_ok=False)
    with (path / "state.npz").open("xb") as handle:
        np.savez_compressed(
            handle,
            rna_components=state.rna_components,
            atac_components=state.atac_components,
            target_mean=state.target_mean,
            target_scale=state.target_scale,
            target_components=state.target_components,
        )
    payload = {
        "schema_version": "masld-bench-observed-multiome-surface-preprocess-v1",
        "components": state.components,
        "seed": state.seed,
        "pickle_used": False,
        **metadata,
    }
    with (path / "metadata.json").open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, sort_keys=True, separators=(",", ":"))
        handle.write("\n")


def load_preprocess(path: Path) -> tuple[SurfacePreprocessState, dict[str, Any]]:
    try:
        metadata = json.loads((path / "metadata.json").read_text(encoding="utf-8"))
        with np.load(path / "state.npz", allow_pickle=False) as archive:
            if set(archive.files) != {"rna_components", "atac_components", "target_mean", "target_scale", "target_components"}:
                raise ObservedMultiomeSurfaceError("preprocess state roster differs")
            arrays = {name: np.asarray(archive[name], dtype=np.float64).copy() for name in archive.files}
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise ObservedMultiomeSurfaceError("preprocess state is invalid") from error
    if not isinstance(metadata, dict) or metadata.get("schema_version") != "masld-bench-observed-multiome-surface-preprocess-v1" or metadata.get("pickle_used") is not False:
        raise ObservedMultiomeSurfaceError("preprocess metadata differs")
    components = metadata.get("components")
    seed = metadata.get("seed")
    if not isinstance(components, int) or not isinstance(seed, int) or components < 1:
        raise ObservedMultiomeSurfaceError("preprocess identity differs")
    if arrays["rna_components"].shape[0] != components or arrays["atac_components"].shape[0] != components or arrays["target_components"].shape[0] != components or arrays["target_components"].shape[1] != arrays["target_mean"].size or arrays["target_mean"].shape != arrays["target_scale"].shape or np.any(arrays["target_scale"] <= 0) or not all(np.isfinite(value).all() for value in arrays.values()):
        raise ObservedMultiomeSurfaceError("preprocess dimensions differ")
    return SurfacePreprocessState(components=components, seed=seed, **arrays), metadata


def profile_deviance_skill(observed: np.ndarray, predicted: np.ndarray, *, pseudocount: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    counts = np.asarray(observed, dtype=np.float64)
    estimate = np.asarray(predicted, dtype=np.float64)
    if counts.ndim != 2 or counts.shape != estimate.shape or counts.shape[1] < 2 or np.any(counts < 0) or np.any(estimate < 0) or not np.isfinite(counts).all() or not np.isfinite(estimate).all() or pseudocount <= 0:
        raise ObservedMultiomeSurfaceError("profile deviance input differs")
    totals = counts.sum(axis=1)
    evaluable = totals > 0
    probabilities = (estimate + pseudocount) / (estimate.sum(axis=1, keepdims=True) + pseudocount * estimate.shape[1])
    uniform = 1.0 / estimate.shape[1]
    model = np.zeros(counts.shape[0], dtype=np.float64)
    null = np.zeros(counts.shape[0], dtype=np.float64)
    for index in np.flatnonzero(evaluable):
        positive = counts[index] > 0
        model[index] = 2.0 * np.sum(counts[index, positive] * np.log(counts[index, positive] / (totals[index] * probabilities[index, positive])))
        null[index] = 2.0 * np.sum(counts[index, positive] * np.log(counts[index, positive] / (totals[index] * uniform)))
    valid = evaluable & (null > 0)
    skill = np.full(counts.shape[0], np.nan, dtype=np.float64)
    skill[valid] = 1.0 - model[valid] / null[valid]
    return skill, model, null
