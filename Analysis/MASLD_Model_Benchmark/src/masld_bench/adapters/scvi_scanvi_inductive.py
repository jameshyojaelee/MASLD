"""Fail-closed helpers for inductive scVI/scANVI label transfer.

The held-study object may be registered by scvi-tools for a frozen forward
pass, but it may never be used by ``train`` or a query-adaptation API.  These
helpers hash the learned module before and after every held-study inference.
"""

from __future__ import annotations

from hashlib import sha256
import json
from typing import Any, Sequence

import numpy as np


class InductiveSCVIError(ValueError):
    """Raised when a task-native scVI/scANVI requirement is not met."""


ROSTER = (
    "cholangiocyte",
    "endothelial",
    "hepatocyte",
    "immune",
    "mesenchymal_stromal",
)
UNLABELED_CATEGORY = "__unlabeled__"
KNN_NEIGHBORS = 15


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return sha256(payload).hexdigest()


def validate_integer_like_counts(matrix: Any, *, chunk_size: int = 1_000_000) -> None:
    """Reject negative, nonfinite, or noninteger values before scVI setup."""

    values = matrix.data if hasattr(matrix, "data") else np.asarray(matrix).reshape(-1)
    if chunk_size < 1:
        raise InductiveSCVIError("raw-count validation chunk size differs")
    for start in range(0, len(values), chunk_size):
        block = np.asarray(values[start : start + chunk_size])
        if (
            not np.all(np.isfinite(block))
            or np.any(block < 0)
            or np.any(block != np.floor(block))
        ):
            raise InductiveSCVIError(
                "input is not nonnegative integer-valued raw UMI counts"
            )


def donor_class_weights(
    donors: Sequence[str], labels: Sequence[str]
) -> np.ndarray:
    """Give every observed donor-class combination equal total reference mass."""

    if len(donors) != len(labels) or len(donors) == 0:
        raise InductiveSCVIError("donor and label vectors differ")
    counts: dict[tuple[str, str], int] = {}
    for donor, label in zip(donors, labels, strict=True):
        if label not in ROSTER:
            raise InductiveSCVIError("reference label is outside the frozen roster")
        key = (str(donor), str(label))
        counts[key] = counts.get(key, 0) + 1
    weights = np.asarray(
        [1.0 / counts[(str(donor), str(label))] for donor, label in zip(donors, labels, strict=True)],
        dtype=np.float64,
    )
    weights /= weights.mean()
    return weights


def weighted_knn_probabilities(
    reference_latent: np.ndarray,
    reference_labels: Sequence[str],
    reference_weights: np.ndarray,
    query_latent: np.ndarray,
    *,
    neighbors: int = KNN_NEIGHBORS,
    n_jobs: int = 1,
) -> np.ndarray:
    """Apply the frozen, nonparametric scVI label-transfer rule."""

    from sklearn.neighbors import NearestNeighbors

    reference = np.asarray(reference_latent, dtype=np.float32)
    query = np.asarray(query_latent, dtype=np.float32)
    weights = np.asarray(reference_weights, dtype=np.float64)
    if (
        reference.ndim != 2
        or query.ndim != 2
        or reference.shape[1] != query.shape[1]
        or len(reference) != len(reference_labels)
        or weights.shape != (len(reference),)
        or not np.all(np.isfinite(reference))
        or not np.all(np.isfinite(query))
        or not np.all(np.isfinite(weights))
        or np.any(weights <= 0)
        or not 0 < neighbors <= len(reference)
    ):
        raise InductiveSCVIError("latent label-transfer inputs differ")
    label_to_index = {label: index for index, label in enumerate(ROSTER)}
    try:
        encoded = np.asarray(
            [label_to_index[str(label)] for label in reference_labels],
            dtype=np.int64,
        )
    except KeyError as error:
        raise InductiveSCVIError("reference label is outside the frozen roster") from error
    order = np.lexsort((np.arange(len(reference)), encoded))
    search = NearestNeighbors(
        n_neighbors=neighbors,
        algorithm="auto",
        n_jobs=max(1, int(n_jobs)),
    )
    search.fit(reference[order])
    _, indices = search.kneighbors(query, return_distance=True)
    original = order[indices]
    probabilities = np.zeros((len(query), len(ROSTER)), dtype=np.float64)
    for row_index, nearest in enumerate(original):
        votes = np.bincount(
            encoded[nearest], weights=weights[nearest], minlength=len(ROSTER)
        )
        if not np.isfinite(votes).all() or votes.sum() <= 0:
            raise InductiveSCVIError("kNN votes are invalid")
        probabilities[row_index] = votes / votes.sum()
    return validate_probabilities(probabilities, rows=len(query))


def validate_probabilities(values: Any, *, rows: int) -> np.ndarray:
    result = np.asarray(values, dtype=np.float64)
    if (
        result.shape != (rows, len(ROSTER))
        or not np.all(np.isfinite(result))
        or np.any(result < 0)
    ):
        raise InductiveSCVIError("classifier probabilities differ")
    totals = result.sum(axis=1)
    if np.any(totals <= 0) or not np.allclose(
        totals, 1.0, rtol=1.0e-5, atol=1.0e-7
    ):
        raise InductiveSCVIError("classifier probabilities do not sum to one")
    result /= totals[:, np.newaxis]
    return result


def module_state_sha256(model: Any) -> str:
    """Hash every learned tensor without serializing executable objects."""

    digest = sha256()
    state = model.module.state_dict()
    for key in sorted(state):
        tensor = state[key].detach().cpu().contiguous()
        array = tensor.numpy()
        digest.update(key.encode("utf-8"))
        digest.update(b"\0")
        digest.update(str(array.dtype).encode("ascii"))
        digest.update(b"\0")
        digest.update(json.dumps(list(array.shape), separators=(",", ":")).encode("ascii"))
        digest.update(b"\0")
        digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def direct_latent_inference(model: Any, query: Any) -> tuple[np.ndarray, str]:
    """Run a direct frozen SCVI encoding and prove learned state did not change."""

    before = module_state_sha256(model)
    latent = np.asarray(model.get_latent_representation(query), dtype=np.float32)
    after = module_state_sha256(model)
    if before != after:
        raise InductiveSCVIError("held-study latent inference changed learned state")
    if latent.ndim != 2 or len(latent) != query.n_obs or not np.all(np.isfinite(latent)):
        raise InductiveSCVIError("held-study latent representation differs")
    return latent, before


def direct_scanvi_probabilities(model: Any, query: Any) -> tuple[np.ndarray, str]:
    """Run the native frozen scANVI classifier with no query training."""

    before = module_state_sha256(model)
    predicted = model.predict(query, soft=True)
    after = module_state_sha256(model)
    if before != after:
        raise InductiveSCVIError("held-study scANVI inference changed learned state")
    if not hasattr(predicted, "columns"):
        raise InductiveSCVIError("scANVI soft prediction lacks named classes")
    columns = tuple(map(str, predicted.columns))
    if set(columns) != set(ROSTER):
        raise InductiveSCVIError(
            f"scANVI class roster differs: observed {columns!r}"
        )
    values = np.asarray(predicted.loc[:, list(ROSTER)], dtype=np.float64)
    return validate_probabilities(values, rows=query.n_obs), before


__all__ = [
    "InductiveSCVIError",
    "KNN_NEIGHBORS",
    "ROSTER",
    "UNLABELED_CATEGORY",
    "canonical_sha256",
    "direct_latent_inference",
    "direct_scanvi_probabilities",
    "donor_class_weights",
    "module_state_sha256",
    "validate_integer_like_counts",
    "validate_probabilities",
    "weighted_knn_probabilities",
]
