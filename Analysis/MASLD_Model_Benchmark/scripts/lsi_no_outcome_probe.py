#!/usr/bin/env python3
"""Synthetic training-only TF-IDF/LSI runtime and leakage probe."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


class LSIProbeError(ValueError):
    """Raised when the observed-ATAC LSI requirement differs."""


def _sha(value) -> str:
    import numpy as np

    array = np.ascontiguousarray(value)
    digest = sha256()
    digest.update(str(array.dtype).encode("ascii"))
    digest.update(json.dumps(list(array.shape)).encode("ascii"))
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def method1_tfidf(counts, training_idf=None):
    """Signac method-1 orientation: log1p(10,000 * TF * training IDF)."""
    import numpy as np

    values = np.asarray(counts, dtype=np.float64)
    if values.ndim != 2 or values.shape[0] == 0 or values.shape[1] == 0:
        raise LSIProbeError("ATAC matrix axes differ")
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise LSIProbeError("ATAC matrix values differ")
    cell_depth = values.sum(axis=1)
    if np.any(cell_depth <= 0):
        raise LSIProbeError("observed all-zero ATAC rows must fail admission")
    if training_idf is None:
        peak_sum = values.sum(axis=0)
        if np.any(peak_sum <= 0):
            raise LSIProbeError("training-only peak inventory contains an empty peak")
        idf = values.shape[0] / peak_sum
    else:
        idf = np.asarray(training_idf, dtype=np.float64)
        if idf.shape != (values.shape[1],) or np.any(idf <= 0):
            raise LSIProbeError("training IDF axis differs")
    tf = values / cell_depth[:, None]
    return np.log1p(10_000.0 * tf * idf[None, :]), idf


def run(output: Path) -> dict:
    import numpy as np
    import sklearn
    from sklearn.decomposition import TruncatedSVD

    if output.exists():
        raise LSIProbeError("output already exists")
    if sklearn.__version__ != "1.9.0":
        raise LSIProbeError("scikit-learn version differs")
    rng = np.random.default_rng(20260824)
    training = rng.poisson(0.7, size=(24, 40)).astype(np.float64)
    training[:, 0] += 1
    training[0, :] += 1
    held = rng.poisson(0.7, size=(7, 40)).astype(np.float64)
    held[:, 0] += 1
    train_tfidf, idf = method1_tfidf(training)
    held_tfidf, reused = method1_tfidf(held, idf)
    if not np.array_equal(idf, reused):
        raise LSIProbeError("held transform changed training IDF")
    model = TruncatedSVD(
        n_components=8,
        algorithm="randomized",
        n_iter=5,
        n_oversamples=10,
        power_iteration_normalizer="auto",
        random_state=20260824,
        tol=0.0,
    )
    train_embedding = model.fit_transform(train_tfidf)
    held_embedding = model.transform(held_tfidf)
    held_individual = np.vstack(
        [model.transform(held_tfidf[index : index + 1]) for index in range(held.shape[0])]
    )
    held_repeat = model.transform(held_tfidf)
    composition_max_abs_difference = float(
        np.max(np.abs(held_embedding - held_individual))
    )
    composition_tolerance = 1.0e-12
    if (
        train_embedding.shape != (24, 8)
        or held_embedding.shape != (7, 8)
        or not np.all(np.isfinite(train_embedding))
        or not np.all(np.isfinite(held_embedding))
        or composition_max_abs_difference > composition_tolerance
        or not np.array_equal(held_embedding, held_repeat)
    ):
        raise LSIProbeError("LSI held transform or determinism differs")
    try:
        method1_tfidf(np.zeros((1, 40)), idf)
    except LSIProbeError:
        observed_zero_rejected = True
    else:
        observed_zero_rejected = False
    if not observed_zero_rejected:
        raise LSIProbeError("observed zero-row rejection fixture differs")
    output.mkdir(parents=True, mode=0o750)
    receipt = {
        "schema_version": "masld-bench-lsi-no-outcome-probe-v1",
        "status": "pass",
        "scikit_learn": sklearn.__version__,
        "algorithm": "Signac_method_1_plus_sklearn_TruncatedSVD",
        "idf_formula": "number_of_outer_training_cells_divided_by_outer_training_peak_sum",
        "scale_factor": 10_000,
        "svd_parameters": {
            "n_components": 8,
            "algorithm": "randomized",
            "n_iter": 5,
            "n_oversamples": 10,
            "power_iteration_normalizer": "auto",
            "random_state": 20260824,
            "tol": 0.0,
        },
        "training_idf_sha256": _sha(idf),
        "svd_components_sha256": _sha(model.components_),
        "held_embedding_sha256": _sha(held_embedding),
        "held_transform_reused_training_idf": True,
        "held_transform_did_not_refit_svd": True,
        "query_composition_invariant": True,
        "query_composition_max_abs_difference": composition_max_abs_difference,
        "query_composition_tolerance": composition_tolerance,
        "repeat_transform_bit_identical": True,
        "observed_all_zero_atac_rejected": True,
        "structurally_missing_atac_mapped_to_zero": False,
        "observed_atac_required_at_query": True,
        "rna_conditioned_atac_eligible": False,
        "sealed_rna_only_eligible": False,
        "synthetic_inputs": True,
        "project_data_read": False,
        "outcomes_read": False,
        "champion_claim_allowed": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    run(arguments.output)


if __name__ == "__main__":
    main()
