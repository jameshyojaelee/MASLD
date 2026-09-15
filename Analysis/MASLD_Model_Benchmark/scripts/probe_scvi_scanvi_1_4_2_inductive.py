#!/usr/bin/env python3
"""Empirically admit direct, zero-fit held-study inference in scvi-tools 1.4.2."""

from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
from scipy import sparse
import scvi
import torch

from masld_bench.adapters.scvi_scanvi_inductive import (
    ROSTER,
    UNLABELED_CATEGORY,
    direct_latent_inference,
    direct_scanvi_probabilities,
    donor_class_weights,
    module_state_sha256,
    weighted_knn_probabilities,
)


def _synthetic_objects() -> tuple[anndata.AnnData, anndata.AnnData, np.ndarray]:
    generator = np.random.default_rng(20260824)
    genes = np.asarray([f"ENSG{index:011d}" for index in range(40)], dtype=str)
    labels = np.repeat(np.asarray(ROSTER, dtype=str), 24)
    rates = np.ones((len(ROSTER), len(genes)), dtype=np.float64)
    for class_index in range(len(ROSTER)):
        rates[class_index, class_index * 6 : class_index * 6 + 6] = 8.0
    encoded = np.asarray([ROSTER.index(label) for label in labels], dtype=np.int64)
    train_x = generator.poisson(rates[encoded]).astype(np.float32)
    query_encoded = np.repeat(np.arange(len(ROSTER), dtype=np.int64), 4)
    query_x = generator.poisson(rates[query_encoded]).astype(np.float32)

    train_obs = pd.DataFrame(
        {
            "row_id": [f"reference-{index}" for index in range(len(labels))],
            "donor_id": [f"reference-donor-{index // 10}" for index in range(len(labels))],
            "study": "reference-study",
            "broad_label": pd.Categorical(
                labels, categories=[*ROSTER, UNLABELED_CATEGORY]
            ),
        }
    )
    query_obs = pd.DataFrame(
        {
            "row_id": [f"query-{index}" for index in range(len(query_encoded))],
            "donor_id": [f"query-donor-{index // 4}" for index in range(len(query_encoded))],
            "study": "entirely-unseen-held-study",
            "broad_label": pd.Categorical(
                [UNLABELED_CATEGORY] * len(query_encoded),
                categories=[*ROSTER, UNLABELED_CATEGORY],
            ),
        }
    )
    train = anndata.AnnData(
        X=sparse.csr_matrix(train_x),
        obs=train_obs,
        var=pd.DataFrame(index=genes),
    )
    query = anndata.AnnData(
        X=sparse.csr_matrix(query_x),
        obs=query_obs,
        var=pd.DataFrame(index=genes),
    )
    return train, query, labels


def run(output: Path) -> None:
    if output.exists():
        raise ValueError("refusing to overwrite inductive runtime probe")
    if scvi.__version__ != "1.4.2":
        raise ValueError(f"scvi-tools version differs: {scvi.__version__}")
    output.mkdir(mode=0o750)
    scvi.settings.seed = 20260824
    train, query, labels = _synthetic_objects()

    scvi.model.SCVI.setup_anndata(train)
    scvi_model = scvi.model.SCVI(
        train,
        n_hidden=32,
        n_latent=8,
        n_layers=1,
        dropout_rate=0.0,
        dispersion="gene",
        gene_likelihood="nb",
        latent_distribution="normal",
    )
    scvi_model.train(
        max_epochs=2,
        train_size=1.0,
        batch_size=32,
        early_stopping=False,
        check_val_every_n_epoch=None,
        enable_progress_bar=False,
    )
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        labels_key="broad_label",
        unlabeled_category=UNLABELED_CATEGORY,
    )
    scanvi_model.train(
        max_epochs=2,
        train_size=1.0,
        batch_size=32,
        early_stopping=False,
        check_val_every_n_epoch=None,
        enable_progress_bar=False,
    )

    reference_latent = np.asarray(
        scvi_model.get_latent_representation(), dtype=np.float32
    )
    scvi_before = module_state_sha256(scvi_model)
    query_latent, scvi_inference_hash = direct_latent_inference(scvi_model, query)
    weights = donor_class_weights(
        list(map(str, train.obs["donor_id"])), list(map(str, labels))
    )
    knn = weighted_knn_probabilities(
        reference_latent,
        labels,
        weights,
        query_latent,
        neighbors=15,
        n_jobs=1,
    )
    scanvi_before = module_state_sha256(scanvi_model)
    scanvi_probabilities, scanvi_inference_hash = direct_scanvi_probabilities(
        scanvi_model, query
    )

    if (
        scvi_before != scvi_inference_hash
        or scanvi_before != scanvi_inference_hash
        or set(map(str, query.obs["study"])) != {"entirely-unseen-held-study"}
        or knn.shape != (query.n_obs, len(ROSTER))
        or scanvi_probabilities.shape != (query.n_obs, len(ROSTER))
    ):
        raise ValueError("direct inductive inference contract differs")

    payload = {
        "schema_version": "masld-bench-scvi-scanvi-inductive-probe-v1",
        "status": "pass",
        "python": platform.python_version(),
        "scvi_tools": scvi.__version__,
        "torch": torch.__version__,
        "torch_cuda_build": torch.version.cuda,
        "cuda_available": torch.cuda.is_available(),
        "train_rows": train.n_obs,
        "query_rows": query.n_obs,
        "genes": train.n_vars,
        "query_study_unseen": True,
        "query_labels_exposed_to_model": False,
        "batch_key": None,
        "categorical_covariates": [],
        "continuous_covariates": [],
        "query_training_executed": False,
        "query_adaptation_executed": False,
        "direct_scvi_inference": True,
        "direct_scanvi_inference": True,
        "scvi_state_unchanged": module_state_sha256(scvi_model) == scvi_before,
        "scanvi_state_unchanged": module_state_sha256(scanvi_model) == scanvi_before,
        "scvi_state_sha256": scvi_before,
        "scanvi_state_sha256": scanvi_before,
        "knn_neighbors": 15,
        "knn_probabilities_normalized": bool(
            np.allclose(knn.sum(axis=1), 1.0)
        ),
        "scanvi_probabilities_normalized": bool(
            np.allclose(scanvi_probabilities.sum(axis=1), 1.0)
        ),
        "project_data_read": False,
        "sealed_outcomes_read": False,
    }
    if not all(
        payload[key]
        for key in (
            "query_study_unseen",
            "direct_scvi_inference",
            "direct_scanvi_inference",
            "scvi_state_unchanged",
            "scanvi_state_unchanged",
            "knn_probabilities_normalized",
            "scanvi_probabilities_normalized",
        )
    ):
        raise ValueError("inductive probe did not pass")
    with (output / "probe.json").open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, sort_keys=True, separators=(",", ":"))
        handle.write("\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    run(arguments.output)
    print(json.dumps({"output": str(arguments.output.resolve()), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
