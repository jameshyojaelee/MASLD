#!/usr/bin/env python3
"""Fit one seed of donor-safe, study-held scVI and scANVI baselines.

scVI is evaluated with a frozen weighted 15-nearest-neighbor rule in its
training-reference latent space.  scANVI is evaluated with its native
classifier.  The held study receives only direct frozen inference: there is no
query-model construction, query training, adaptation, or held-label stopping.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from masld_bench.adapters.hvg_pca_logistic import _log_normalize, _select_hvgs
from masld_bench.adapters.scvi_scanvi_inductive import (
    KNN_NEIGHBORS,
    ROSTER,
    UNLABELED_CATEGORY,
    canonical_sha256,
    direct_latent_inference,
    direct_scanvi_probabilities,
    donor_class_weights,
    module_state_sha256,
    validate_integer_like_counts,
    validate_probabilities,
    weighted_knn_probabilities,
)


DATASET_VIEW_ID = "resource_atlas_frozen_screen_50000_v1"
SPLIT_ID = "resource_atlas_study_outer_5fold_v1"
MODEL_IDS = ("scvi_baseline", "scanvi_baseline")
SEEDS = (20260824, 20260825, 20260826)
N_TOP_HVG = 4_000
HVG_MEAN_BINS = 20
HVG_TARGET_SUM = 10_000.0
N_HIDDEN = 128
N_LATENT = 30
N_LAYERS = 2
DROPOUT_RATE = 0.1
DISPERSION = "gene"
GENE_LIKELIHOOD = "nb"
LATENT_DISTRIBUTION = "normal"
SCVI_EPOCHS = 200
SCANVI_EPOCHS = 100
BATCH_SIZE = 512


class SCVIFrozenScreenError(ValueError):
    """Raised when a frozen-screen input, fit, or prediction differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json_exclusive(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(
            value,
            handle,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
        handle.write("\n")


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SCVIFrozenScreenError(f"TSV lacks a header: {path}")
        return list(reader)


def _write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]
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


def _save_array(path: Path, value: Any) -> None:
    with path.open("xb") as handle:
        np.save(handle, np.asarray(value), allow_pickle=False)


def _training_and_query_objects(
    *,
    raw_train: Any,
    raw_query: Any,
    gene_ids: Sequence[str],
    train_row_ids: Sequence[str],
    train_donors: Sequence[str],
    train_studies: Sequence[str],
    train_labels: Sequence[str],
    query_row_ids: Sequence[str],
    query_donors: Sequence[str],
    query_studies: Sequence[str],
) -> tuple[Any, Any]:
    import anndata
    import pandas as pd
    from scipy import sparse

    if set(train_labels) != set(ROSTER):
        raise SCVIFrozenScreenError("outer training labels lack the frozen roster")
    if set(train_studies) & set(query_studies):
        raise SCVIFrozenScreenError("one study crosses the outer firewall")
    categories = [*ROSTER, UNLABELED_CATEGORY]
    variables = pd.DataFrame(
        {"ensembl_id": list(map(str, gene_ids))},
        index=pd.Index(list(map(str, gene_ids)), name="ensembl_id"),
    )
    training = anndata.AnnData(
        X=sparse.csr_matrix(raw_train, dtype=np.float32),
        obs=pd.DataFrame(
            {
                "row_id": list(map(str, train_row_ids)),
                "donor_id": list(map(str, train_donors)),
                "study": list(map(str, train_studies)),
                "broad_label": pd.Categorical(train_labels, categories=categories),
            },
            index=pd.Index(list(map(str, train_row_ids)), name="row_id"),
        ),
        var=variables.copy(),
    )
    query = anndata.AnnData(
        X=sparse.csr_matrix(raw_query, dtype=np.float32),
        obs=pd.DataFrame(
            {
                "row_id": list(map(str, query_row_ids)),
                "donor_id": list(map(str, query_donors)),
                "study": list(map(str, query_studies)),
                "broad_label": pd.Categorical(
                    [UNLABELED_CATEGORY] * len(query_row_ids), categories=categories
                ),
            },
            index=pd.Index(list(map(str, query_row_ids)), name="row_id"),
        ),
        var=variables.copy(),
    )
    validate_integer_like_counts(training.X)
    validate_integer_like_counts(query.X)
    if set(map(str, query.obs["broad_label"])) != {UNLABELED_CATEGORY}:
        raise SCVIFrozenScreenError("held-study truth entered the model-facing object")
    return training, query


def _train_one_fold(
    *,
    raw_train: Any,
    raw_query: Any,
    gene_ids: Sequence[str],
    train_row_ids: np.ndarray,
    train_donors: np.ndarray,
    train_studies: np.ndarray,
    train_labels: np.ndarray,
    query_row_ids: np.ndarray,
    query_donors: np.ndarray,
    query_studies: np.ndarray,
    seed: int,
    fold: int,
    output: Path,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    import scvi
    import torch

    if scvi.__version__ != "1.4.2":
        raise SCVIFrozenScreenError(
            f"registered scvi-tools 1.4.2 differs: {scvi.__version__}"
        )
    require_cuda = os.environ.get("MASLD_REQUIRE_CUDA", "0") == "1"
    if require_cuda and not torch.cuda.is_available():
        raise SCVIFrozenScreenError("production wrapper requires one visible CUDA GPU")
    accelerator = "gpu" if torch.cuda.is_available() else "cpu"
    training, query = _training_and_query_objects(
        raw_train=raw_train,
        raw_query=raw_query,
        gene_ids=gene_ids,
        train_row_ids=train_row_ids,
        train_donors=train_donors,
        train_studies=train_studies,
        train_labels=train_labels,
        query_row_ids=query_row_ids,
        query_donors=query_donors,
        query_studies=query_studies,
    )
    model_seed = seed + fold
    scvi.settings.seed = model_seed
    scvi.model.SCVI.setup_anndata(training)
    scvi_model = scvi.model.SCVI(
        training,
        n_hidden=N_HIDDEN,
        n_latent=N_LATENT,
        n_layers=N_LAYERS,
        dropout_rate=DROPOUT_RATE,
        dispersion=DISPERSION,
        gene_likelihood=GENE_LIKELIHOOD,
        latent_distribution=LATENT_DISTRIBUTION,
    )
    scvi_model.train(
        max_epochs=SCVI_EPOCHS,
        train_size=1.0,
        batch_size=BATCH_SIZE,
        early_stopping=False,
        check_val_every_n_epoch=None,
        enable_progress_bar=False,
        accelerator=accelerator,
        devices=1,
    )
    scvi_training_hash = module_state_sha256(scvi_model)
    reference_latent = np.asarray(
        scvi_model.get_latent_representation(), dtype=np.float32
    )
    if (
        reference_latent.shape != (training.n_obs, N_LATENT)
        or not np.all(np.isfinite(reference_latent))
    ):
        raise SCVIFrozenScreenError("training-reference latent representation differs")

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        labels_key="broad_label",
        unlabeled_category=UNLABELED_CATEGORY,
    )
    scanvi_model.train(
        max_epochs=SCANVI_EPOCHS,
        train_size=1.0,
        batch_size=BATCH_SIZE,
        early_stopping=False,
        check_val_every_n_epoch=None,
        enable_progress_bar=False,
        accelerator=accelerator,
        devices=1,
    )
    scanvi_training_hash = module_state_sha256(scanvi_model)

    scvi_model.save(output / "scvi_model", overwrite=False, save_anndata=False)
    scanvi_model.save(output / "scanvi_model", overwrite=False, save_anndata=False)
    query_latent, scvi_inference_hash = direct_latent_inference(scvi_model, query)
    reference_weights = donor_class_weights(train_donors, train_labels)
    scvi_probabilities = weighted_knn_probabilities(
        reference_latent,
        train_labels,
        reference_weights,
        query_latent,
        neighbors=KNN_NEIGHBORS,
        n_jobs=max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))),
    )
    scanvi_probabilities, scanvi_inference_hash = direct_scanvi_probabilities(
        scanvi_model, query
    )
    if (
        scvi_training_hash != scvi_inference_hash
        or scanvi_training_hash != scanvi_inference_hash
        or module_state_sha256(scvi_model) != scvi_training_hash
        or module_state_sha256(scanvi_model) != scanvi_training_hash
    ):
        raise SCVIFrozenScreenError("held-study inference changed learned state")

    _save_array(output / "reference_latent.npy", reference_latent)
    _save_array(output / "query_latent.npy", query_latent)
    receipt = {
        "seed": seed,
        "model_seed": model_seed,
        "outer_fold": fold,
        "training_rows": training.n_obs,
        "training_donors": len(set(map(str, train_donors))),
        "training_studies": sorted(set(map(str, train_studies))),
        "query_rows": query.n_obs,
        "query_donors": len(set(map(str, query_donors))),
        "query_studies": sorted(set(map(str, query_studies))),
        "selected_feature_ids_sha256": canonical_sha256(list(map(str, gene_ids))),
        "scvi_state_sha256": scvi_training_hash,
        "scanvi_state_sha256": scanvi_training_hash,
        "scvi_epochs_completed": SCVI_EPOCHS,
        "scanvi_epochs_completed": SCANVI_EPOCHS,
        "fixed_epochs": True,
        "early_stopping": False,
        "train_size": 1.0,
        "held_labels_in_model_object": False,
        "held_query_training": False,
        "held_query_adaptation": False,
        "held_query_parameter_change": False,
        "batch_key": None,
        "categorical_covariates": [],
        "continuous_covariates": [],
        "observed_library_size_only": True,
        "biological_unit": "donor",
        "accelerator": accelerator,
    }
    return {
        "scvi_baseline": validate_probabilities(
            scvi_probabilities, rows=query.n_obs
        ),
        "scanvi_baseline": validate_probabilities(
            scanvi_probabilities, rows=query.n_obs
        ),
    }, receipt


def run(
    *,
    source: Path,
    split: Path,
    output: Path,
    expected_source_artifacts_sha256: str,
    expected_split_artifacts_sha256: str,
    seed: int,
) -> None:
    import anndata

    if output.exists():
        raise SCVIFrozenScreenError("refusing to overwrite scVI/scANVI output")
    if seed not in SEEDS:
        raise SCVIFrozenScreenError("seed is outside the frozen roster")
    if sha256_file(source / "ARTIFACTS.json") != expected_source_artifacts_sha256:
        raise SCVIFrozenScreenError("source ARTIFACTS SHA-256 differs")
    if sha256_file(split / "ARTIFACTS.json") != expected_split_artifacts_sha256:
        raise SCVIFrozenScreenError("split ARTIFACTS SHA-256 differs")
    source_manifest = json.loads(
        (source / "ARTIFACTS.json").read_text(encoding="utf-8")
    )
    split_manifest = json.loads(
        (split / "ARTIFACTS.json").read_text(encoding="utf-8")
    )
    if (
        source_manifest["metadata"].get("subset_id") != DATASET_VIEW_ID
        or source_manifest["metadata"].get("row_count") != 50_000
        or split_manifest["metadata"].get("split_id") != SPLIT_ID
        or split_manifest["metadata"].get("target_labels_used_for_assignment") is not False
        or split_manifest["metadata"].get("sealed_outcomes_read") is not False
    ):
        raise SCVIFrozenScreenError("source or split metadata differs")

    selection = _read_tsv(source / "selection.tsv")
    split_rows = _read_tsv(split / "row_outer_folds.tsv")
    inner_rows = _read_tsv(split / "inner_donor_folds.tsv")
    adata = anndata.read_h5ad(source / "resource_atlas_frozen_screen_50000.h5ad")
    row_ids = np.asarray(list(map(str, adata.obs["row_id"])), dtype=str)
    donors = np.asarray(list(map(str, adata.obs["donor_id"])), dtype=str)
    studies = np.asarray(list(map(str, adata.obs["dataset"])), dtype=str)
    labels = np.asarray(list(map(str, adata.obs["broad_label"])), dtype=str)
    expected_row_ids = [row["row_id"] for row in selection]
    if (
        adata.shape != (50_000, 37_533)
        or row_ids.tolist() != expected_row_ids
        or [row["row_id"] for row in split_rows] != expected_row_ids
        or donors.tolist() != [row["donor_id"] for row in split_rows]
        or studies.tolist() != [row["dataset"] for row in split_rows]
        or set(labels) != set(ROSTER)
    ):
        raise SCVIFrozenScreenError("Atlas, selection, and split row contracts differ")
    outer = np.asarray([int(row["outer_fold"]) for row in split_rows], dtype=np.int8)
    if (
        set(outer) != set(range(5))
        or any(len(set(outer[donors == donor])) != 1 for donor in set(donors))
        or any(len(set(outer[studies == study])) != 1 for study in set(studies))
    ):
        raise SCVIFrozenScreenError("one donor or study crosses outer folds")
    validate_integer_like_counts(adata.X)
    gene_ids = list(map(str, adata.var["ensembl_id"]))
    if len(gene_ids) != len(set(gene_ids)):
        raise SCVIFrozenScreenError("Atlas Ensembl identifiers are not unique")

    inner_by_outer: dict[int, dict[str, int]] = defaultdict(dict)
    for row in inner_rows:
        held = int(row["held_outer_fold"])
        donor = row["donor_id"]
        inner_fold = int(row["inner_fold"])
        if held not in range(5) or inner_fold not in range(5):
            raise SCVIFrozenScreenError("inner or held fold is outside the roster")
        if donor in inner_by_outer[held]:
            raise SCVIFrozenScreenError("inner donor assignment repeats")
        inner_by_outer[held][donor] = inner_fold

    output.mkdir(mode=0o750)
    folds_root = output / "folds"
    folds_root.mkdir()
    probabilities = {
        model_id: np.full((len(row_ids), len(ROSTER)), np.nan, dtype=np.float64)
        for model_id in MODEL_IDS
    }
    fold_receipts: list[dict[str, Any]] = []
    for fold in range(5):
        fold_root = folds_root / f"fold{fold}"
        fold_root.mkdir()
        train = np.flatnonzero(outer != fold)
        query = np.flatnonzero(outer == fold)
        if (
            set(inner_by_outer[fold]) != set(donors[train])
            or set(inner_by_outer[fold].values()) != set(range(5))
            or set(studies[train]) & set(studies[query])
        ):
            raise SCVIFrozenScreenError(
                "frozen inner assignment or study-held firewall differs"
            )
        normalized_training = _log_normalize(adata.X[train], HVG_TARGET_SUM)
        selected, feature_records = _select_hvgs(
            normalized_training,
            gene_ids,
            n_top=N_TOP_HVG,
            n_bins=HVG_MEAN_BINS,
        )
        selected_gene_ids = [gene_ids[index] for index in selected]
        _write_tsv(
            fold_root / "selected_features.tsv",
            (
                "feature_index",
                "ensembl_id",
                "training_mean",
                "training_variance",
                "normalized_dispersion",
            ),
            feature_records,
        )
        fold_predictions, fold_receipt = _train_one_fold(
            raw_train=adata.X[train][:, selected],
            raw_query=adata.X[query][:, selected],
            gene_ids=selected_gene_ids,
            train_row_ids=row_ids[train],
            train_donors=donors[train],
            train_studies=studies[train],
            train_labels=labels[train],
            query_row_ids=row_ids[query],
            query_donors=donors[query],
            query_studies=studies[query],
            seed=seed,
            fold=fold,
            output=fold_root,
        )
        fold_receipt.update(
            {
                "hvg_fit_rows": len(train),
                "hvg_fit_donors": len(set(donors[train])),
                "hvg_fit_studies": sorted(set(studies[train])),
                "hvg_fit_on_outer_training_only": True,
                "inner_donor_assignment_sha256": canonical_sha256(
                    sorted(inner_by_outer[fold].items())
                ),
                "inner_folds_used_for_model_selection": False,
                "inner_folds_used_for_early_stopping": False,
            }
        )
        write_json_exclusive(fold_root / "fold_receipt.json", fold_receipt)
        fold_receipts.append(fold_receipt)
        for model_id, values in fold_predictions.items():
            probabilities[model_id][query] = values

    predictions_root = output / "predictions"
    predictions_root.mkdir()
    fields = (
        "row_id",
        "donor_id",
        "dataset",
        "outer_fold",
        "predicted_class",
        *(f"probability::{label}" for label in ROSTER),
    )
    for model_id, values in probabilities.items():
        values = validate_probabilities(values, rows=len(row_ids))
        _write_tsv(
            predictions_root / f"{model_id}--seed-{seed}.tsv",
            fields,
            (
                {
                    "row_id": row_ids[index],
                    "donor_id": donors[index],
                    "dataset": studies[index],
                    "outer_fold": int(outer[index]),
                    "predicted_class": ROSTER[int(np.argmax(values[index]))],
                    **{
                        f"probability::{label}": format(
                            float(values[index, class_index]), ".17g"
                        )
                        for class_index, label in enumerate(ROSTER)
                    },
                }
                for index in range(len(row_ids))
            ),
        )

    receipt = {
        "schema_version": "masld-bench-scvi-scanvi-study-screen-v1",
        "status": "pass_development_predictions",
        "dataset_view_id": DATASET_VIEW_ID,
        "split_id": SPLIT_ID,
        "rows": len(row_ids),
        "donors": len(set(donors)),
        "studies": len(set(studies)),
        "outer_folds": 5,
        "inner_folds": 5,
        "seed": seed,
        "models": list(MODEL_IDS),
        "biological_unit": "donor",
        "outer_split_unit": "study",
        "inner_split_unit": "donor",
        "parameters": {
            "input": "raw_nonnegative_integer_umi_counts",
            "n_top_hvg": N_TOP_HVG,
            "hvg_selection": "training_only_log1p_library_size_dispersion",
            "hvg_mean_bins": HVG_MEAN_BINS,
            "n_hidden": N_HIDDEN,
            "n_latent": N_LATENT,
            "n_layers": N_LAYERS,
            "dropout_rate": DROPOUT_RATE,
            "dispersion": DISPERSION,
            "gene_likelihood": GENE_LIKELIHOOD,
            "latent_distribution": LATENT_DISTRIBUTION,
            "scvi_epochs": SCVI_EPOCHS,
            "scanvi_epochs": SCANVI_EPOCHS,
            "batch_size": BATCH_SIZE,
            "early_stopping": False,
            "train_size": 1.0,
            "scvi_label_transfer": "donor_class_balanced_weighted_15nn",
            "scanvi_label_transfer": "native_classifier",
            "batch_key": None,
            "categorical_covariates": [],
            "continuous_covariates": [],
        },
        "folds": fold_receipts,
        "metrics_calculated": False,
        "common_head_fit": False,
        "held_query_training": False,
        "held_query_adaptation": False,
        "development_labels_read": ["broad_label_training_rows_only_for_fit"],
        "prediction_tables_contain_observed_labels": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
    }
    write_json_exclusive(output / "prediction_receipt.json", receipt)


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--source", required=True, type=Path)
    value.add_argument("--split", required=True, type=Path)
    value.add_argument("--output", required=True, type=Path)
    value.add_argument("--expected-source-artifacts-sha256", required=True)
    value.add_argument("--expected-split-artifacts-sha256", required=True)
    value.add_argument("--seed", required=True, type=int, choices=SEEDS)
    return value


def main() -> int:
    arguments = parser().parse_args()
    run(
        source=arguments.source,
        split=arguments.split,
        output=arguments.output,
        expected_source_artifacts_sha256=arguments.expected_source_artifacts_sha256,
        expected_split_artifacts_sha256=arguments.expected_split_artifacts_sha256,
        seed=arguments.seed,
    )
    print(json.dumps({"output": str(arguments.output.resolve()), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
