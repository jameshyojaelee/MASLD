#!/usr/bin/env python3
"""Fit one identical common head on every encoder of the 50,000-cell atlas field.

Panel 7F asks a single question: how does a liver-specific encoder trained here
compare with pretrained foundation encoders when the task, the rows, the split
and the head are held identical.  This producer answers the head half of it.

Every feature block below is scored through the same head implementation that
already produced the frozen Geneformer and TranscriptFormer shards; the head
functions are imported from ``fit_predict_common_cell_heads_study_50000`` rather
than restated, so "identical head" is a property of the import and not of a
comment.  Blocks differ only in how their features are produced:

* ``scvi_liver_latent`` is refit on the outer-training studies of each fold, so
  its features are fold-specific by construction (inductive liver encoder).
* the five pretrained blocks are frozen 50,000-row matrices reused unchanged
  across folds.
* ``hvg_pca_task_native`` is refit per fold from raw counts, giving a task-native
  representation that never saw any pretraining corpus.

Shards are written the moment they finish.  A timeout therefore costs the
unfinished shards only, never the finished ones.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from masld_bench.artifacts import freeze_tree, write_json_exclusive  # noqa: E402
from masld_bench.hashing import canonical_sha256, sha256_file  # noqa: E402
from scripts.fit_predict_common_cell_heads_study_50000 import (  # noqa: E402
    BATCH_SIZE,
    DATASET_VIEW_ID,
    FIXED_EPOCHS,
    HEADS,
    LEARNING_RATE,
    OUTER_FOLDS,
    PREDICTION_FIELDS,
    ROSTER,
    SCREEN_SEEDS,
    SPLIT_ID,
    WEIGHT_DECAY,
    _fit_fixed_head,
    _fit_temperature,
    _read_tsv,
    _validate_cuda_runtime,
    _write_tsv,
    donor_class_weights,
    validate_inner_assignment,
    validate_study_outer_split,
)

# The liver encoder was fit under its own seed bundle.  Screen seed j is paired
# with encoder seed j so the three-seed ensemble spans encoder and head noise
# together rather than head noise alone.
SCVI_ENCODER_SEEDS = (20260824, 20260825, 20260826)
N_TOP_HVG = 2_000
N_HVG_BINS = 20
N_PCA = 50
TARGET_SUM = 10_000.0

FAMILY_LIVER = "liver_specific_trained_here"
FAMILY_PRETRAINED = "pretrained_foundation_encoder"
FAMILY_TASK_NATIVE = "task_native_no_pretraining"


class EncoderFieldError(ValueError):
    """Raised when a field input, firewall, or fit does not meet the contract."""


# --------------------------------------------------------------------------
# feature blocks
# --------------------------------------------------------------------------


class FeatureBlock:
    """A named source of 50,000 x d features for one (screen seed, outer fold)."""

    block_id: str
    family: str
    width: int

    def identity(self) -> dict[str, Any]:
        raise NotImplementedError

    def features(self, *, screen_seed: int, held_outer: int) -> np.ndarray:
        raise NotImplementedError


class FrozenMatrixBlock(FeatureBlock):
    """A pretrained encoder's frozen 50,000-row matrix, reused across folds."""

    def __init__(
        self,
        *,
        block_id: str,
        bundle_path: Path,
        receipt_path: Path,
        row_ids: np.ndarray,
        exposure_status: str,
    ) -> None:
        self.block_id = block_id
        self.family = FAMILY_PRETRAINED
        self._bundle_path = Path(bundle_path)
        self._receipt_path = Path(receipt_path)
        self._exposure_status = exposure_status
        with np.load(self._bundle_path, allow_pickle=False) as bundle:
            if set(bundle.files) != {"embeddings", "outer_folds", "row_ids"}:
                raise EncoderFieldError(
                    f"{block_id}: embedding bundle fields differ"
                )
            matrix = np.asarray(bundle["embeddings"], dtype=np.float32)
            bundle_row_ids = np.asarray(bundle["row_ids"]).astype(str)
        if (
            matrix.ndim != 2
            or matrix.shape[0] != len(row_ids)
            or not np.all(np.isfinite(matrix))
        ):
            raise EncoderFieldError(f"{block_id}: embedding matrix shape or values differ")
        # Row-order guard.  A silently permuted matrix would still fit a head and
        # still produce plausible metrics, so the identity is asserted, not assumed.
        if bundle_row_ids.tolist() != row_ids.tolist():
            raise EncoderFieldError(f"{block_id}: embedding row order differs from fixture")
        self._matrix = matrix
        self.width = int(matrix.shape[1])
        self._receipt = json.loads(self._receipt_path.read_text(encoding="utf-8"))
        if self._receipt.get("rows") != len(row_ids):
            raise EncoderFieldError(f"{block_id}: embedding receipt row count differs")
        if self._receipt.get("sealed_outcomes_read") is not False:
            raise EncoderFieldError(f"{block_id}: embedding receipt firewall differs")

    def identity(self) -> dict[str, Any]:
        return {
            "block_id": self.block_id,
            "family": self.family,
            "width": self.width,
            "features_refit_per_outer_fold": False,
            "exposure_status": self._exposure_status,
            "embedding_bundle_sha256": sha256_file(self._bundle_path),
            "embedding_receipt_sha256": sha256_file(self._receipt_path),
            "embedding_bundle_path": self._bundle_path.as_posix(),
            "checkpoint_sha256": self._receipt.get("checkpoint_sha256"),
            "representation_policy": self._receipt.get("policies"),
        }

    def features(self, *, screen_seed: int, held_outer: int) -> np.ndarray:
        return self._matrix


class ScviLatentBlock(FeatureBlock):
    """The liver encoder's latent, refit on the outer-training studies per fold."""

    def __init__(
        self,
        *,
        bundle_root: Path,
        row_ids: np.ndarray,
        outer: np.ndarray,
        datasets: np.ndarray,
    ) -> None:
        self.block_id = "scvi_liver_latent"
        self.family = FAMILY_LIVER
        self._root = Path(bundle_root)
        self._row_count = len(row_ids)
        self._outer = outer
        self._datasets = datasets
        self._seed_map = dict(zip(SCREEN_SEEDS, SCVI_ENCODER_SEEDS, strict=True))
        widths: set[int] = set()
        self._fold_receipts: dict[tuple[int, int], dict[str, Any]] = {}
        for screen_seed, encoder_seed in self._seed_map.items():
            for fold in OUTER_FOLDS:
                fold_root = self._root / f"seed-{encoder_seed}" / "folds" / f"fold{fold}"
                receipt = json.loads(
                    (fold_root / "fold_receipt.json").read_text(encoding="utf-8")
                )
                held = sorted(set(datasets[outer == fold].tolist()))
                trained = sorted(set(datasets[outer != fold].tolist()))
                # Exposure firewall.  The liver encoder must have been fit on the
                # outer-training studies only; if a held-out study appears in its
                # training roster the whole panel claim collapses.
                if (
                    sorted(map(str, receipt.get("query_studies", []))) != held
                    or sorted(map(str, receipt.get("training_studies", []))) != trained
                    or set(receipt.get("query_studies", [])) & set(receipt.get("training_studies", []))
                    or int(receipt.get("outer_fold", -1)) != fold
                    or int(receipt.get("seed", -1)) != encoder_seed
                    or receipt.get("held_query_training") is not False
                    or receipt.get("held_query_adaptation") is not False
                    or receipt.get("held_query_parameter_change") is not False
                    or receipt.get("hvg_fit_on_outer_training_only") is not True
                ):
                    raise EncoderFieldError(
                        f"scvi_liver_latent: fold {fold} seed {encoder_seed} exposure firewall differs"
                    )
                reference = np.load(fold_root / "reference_latent.npy", mmap_mode="r")
                query = np.load(fold_root / "query_latent.npy", mmap_mode="r")
                if (
                    reference.shape[0] != int((outer != fold).sum())
                    or query.shape[0] != int((outer == fold).sum())
                    or reference.shape[1] != query.shape[1]
                ):
                    raise EncoderFieldError(
                        f"scvi_liver_latent: fold {fold} latent row counts differ from the split"
                    )
                widths.add(int(reference.shape[1]))
                self._fold_receipts[(screen_seed, fold)] = receipt
        if len(widths) != 1:
            raise EncoderFieldError("scvi_liver_latent: latent width differs across folds")
        self.width = widths.pop()

    def identity(self) -> dict[str, Any]:
        return {
            "block_id": self.block_id,
            "family": self.family,
            "width": self.width,
            "features_refit_per_outer_fold": True,
            "exposure_status": "trained_here_on_outer_training_studies_only",
            "bundle_root": self._root.as_posix(),
            "screen_seed_to_encoder_seed": {
                str(key): value for key, value in self._seed_map.items()
            },
            "fold_receipt_sha256": canonical_sha256(
                {
                    f"{seed}|{fold}": receipt
                    for (seed, fold), receipt in sorted(self._fold_receipts.items())
                }
            ),
        }

    def features(self, *, screen_seed: int, held_outer: int) -> np.ndarray:
        encoder_seed = self._seed_map[screen_seed]
        fold_root = self._root / f"seed-{encoder_seed}" / "folds" / f"fold{held_outer}"
        reference = np.asarray(
            np.load(fold_root / "reference_latent.npy"), dtype=np.float32
        )
        query = np.asarray(np.load(fold_root / "query_latent.npy"), dtype=np.float32)
        matrix = np.empty((self._row_count, self.width), dtype=np.float32)
        # The producer wrote the reference block in fixture order over the
        # outer-training rows and the query block in fixture order over the
        # held-out rows.  reconstruct_scvi_knn_probabilities() proves this
        # placement by reproducing the frozen scvi_baseline predictions.
        matrix[self._outer != held_outer] = reference
        matrix[self._outer == held_outer] = query
        if not np.all(np.isfinite(matrix)):
            raise EncoderFieldError("scvi_liver_latent: latent contains nonfinite values")
        return matrix


class HvgPcaBlock(FeatureBlock):
    """HVG + PCA refit per outer fold from raw counts: the task-native block."""

    def __init__(self, *, fixture: Path, row_ids: np.ndarray, outer: np.ndarray) -> None:
        self.block_id = "hvg_pca_task_native"
        self.family = FAMILY_TASK_NATIVE
        self.width = N_PCA
        self._fixture = Path(fixture)
        self._row_ids = row_ids
        self._outer = outer
        self._adata: Any = None
        self._gene_ids: list[str] | None = None
        self._cache: dict[tuple[int, int], np.ndarray] = {}

    def _load(self) -> None:
        if self._adata is not None:
            return
        import anndata

        adata = anndata.read_h5ad(
            self._fixture / "resource_atlas_frozen_screen_50000.h5ad"
        )
        if (
            adata.n_obs != len(self._row_ids)
            or list(map(str, adata.obs["row_id"])) != self._row_ids.tolist()
        ):
            raise EncoderFieldError("hvg_pca_task_native: fixture row order differs")
        self._adata = adata
        self._gene_ids = list(map(str, adata.var["ensembl_id"]))

    def identity(self) -> dict[str, Any]:
        return {
            "block_id": self.block_id,
            "family": self.family,
            "width": self.width,
            "features_refit_per_outer_fold": True,
            "exposure_status": "no_pretraining_corpus",
            "recipe": {
                "normalization": f"library_size_{int(TARGET_SUM)}_log1p",
                "n_top_hvg": N_TOP_HVG,
                "n_hvg_bins": N_HVG_BINS,
                "n_pca": N_PCA,
                "fit_scope": "outer_training_rows_only",
            },
        }

    def features(self, *, screen_seed: int, held_outer: int) -> np.ndarray:
        key = (screen_seed, held_outer)
        if key in self._cache:
            return self._cache[key]
        from sklearn.decomposition import PCA

        from masld_bench.adapters.hvg_pca_logistic import _log_normalize, _select_hvgs

        self._load()
        assert self._adata is not None and self._gene_ids is not None
        train = np.flatnonzero(self._outer != held_outer)
        normalized_train = _log_normalize(self._adata.X[train], TARGET_SUM)
        selected, _ = _select_hvgs(
            normalized_train, self._gene_ids, n_top=N_TOP_HVG, n_bins=N_HVG_BINS
        )
        train_dense = (
            normalized_train[:, selected].toarray().astype(np.float32, copy=False)
        )
        pca = PCA(
            n_components=N_PCA,
            svd_solver="randomized",
            whiten=False,
            random_state=screen_seed + held_outer,
        )
        pca.fit(train_dense)
        del train_dense
        normalized_all = _log_normalize(self._adata.X, TARGET_SUM)
        all_dense = normalized_all[:, selected].toarray().astype(np.float32, copy=False)
        matrix = pca.transform(all_dense).astype(np.float32, copy=False)
        del all_dense, normalized_all, normalized_train
        if matrix.shape != (len(self._row_ids), N_PCA) or not np.all(np.isfinite(matrix)):
            raise EncoderFieldError("hvg_pca_task_native: projection shape or values differ")
        self._cache = {key: matrix}
        return matrix


# --------------------------------------------------------------------------
# shared fixture / split load
# --------------------------------------------------------------------------


def load_split(source: Path, split: Path) -> dict[str, Any]:
    selection = _read_tsv(source / "selection.tsv")
    split_rows = _read_tsv(split / "row_outer_folds.tsv")
    inner_rows = _read_tsv(split / "inner_donor_folds.tsv")
    row_ids = np.asarray([row["row_id"] for row in selection], dtype=str)
    donors = np.asarray([row["donor_id"] for row in selection], dtype=str)
    datasets = np.asarray([row["dataset"] for row in selection], dtype=str)
    labels = np.asarray([row["broad_label"] for row in selection], dtype=str)
    if (
        [row["row_id"] for row in split_rows] != row_ids.tolist()
        or [row["donor_id"] for row in split_rows] != donors.tolist()
        or [row["dataset"] for row in split_rows] != datasets.tolist()
        or set(labels.tolist()) != set(ROSTER)
        or len(set(row_ids.tolist())) != len(row_ids)
    ):
        raise EncoderFieldError("source and split rows differ")
    outer = np.asarray([int(row["outer_fold"]) for row in split_rows], dtype=np.int8)
    validate_study_outer_split(donors, datasets, outer)
    targets = np.asarray([ROSTER.index(label) for label in labels], dtype=np.int64)
    inner_by_outer: dict[int, dict[str, int]] = {fold: {} for fold in OUTER_FOLDS}
    for row in inner_rows:
        held = int(row["held_outer_fold"])
        if held not in OUTER_FOLDS or row["donor_id"] in inner_by_outer[held]:
            raise EncoderFieldError("inner donor or held fold differs")
        inner_by_outer[held][row["donor_id"]] = int(row["inner_fold"])
    for held_outer in OUTER_FOLDS:
        validate_inner_assignment(donors[outer != held_outer], inner_by_outer[held_outer])
    return {
        "row_ids": row_ids,
        "donors": donors,
        "datasets": datasets,
        "targets": targets,
        "outer": outer,
        "inner_by_outer": inner_by_outer,
    }


# --------------------------------------------------------------------------
# the identical head, applied to one block / head / seed / fold
# --------------------------------------------------------------------------


def fit_shard(
    *,
    block: FeatureBlock,
    head_id: str,
    screen_seed: int,
    held_outer: int,
    context: Mapping[str, Any],
) -> dict[str, Any]:
    """Reproduce the frozen common-head shard protocol on one feature block.

    Seeds, epoch count, optimizer, inner-fold temperature cross-fitting and the
    final outer-training refit follow run_shard() in the common-head producer
    exactly; only the feature source is parameterised.
    """

    import torch

    features = block.features(screen_seed=screen_seed, held_outer=held_outer)
    targets = context["targets"]
    donors = context["donors"]
    datasets = context["datasets"]
    outer = context["outer"]
    training = np.flatnonzero(outer != held_outer)
    test = np.flatnonzero(outer == held_outer)
    if set(datasets[training].tolist()) & set(datasets[test].tolist()):
        raise EncoderFieldError("shard train/test studies overlap")
    if features.shape != (len(outer), block.width):
        raise EncoderFieldError("block feature matrix shape differs")
    row_inner_folds = validate_inner_assignment(
        donors[training], context["inner_by_outer"][held_outer]
    )
    oof_logits = np.full((len(training), len(ROSTER)), np.nan, dtype=np.float64)
    oof_coverage = np.zeros(len(training), dtype=np.int8)
    inner_records: list[dict[str, Any]] = []
    for inner_fold in OUTER_FOLDS:
        validation_local = row_inner_folds == inner_fold
        seed = screen_seed + held_outer * 101 + inner_fold * 17 + HEADS.index(head_id)
        model, _mean, _scale, logits, _state, record = _fit_fixed_head(
            head_id=head_id,
            features=features,
            targets=targets,
            donors=donors,
            fit_indices=training[~validation_local],
            prediction_indices=training[validation_local],
            seed=seed,
        )
        oof_logits[validation_local] = logits
        oof_coverage[validation_local] += 1
        inner_records.append({"inner_fold": inner_fold, "seed": seed, **record})
        del model
        torch.cuda.empty_cache()
    if not np.all(oof_coverage == 1) or not np.all(np.isfinite(oof_logits)):
        raise EncoderFieldError("OOF logits lack exact one-time coverage")
    temperature = _fit_temperature(
        oof_logits,
        targets[training],
        donor_class_weights(donors[training], targets[training]),
    )
    final_seed = screen_seed + held_outer * 101 + 1009 + HEADS.index(head_id)
    final_model, _fmean, _fscale, final_logits, _fstate, final_record = _fit_fixed_head(
        head_id=head_id,
        features=features,
        targets=targets,
        donors=donors,
        fit_indices=training,
        prediction_indices=test,
        seed=final_seed,
    )
    final_model.eval()
    with torch.no_grad():
        probabilities = (
            torch.softmax(torch.from_numpy(final_logits).to("cuda") / temperature, dim=1)
            .double()
            .cpu()
            .numpy()
        )
    probabilities /= probabilities.sum(axis=1, keepdims=True)
    if (
        not np.all(np.isfinite(probabilities))
        or np.any(probabilities < 0.0)
        or not np.allclose(probabilities.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-12)
    ):
        raise EncoderFieldError("common-head probabilities differ")
    del final_model
    torch.cuda.empty_cache()
    return {
        "test_indices": test,
        "probabilities": probabilities,
        "temperature": temperature,
        "inner_models": inner_records,
        "final_record": final_record,
        "final_seed": final_seed,
        "training_rows": int(len(training)),
        "training_donors": len(set(donors[training].tolist())),
        "training_studies": sorted(set(datasets[training].tolist())),
        "test_rows": int(len(test)),
        "test_donors": len(set(donors[test].tolist())),
        "test_studies": sorted(set(datasets[test].tolist())),
    }


def write_shard(
    *,
    output: Path,
    block: FeatureBlock,
    head_id: str,
    screen_seed: int,
    held_outer: int,
    result: Mapping[str, Any],
    context: Mapping[str, Any],
) -> None:
    output.mkdir(mode=0o750, parents=True)
    probabilities = result["probabilities"]
    test = result["test_indices"]
    donors = context["donors"]
    datasets = context["datasets"]
    outer = context["outer"]
    row_ids = context["row_ids"]
    predicted = [ROSTER[int(np.argmax(row))] for row in probabilities]
    _write_tsv(
        output / "predictions.tsv",
        PREDICTION_FIELDS,
        (
            {
                "row_id": row_ids[index],
                "donor_id": donors[index],
                "dataset": datasets[index],
                "outer_fold": int(outer[index]),
                "predicted_class": predicted[position],
                **{
                    f"probability::{label}": format(
                        float(probabilities[position, class_index]), ".17g"
                    )
                    for class_index, label in enumerate(ROSTER)
                },
            }
            for position, index in enumerate(test)
        ),
    )
    receipt = {
        "schema_version": "masld-bench-encoder-field-common-head-shard-v1",
        "status": "pass_development_prediction_shard",
        **block.identity(),
        "dataset_view_id": DATASET_VIEW_ID,
        "split_id": SPLIT_ID,
        "head_id": head_id,
        "head_source_module": "scripts.fit_predict_common_cell_heads_study_50000",
        "head_hyperparameters": {
            "fixed_epochs": FIXED_EPOCHS,
            "batch_size": BATCH_SIZE,
            "learning_rate": LEARNING_RATE,
            "weight_decay": WEIGHT_DECAY,
        },
        "screen_seed": screen_seed,
        "outer_fold": held_outer,
        "training_rows": result["training_rows"],
        "training_donors": result["training_donors"],
        "training_studies": result["training_studies"],
        "test_rows": result["test_rows"],
        "test_donors": result["test_donors"],
        "test_studies": result["test_studies"],
        "calibration": "fixed_epoch_donor_cross_fitted_temperature",
        "temperature": result["temperature"],
        "inner_models": result["inner_models"],
        "inner_models_used_for_outer_prediction": False,
        "checkpoint_selection_uses_validation_labels": False,
        "outer_prediction_model": "fixed_epoch_final_outer_training_refit",
        "final_fit_seed": result["final_seed"],
        "metrics_calculated": False,
        "development_labels_read": ["broad_label"],
        "prediction_tables_contain_observed_labels": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
    }
    write_json_exclusive(output / "prediction_receipt.json", receipt)


def run(
    *,
    source: Path,
    split: Path,
    output: Path,
    blocks: Sequence[FeatureBlock],
    heads: Sequence[str],
    context: Mapping[str, Any],
) -> None:
    shard_root = output / "shards"
    shard_root.mkdir(mode=0o750, parents=True, exist_ok=True)
    # Head is innermost on purpose: the refit-per-fold blocks recompute their
    # features for each (seed, fold), and both heads then share one computation.
    for block in blocks:
        for screen_seed in SCREEN_SEEDS:
            for held_outer in OUTER_FOLDS:
                for head_id in heads:
                    name = (
                        f"{block.block_id}__{head_id}"
                        f"__seed{screen_seed}__fold{held_outer}"
                    )
                    target = shard_root / name
                    if (target / "prediction_receipt.json").is_file():
                        continue
                    result = fit_shard(
                        block=block,
                        head_id=head_id,
                        screen_seed=screen_seed,
                        held_outer=held_outer,
                        context=context,
                    )
                    write_shard(
                        output=target,
                        block=block,
                        head_id=head_id,
                        screen_seed=screen_seed,
                        held_outer=held_outer,
                        result=result,
                        context=context,
                    )
                    # Written the moment it exists: a wall-clock timeout costs
                    # only the shard in flight.
                    print(
                        json.dumps(
                            {
                                "shard": name,
                                "test_rows": result["test_rows"],
                                "temperature": result["temperature"],
                            },
                            sort_keys=True,
                        ),
                        flush=True,
                    )


def build_blocks(
    *,
    registry: Mapping[str, Any],
    context: Mapping[str, Any],
    fixture: Path,
    only: Sequence[str] | None = None,
) -> list[FeatureBlock]:
    blocks: list[FeatureBlock] = []
    for entry in registry["frozen_matrix_blocks"]:
        if only and entry["block_id"] not in only:
            continue
        blocks.append(
            FrozenMatrixBlock(
                block_id=entry["block_id"],
                bundle_path=ROOT / entry["bundle_relative"],
                receipt_path=ROOT / entry["receipt_relative"],
                row_ids=context["row_ids"],
                exposure_status=entry["exposure_status"],
            )
        )
    if not only or "scvi_liver_latent" in only:
        blocks.append(
            ScviLatentBlock(
                bundle_root=ROOT / registry["scvi_bundle_relative"],
                row_ids=context["row_ids"],
                outer=context["outer"],
                datasets=context["datasets"],
            )
        )
    if not only or "hvg_pca_task_native" in only:
        blocks.append(
            HvgPcaBlock(
                fixture=fixture, row_ids=context["row_ids"], outer=context["outer"]
            )
        )
    return blocks


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--source", required=True, type=Path)
    value.add_argument("--split", required=True, type=Path)
    value.add_argument("--registry", required=True, type=Path)
    value.add_argument("--output", required=True, type=Path)
    value.add_argument("--only-block", action="append", default=None)
    value.add_argument("--head", action="append", default=None, choices=list(HEADS))
    value.add_argument("--freeze", action="store_true")
    return value


def main() -> int:
    arguments = parser().parse_args()
    for variable in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        if not os.environ.get(variable):
            raise EncoderFieldError(f"{variable} must be pinned to --cpus-per-task")
    # The frozen Geneformer and TranscriptFormer shards this experiment must
    # reproduce were fit on one L40S under a deterministic CUBLAS workspace.
    # A different device silently produces a different head, so refuse early.
    _validate_cuda_runtime()
    registry = json.loads(arguments.registry.read_text(encoding="utf-8"))
    context = load_split(arguments.source, arguments.split)
    blocks = build_blocks(
        registry=registry,
        context=context,
        fixture=arguments.source,
        only=arguments.only_block,
    )
    run(
        source=arguments.source,
        split=arguments.split,
        output=arguments.output,
        blocks=blocks,
        heads=tuple(arguments.head) if arguments.head else HEADS,
        context=context,
    )
    if arguments.freeze:
        freeze_tree(
            arguments.output,
            {
                "artifact_class": "encoder_field_common_head_predictions",
                "dataset_view_id": DATASET_VIEW_ID,
                "split_id": SPLIT_ID,
                "blocks": [block.block_id for block in blocks],
                "heads": list(arguments.head) if arguments.head else list(HEADS),
                "metrics_calculated": False,
                "sealed_outcomes_read": False,
                "status": "passed",
            },
        )
    print(json.dumps({"output": arguments.output.resolve().as_posix(), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
