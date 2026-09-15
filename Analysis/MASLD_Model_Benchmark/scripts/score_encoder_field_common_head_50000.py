#!/usr/bin/env python3
"""Score the encoder field under one shared donor bootstrap and one exposure audit.

Every model in the field is resampled with the *same* donor multiplicities, so
the deltas below are paired: draw b of model A and draw b of model B resample
the identical donors.  Marginal intervals on two models would not license a
comparison; these do.

No score is emitted without the exposure status of the encoder that produced it,
and exposure is resolved from the audit and crosswalk files at write time rather
than transcribed.  Because every outer fold holds out whole studies, exposure is
reported per fold: a fold whose held-out study appears in an encoder's
pretraining corpus is labelled confounded for that encoder and is never averaged
together with its clean folds.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for entry in (str(ROOT), str(ROOT / "src")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from masld_bench.artifacts import freeze_tree, write_json_exclusive  # noqa: E402
from masld_bench.hashing import canonical_sha256, sha256_file  # noqa: E402
from scripts.fit_predict_common_cell_heads_study_50000 import (  # noqa: E402
    HEADS,
    OUTER_FOLDS,
    SCREEN_SEEDS,
)
from scripts.fit_predict_encoder_field_common_head_50000 import (  # noqa: E402
    FAMILY_LIVER,
    FAMILY_PRETRAINED,
    FAMILY_TASK_NATIVE,
    EncoderFieldError,
    load_split,
)
from scripts.score_cell_baselines_study_50000 import (  # noqa: E402
    BOOTSTRAP_RESAMPLES,
    BOOTSTRAP_SEED,
    DATASET_VIEW_ID,
    ROSTER,
    SPLIT_ID,
    build_multiplicities,
    endpoints_from_multiplicities,
    interval,
    score_subset,
    sufficient_statistics,
)

# Exposure states that leave a held-out study usable as a clean read on an
# encoder.  "unknown" is deliberately not clean: an unresolved corpus is not
# evidence of absence.
CLEAN_EXPOSURE = frozenset({"clean_declared", "target_label_unexposed"})
CONFOUNDED_EXPOSURE = frozenset({"encoder_seen", "continual_seen", "reference_only"})


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise EncoderFieldError(f"TSV lacks a header: {path}")
        return list(reader.fieldnames), list(reader)


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def probability_block(rows: Sequence[Mapping[str, str]]) -> np.ndarray:
    fields = [f"probability::{label}" for label in ROSTER]
    values = np.asarray(
        [[float(row[field]) for field in fields] for row in rows], dtype=np.float64
    )
    if (
        values.shape != (len(rows), len(ROSTER))
        or not np.all(np.isfinite(values))
        or np.any(values < 0.0)
        or not np.allclose(values.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-7)
    ):
        raise EncoderFieldError("prediction probabilities differ")
    return values


def assemble_common_head_model(
    *, shard_root: Path, block_id: str, head_id: str, row_ids: np.ndarray, outer: np.ndarray
) -> tuple[np.ndarray, dict[str, Any]]:
    """Build the three-seed ensemble for one (block, head) from its 15 shards."""

    order = {row_id: index for index, row_id in enumerate(row_ids.tolist())}
    per_seed: list[np.ndarray] = []
    receipts: dict[str, Any] = {}
    for screen_seed in SCREEN_SEEDS:
        matrix = np.full((len(row_ids), len(ROSTER)), np.nan, dtype=np.float64)
        coverage = np.zeros(len(row_ids), dtype=np.int8)
        for fold in OUTER_FOLDS:
            name = f"{block_id}__{head_id}__seed{screen_seed}__fold{fold}"
            shard = shard_root / name
            _fields, rows = read_tsv(shard / "predictions.tsv")
            positions = np.asarray([order[row["row_id"]] for row in rows], dtype=np.int64)
            if set(int(row["outer_fold"]) for row in rows) != {fold} or not np.all(
                outer[positions] == fold
            ):
                raise EncoderFieldError(f"{name}: shard rows are not the held-out fold")
            matrix[positions] = probability_block(rows)
            coverage[positions] += 1
            receipts[name] = json.loads(
                (shard / "prediction_receipt.json").read_text(encoding="utf-8")
            )
        if not np.all(coverage == 1) or not np.all(np.isfinite(matrix)):
            raise EncoderFieldError(
                f"{block_id}/{head_id}/seed{screen_seed}: rows lack exact one-time coverage"
            )
        per_seed.append(matrix)
    ensemble = np.mean(np.stack(per_seed, axis=0), axis=0)
    ensemble /= ensemble.sum(axis=1, keepdims=True)
    identity = {
        key: receipts[next(iter(receipts))][key]
        for key in ("family", "width", "exposure_status", "features_refit_per_outer_fold")
        if key in receipts[next(iter(receipts))]
    }
    identity["shard_receipts_sha256"] = canonical_sha256(receipts)
    identity["shards"] = len(receipts)
    return ensemble, identity


def load_prior_model(
    *, path: Path, row_ids: np.ndarray
) -> np.ndarray:
    _fields, rows = read_tsv(path)
    if len(rows) != len(row_ids):
        raise EncoderFieldError(f"{path}: prior prediction table row count differs")
    if [row["row_id"] for row in rows] != row_ids.tolist():
        raise EncoderFieldError(f"{path}: prior prediction table row order differs")
    return probability_block(rows)


def resolve_exposure(entry: Mapping[str, Any], studies: Sequence[str]) -> dict[str, Any]:
    """Read per-study exposure out of the audit and crosswalk, never a literal."""

    audit_path = ROOT / entry["audit_relative"]
    crosswalk_path = ROOT / entry["crosswalk_relative"]
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
    checkpoint = audit["checkpoint_findings"][entry["checkpoint_key"]]
    per_study: dict[str, str] = {}
    for study in studies:
        finding = crosswalk["findings"][study]
        if entry.get("crosswalk_corpus_indirection"):
            corpus = checkpoint["corpus"]
            per_study[study] = finding["corpus_findings"][corpus]["exposure_state"]
        else:
            per_study[study] = finding["exposure_state"]
    return {
        "checkpoint_exposure_state": checkpoint.get("exposure_state"),
        "checkpoint_corpus": checkpoint.get("corpus"),
        "corpus_annotation_state": checkpoint.get("corpus_annotation_state"),
        "sealed_champion_eligibility": audit.get("sealed_champion_eligibility"),
        "per_study_exposure": per_study,
        "audit_sha256": sha256_file(audit_path),
        "crosswalk_sha256": sha256_file(crosswalk_path),
    }


def fold_disposition(per_study: Mapping[str, str], fold_studies: Sequence[str]) -> str:
    states = {per_study[study] for study in fold_studies}
    if states & CONFOUNDED_EXPOSURE:
        return "confounded_by_exposure"
    if states - CLEAN_EXPOSURE:
        return "exposure_unresolved"
    return "clean"


def run(
    *,
    source: Path,
    split: Path,
    registry_path: Path,
    shards: Path,
    output: Path,
) -> None:
    if output.exists():
        raise EncoderFieldError("refusing to overwrite scorer output")
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    context = load_split(source, split)
    row_ids = context["row_ids"]
    donors = context["donors"]
    datasets = context["datasets"]
    truth = context["targets"].astype(np.int8)
    outer = context["outer"]
    unique_donors = sorted(set(donors.tolist()))
    unique_studies = sorted(set(datasets.tolist()))
    donor_to_index = {donor: index for index, donor in enumerate(unique_donors)}
    donor_indices = np.asarray(
        [donor_to_index[donor] for donor in donors], dtype=np.int16
    )
    donor_count = len(unique_donors)
    donor_studies = np.asarray(
        [next(iter(set(datasets[donors == donor].tolist()))) for donor in unique_donors],
        dtype=str,
    )
    fold_studies = {
        fold: sorted(set(datasets[outer == fold].tolist())) for fold in OUTER_FOLDS
    }
    study_to_fold = {
        study: int(outer[datasets == study][0]) for study in unique_studies
    }

    exposure_by_block: dict[str, dict[str, Any]] = {}
    for entry in registry["frozen_matrix_blocks"]:
        exposure_by_block[entry["block_id"]] = resolve_exposure(entry, unique_studies)
    # The liver encoder has no pretraining corpus to audit: its exposure is the
    # split itself, and the split holds out whole studies, so every fold is clean
    # by construction.  This is recorded, not assumed, by the fold receipts that
    # ScviLatentBlock already validated.
    exposure_by_block["scvi_liver_latent"] = {
        "checkpoint_exposure_state": "trained_here_on_outer_training_studies_only",
        "checkpoint_corpus": "resource atlas outer-training studies of each fold",
        "corpus_annotation_state": "scvi_unsupervised_plus_scanvi_labels_on_training_studies",
        "sealed_champion_eligibility": "not_applicable_no_pretraining_checkpoint",
        "per_study_exposure": {study: "clean_declared" for study in unique_studies},
        "audit_sha256": None,
        "crosswalk_sha256": None,
    }
    exposure_by_block["hvg_pca_task_native"] = {
        "checkpoint_exposure_state": "no_pretraining_corpus",
        "checkpoint_corpus": None,
        "corpus_annotation_state": "none",
        "sealed_champion_eligibility": "not_applicable_no_pretraining_checkpoint",
        "per_study_exposure": {study: "clean_declared" for study in unique_studies},
        "audit_sha256": None,
        "crosswalk_sha256": None,
    }

    models: dict[str, np.ndarray] = {}
    model_meta: dict[str, dict[str, Any]] = {}
    block_ids = [entry["block_id"] for entry in registry["frozen_matrix_blocks"]]
    block_ids += ["scvi_liver_latent", "hvg_pca_task_native"]
    for block_id in block_ids:
        for head_id in HEADS:
            model_id = f"{block_id}::{head_id}"
            probabilities, identity = assemble_common_head_model(
                shard_root=shards,
                block_id=block_id,
                head_id=head_id,
                row_ids=row_ids,
                outer=outer,
            )
            models[model_id] = probabilities
            family = identity.get("family") or (
                FAMILY_LIVER
                if block_id == "scvi_liver_latent"
                else FAMILY_TASK_NATIVE
                if block_id == "hvg_pca_task_native"
                else FAMILY_PRETRAINED
            )
            model_meta[model_id] = {
                "block_id": block_id,
                "head_id": head_id,
                "head_is_the_shared_common_head": True,
                "family": family,
                **{k: v for k, v in identity.items() if k != "family"},
                **exposure_by_block[block_id],
            }
    for entry in registry["prior_frozen_prediction_sets"]:
        model_id = entry["model_id"]
        models[model_id] = load_prior_model(
            path=ROOT / entry["predictions_relative"], row_ids=row_ids
        )
        block_id = (
            "scvi_liver_latent"
            if entry["family"] == FAMILY_LIVER
            else "hvg_pca_task_native"
        )
        model_meta[model_id] = {
            "block_id": block_id,
            "head_id": entry["head"],
            "head_is_the_shared_common_head": False,
            "family": entry["family"],
            "predictions_sha256": sha256_file(ROOT / entry["predictions_relative"]),
            **exposure_by_block[block_id],
        }

    global_present = np.zeros((donor_count, len(ROSTER)), dtype=bool)
    for donor in range(donor_count):
        global_present[donor, np.unique(truth[donor_indices == donor])] = True
    multiplicities = build_multiplicities(
        donor_studies, global_present, n_resamples=BOOTSTRAP_RESAMPLES, seed=BOOTSTRAP_SEED
    )
    if multiplicities.shape != (BOOTSTRAP_RESAMPLES, donor_count) or np.any(
        multiplicities.sum(axis=1) != donor_count
    ):
        raise EncoderFieldError("donor bootstrap multiplicities differ")

    output.mkdir(mode=0o750, parents=True)
    all_donor_indices = np.arange(donor_count)
    all_classes = np.ones(len(ROSTER), dtype=bool)
    point_metrics: dict[str, dict[str, float]] = {}
    distributions: dict[str, np.ndarray] = {}
    metrics_rows: list[dict[str, Any]] = []
    study_rows: list[dict[str, Any]] = []
    interval_rows: list[dict[str, Any]] = []

    for model_id, probabilities in models.items():
        overall = score_subset(truth, donors, probabilities)
        study_values: dict[str, float] = {}
        stats = sufficient_statistics(
            truth,
            np.argmax(probabilities, axis=1),
            probabilities,
            donor_indices,
            donors=donor_count,
        )
        for study in unique_studies:
            selected = datasets == study
            value = score_subset(truth[selected], donors[selected], probabilities[selected])
            study_values[study] = float(value["macro_f1"])
            selected_donors = np.flatnonzero(donor_studies == study)
            study_classes = stats["present"][selected_donors].any(axis=0)
            study_draws = endpoints_from_multiplicities(
                stats, multiplicities, selected_donors, study_classes
            )["macro_f1"]
            distributions[f"{model_id}::study_macro_f1::{study}"] = study_draws
            disposition = fold_disposition(
                model_meta[model_id]["per_study_exposure"], [study]
            )
            study_rows.append(
                {
                    "model_id": model_id,
                    "block_id": model_meta[model_id]["block_id"],
                    "head_id": model_meta[model_id]["head_id"],
                    "family": model_meta[model_id]["family"],
                    "held_out_study": study,
                    "outer_fold": study_to_fold[study],
                    "rows": value["rows"],
                    "donors": value["donors"],
                    "study_exposure_state": model_meta[model_id]["per_study_exposure"][
                        study
                    ],
                    "exposure_disposition": disposition,
                    "donor_class_balanced_macro_f1": study_values[study],
                    **{
                        key: val
                        for key, val in interval(study_draws).items()
                        if key in {"lower", "median", "upper"}
                    },
                }
            )
        # Internal consistency: the bootstrap machinery and the point-estimate
        # machinery must agree at multiplicity one.  A wrong donor index map
        # would leave both self-consistent but disagree here.
        identity_draw = endpoints_from_multiplicities(
            stats,
            np.ones((1, donor_count), dtype=np.int16),
            all_donor_indices,
            all_classes,
        )
        if (
            abs(float(identity_draw["macro_f1"][0]) - float(overall["macro_f1"])) > 1.0e-12
            or abs(float(identity_draw["brier"][0]) - float(overall["multiclass_brier"]))
            > 1.0e-12
        ):
            raise EncoderFieldError(
                f"{model_id}: sufficient-statistic endpoint rederivation differs"
            )
        study_balanced = float(np.mean(list(study_values.values())))
        point_metrics[model_id] = {
            "donor_class_balanced_macro_f1": float(overall["macro_f1"]),
            "study_balanced_macro_f1": study_balanced,
            "multiclass_brier": float(overall["multiclass_brier"]),
            "top_label_ece": float(overall["calibration"]["top_label_ece"]),
        }
        donor_bootstrap = endpoints_from_multiplicities(
            stats, multiplicities, all_donor_indices, all_classes
        )
        study_parts = [
            distributions[f"{model_id}::study_macro_f1::{study}"]
            for study in unique_studies
        ]
        metric_draws = {
            "donor_class_balanced_macro_f1": donor_bootstrap["macro_f1"],
            "study_balanced_macro_f1": np.mean(np.stack(study_parts, axis=1), axis=1),
            "multiclass_brier": donor_bootstrap["brier"],
            "top_label_ece": donor_bootstrap["top_label_ece"],
        }
        for metric, values in metric_draws.items():
            distributions[f"{model_id}::{metric}"] = values.astype(np.float64, copy=False)
            interval_rows.append(
                {
                    "model_id": model_id,
                    "block_id": model_meta[model_id]["block_id"],
                    "head_id": model_meta[model_id]["head_id"],
                    "family": model_meta[model_id]["family"],
                    "metric": metric,
                    "estimate": point_metrics[model_id][metric],
                    **interval(values),
                }
            )
        metrics_rows.append(
            {
                "model_id": model_id,
                "block_id": model_meta[model_id]["block_id"],
                "head_id": model_meta[model_id]["head_id"],
                "family": model_meta[model_id]["family"],
                "head_is_the_shared_common_head": model_meta[model_id][
                    "head_is_the_shared_common_head"
                ],
                "checkpoint_exposure_state": model_meta[model_id][
                    "checkpoint_exposure_state"
                ],
                "clean_folds": sum(
                    1
                    for fold in OUTER_FOLDS
                    if fold_disposition(
                        model_meta[model_id]["per_study_exposure"], fold_studies[fold]
                    )
                    == "clean"
                ),
                "confounded_folds": sum(
                    1
                    for fold in OUTER_FOLDS
                    if fold_disposition(
                        model_meta[model_id]["per_study_exposure"], fold_studies[fold]
                    )
                    == "confounded_by_exposure"
                ),
                "unresolved_folds": sum(
                    1
                    for fold in OUTER_FOLDS
                    if fold_disposition(
                        model_meta[model_id]["per_study_exposure"], fold_studies[fold]
                    )
                    == "exposure_unresolved"
                ),
                **point_metrics[model_id],
            }
        )

    # ------------------------------------------------------------------
    # paired deltas under the shared multiplicities
    # ------------------------------------------------------------------
    common_head_models = [
        model_id
        for model_id in models
        if model_meta[model_id]["head_is_the_shared_common_head"]
    ]
    delta_rows: list[dict[str, Any]] = []

    def add_delta(left: str, right: str, comparison: str) -> None:
        for metric, higher_is_better in (
            ("donor_class_balanced_macro_f1", True),
            ("study_balanced_macro_f1", True),
            ("multiclass_brier", False),
        ):
            left_draws = distributions[f"{left}::{metric}"]
            right_draws = distributions[f"{right}::{metric}"]
            values = (
                left_draws - right_draws if higher_is_better else right_draws - left_draws
            )
            observed = (
                point_metrics[left][metric] - point_metrics[right][metric]
                if higher_is_better
                else point_metrics[right][metric] - point_metrics[left][metric]
            )
            bounds = interval(values)
            delta_rows.append(
                {
                    "comparison": comparison,
                    "model_id": left,
                    "reference_id": right,
                    "model_exposure": model_meta[left]["checkpoint_exposure_state"],
                    "reference_exposure": model_meta[right]["checkpoint_exposure_state"],
                    "paired": True,
                    "metric": metric,
                    "oriented_delta": observed,
                    "lower": bounds["lower"],
                    "median": bounds["median"],
                    "upper": bounds["upper"],
                    "interval_crosses_zero": bool(
                        bounds["lower"] <= 0.0 <= bounds["upper"]
                    ),
                    "probability_model_better": float(np.mean(values > 0.0)),
                    "n_resamples": BOOTSTRAP_RESAMPLES,
                    "confidence_level": 0.95,
                    "seed": BOOTSTRAP_SEED,
                }
            )

    for head_id in HEADS:
        liver = f"scvi_liver_latent::{head_id}"
        task_native = f"hvg_pca_task_native::{head_id}"
        for model_id in common_head_models:
            if model_meta[model_id]["head_id"] != head_id:
                continue
            if model_id != task_native:
                add_delta(model_id, task_native, "vs_task_native_same_head")
            if model_id != liver:
                add_delta(model_id, liver, "vs_liver_encoder_same_head")

    # Per-fold paired deltas against the liver encoder, so a confounded fold can
    # never be silently averaged into a clean one.
    fold_delta_rows: list[dict[str, Any]] = []
    for head_id in HEADS:
        liver = f"scvi_liver_latent::{head_id}"
        for model_id in common_head_models:
            if model_meta[model_id]["head_id"] != head_id or model_id == liver:
                continue
            for study in unique_studies:
                left = distributions[f"{model_id}::study_macro_f1::{study}"]
                right = distributions[f"{liver}::study_macro_f1::{study}"]
                values = left - right
                bounds = interval(values)
                observed = float(
                    next(
                        row["donor_class_balanced_macro_f1"]
                        for row in study_rows
                        if row["model_id"] == model_id and row["held_out_study"] == study
                    )
                    - next(
                        row["donor_class_balanced_macro_f1"]
                        for row in study_rows
                        if row["model_id"] == liver and row["held_out_study"] == study
                    )
                )
                fold_delta_rows.append(
                    {
                        "model_id": model_id,
                        "reference_id": liver,
                        "held_out_study": study,
                        "outer_fold": study_to_fold[study],
                        "model_study_exposure": model_meta[model_id][
                            "per_study_exposure"
                        ][study],
                        "exposure_disposition": fold_disposition(
                            model_meta[model_id]["per_study_exposure"], [study]
                        ),
                        "metric": "donor_class_balanced_macro_f1",
                        "oriented_delta": observed,
                        "lower": bounds["lower"],
                        "median": bounds["median"],
                        "upper": bounds["upper"],
                        "interval_crosses_zero": bool(
                            bounds["lower"] <= 0.0 <= bounds["upper"]
                        ),
                        "probability_model_better": float(np.mean(values > 0.0)),
                    }
                )

    # ------------------------------------------------------------------
    # the comparison that is actually clean
    #
    # Each encoder has its own set of held-out studies that its pretraining
    # corpus demonstrably did not contain.  Restricting both sides of a pair to
    # that encoder's clean studies is the only contrast in this panel that is
    # not confounded by exposure.  An encoder with no clean study gets a row
    # saying so and no number, because averaging its confounded folds would
    # manufacture a clean-looking score.
    # ------------------------------------------------------------------
    clean_rows: list[dict[str, Any]] = []
    for head_id in HEADS:
        liver = f"scvi_liver_latent::{head_id}"
        for model_id in common_head_models:
            if model_meta[model_id]["head_id"] != head_id or model_id == liver:
                continue
            per_study = model_meta[model_id]["per_study_exposure"]
            clean_studies = [
                study
                for study in unique_studies
                if fold_disposition(per_study, [study]) == "clean"
            ]
            row: dict[str, Any] = {
                "model_id": model_id,
                "reference_id": liver,
                "head_id": head_id,
                "clean_studies": ",".join(clean_studies),
                "clean_study_count": len(clean_studies),
                "confounded_studies": ",".join(
                    study
                    for study in unique_studies
                    if fold_disposition(per_study, [study]) == "confounded_by_exposure"
                ),
                "unresolved_studies": ",".join(
                    study
                    for study in unique_studies
                    if fold_disposition(per_study, [study]) == "exposure_unresolved"
                ),
            }
            if not clean_studies:
                row.update(
                    {
                        "model_clean_macro_f1": "",
                        "reference_clean_macro_f1": "",
                        "oriented_delta": "",
                        "lower": "",
                        "median": "",
                        "upper": "",
                        "interval_crosses_zero": "",
                        "probability_model_better": "",
                        "verdict": "no_clean_held_out_study_exists_for_this_encoder",
                    }
                )
                clean_rows.append(row)
                continue
            model_draws = np.mean(
                np.stack(
                    [
                        distributions[f"{model_id}::study_macro_f1::{study}"]
                        for study in clean_studies
                    ],
                    axis=1,
                ),
                axis=1,
            )
            liver_draws = np.mean(
                np.stack(
                    [
                        distributions[f"{liver}::study_macro_f1::{study}"]
                        for study in clean_studies
                    ],
                    axis=1,
                ),
                axis=1,
            )
            by_study = {
                (r["model_id"], r["held_out_study"]): r[
                    "donor_class_balanced_macro_f1"
                ]
                for r in study_rows
            }
            model_point = float(
                np.mean([by_study[(model_id, study)] for study in clean_studies])
            )
            liver_point = float(
                np.mean([by_study[(liver, study)] for study in clean_studies])
            )
            values = model_draws - liver_draws
            bounds = interval(values)
            crosses = bool(bounds["lower"] <= 0.0 <= bounds["upper"])
            row.update(
                {
                    "model_clean_macro_f1": model_point,
                    "reference_clean_macro_f1": liver_point,
                    "oriented_delta": model_point - liver_point,
                    "lower": bounds["lower"],
                    "median": bounds["median"],
                    "upper": bounds["upper"],
                    "interval_crosses_zero": crosses,
                    "probability_model_better": float(np.mean(values > 0.0)),
                    "verdict": "interval_crosses_zero"
                    if crosses
                    else (
                        "model_above_liver_encoder"
                        if model_point > liver_point
                        else "model_below_liver_encoder"
                    ),
                }
            )
            clean_rows.append(row)
    write_tsv(
        output / "clean_fold_paired_deltas_vs_liver_encoder.tsv",
        list(clean_rows[0]),
        clean_rows,
    )

    write_tsv(output / "model_metrics.tsv", list(metrics_rows[0]), metrics_rows)
    write_tsv(output / "study_metrics.tsv", list(study_rows[0]), study_rows)
    write_tsv(output / "bootstrap_intervals.tsv", list(interval_rows[0]), interval_rows)
    write_tsv(output / "paired_deltas.tsv", list(delta_rows[0]), delta_rows)
    write_tsv(
        output / "paired_deltas_by_held_out_study.tsv",
        list(fold_delta_rows[0]),
        fold_delta_rows,
    )
    with (output / "bootstrap_distributions.npz").open("xb") as handle:
        np.savez_compressed(handle, **distributions)
    write_json_exclusive(
        output / "model_exposure.json",
        {model_id: model_meta[model_id] for model_id in sorted(model_meta)},
    )
    receipt = {
        "schema_version": "masld-bench-encoder-field-evaluator-v1",
        "status": "pass_development_evaluation_only",
        "evaluator_id": "encoder_field_common_head_50000_evaluation_only_v1",
        "dataset_view_id": DATASET_VIEW_ID,
        "split_id": SPLIT_ID,
        "rows": int(len(row_ids)),
        "donors": donor_count,
        "studies": len(unique_studies),
        "classes": list(ROSTER),
        "models": sorted(models),
        "shared_common_head_models": sorted(common_head_models),
        "heads": list(HEADS),
        "screen_seeds": list(SCREEN_SEEDS),
        "fold_held_out_studies": {str(k): v for k, v in fold_studies.items()},
        "models_fit": False,
        "preprocessing_fit": False,
        "calibration_fit": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
        "clinical_claim_allowed": False,
        "exposure_resolution": "read_from_config/artifacts/models/*/exposure_audit.json_and_development_crosswalk.json_at_write_time",
        "bootstrap": {
            "method": "study_stratified_class_coverage_preserving_donor_cluster_bootstrap",
            "resamples": BOOTSTRAP_RESAMPLES,
            "seed": BOOTSTRAP_SEED,
            "donors_are_independent_units": True,
            "cells_are_independent_units": False,
            "shared_multiplicities_across_all_models": True,
            "deltas_are_paired_on_the_same_draws": True,
            "multiplicities_sha256": canonical_sha256(multiplicities.tolist()),
            "distribution_sha256": canonical_sha256(
                {key: values.tolist() for key, values in sorted(distributions.items())}
            ),
        },
        "metric_definitions": {
            "donor_class_balanced_macro_f1": "Every class receives equal mass, every donor containing that class receives equal within-class mass, and cells divide their donor-class share.",
            "study_balanced_macro_f1": "Arithmetic mean of the seven study-specific donor-class-balanced macro-F1 values.",
            "multiclass_brier": "Donor-class-balanced sum of squared error over the frozen five-class probability vector.",
            "oriented_delta": "Positive always favours model_id over reference_id, for both higher-is-better and lower-is-better metrics.",
        },
        "exposure_policy": {
            "clean": sorted(CLEAN_EXPOSURE),
            "confounded": sorted(CONFOUNDED_EXPOSURE),
            "unresolved": "any other state, including unknown; an unresolved corpus is not evidence of no exposure",
        },
        "registry_sha256": sha256_file(registry_path),
    }
    write_json_exclusive(output / "evaluation_receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "encoder_field_common_head_evaluation",
            "dataset_view_id": DATASET_VIEW_ID,
            "split_id": SPLIT_ID,
            "rows": int(len(row_ids)),
            "donors": donor_count,
            "studies": len(unique_studies),
            "models_fit": False,
            "histology_read": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )


def main() -> int:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--source", required=True, type=Path)
    value.add_argument("--split", required=True, type=Path)
    value.add_argument("--registry", required=True, type=Path)
    value.add_argument("--shards", required=True, type=Path)
    value.add_argument("--output", required=True, type=Path)
    arguments = value.parse_args()
    run(
        source=arguments.source,
        split=arguments.split,
        registry_path=arguments.registry,
        shards=arguments.shards,
        output=arguments.output,
    )
    print(json.dumps({"output": arguments.output.resolve().as_posix(), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
