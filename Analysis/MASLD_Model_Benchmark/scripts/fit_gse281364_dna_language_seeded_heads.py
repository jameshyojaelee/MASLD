#!/usr/bin/env python3
"""Fit five block-bootstrap DNA-language heads under 524-kb-safe folds."""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation
import gzip
import io
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from sklearn.linear_model import Ridge

from masld_bench.artifacts import verify_frozen_tree
from scripts.evaluate_gse281364_dna_lm_common_lane import load_outcomes
from scripts.fit_gse281364_enformer_sei_heads import (
    _bootstrap_indices,
    _training_transform,
)
from scripts.fit_gse281364_open_sequence_heads import fit_inner_projection
from scripts.gse281364_dna_lm_native_contract import apply_projection, head_features
from scripts.reconcile_gse281364_dna_language_seeded_features import (
    CANDIDATES,
    CONTEXTS,
    HEADS,
    INNER_BOOTSTRAP_SALT_OFFSET,
    MISSING_BASELINES,
    MODELS,
    OUTER_BOOTSTRAP_SALT_OFFSET,
    PREDICTION_FIELDS,
    SEEDS,
    TRAINING_MEAN_BOOTSTRAP_SALT_OFFSET,
    file_sha256,
    load_config,
)


SCHEMA = "masld-bench-gse281364-dna-language-seeded-head-fit-v1"
ROW_FIELDS = (
    "seed",
    "row_hash",
    "unit_hash",
    "block_hash",
    "stratum",
    "outer_fold",
    "study_id",
    "assay_context_id",
    "element_id",
    "source_locus_group_id",
    "long_range_block_id",
)


class DNASeededHeadFitError(RuntimeError):
    """Raised when a feature or held-block fit boundary differs."""


def _load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise DNASeededHeadFitError(f"invalid JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise DNASeededHeadFitError(f"JSON object required: {path}")
    return value


def _tree(root: Path, binding: Mapping[str, Any], *, label: str) -> Path:
    value = binding.get("tree_path")
    if not isinstance(value, str) or not value:
        raise DNASeededHeadFitError(f"{label} tree path missing")
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise DNASeededHeadFitError(f"unsafe {label} tree path")
    path = (root / relative).resolve(strict=True)
    path.relative_to(root)
    if file_sha256(path / "ARTIFACTS.json") != binding.get("artifacts_sha256"):
        raise DNASeededHeadFitError(f"{label} ARTIFACTS identity differs")
    verify_frozen_tree(path)
    return path


def _read_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise DNASeededHeadFitError(f"TSV schema differs: {path}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise DNASeededHeadFitError(f"empty TSV: {path}")
    return rows


def _write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    if not fields:
        raise DNASeededHeadFitError(f"cannot write empty TSV: {path}")
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_gzip_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=PREDICTION_FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)


def _load_alleles(path: Path, elements: Sequence[str]) -> np.ndarray:
    wanted = set(elements)
    pairs: dict[str, tuple[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["element_id"] in wanted:
                pair = (row["ref"], row["alt"])
                if row["element_id"] in pairs or pair[0] not in "ACGT" or pair[1] not in "ACGT" or pair[0] == pair[1]:
                    raise DNASeededHeadFitError("allele identity differs")
                pairs[row["element_id"]] = pair
    if set(pairs) != wanted:
        raise DNASeededHeadFitError("allele coverage differs")
    bases = {base: index for index, base in enumerate("ACGT")}
    features = np.zeros((len(elements), 16), dtype=np.float64)
    for index, element in enumerate(elements):
        ref, alt = pairs[element]
        features[index, 4 * bases[ref] + bases[alt]] = 1.0
    return features


def _load_frozen_shared_control(
    root: Path,
    config: Mapping[str, Any],
    rows: Sequence[Mapping[str, str]],
) -> tuple[dict[tuple[int, str], dict[str, str]], list[dict[str, Any]]]:
    binding = config["shared_control_reference"]
    tree = _tree(root, binding, label="shared_control_reference")
    member = tree / binding["member"]
    if file_sha256(member) != binding["member_sha256"]:
        raise DNASeededHeadFitError("shared-control prediction identity differs")
    authority = {(int(row["seed"]), row["row_hash"]): row for row in rows}
    frozen_rows: dict[tuple[int, str], dict[str, str]] = {}
    with gzip.open(member, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != PREDICTION_FIELDS:
            raise DNASeededHeadFitError("shared-control prediction schema differs")
        for row in reader:
            if row["model_id"] != binding["model_id"] or row["head_id"] != binding["head_id"]:
                continue
            seed = int(row["seed"])
            identity = (seed, row["row_hash"])
            expected = authority.get(identity)
            if seed not in SEEDS or identity in frozen_rows or expected is None:
                raise DNASeededHeadFitError("shared-control row identity differs")
            if any(row[field] != expected[field] for field in ROW_FIELDS):
                raise DNASeededHeadFitError("shared-control row metadata differs")
            if (
                row["experimental_replicates"] != "4"
                or row["biological_donors"] != "0"
                or row["outcome_role"] != "exposed_development_MPRA_only"
            ):
                raise DNASeededHeadFitError("shared-control topology differs")
            try:
                value = Decimal(row["prediction"])
            except InvalidOperation as error:
                raise DNASeededHeadFitError("shared-control prediction is not numeric") from error
            if not value.is_finite():
                raise DNASeededHeadFitError("shared-control prediction is non-finite")
            frozen_rows[identity] = dict(row)
    if len(frozen_rows) != int(binding["expected_rows"]):
        raise DNASeededHeadFitError("shared-control prediction coverage differs")
    selections = []
    with (tree / "head_selection.tsv").open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["model_id"] != binding["model_id"] or row["head_id"] != binding["head_id"]:
                continue
            if any(
                row[field] != "False"
                for field in (
                    "held_outcomes_used_for_fit_or_tuning",
                    "held_features_used_for_preprocessing_selection_or_tuning",
                    "held_blocks_overlap_training",
                    "posthoc_calibration",
                )
            ):
                raise DNASeededHeadFitError("shared-control selection firewall differs")
            selections.append({
                **row,
                "shared_control_prediction_policy": "reuse_bit_exact_frozen_enformer_sei_oof",
                "shared_control_source_artifacts_sha256": binding["artifacts_sha256"],
            })
    if len(selections) != 50:
        raise DNASeededHeadFitError("shared-control selection denominator differs")
    return frozen_rows, selections


def _raw_tree_and_member(config: Mapping[str, Any], model: str) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    record = next(row for row in config["models"] if row["model_id"] == model)
    source = config[record["raw_source"]]
    return source, source["models"][model]


def _load_raw_embeddings(
    root: Path,
    config: Mapping[str, Any],
    model: str,
    projection_root: Path,
    elements: Sequence[str],
) -> np.ndarray:
    source, member = _raw_tree_and_member(config, model)
    raw_tree = _tree(root, source, label=f"{model}_raw")
    raw_path = raw_tree / f"raw/{model}/allele_embeddings.npz"
    if file_sha256(raw_path) != member["matrix_sha256"]:
        raise DNASeededHeadFitError(f"{model} raw matrix identity differs")
    with np.load(raw_path, allow_pickle=False) as data:
        fixture_ids = data["fixture_ids"].astype(str)
        embeddings = data["embeddings"]
        allele_order = tuple(data["allele_order"].astype(str).tolist())
    with np.load(projection_root / f"{model}/heldout_fold0/head_features.npz", allow_pickle=False) as data:
        projected_fixture_ids = data["fixture_ids"].astype(str)
        projected_elements = data["element_ids"].astype(str)
    mapping = dict(zip(projected_fixture_ids.tolist(), projected_elements.tolist()))
    element_index = {mapping[fixture]: index for index, fixture in enumerate(fixture_ids)}
    record = next(row for row in config["models"] if row["model_id"] == model)
    if (
        allele_order != ("REF", "ALT", "REF_RC", "ALT_RC")
        or embeddings.shape != (1033, 4, int(record["hidden_width"]))
        or len(mapping) != 1033
        or set(element_index) != set(elements)
        or not np.isfinite(embeddings).all()
    ):
        raise DNASeededHeadFitError(f"{model} raw alignment differs")
    return embeddings[np.asarray([element_index[element] for element in elements])]


def _load_outer_features(
    projection_root: Path,
    model: str,
    held_fold: int,
    elements: Sequence[str],
    folds: np.ndarray,
) -> np.ndarray:
    path = projection_root / f"{model}/heldout_fold{held_fold}/head_features.npz"
    with np.load(path, allow_pickle=False) as data:
        ids = data["element_ids"].astype(str)
        saved_folds = data["outer_folds"].astype(np.int64)
        order = tuple(data["feature_block_order"].astype(str).tolist())
        features = data["features"]
    lookup = {element: index for index, element in enumerate(ids)}
    reorder = np.asarray([lookup[element] for element in elements])
    if (
        len(lookup) != 1033
        or order != ("REF", "ALT", "ALT_minus_REF", "absolute_ALT_minus_REF")
        or features.shape != (1033, 1024)
        or not np.array_equal(saved_folds[reorder], folds)
        or not np.isfinite(features).all()
    ):
        raise DNASeededHeadFitError(f"{model} outer projection differs")
    return features[reorder].astype(np.float64)


def fit_seeded_dual_projection_fold(
    *,
    inner_features: np.ndarray,
    outer_features: np.ndarray,
    outcomes: np.ndarray,
    folds: np.ndarray,
    block_ids: np.ndarray,
    held_fold: int,
    seed: int,
    top_k_grid: Sequence[int],
    alphas: Sequence[float],
) -> tuple[np.ndarray, dict[str, Any], dict[str, np.ndarray]]:
    """Select on inner-training projections and refit on outer-training projections."""
    if inner_features.shape != outer_features.shape or outcomes.shape != (inner_features.shape[0],):
        raise DNASeededHeadFitError("dual projection geometry differs")
    inner_fold = (held_fold + 1) % 5
    outer_test = folds == held_fold
    outer_train = ~outer_test
    inner_valid = folds == inner_fold
    inner_train = outer_train & ~inner_valid
    if set(block_ids[outer_test]) & set(block_ids[outer_train]):
        raise DNASeededHeadFitError("held long-range block overlaps training")
    # Reuse the frozen Enformer/Sei draws so shared controls and seed-level
    # comparisons are identical across family-native campaigns.
    inner_boot = _bootstrap_indices(
        inner_train,
        block_ids,
        seed=seed,
        salt=INNER_BOOTSTRAP_SALT_OFFSET + held_fold,
    )
    mean, scale, rank = _training_transform(inner_features, outcomes, inner_boot)
    choices: list[tuple[float, int, float]] = []
    for k in sorted({min(int(value), int(rank.size)) for value in top_k_grid if int(value) > 0}):
        if k == 0:
            continue
        selected = rank[:k]
        x_train = (inner_features[inner_boot][:, selected] - mean[selected]) / scale[selected]
        x_valid = (inner_features[inner_valid][:, selected] - mean[selected]) / scale[selected]
        for alpha in alphas:
            model = Ridge(alpha=float(alpha), fit_intercept=True, solver="lsqr", tol=1.0e-7)
            model.fit(x_train, outcomes[inner_boot])
            prediction = model.predict(x_valid)
            choices.append((float(np.sqrt(np.mean((outcomes[inner_valid] - prediction) ** 2))), k, float(alpha)))
    constant_rmse = float(np.sqrt(np.mean((outcomes[inner_valid] - outcomes[inner_boot].mean()) ** 2)))
    best_rmse, best_k, best_alpha = min([*choices, (constant_rmse, 0, math.inf)], key=lambda value: (value[0], value[1], value[2]))
    outer_boot = _bootstrap_indices(
        outer_train,
        block_ids,
        seed=seed,
        salt=OUTER_BOOTSTRAP_SALT_OFFSET + held_fold,
    )
    outer_mean, outer_scale, outer_rank = _training_transform(outer_features, outcomes, outer_boot)
    if best_k == 0 or outer_rank.size == 0:
        selected = np.empty(0, dtype=np.int64)
        coefficient = np.empty(0, dtype=np.float64)
        intercept = float(outcomes[outer_boot].mean())
        prediction = np.full(int(outer_test.sum()), intercept)
        selected_alpha: float | str = "constant_training_mean"
    else:
        selected = outer_rank[: min(best_k, int(outer_rank.size))]
        x_train = (outer_features[outer_boot][:, selected] - outer_mean[selected]) / outer_scale[selected]
        x_test = (outer_features[outer_test][:, selected] - outer_mean[selected]) / outer_scale[selected]
        model = Ridge(alpha=best_alpha, fit_intercept=True, solver="lsqr", tol=1.0e-7)
        model.fit(x_train, outcomes[outer_boot])
        prediction = np.asarray(model.predict(x_test), dtype=np.float64)
        coefficient = np.asarray(model.coef_, dtype=np.float64)
        intercept = float(model.intercept_)
        selected_alpha = best_alpha
    receipt = {
        "held_out_fold": held_fold,
        "inner_validation_fold": inner_fold,
        "seed": seed,
        "inner_training_elements": int(inner_train.sum()),
        "inner_validation_elements": int(inner_valid.sum()),
        "outer_training_elements": int(outer_train.sum()),
        "outer_test_elements": int(outer_test.sum()),
        "inner_training_blocks": len(set(block_ids[inner_train].tolist())),
        "inner_validation_blocks": len(set(block_ids[inner_valid].tolist())),
        "outer_training_blocks": len(set(block_ids[outer_train].tolist())),
        "outer_test_blocks": len(set(block_ids[outer_test].tolist())),
        "inner_bootstrap_rows": int(inner_boot.size),
        "outer_bootstrap_rows": int(outer_boot.size),
        "selected_feature_count": int(selected.size),
        "selected_alpha": selected_alpha,
        "selected_inner_rmse": best_rmse,
        "constant_inner_rmse": constant_rmse,
        "inner_projection_fit_on_inner_training_only": True,
        "outer_projection_fit_on_outer_training_only": True,
        "held_outcomes_used_for_fit_or_tuning": False,
        "held_features_used_for_preprocessing_selection_or_tuning": False,
        "held_blocks_overlap_training": False,
        "posthoc_calibration": False,
    }
    state = {
        "selected_feature_indices0": selected,
        "selected_feature_mean": outer_mean[selected],
        "selected_feature_scale": outer_scale[selected],
        "coefficient": coefficient,
        "intercept": np.asarray(intercept),
        "selected_alpha": np.asarray(best_alpha),
    }
    return prediction, receipt, state


def fit_campaign(arguments: argparse.Namespace) -> dict[str, Any]:
    root = arguments.root.resolve(strict=True)
    config = load_config(arguments.config.resolve(strict=True))
    if arguments.output.exists():
        raise DNASeededHeadFitError("fit output exists")
    task_tree = arguments.task_spec_tree.resolve(strict=True)
    if file_sha256(task_tree / "ARTIFACTS.json") != arguments.task_spec_sha256:
        raise DNASeededHeadFitError("TaskSpec identity differs")
    verify_frozen_tree(task_tree)
    task = _load_json(task_tree / "receipt.json")
    if task.get("status") != "pass_complete_feature_reconciliation_and_prespecified_taskspec" or task.get("outcomes_read") or task.get("feature_incomplete_models"):
        raise DNASeededHeadFitError("TaskSpec feature or outcome boundary differs")
    open_tree = _tree(root, config["open_projection_authority"], label="open_projection")
    open_projection = open_tree / "open_sequence_projection"
    nt_projection = arguments.nt_projection_tree.resolve(strict=True)
    if file_sha256(nt_projection / "ARTIFACTS.json") != arguments.nt_projection_sha256:
        raise DNASeededHeadFitError("restricted projection identity differs")
    verify_frozen_tree(nt_projection)
    row_tree = _tree(root, config["row_authority"], label="row_authority")
    outcome_tree = _tree(root, config["outcome_authority"], label="outcome_authority")
    allele_tree = _tree(root, config["allele_authority"], label="allele_authority")
    row_path = row_tree / config["row_authority"]["member"]
    outcome_path = outcome_tree / config["outcome_authority"]["member"]
    allele_path = allele_tree / config["allele_authority"]["member"]
    rows = _read_tsv(row_path, ROW_FIELDS)
    if len(rows) != 10330:
        raise DNASeededHeadFitError("row denominator differs")
    metadata: dict[str, tuple[str, str, str]] = {}
    for row in rows:
        value = (row["source_locus_group_id"], row["long_range_block_id"], row["outer_fold"])
        if metadata.setdefault(row["element_id"], value) != value:
            raise DNASeededHeadFitError("element metadata differs")
    # Preserve the frozen row-universe insertion order.  Re-sorting the same
    # rows changes floating-point reduction order and breaks bit identity of
    # controls shared with the Enformer/Sei campaign.
    elements = list(metadata)
    element_index = {element: index for index, element in enumerate(elements)}
    folds = np.asarray([int(metadata[element][2].removeprefix("fold-")) for element in elements], dtype=np.int64)
    blocks = np.asarray([metadata[element][1] for element in elements])
    if len(elements) != 1033 or len(set(blocks.tolist())) != 239:
        raise DNASeededHeadFitError("element or block census differs")
    outcomes = load_outcomes(outcome_path, set(elements))
    targets = {context: np.asarray([outcomes[(element, context)]["mean"] for element in elements]) for context in CONTEXTS}
    # Verify the allele authority even though the shared prediction control is
    # reused verbatim to avoid cross-CPU last-bit drift.
    _load_alleles(allele_path, elements)

    projection_roots = {model: (open_projection if model != "nucleotide_transformer" else nt_projection) for model in MODELS}
    raw = {model: _load_raw_embeddings(root, config, model, projection_roots[model], elements) for model in MODELS}
    inner_features: dict[tuple[str, int], np.ndarray] = {}
    outer_features: dict[tuple[str, int], np.ndarray] = {}
    projection_rows = []
    arguments.output.mkdir(parents=True, mode=0o750)
    (arguments.output / "inner_projections").mkdir()
    (arguments.output / "heads").mkdir()
    for model in MODELS:
        model_root = arguments.output / "inner_projections" / model
        model_root.mkdir()
        for held_fold in range(5):
            inner_fold = (held_fold + 1) % 5
            inner_train = (folds != held_fold) & (folds != inner_fold)
            parameters = fit_inner_projection(raw[model], inner_train)
            inner = head_features(apply_projection(raw[model], parameters)).astype(np.float64)
            outer = _load_outer_features(projection_roots[model], model, held_fold, elements, folds)
            inner_features[(model, held_fold)] = inner
            outer_features[(model, held_fold)] = outer
            np.savez_compressed(model_root / f"outer{held_fold}_inner_projection.npz", **parameters)
            projection_rows.append({
                "model_id": model,
                "held_out_fold": held_fold,
                "inner_validation_fold": inner_fold,
                "inner_projection_training_elements": int(inner_train.sum()),
                "inner_projection_fit_on_inner_training_only": "true",
                "outer_projection_fit_on_outer_training_only": "true",
                "held_outcomes_used": "false",
                "obsolete_6kb_folds_reused": "false",
            })

    shared_control_rows, selections = _load_frozen_shared_control(root, config, rows)
    predictions: dict[tuple[str, str, int, str], np.ndarray] = {}
    alphas = tuple(float(value) for value in config["fit"]["alpha_grid"])
    for context in CONTEXTS:
        y = targets[context]
        for seed in SEEDS:
            predictions[("available_simple_controls", "zero", seed, context)] = np.zeros(1033)
            mean_oof = np.empty(1033)
            for held_fold in range(5):
                test = folds == held_fold
                train = ~test
                sampled = _bootstrap_indices(
                    train,
                    blocks,
                    seed=seed,
                    salt=TRAINING_MEAN_BOOTSTRAP_SALT_OFFSET + held_fold,
                )
                mean_oof[test] = float(y[sampled].mean())
                selections.append({
                    "model_id": "available_simple_controls", "head_id": "outer_training_mean", "assay_context_id": context,
                    "seed": seed, "held_out_fold": held_fold, "outer_training_elements": int(train.sum()),
                    "outer_test_elements": int(test.sum()), "selected_feature_count": 0, "selected_alpha": "not_applicable",
                    "held_outcomes_used_for_fit_or_tuning": False, "held_features_used_for_preprocessing_selection_or_tuning": False,
                    "held_blocks_overlap_training": False, "posthoc_calibration": False,
                })
            predictions[("available_simple_controls", "outer_training_mean", seed, context)] = mean_oof
        for model in MODELS:
            for head in HEADS:
                feature_slice = slice(512, 768) if head == "delta_ridge" else slice(0, 1024)
                top_k = config["fit"]["delta_top_k_grid"] if head == "delta_ridge" else config["fit"]["full_top_k_grid"]
                for seed in SEEDS:
                    oof = np.empty(1033)
                    for held_fold in range(5):
                        test = folds == held_fold
                        prediction, receipt, state = fit_seeded_dual_projection_fold(
                            inner_features=inner_features[(model, held_fold)][:, feature_slice],
                            outer_features=outer_features[(model, held_fold)][:, feature_slice],
                            outcomes=y, folds=folds, block_ids=blocks, held_fold=held_fold, seed=seed,
                            top_k_grid=top_k, alphas=alphas,
                        )
                        oof[test] = prediction
                        state_root = arguments.output / "heads" / model / head / context
                        state_root.mkdir(parents=True, exist_ok=True)
                        np.savez_compressed(state_root / f"seed{seed}_fold{held_fold}.npz", **state)
                        selections.append({"model_id": model, "head_id": head, "assay_context_id": context, **receipt})
                    if not np.isfinite(oof).all():
                        raise DNASeededHeadFitError("OOF prediction coverage differs")
                    predictions[(model, head, seed, context)] = oof

    output_rows = []
    for model, head in CANDIDATES:
        for row in rows:
            seed = int(row["seed"])
            if model == "available_simple_controls" and head == "allele_identity_ridge":
                output_rows.append(dict(shared_control_rows[(seed, row["row_hash"])]))
                continue
            index = element_index[row["element_id"]]
            value = float(predictions[(model, head, seed, row["assay_context_id"])][index])
            if not math.isfinite(value):
                raise DNASeededHeadFitError("non-finite prediction")
            output_rows.append({
                **row,
                "model_id": model,
                "head_id": head,
                "prediction": format(value, ".17g"),
                "experimental_replicates": 4,
                "biological_donors": 0,
                "outcome_role": "exposed_development_MPRA_only",
            })
    expected = len(CANDIDATES) * 10330
    if len(output_rows) != expected:
        raise DNASeededHeadFitError("prediction row count differs")
    _write_gzip_tsv(arguments.output / "oof_predictions.tsv.gz", output_rows)
    _write_tsv(arguments.output / "head_selection.tsv", selections)
    _write_tsv(arguments.output / "projection_audit.tsv", projection_rows)
    receipt = {
        "schema_version": SCHEMA,
        "status": "pass_complete_dna_language_seeded_oof_fit",
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "models": list(MODELS),
        "candidate_count": len(CANDIDATES),
        "prediction_rows": expected,
        "rows_per_candidate": 10330,
        "elements": 1033,
        "long_range_blocks": 239,
        "fixed_seeds": list(SEEDS),
        "five_seed_strategy": "outer_training_long_range_block_bootstrap",
        "bootstrap_draw_contract": "exact_frozen_enformer_sei_offsets",
        "shared_control_prediction_policy": "reuse_bit_exact_frozen_enformer_sei_oof",
        "shared_control_source_artifacts_sha256": config["shared_control_reference"]["artifacts_sha256"],
        "shared_control_oof_rows": int(config["shared_control_reference"]["expected_rows"]),
        "shared_control_rows_copied_without_numerical_transformation": True,
        "shared_control_refit_performed": False,
        "shared_control_counted_as_new_evidence": False,
        "shared_control_seeds_counted_as_independent_evidence": False,
        "bootstrap_salt_offsets": {
            "inner": INNER_BOOTSTRAP_SALT_OFFSET,
            "outer": OUTER_BOOTSTRAP_SALT_OFFSET,
            "training_mean": TRAINING_MEAN_BOOTSTRAP_SALT_OFFSET,
        },
        "inner_projection_fit_on_inner_training_only": True,
        "outer_projection_fit_on_outer_training_only": True,
        "features_imputed": False,
        "obsolete_6kb_folds_reused": False,
        "task_spec_artifacts_sha256": arguments.task_spec_sha256,
        "nucleotide_transformer_projection_artifacts_sha256": arguments.nt_projection_sha256,
        "biological_donors": 0,
        "experimental_replicates_per_context": 4,
        "mandatory_task_native_baselines_complete": False,
        "missing_baselines": list(MISSING_BASELINES),
        "metrics_calculated": False,
        "models_ranked": False,
        "shortlist_created": False,
        "champion_claim": False,
        "sealed_assets_read": False,
        "base_checkpoints_tuned": False,
        "posthoc_calibration": False,
        "stack_fit": False,
    }
    (arguments.output / "receipt.json").write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--task-spec-tree", type=Path, required=True)
    parser.add_argument("--task-spec-sha256", required=True)
    parser.add_argument("--nt-projection-tree", type=Path, required=True)
    parser.add_argument("--nt-projection-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    fit_campaign(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
