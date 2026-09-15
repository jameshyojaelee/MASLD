#!/usr/bin/env python3
"""Fit five block-bootstrap, outer-training-only Enformer/Sei MPRA heads."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import io
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from sklearn.linear_model import Ridge

from masld_bench.artifacts import ArtifactError, verify_frozen_tree
from scripts.evaluate_gse281364_dna_lm_common_lane import load_outcomes
from scripts.freeze_gse281364_enformer_sei_head_taskspec import (
    CANDIDATES,
    CONTEXTS,
    MISSING_BASELINES,
    PREDICTION_FIELDS,
    SEEDS,
    file_sha256,
    load_config,
    validate_config,
)


SCHEMA = "masld-bench-gse281364-enformer-sei-head-fit-v1"
BASE_FIELDS = (
    "element_id",
    "source_locus_group_id",
    "long_range_block_id",
    "outer_fold",
    "enformer_feature_index0",
    "sei_feature_index0",
    "legacy_outer_fold",
    "fold_reassigned",
    "enformer_hepg2_accessibility_sad",
    "enformer_hepg2_accessibility_sar",
    "enformer_liver_accessibility_sad",
    "enformer_liver_accessibility_sar",
    "enformer_all_accessibility_mean_sad",
    "enformer_all_accessibility_mean_sar",
    "sei_max_abs_sequence_class_index0",
    "sei_signed_max_abs_score",
    "sei_max_abs_score",
    "sei_mean_signed_score",
)
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
    "enformer_feature_index0",
    "sei_feature_index0",
    "legacy_outer_fold",
    "fold_reassigned",
)


class EnformerSeiHeadFitError(RuntimeError):
    """Raised when a fit input or isolation invariant differs."""


def _load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise EnformerSeiHeadFitError(f"invalid JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise EnformerSeiHeadFitError(f"JSON object required: {path}")
    return value


def _safe_tree(root: Path, binding: Mapping[str, Any], *, label: str) -> Path:
    value = binding.get("tree_path")
    if not isinstance(value, str) or not value:
        raise EnformerSeiHeadFitError(f"{label} tree path missing")
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise EnformerSeiHeadFitError(f"unsafe {label} tree path")
    try:
        tree = (root / relative).resolve(strict=True)
        tree.relative_to(root)
    except (OSError, ValueError) as error:
        raise EnformerSeiHeadFitError(f"{label} tree missing or escapes root") from error
    if file_sha256(tree / "ARTIFACTS.json") != binding.get("artifacts_sha256"):
        raise EnformerSeiHeadFitError(f"{label} ARTIFACTS identity differs")
    try:
        verify_frozen_tree(tree)
    except ArtifactError as error:
        raise EnformerSeiHeadFitError(f"{label} frozen tree differs: {error}") from error
    return tree


def _read_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise EnformerSeiHeadFitError(f"TSV schema differs: {path}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise EnformerSeiHeadFitError(f"empty TSV: {path}")
    return rows


def _read_gzip_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise EnformerSeiHeadFitError(f"gzip TSV schema differs: {path}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise EnformerSeiHeadFitError(f"empty gzip TSV: {path}")
    return rows


def _write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise EnformerSeiHeadFitError(f"cannot write empty TSV: {path}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_gzip_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerows(rows)


def _array_sha256(values: np.ndarray) -> str:
    contiguous = np.ascontiguousarray(values)
    digest = sha256()
    digest.update(str(contiguous.dtype).encode())
    digest.update(str(contiguous.shape).encode())
    digest.update(contiguous.tobytes())
    return digest.hexdigest()


def _bootstrap_indices(
    eligible: np.ndarray, block_ids: np.ndarray, *, seed: int, salt: int
) -> np.ndarray:
    indices = np.flatnonzero(eligible)
    blocks = sorted(set(block_ids[indices].tolist()))
    if len(blocks) < 2:
        raise EnformerSeiHeadFitError("training block bootstrap has fewer than two blocks")
    by_block = {block: indices[block_ids[indices] == block] for block in blocks}
    rng = np.random.default_rng(np.random.SeedSequence([seed, salt]))
    sampled = rng.integers(0, len(blocks), size=len(blocks))
    return np.concatenate([by_block[blocks[index]] for index in sampled])


def _training_transform(
    features: np.ndarray, outcomes: np.ndarray, training_indices: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    training = features[training_indices]
    mean = training.mean(axis=0)
    scale = training.std(axis=0)
    usable = np.flatnonzero(scale > 1.0e-8)
    if usable.size == 0:
        return mean, np.maximum(scale, 1.0e-8), usable
    normalized = (training[:, usable] - mean[usable]) / scale[usable]
    centered_y = outcomes[training_indices] - outcomes[training_indices].mean()
    denominator = np.sqrt(np.sum(normalized**2, axis=0) * np.sum(centered_y**2))
    correlation = np.divide(
        np.abs(normalized.T @ centered_y),
        denominator,
        out=np.zeros(usable.size, dtype=np.float64),
        where=denominator > 0,
    )
    order = np.lexsort((usable, -correlation))
    return mean, np.maximum(scale, 1.0e-8), usable[order]


def fit_seeded_outer_fold(
    *,
    features: np.ndarray,
    outcomes: np.ndarray,
    folds: np.ndarray,
    block_ids: np.ndarray,
    held_fold: int,
    seed: int,
    top_k_grid: Sequence[int],
    alphas: Sequence[float],
) -> tuple[np.ndarray, dict[str, Any], dict[str, np.ndarray]]:
    """Fit one held-fold head without reading held outcomes for any fit choice."""
    if features.ndim != 2 or outcomes.shape != (features.shape[0],):
        raise EnformerSeiHeadFitError("feature/outcome shapes differ")
    inner_fold = (held_fold + 1) % 5
    outer_test = folds == held_fold
    outer_train = ~outer_test
    inner_valid = folds == inner_fold
    inner_train = outer_train & ~inner_valid
    if (
        np.any(inner_train & inner_valid)
        or np.any(inner_train & outer_test)
        or np.any(inner_valid & outer_test)
        or set(block_ids[outer_test]) & set(block_ids[outer_train])
    ):
        raise EnformerSeiHeadFitError("held-block isolation differs")
    inner_boot = _bootstrap_indices(
        inner_train, block_ids, seed=seed, salt=10000 + held_fold
    )
    inner_mean, inner_scale, inner_rank = _training_transform(
        features, outcomes, inner_boot
    )
    candidates: list[tuple[float, int, float]] = []
    effective_k = sorted({min(int(k), int(inner_rank.size)) for k in top_k_grid if int(k) > 0})
    for k in effective_k:
        if k == 0:
            continue
        selected = inner_rank[:k]
        x_train = (features[inner_boot][:, selected] - inner_mean[selected]) / inner_scale[selected]
        x_valid = (features[inner_valid][:, selected] - inner_mean[selected]) / inner_scale[selected]
        for alpha in alphas:
            model = Ridge(alpha=float(alpha), fit_intercept=True, solver="lsqr", tol=1.0e-7)
            model.fit(x_train, outcomes[inner_boot])
            prediction = model.predict(x_valid)
            rmse = float(np.sqrt(np.mean((outcomes[inner_valid] - prediction) ** 2)))
            candidates.append((rmse, k, float(alpha)))
    constant_rmse = float(
        np.sqrt(
            np.mean(
                (outcomes[inner_valid] - float(outcomes[inner_boot].mean())) ** 2
            )
        )
    )
    best_rmse, best_k, best_alpha = min(
        [*candidates, (constant_rmse, 0, math.inf)],
        key=lambda value: (value[0], value[1], value[2]),
    )
    outer_boot = _bootstrap_indices(
        outer_train, block_ids, seed=seed, salt=20000 + held_fold
    )
    outer_mean, outer_scale, outer_rank = _training_transform(
        features, outcomes, outer_boot
    )
    if best_k == 0 or outer_rank.size == 0:
        selected = np.empty(0, dtype=np.int64)
        coefficient = np.empty(0, dtype=np.float64)
        intercept = float(outcomes[outer_boot].mean())
        prediction = np.full(int(outer_test.sum()), intercept, dtype=np.float64)
        selected_alpha: float | str = "constant_training_mean"
    else:
        selected = outer_rank[: min(best_k, int(outer_rank.size))]
        x_train = (features[outer_boot][:, selected] - outer_mean[selected]) / outer_scale[selected]
        x_test = (features[outer_test][:, selected] - outer_mean[selected]) / outer_scale[selected]
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
        "inner_bootstrap_unique_blocks": len(set(block_ids[inner_boot].tolist())),
        "outer_bootstrap_rows": int(outer_boot.size),
        "outer_bootstrap_unique_blocks": len(set(block_ids[outer_boot].tolist())),
        "selected_feature_count": int(selected.size),
        "selected_alpha": selected_alpha,
        "selected_inner_rmse": best_rmse,
        "constant_inner_rmse": constant_rmse,
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
        "intercept": np.asarray(intercept, dtype=np.float64),
        "selected_alpha": np.asarray(best_alpha, dtype=np.float64),
    }
    return prediction, receipt, state


def _load_feature_matrix(
    root: Path,
    contract: Mapping[str, Any],
    model: str,
    feature: str,
    expected_shape: tuple[int, int],
) -> np.ndarray:
    binding = (
        contract["enformer"]["feature_matrices"][feature]
        if model == "enformer"
        else contract["sei"]["class_matrix"]
    )
    tree = _safe_tree(root, binding, label=f"{model}_{feature}")
    member = tree / binding["member"]
    if file_sha256(member) != binding["sha256"]:
        raise EnformerSeiHeadFitError(f"{model} {feature} member differs")
    values = np.load(member, allow_pickle=False)
    if values.shape != expected_shape or str(values.dtype) != "float32" or not np.isfinite(values).all():
        raise EnformerSeiHeadFitError(f"{model} {feature} matrix contract differs")
    return values.astype(np.float64)


def _load_allele_features(path: Path, elements: Sequence[str]) -> np.ndarray:
    wanted = set(elements)
    pairs: dict[str, tuple[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"element_id", "ref", "alt"}
        if not required.issubset(reader.fieldnames or ()):
            raise EnformerSeiHeadFitError("allele manifest schema differs")
        for row in reader:
            if row["element_id"] in wanted:
                value = (row["ref"], row["alt"])
                if row["element_id"] in pairs or value[0] not in "ACGT" or value[1] not in "ACGT" or value[0] == value[1]:
                    raise EnformerSeiHeadFitError("selected allele identity differs")
                pairs[row["element_id"]] = value
    if set(pairs) != wanted:
        raise EnformerSeiHeadFitError("selected allele coverage differs")
    bases = {base: index for index, base in enumerate("ACGT")}
    features = np.zeros((len(elements), 16), dtype=np.float64)
    for index, element in enumerate(elements):
        ref, alt = pairs[element]
        features[index, 4 * bases[ref] + bases[alt]] = 1.0
    return features


def fit_campaign(arguments: argparse.Namespace) -> dict[str, Any]:
    root = arguments.root.resolve(strict=True)
    config = load_config(arguments.config.resolve(strict=True))
    validate_config(config)
    if arguments.output.exists():
        raise EnformerSeiHeadFitError("fit output already exists")
    task_tree = arguments.task_spec_tree.resolve(strict=True)
    if file_sha256(task_tree / "ARTIFACTS.json") != arguments.task_spec_sha256:
        raise EnformerSeiHeadFitError("TaskSpec ARTIFACTS identity differs")
    verify_frozen_tree(task_tree)
    task_receipt = _load_json(task_tree / "task_spec.json")
    if task_receipt.get("status") != "pass_prespecified_outcome_blind_taskspec" or task_receipt.get("outcomes_read"):
        raise EnformerSeiHeadFitError("TaskSpec outcome firewall differs")

    static_tree = _safe_tree(root, config["static_authority"], label="static_authority")
    outcome_tree = _safe_tree(root, config["outcome_authority"], label="outcome_authority")
    allele_tree = _safe_tree(root, config["allele_authority"], label="allele_authority")
    static_binding = config["static_authority"]
    for member_key, hash_key in (
        ("base_member", "base_sha256"),
        ("row_member", "row_sha256"),
        ("contract_member", "contract_sha256"),
    ):
        if file_sha256(static_tree / static_binding[member_key]) != static_binding[hash_key]:
            raise EnformerSeiHeadFitError("static member identity differs")
    outcome_path = outcome_tree / config["outcome_authority"]["member"]
    allele_path = allele_tree / config["allele_authority"]["member"]
    if file_sha256(outcome_path) != config["outcome_authority"]["member_sha256"]:
        raise EnformerSeiHeadFitError("outcome member identity differs")
    if file_sha256(allele_path) != config["allele_authority"]["member_sha256"]:
        raise EnformerSeiHeadFitError("allele member identity differs")

    base_rows = _read_tsv(static_tree / static_binding["base_member"], BASE_FIELDS)
    row_rows = _read_gzip_tsv(static_tree / static_binding["row_member"], ROW_FIELDS)
    if len(base_rows) != 1033 or len(row_rows) != 10330:
        raise EnformerSeiHeadFitError("static row denominator differs")
    base_rows.sort(key=lambda row: int(row["enformer_feature_index0"]))
    if [int(row["enformer_feature_index0"]) for row in base_rows] != list(range(1033)) or [int(row["sei_feature_index0"]) for row in base_rows] != list(range(1033)):
        raise EnformerSeiHeadFitError("feature index map differs")
    elements = [row["element_id"] for row in base_rows]
    if len(set(elements)) != 1033:
        raise EnformerSeiHeadFitError("element identity differs")
    folds = np.asarray([int(row["outer_fold"].removeprefix("fold-")) for row in base_rows], dtype=np.int64)
    block_ids = np.asarray([row["long_range_block_id"] for row in base_rows])
    if set(folds.tolist()) != set(range(5)) or len(set(block_ids.tolist())) != 239:
        raise EnformerSeiHeadFitError("fold or block census differs")
    for block in set(block_ids.tolist()):
        if len(set(folds[block_ids == block].tolist())) != 1:
            raise EnformerSeiHeadFitError("long-range block crosses outer folds")
    element_index = {element: index for index, element in enumerate(elements)}
    row_lookup = {(int(row["seed"]), row["row_hash"]): row for row in row_rows}
    if len(row_lookup) != 10330:
        raise EnformerSeiHeadFitError("row authority identities differ")

    outcomes = load_outcomes(outcome_path, set(elements))
    targets = {
        context: np.asarray([outcomes[(element, context)]["mean"] for element in elements], dtype=np.float64)
        for context in CONTEXTS
    }
    if any(not np.isfinite(values).all() for values in targets.values()):
        raise EnformerSeiHeadFitError("aggregated development outcome is non-finite")
    contract = _load_json(static_tree / static_binding["contract_member"])
    feature_sets = {
        ("available_simple_controls", "allele_identity_ridge"): _load_allele_features(allele_path, elements),
        ("enformer_crested_restricted_port", "rc_ensemble_sad_ridge"): _load_feature_matrix(root, contract, "enformer", "rc_ensemble_sad", (1033, 5313)),
        ("enformer_crested_restricted_port", "rc_ensemble_sar_ridge"): _load_feature_matrix(root, contract, "enformer", "rc_ensemble_sar", (1033, 5313)),
        ("sei", "native_40class_ridge"): _load_feature_matrix(root, contract, "sei", "native_40class", (1033, 40)),
    }
    grids = {
        ("available_simple_controls", "allele_identity_ridge"): config["fit"]["allele_top_k_grid"],
        ("enformer_crested_restricted_port", "rc_ensemble_sad_ridge"): config["fit"]["enformer_top_k_grid"],
        ("enformer_crested_restricted_port", "rc_ensemble_sar_ridge"): config["fit"]["enformer_top_k_grid"],
        ("sei", "native_40class_ridge"): config["fit"]["sei_top_k_grid"],
    }
    alphas = tuple(float(value) for value in config["fit"]["alpha_grid"])
    arguments.output.mkdir(parents=True, mode=0o750)
    (arguments.output / "heads").mkdir()
    predictions: dict[tuple[str, str, int, str], np.ndarray] = {}
    selection_rows: list[dict[str, Any]] = []
    for context in CONTEXTS:
        y = targets[context]
        for seed in SEEDS:
            zero = np.zeros(1033, dtype=np.float64)
            training_mean = np.empty(1033, dtype=np.float64)
            for held_fold in range(5):
                test = folds == held_fold
                train = ~test
                bootstrap = _bootstrap_indices(train, block_ids, seed=seed, salt=30000 + held_fold)
                training_mean[test] = float(y[bootstrap].mean())
                selection_rows.append({
                    "model_id": "available_simple_controls",
                    "head_id": "outer_training_mean",
                    "assay_context_id": context,
                    "seed": seed,
                    "held_out_fold": held_fold,
                    "outer_training_elements": int(train.sum()),
                    "outer_test_elements": int(test.sum()),
                    "outer_training_blocks": len(set(block_ids[train].tolist())),
                    "outer_test_blocks": len(set(block_ids[test].tolist())),
                    "outer_bootstrap_rows": int(bootstrap.size),
                    "outer_bootstrap_unique_blocks": len(set(block_ids[bootstrap].tolist())),
                    "selected_feature_count": 0,
                    "selected_alpha": "not_applicable",
                    "held_outcomes_used_for_fit_or_tuning": False,
                    "held_features_used_for_preprocessing_selection_or_tuning": False,
                    "held_blocks_overlap_training": False,
                    "posthoc_calibration": False,
                })
            predictions[("available_simple_controls", "zero", seed, context)] = zero
            predictions[("available_simple_controls", "outer_training_mean", seed, context)] = training_mean

        for candidate, features in feature_sets.items():
            model_id, head_id = candidate
            for context_seed in SEEDS:
                oof = np.empty(1033, dtype=np.float64)
                for held_fold in range(5):
                    test = folds == held_fold
                    prediction, selection, state = fit_seeded_outer_fold(
                        features=features,
                        outcomes=y,
                        folds=folds,
                        block_ids=block_ids,
                        held_fold=held_fold,
                        seed=context_seed,
                        top_k_grid=grids[candidate],
                        alphas=alphas,
                    )
                    oof[test] = prediction
                    state_root = arguments.output / "heads" / model_id / head_id / context
                    state_root.mkdir(parents=True, exist_ok=True)
                    np.savez_compressed(state_root / f"seed{context_seed}_fold{held_fold}.npz", **state)
                    selection_rows.append({
                        "model_id": model_id,
                        "head_id": head_id,
                        "assay_context_id": context,
                        **selection,
                    })
                if not np.isfinite(oof).all():
                    raise EnformerSeiHeadFitError("OOF prediction coverage differs")
                predictions[(model_id, head_id, context_seed, context)] = oof

    prediction_rows = []
    for model_id, head_id in CANDIDATES:
        for row in row_rows:
            seed = int(row["seed"])
            context = row["assay_context_id"]
            index = element_index[row["element_id"]]
            value = float(predictions[(model_id, head_id, seed, context)][index])
            if not math.isfinite(value):
                raise EnformerSeiHeadFitError("non-finite standardized prediction")
            prediction_rows.append({
                **{field: row[field] for field in PREDICTION_FIELDS[:11]},
                "model_id": model_id,
                "head_id": head_id,
                "prediction": format(value, ".17g"),
                "experimental_replicates": 4,
                "biological_donors": 0,
                "outcome_role": "exposed_development_MPRA_only",
            })
    expected_rows = len(CANDIDATES) * 10330
    if len(prediction_rows) != expected_rows:
        raise EnformerSeiHeadFitError("prediction denominator differs")
    _write_gzip_tsv(arguments.output / "oof_predictions.tsv.gz", PREDICTION_FIELDS, prediction_rows)
    _write_tsv(arguments.output / "head_selection.tsv", selection_rows)
    feature_receipt = {
        f"{model}/{head}": {"shape": list(values.shape), "array_sha256": _array_sha256(values)}
        for (model, head), values in feature_sets.items()
    }
    receipt = {
        "schema_version": SCHEMA,
        "status": "pass_restricted_development_oof_fit",
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "candidate_count": len(CANDIDATES),
        "prediction_rows": expected_rows,
        "rows_per_candidate": 10330,
        "elements": 1033,
        "long_range_blocks": 239,
        "contexts": list(CONTEXTS),
        "fixed_seeds": list(SEEDS),
        "five_seed_strategy": "outer_training_long_range_block_bootstrap",
        "seeds_are_biological_replicates": False,
        "feature_receipt": feature_receipt,
        "task_spec_artifacts_sha256": arguments.task_spec_sha256,
        "outcome_artifacts_sha256": config["outcome_authority"]["artifacts_sha256"],
        "biological_donors": 0,
        "experimental_replicates_per_context": 4,
        "experimental_replicates_used_as_independent_donors": False,
        "mandatory_task_native_baselines_complete": False,
        "missing_baselines": list(MISSING_BASELINES),
        "family_native_feature_spaces_combined": False,
        "metrics_calculated": False,
        "models_ranked_as_interchangeable": False,
        "shortlist_created": False,
        "champion_claim": False,
        "native_equivalence_claim": False,
        "external_evaluation": False,
        "sealed_assets_read": False,
        "base_checkpoints_tuned": False,
        "posthoc_calibration": False,
        "stack_fit": False,
        "residual_complementarity_calculated": False,
        "global_frozen_census_changed": False,
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
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    fit_campaign(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
