#!/usr/bin/env python3
"""Fit GSE135251 NAS-score regressions in the registered per-sample representations.

The endpoint is the continuous 0-to-8 activity score and is never binarised.
Both cohorts are RNA-seq, so both sides go to log2(1 + CPM) with the library
total taken over that sample's complete deposited axis, and only then onto the
shared gene axis and into one of the two registered per-sample transforms.  No
statistic is ever pooled across samples.

Ridge rather than a sparse penalty, deliberately.  A previous lane in this
campaign shipped three models whose coefficients had all been shrunk to exactly
zero and which still scored; ridge has a finite solution at every alpha and
keeps its coefficients, so that failure cannot recur here.  The degeneracy guard
is retained anyway so a sparse variant added later inherits it.

Configurations are ranked by Spearman in excess of their own permutation null,
never by raw Spearman, because two scorers with different tie structures do not
share a null.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from scripts.fit_gse267145_fibrosis_transfer_source_models import (
    REPRESENTATION_FUNCTIONS,
    assert_row_independent,
    build_representation,
)


PRIOR_MODEL_ID = "training_mean_nas_sum"
LEARNED_MODEL_IDS = (
    "gene_median_pca_ridge",
    "gene_median_ridge",
    "per_array_rank_pca_ridge",
    "per_array_rank_ridge",
)
MODEL_IDS = (PRIOR_MODEL_ID,) + LEARNED_MODEL_IDS


class NasSourceFitError(RuntimeError):
    """Raised when the source fit would violate its frozen recipe."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise NasSourceFitError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def log2_cpm_complete_axis(matrix: np.ndarray) -> np.ndarray:
    """log2(1 + CPM), library total over that sample's complete deposited axis."""

    values = np.asarray(matrix, dtype=np.float64)
    if values.ndim != 2 or not np.all(np.isfinite(values)) or np.any(values < 0):
        raise NasSourceFitError("RNA matrix is not finite and nonnegative")
    totals = values.sum(axis=1)
    if np.any(totals <= 0):
        raise NasSourceFitError("a participant RNA library is empty")
    return np.log2(1.0 + values * (1_000_000.0 / totals[:, None]))


def rank_vector(values: np.ndarray) -> np.ndarray:
    from scipy.stats import rankdata

    return rankdata(np.asarray(values, dtype=np.float64), method="average")


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman as Pearson on average ranks, computed directly for speed."""

    x = rank_vector(a)
    y = rank_vector(b)
    x = x - x.mean()
    y = y - y.mean()
    nx = float(np.linalg.norm(x))
    ny = float(np.linalg.norm(y))
    if nx <= 0 or ny <= 0:
        return 0.0
    return float(np.dot(x, y) / (nx * ny))


def spearman_null(
    observed: np.ndarray, scores: np.ndarray, *, replicates: int, seed: int
) -> dict[str, float]:
    """Permutation null for this exact score vector.

    The scores are held fixed and the outcome is permuted, so the null inherits
    the scorer's own tie structure.  Outcome ties do not narrow a rank-
    correlation null, but score ties do, and that is what makes raw Spearman
    non-comparable across configurations.
    """

    rng = np.random.default_rng(seed)
    score_ranks = rank_vector(scores)
    score_ranks = score_ranks - score_ranks.mean()
    score_norm = float(np.linalg.norm(score_ranks))
    labels = rank_vector(observed)
    labels = labels - labels.mean()
    label_norm = float(np.linalg.norm(labels))
    if score_norm <= 0 or label_norm <= 0:
        return {
            "spearman": 0.0, "null_mean": 0.0, "null_p95": 0.0,
            "excess_over_null_mean": 0.0, "permutation_p_value": 1.0,
        }
    working = labels.copy()
    draws = np.empty(replicates, dtype=np.float64)
    for index in range(replicates):
        rng.shuffle(working)
        draws[index] = float(np.dot(working, score_ranks) / (score_norm * label_norm))
    point = float(np.dot(labels, score_ranks) / (score_norm * label_norm))
    return {
        "spearman": point,
        "null_mean": float(np.mean(draws)),
        "null_p95": float(np.quantile(draws, 0.95)),
        "excess_over_null_mean": float(point - np.mean(draws)),
        "permutation_p_value": float(
            (1.0 + float(np.sum(draws >= point))) / (replicates + 1.0)
        ),
    }


def select_feature_indices(
    *, eligibility: np.ndarray, representation: np.ndarray,
    gene_ids: Sequence[str], feature_count: int,
) -> np.ndarray:
    ids = np.asarray(gene_ids, dtype=str)
    if eligibility.shape != representation.shape or representation.shape[1] != len(ids):
        raise NasSourceFitError("feature-selection axes differ")
    minimum_positive = max(3, math.ceil(0.10 * representation.shape[0]))
    eligible = np.flatnonzero(np.sum(eligibility > 0, axis=0) >= minimum_positive)
    if len(eligible) < feature_count:
        raise NasSourceFitError("shared axis has too few eligible genes")
    variances = np.var(representation[:, eligible], axis=0, ddof=0)
    order = np.lexsort((ids[eligible], -variances))
    return np.asarray(np.sort(eligible[order[:feature_count]]), dtype=np.int64)


def fit_pipeline(
    *, config: Mapping[str, Any], representation: np.ndarray, eligibility: np.ndarray,
    gene_ids: Sequence[str], fitting: np.ndarray, targets: np.ndarray,
    feature_count: int, hyperparameters: Mapping[str, Any], seed: int,
) -> dict[str, Any]:
    from sklearn.decomposition import PCA
    from sklearn.linear_model import Ridge

    selected = select_feature_indices(
        eligibility=eligibility[fitting], representation=representation[fitting],
        gene_ids=gene_ids, feature_count=feature_count,
    )
    dense = representation[fitting][:, selected]
    center = dense.mean(axis=0)
    scale = dense.std(axis=0)
    scale[scale <= np.finfo(np.float64).eps] = 1.0
    scaled = (dense - center) / scale
    state: dict[str, Any] = {"selected_indices": selected, "center": center, "scale": scale}
    if config["reducer"] == "pca":
        components = int(hyperparameters["pca_components"])
        if components >= min(scaled.shape):
            raise NasSourceFitError("partition cannot support the requested PCA size")
        pca = PCA(n_components=components, svd_solver="randomized",
                  whiten=False, random_state=seed + 90_000)
        design = pca.fit_transform(scaled)
        state["pca_mean"] = np.asarray(pca.mean_, dtype=np.float64)
        state["pca_components"] = np.asarray(pca.components_, dtype=np.float64)
    elif config["reducer"] == "none":
        design = scaled
    else:
        raise NasSourceFitError(f"unregistered reducer: {config['reducer']}")
    model = Ridge(alpha=float(hyperparameters["alpha"]), fit_intercept=True,
                  solver="auto", random_state=seed + 91_000)
    model.fit(design, targets[fitting])
    state["coef"] = np.asarray(model.coef_, dtype=np.float64).reshape(1, -1)
    state["intercept"] = np.asarray(
        np.atleast_1d(model.intercept_), dtype=np.float64
    )
    state["nonzero_coefficients"] = np.asarray(
        [int(np.sum(np.abs(state["coef"]) > 0))], dtype=np.int64
    )
    state["_design_columns"] = design.shape[1]
    return state


def apply_pipeline(
    *, state: Mapping[str, Any], representation: np.ndarray
) -> np.ndarray:
    """Apply a frozen state and return the predicted activity score."""

    selected = np.asarray(state["selected_indices"], dtype=np.int64)
    center = np.asarray(state["center"], dtype=np.float64)
    scale = np.asarray(state["scale"], dtype=np.float64)
    scaled = (np.asarray(representation, dtype=np.float64)[:, selected] - center) / scale
    if "pca_components" in state:
        design = (scaled - np.asarray(state["pca_mean"], dtype=np.float64)) @ np.asarray(
            state["pca_components"], dtype=np.float64
        ).T
    else:
        design = scaled
    return (
        design @ np.asarray(state["coef"], dtype=np.float64).T
        + np.asarray(state["intercept"], dtype=np.float64)
    ).ravel()


def out_of_fold_scores(
    *, config: Mapping[str, Any], representation: np.ndarray, eligibility: np.ndarray,
    gene_ids: Sequence[str], outer_folds: np.ndarray, targets: np.ndarray,
    feature_count: int, hyperparameters: Mapping[str, Any], seed: int,
) -> tuple[np.ndarray, bool]:
    scores = np.full(len(targets), np.nan, dtype=np.float64)
    uses_data = True
    for fold in sorted({int(v) for v in outer_folds}):
        fitting = np.flatnonzero(outer_folds != fold)
        held = np.flatnonzero(outer_folds == fold)
        state = fit_pipeline(
            config=config, representation=representation, eligibility=eligibility,
            gene_ids=gene_ids, fitting=fitting, targets=targets,
            feature_count=feature_count, hyperparameters=hyperparameters,
            seed=seed + 10_000 * fold,
        )
        state.pop("_design_columns", None)
        uses_data = uses_data and int(state["nonzero_coefficients"][0]) > 0
        scores[held] = apply_pipeline(state=state, representation=representation[held])
    if not np.all(np.isfinite(scores)):
        raise NasSourceFitError("source out-of-fold scores are incomplete")
    return scores, uses_data


def select_hyperparameters(
    *, model_id: str, config: Mapping[str, Any], grid: Sequence[Mapping[str, Any]],
    representation: np.ndarray, eligibility: np.ndarray, gene_ids: Sequence[str],
    outer_folds: np.ndarray, targets: np.ndarray, feature_count: int,
    seed: int, null_replicates: int, full_indices: np.ndarray,
) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    """Rank configurations by Spearman in excess of their own permutation null."""

    rows: list[dict[str, Any]] = []
    best: tuple[float, float, float] | None = None
    chosen: dict[str, Any] | None = None
    for candidate in grid:
        scores, folds_use_data = out_of_fold_scores(
            config=config, representation=representation, eligibility=eligibility,
            gene_ids=gene_ids, outer_folds=outer_folds, targets=targets,
            feature_count=feature_count, hyperparameters=candidate, seed=seed,
        )
        full_state = fit_pipeline(
            config=config, representation=representation, eligibility=eligibility,
            gene_ids=gene_ids, fitting=full_indices, targets=targets,
            feature_count=feature_count, hyperparameters=candidate, seed=seed,
        )
        full_state.pop("_design_columns", None)
        full_nonzero = int(full_state["nonzero_coefficients"][0])
        full_scores = apply_pipeline(state=full_state, representation=representation)
        full_distinct = int(len(np.unique(np.round(full_scores, 12))))
        null = spearman_null(
            targets, scores, replicates=null_replicates, seed=seed + 55_000
        )
        distinct = int(len(np.unique(np.round(scores, 12))))
        eligible = bool(
            folds_use_data and full_nonzero > 0 and full_distinct > 1 and distinct > 1
        )
        rows.append({
            "model_id": model_id,
            "alpha": float(candidate["alpha"]),
            "pca_components": int(candidate.get("pca_components", 0)),
            "source_pooled_oof_spearman": null["spearman"],
            "source_pooled_oof_null_mean": null["null_mean"],
            "source_pooled_oof_null_p95": null["null_p95"],
            "source_pooled_oof_excess_over_null_mean": null["excess_over_null_mean"],
            "source_pooled_oof_permutation_p_value": null["permutation_p_value"],
            "distinct_oof_scores": distinct,
            "distinct_full_fit_scores": full_distinct,
            "full_fit_nonzero_coefficients": full_nonzero,
            "every_outer_fold_uses_the_data": folds_use_data,
            "eligible": eligible,
        })
        if not eligible:
            continue
        key = (
            -null["excess_over_null_mean"],
            float(candidate["alpha"]),
            float(candidate.get("pca_components", 0)),
        )
        if best is None or key < best:
            best = key
            chosen = dict(candidate)
    eligible_count = int(sum(1 for row in rows if row["eligible"]))
    for row in rows:
        row["eligible_configurations_for_this_model"] = eligible_count
    return chosen, rows


def _check_hash(path: Path, expected: str, label: str) -> None:
    if sha256_file(path) != expected:
        raise NasSourceFitError(f"{label} SHA-256 differs")


def run_fit(
    *, benchmark_root: Path, contract_path: Path, contract_sha256: str, output: Path
) -> dict[str, Any]:
    from masld_bench.artifacts import verify_frozen_tree

    if output.exists():
        raise NasSourceFitError(f"refusing to overwrite source fit: {output}")
    _check_hash(contract_path, contract_sha256, "source fit contract")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if (
        contract.get("status") != "registered_source_fit_pending"
        or contract.get("firewall", {}).get("target_outcomes_read") is not False
    ):
        raise NasSourceFitError("source fit contract differs")

    activation_spec = contract["activation"]
    activation = benchmark_root / activation_spec["path"]
    _check_hash(activation / "ARTIFACTS.json", activation_spec["artifacts_sha256"], "activation")
    verify_frozen_tree(activation)
    source_spec = contract["source_inputs"]
    source = benchmark_root / source_spec["path"]
    _check_hash(source / "ARTIFACTS.json", source_spec["artifacts_sha256"], "source")
    verify_frozen_tree(source)

    _, axis_rows = read_tsv(activation / "common_stable_gene_axis.tsv")
    shared_ids = [row["stable_gene_id"] for row in axis_rows]
    if len(shared_ids) != activation_spec["shared_genes"]:
        raise NasSourceFitError("shared gene axis census differs")

    _, participants = read_tsv(source / "molecular" / "participant_axis.tsv")
    _, feature_rows = read_tsv(source / "molecular" / "rna_feature_axis.tsv")
    _, endpoint_rows = read_tsv(source / "outcomes" / "participant_endpoints.tsv")
    _, fold_rows = read_tsv(source / "folds" / "participant_outer_folds.tsv")
    ids = [row["participant_id"] for row in participants]
    if (
        len(ids) != source_spec["participants"]
        or ids != [row["participant_id"] for row in endpoint_rows]
        or ids != [row["participant_id"] for row in fold_rows]
    ):
        raise NasSourceFitError("source participant axes differ")

    endpoint_field = contract["endpoint"]["source_field"]
    low, high = contract["endpoint"]["source_scale_range"]
    targets = np.asarray([float(row[endpoint_field]) for row in endpoint_rows])
    if np.any(targets < low) or np.any(targets > high) or not np.all(np.isfinite(targets)):
        raise NasSourceFitError("a source activity score is off its recorded scale")
    census = Counter(int(v) for v in targets)
    if {str(k): v for k, v in sorted(census.items())} != contract["endpoint"]["expected_census"]:
        raise NasSourceFitError("source activity-score census differs")

    fold_field = contract["source_inputs"]["fold_field"]
    if any(not row[fold_field].isdigit() for row in fold_rows):
        raise NasSourceFitError("the chosen fold column does not cover every participant")
    outer_folds = np.asarray([int(row[fold_field]) for row in fold_rows], dtype=np.int64)

    source_ids = [row["stable_gene_id"] for row in feature_rows]
    lookup = {value: index for index, value in enumerate(source_ids)}
    if any(value not in lookup for value in shared_ids):
        raise NasSourceFitError("a shared gene is absent from the source axis")
    raw = np.load(source / "molecular" / "rna_values.npy", mmap_mode="r", allow_pickle=False)
    if raw.shape != (source_spec["participants"], len(source_ids)):
        raise NasSourceFitError("source RNA matrix shape differs")
    columns = [lookup[value] for value in shared_ids]
    log2_shared = log2_cpm_complete_axis(raw)[:, columns]
    eligibility = np.asarray(raw, dtype=np.float64)[:, columns]

    recipe = contract["recipe"]
    feature_count = int(recipe["feature_count"])
    seeds = [int(v) for v in recipe["model_seeds"]]
    selection_seed = int(recipe["hyperparameter_selection_seed"])
    null_replicates = int(recipe["source_null_replicates"])
    representations = {
        name: build_representation(name=name, log2_shared=log2_shared)
        for name in sorted({contract["models"][m]["representation"] for m in LEARNED_MODEL_IDS})
    }
    for name in representations:
        assert_row_independent(name, log2_shared)

    output.mkdir(parents=True)
    model_root = output / "models"
    model_root.mkdir()
    grid_rows: list[dict[str, Any]] = []
    seed_rows: list[dict[str, Any]] = []
    ranking: list[dict[str, Any]] = []
    unfittable: list[dict[str, Any]] = []
    chosen_by_model: dict[str, dict[str, Any]] = {}
    full_indices = np.arange(len(targets), dtype=np.int64)

    for model_id in LEARNED_MODEL_IDS:
        config = contract["models"][model_id]
        base = contract["hyperparameter_grid"]["ridge"]
        if config["reducer"] == "pca":
            grid = [{"alpha": a, "pca_components": k}
                    for a in base["alpha"] for k in base["pca_components"]]
        else:
            grid = [{"alpha": a} for a in base["alpha"]]
        chosen, rows = select_hyperparameters(
            model_id=model_id, config=config, grid=grid,
            representation=representations[config["representation"]],
            eligibility=eligibility, gene_ids=shared_ids, outer_folds=outer_folds,
            targets=targets, feature_count=feature_count, seed=selection_seed,
            null_replicates=null_replicates, full_indices=full_indices,
        )
        grid_rows.extend(rows)
        if chosen is None:
            unfittable.append({
                "model_id": model_id,
                "grid_configurations": len(rows),
                "eligible_grid_configurations": 0,
                "reason": "no registered configuration kept a non-zero coefficient in every outer fold and the full fit",
            })
            continue
        chosen_by_model[model_id] = chosen
        for seed in seeds:
            oof, folds_use_data = out_of_fold_scores(
                config=config, representation=representations[config["representation"]],
                eligibility=eligibility, gene_ids=shared_ids, outer_folds=outer_folds,
                targets=targets, feature_count=feature_count,
                hyperparameters=chosen, seed=seed,
            )
            state = fit_pipeline(
                config=config, representation=representations[config["representation"]],
                eligibility=eligibility, gene_ids=shared_ids, fitting=full_indices,
                targets=targets, feature_count=feature_count,
                hyperparameters=chosen, seed=seed,
            )
            design_columns = int(state.pop("_design_columns"))
            null = spearman_null(targets, oof, replicates=null_replicates, seed=seed + 55_000)
            if int(state["nonzero_coefficients"][0]) == 0 or not folds_use_data:
                raise NasSourceFitError(f"{model_id} is degenerate at seed {seed}")
            current = model_root / model_id / f"seed_{seed}"
            current.mkdir(parents=True)
            np.savez_compressed(
                current / "model_state.npz",
                **{k: np.asarray(v) for k, v in state.items()},
            )
            state_sha = sha256_file(current / "model_state.npz")
            with (current / "receipt.json").open("x", encoding="utf-8") as handle:
                handle.write(json.dumps({
                    "schema_version": "masld-bench-gse135251-nas-source-seed-v1",
                    "model_id": model_id, "model_seed": int(seed),
                    "representation": config["representation"], "reducer": config["reducer"],
                    "hyperparameters": chosen, "design_columns": design_columns,
                    "source_pooled_oof_spearman": null["spearman"],
                    "source_pooled_oof_null_mean": null["null_mean"],
                    "source_pooled_oof_excess_over_null_mean": null["excess_over_null_mean"],
                    "source_pooled_oof_permutation_p_value": null["permutation_p_value"],
                    "model_state_sha256": state_sha,
                    "target_outcomes_read": False,
                    "target_expression_values_read": False,
                }, indent=2, sort_keys=True) + "\n")
            seed_rows.append({
                "model_id": model_id, "model_seed": seed,
                "representation": config["representation"],
                "model_state_path": f"models/{model_id}/seed_{seed}/model_state.npz",
                "model_state_sha256": state_sha,
                "source_pooled_oof_spearman": null["spearman"],
                "source_pooled_oof_null_mean": null["null_mean"],
                "source_pooled_oof_excess_over_null_mean": null["excess_over_null_mean"],
                "source_pooled_oof_permutation_p_value": null["permutation_p_value"],
                "nonzero_coefficients": int(state["nonzero_coefficients"][0]),
            })
        model_rows = [r for r in seed_rows if r["model_id"] == model_id]
        ranking.append({
            "model_id": model_id,
            "mean_source_pooled_oof_spearman": float(np.mean(
                [r["source_pooled_oof_spearman"] for r in model_rows])),
            "mean_source_pooled_oof_null_mean": float(np.mean(
                [r["source_pooled_oof_null_mean"] for r in model_rows])),
            "mean_source_pooled_oof_excess_over_null_mean": float(np.mean(
                [r["source_pooled_oof_excess_over_null_mean"] for r in model_rows])),
            "min_source_pooled_oof_permutation_p_value": float(np.min(
                [r["source_pooled_oof_permutation_p_value"] for r in model_rows])),
            "selected_alpha": float(chosen["alpha"]),
            "selected_pca_components": int(chosen.get("pca_components", 0)),
            "eligible_grid_configurations": int(rows[0]["eligible_configurations_for_this_model"]),
        })

    if not ranking:
        raise NasSourceFitError("every learned model is unfittable on this source")
    ranking.sort(key=lambda r: (-r["mean_source_pooled_oof_excess_over_null_mean"], r["model_id"]))
    selected_source_model_id = ranking[0]["model_id"]

    write_tsv(output / "hyperparameter_selection.tsv", tuple(grid_rows[0]), grid_rows)
    write_tsv(output / "seed_model_index.tsv", tuple(seed_rows[0]), seed_rows)
    write_tsv(output / "source_model_ranking.tsv", tuple(ranking[0]), ranking)

    receipt = {
        "schema_version": "masld-bench-gse135251-nas-sum-source-fit-v1",
        "status": "passed_source_nas_sum_fit_preprocessing_locked",
        "source_series": "GSE135251", "target_series": "GSE267145",
        "endpoint": "continuous_nafld_activity_score_0_to_8",
        "endpoint_binarised": False,
        "source_participants": int(len(targets)),
        "source_activity_mean": float(np.mean(targets)),
        "source_activity_sd": float(np.std(targets, ddof=1)),
        "source_activity_census": {str(k): v for k, v in sorted(census.items())},
        "shared_genes": len(shared_ids),
        "selected_features": feature_count,
        "model_seeds": seeds,
        "seeds_are_biological_replicates": False,
        "fold_field": fold_field,
        "hyperparameter_selection": "source_pooled_oof_spearman_in_excess_of_its_own_permutation_null",
        "source_null_replicates": null_replicates,
        "selected_hyperparameters": chosen_by_model,
        "source_model_ranking": ranking,
        "selected_source_model_id": selected_source_model_id,
        "fitted_model_ids": sorted(chosen_by_model),
        "unfittable_model_ids": sorted(r["model_id"] for r in unfittable),
        "unfittable_models": unfittable,
        "no_source_model_exceeds_its_own_null": not any(
            r["mean_source_pooled_oof_excess_over_null_mean"] > 0.0 for r in ranking),
        "training_mean_nas_sum": float(np.mean(targets)),
        "per_sample_transform_is_row_independent": True,
        "cross_sample_pooling_performed": False,
        "target_expression_values_read": False,
        "target_outcomes_read": False,
        "project_sealed": False,
        "champion_claim_allowed": False,
    }
    with (output / "fit_receipt.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-root", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--contract-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    a = parser.parse_args()
    receipt = run_fit(benchmark_root=a.benchmark_root, contract_path=a.contract,
                      contract_sha256=a.contract_sha256, output=a.output)
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
