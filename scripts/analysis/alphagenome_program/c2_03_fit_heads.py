#!/usr/bin/env python3
"""C2/C3 step 3: delta ridge and allele-identity heads on the Currin caQTL label.

The head recipe is the shipped reporter recipe
(`Analysis/MASLD_Model_Benchmark/scripts/fit_gse281364_dna_language_seeded_heads.py` and the
helpers it imports): fold-safe whitened SVD projection fitted on training alleles only, head
features [REF, ALT, ALT-REF, |ALT-REF|] with reverse-complement forms averaged into each allele,
delta slice [512:768], absolute-training-Pearson top-k selection, ridge alpha chosen by
inner-validation RMSE against a constant-training-mean floor, outer-training 1-Mb block bootstrap
per seed.  Only the fold definition (ChromBPNet chromosome groups) and the label (Currin beta_alt)
differ from the reporter run.

The allele-identity ridge control differs from the delta head in one component: its features are the
16 one-hot (REF, ALT) indicators instead of the projected delta block.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.linear_model import Ridge

STD_FLOOR = 1.0e-6
EIGENVALUE_ABSOLUTE_FLOOR = 1.0e-8
EIGENVALUE_RELATIVE_FLOOR = 1.0e-6
PROJECTION_WIDTH = 256
INNER_SALT = 10_000
OUTER_SALT = 20_000
SEEDS = (1103, 2909, 4721, 6673, 8111)
ALPHAS = (0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0)
DELTA_TOP_K = (32, 128, 256)
FULL_TOP_K = (32, 128, 512, 1024)
ALLELE_TOP_K = (16,)
DELTA_SLICE = slice(512, 768)
FULL_SLICE = slice(0, 1024)
CURVE_SIZES = (800, 2000, 5000, 12500, 0)  # 0 means the whole training pool
BASES = "ACGT"


class HeadFitError(RuntimeError):
    """Raised when a fold, feature or projection boundary differs."""


# ---------------------------------------------------------------- projection


def fit_projection(embeddings: np.ndarray, training_mask: np.ndarray) -> dict[str, np.ndarray]:
    values = np.asarray(embeddings, dtype=np.float64)
    mask = np.asarray(training_mask, dtype=bool)
    if values.ndim != 3 or values.shape[1] != 4 or int(mask.sum()) < 2:
        raise HeadFitError("projection fit request differs")
    training = values[mask].reshape(-1, values.shape[2])
    mean = training.mean(axis=0)
    standard_deviation = np.maximum(training.std(axis=0), STD_FLOOR)
    normalized = (training - mean) / standard_deviation
    width = min(PROJECTION_WIDTH, values.shape[2])
    _, singular_values, right_vectors = np.linalg.svd(normalized, full_matrices=False)
    components = right_vectors[:width].copy()
    for index in range(width):
        pivot = int(np.argmax(np.abs(components[index])))
        if components[index, pivot] < 0:
            components[index] *= -1
    eigenvalues = singular_values[:width] ** 2 / max(1, training.shape[0] - 1)
    floor = max(EIGENVALUE_ABSOLUTE_FLOOR, EIGENVALUE_RELATIVE_FLOOR * float(eigenvalues[0]))
    return {
        "mean": mean,
        "standard_deviation": standard_deviation,
        "components": components,
        "whitening_scale": np.sqrt(np.maximum(eigenvalues, floor)),
    }


def apply_projection(embeddings: np.ndarray, parameters: Mapping[str, np.ndarray]) -> np.ndarray:
    values = np.asarray(embeddings, dtype=np.float64)
    normalized = (values - parameters["mean"]) / parameters["standard_deviation"]
    projected = normalized @ parameters["components"].T
    projected /= parameters["whitening_scale"]
    if not np.isfinite(projected).all():
        raise HeadFitError("projected embeddings are not finite")
    return projected.astype(np.float32)


def head_features(projected: np.ndarray) -> np.ndarray:
    values = np.asarray(projected, dtype=np.float32)
    reference = (values[:, 0] + values[:, 2]) / 2.0
    alternative = (values[:, 1] + values[:, 3]) / 2.0
    difference = alternative - reference
    return np.concatenate((reference, alternative, difference, np.abs(difference)), axis=1)


# ---------------------------------------------------------------- head fit


def bootstrap_indices(eligible: np.ndarray, block_ids: np.ndarray, *, seed: int, salt: int) -> np.ndarray:
    indices = np.flatnonzero(eligible)
    labels = block_ids[indices]
    order = np.argsort(labels, kind="stable")
    sorted_labels, sorted_indices = labels[order], indices[order]
    unique, starts = np.unique(sorted_labels, return_index=True)
    if unique.size < 2:
        raise HeadFitError("training block bootstrap has fewer than two blocks")
    ends = np.append(starts[1:], sorted_labels.size)
    groups = [sorted_indices[start:end] for start, end in zip(starts, ends)]
    rng = np.random.default_rng(np.random.SeedSequence([seed, salt]))
    sampled = rng.integers(0, len(groups), size=len(groups))
    return np.concatenate([groups[index] for index in sampled])


def training_transform(
    features: np.ndarray, outcomes: np.ndarray, training_indices: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    training = features[training_indices]
    mean = training.mean(axis=0)
    scale = training.std(axis=0)
    usable = np.flatnonzero(scale > 1.0e-8)
    if usable.size == 0:
        return mean, np.maximum(scale, 1.0e-8), usable
    normalized = (training[:, usable] - mean[usable]) / scale[usable]
    centered = outcomes[training_indices] - outcomes[training_indices].mean()
    denominator = np.sqrt(np.sum(normalized**2, axis=0) * np.sum(centered**2))
    correlation = np.divide(
        np.abs(normalized.T @ centered),
        denominator,
        out=np.zeros(usable.size, dtype=np.float64),
        where=denominator > 0,
    )
    order = np.lexsort((usable, -correlation))
    return mean, np.maximum(scale, 1.0e-8), usable[order]


def fit_fold(
    *,
    inner_features: np.ndarray,
    outer_features: np.ndarray,
    outcomes: np.ndarray,
    folds: np.ndarray,
    block_ids: np.ndarray,
    held_fold: int,
    seed: int,
    top_k_grid: Sequence[int],
    train_pool: np.ndarray,
) -> tuple[np.ndarray, dict[str, object]]:
    inner_fold = (held_fold + 1) % 5
    outer_test = folds == held_fold
    outer_train = (~outer_test) & train_pool
    inner_valid = (folds == inner_fold) & train_pool
    inner_train = outer_train & ~inner_valid
    if set(block_ids[outer_test]) & set(block_ids[outer_train]):
        raise HeadFitError("held 1-Mb block overlaps training")
    inner_boot = bootstrap_indices(inner_train, block_ids, seed=seed, salt=INNER_SALT + held_fold)
    mean, scale, rank = training_transform(inner_features, outcomes, inner_boot)
    inner_train_matrix = inner_features[inner_boot]
    inner_valid_matrix = inner_features[inner_valid]
    inner_train_y = outcomes[inner_boot]
    inner_valid_y = outcomes[inner_valid]
    choices: list[tuple[float, int, float]] = []
    for k in sorted({min(int(value), int(rank.size)) for value in top_k_grid if int(value) > 0}):
        if k == 0:
            continue
        selected = rank[:k]
        x_train = (inner_train_matrix[:, selected] - mean[selected]) / scale[selected]
        x_valid = (inner_valid_matrix[:, selected] - mean[selected]) / scale[selected]
        for alpha in ALPHAS:
            model = Ridge(alpha=float(alpha), fit_intercept=True, solver="lsqr", tol=1.0e-7)
            model.fit(x_train, inner_train_y)
            prediction = model.predict(x_valid)
            choices.append(
                (float(np.sqrt(np.mean((inner_valid_y - prediction) ** 2))), k, float(alpha))
            )
    constant_rmse = float(np.sqrt(np.mean((inner_valid_y - inner_train_y.mean()) ** 2)))
    best_rmse, best_k, best_alpha = min(
        [*choices, (constant_rmse, 0, math.inf)], key=lambda value: (value[0], value[1], value[2])
    )
    outer_boot = bootstrap_indices(outer_train, block_ids, seed=seed, salt=OUTER_SALT + held_fold)
    outer_mean, outer_scale, outer_rank = training_transform(outer_features, outcomes, outer_boot)
    if best_k == 0 or outer_rank.size == 0:
        intercept = float(outcomes[outer_boot].mean())
        prediction = np.full(int(outer_test.sum()), intercept)
        selected_alpha: float | str = "constant_training_mean"
        selected_size = 0
    else:
        selected = outer_rank[: min(best_k, int(outer_rank.size))]
        x_train = (outer_features[outer_boot][:, selected] - outer_mean[selected]) / outer_scale[selected]
        x_test = (outer_features[outer_test][:, selected] - outer_mean[selected]) / outer_scale[selected]
        model = Ridge(alpha=best_alpha, fit_intercept=True, solver="lsqr", tol=1.0e-7)
        model.fit(x_train, outcomes[outer_boot])
        prediction = np.asarray(model.predict(x_test), dtype=np.float64)
        selected_alpha = best_alpha
        selected_size = int(selected.size)
    receipt = {
        "held_out_fold": held_fold,
        "inner_validation_fold": inner_fold,
        "seed": seed,
        "inner_training_variants": int(inner_train.sum()),
        "inner_validation_variants": int(inner_valid.sum()),
        "outer_training_variants": int(outer_train.sum()),
        "outer_test_variants": int(outer_test.sum()),
        "outer_training_blocks": len(set(block_ids[outer_train].tolist())),
        "outer_test_blocks": len(set(block_ids[outer_test].tolist())),
        "selected_feature_count": selected_size,
        "selected_alpha": selected_alpha,
        "selected_inner_rmse": best_rmse,
        "constant_inner_rmse": constant_rmse,
        "held_outcomes_used_for_fit_or_tuning": False,
        "held_features_used_for_preprocessing_selection_or_tuning": False,
        "held_blocks_overlap_training": False,
        "posthoc_calibration": False,
    }
    return prediction, receipt


# ---------------------------------------------------------------- metrics


def fast_spearman(observed: np.ndarray, predicted: np.ndarray) -> float:
    if observed.size < 3 or np.all(predicted == predicted[0]) or np.all(observed == observed[0]):
        return float("nan")
    left = rankdata(observed)
    right = rankdata(predicted)
    left = left - left.mean()
    right = right - right.mean()
    denominator = float(np.sqrt(np.sum(left**2) * np.sum(right**2)))
    return float("nan") if denominator == 0 else float(np.sum(left * right) / denominator)


def fisher_macro(values: Sequence[float]) -> float:
    """Reporter-native macro: Fisher z-average, undefined if any component is undefined."""
    array = np.asarray(values, dtype=np.float64)
    if not np.isfinite(array).all():
        return float("nan")
    clipped = np.clip(array, -0.999999, 0.999999)
    return float(np.tanh(np.mean(np.arctanh(clipped))))


def macro_over_folds(observed: np.ndarray, predicted: np.ndarray, folds: np.ndarray) -> float:
    return fisher_macro([fast_spearman(observed[folds == f], predicted[folds == f]) for f in range(5)])


# ---------------------------------------------------------------- campaign


def build_fold_features(
    embeddings: np.ndarray,
    folds: np.ndarray,
    held_fold: int,
    train_pool: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    inner_fold = (held_fold + 1) % 5
    outer_train = (folds != held_fold) & train_pool
    inner_train = outer_train & (folds != inner_fold)
    inner = head_features(apply_projection(embeddings, fit_projection(embeddings, inner_train))).astype(np.float64)
    outer = head_features(apply_projection(embeddings, fit_projection(embeddings, outer_train))).astype(np.float64)
    return inner, outer


def allele_identity_features(ref: np.ndarray, alt: np.ndarray) -> np.ndarray:
    index = {base: position for position, base in enumerate(BASES)}
    features = np.zeros((len(ref), 16), dtype=np.float64)
    for row, (r, a) in enumerate(zip(ref, alt)):
        features[row, 4 * index[r] + index[a]] = 1.0
    return features


def run_arm(
    *,
    inner_features_by_fold: Mapping[int, np.ndarray],
    outer_features_by_fold: Mapping[int, np.ndarray],
    outcomes: np.ndarray,
    folds: np.ndarray,
    blocks: np.ndarray,
    feature_slice: slice,
    top_k_grid: Sequence[int],
    train_pool: np.ndarray,
    audit: list[dict[str, object]] | None = None,
    arm_name: str = "",
) -> dict[int, np.ndarray]:
    predictions: dict[int, np.ndarray] = {}
    for seed in SEEDS:
        oof = np.full(len(outcomes), np.nan)
        for held_fold in range(5):
            test = folds == held_fold
            prediction, receipt = fit_fold(
                inner_features=inner_features_by_fold[held_fold][:, feature_slice],
                outer_features=outer_features_by_fold[held_fold][:, feature_slice],
                outcomes=outcomes,
                folds=folds,
                block_ids=blocks,
                held_fold=held_fold,
                seed=seed,
                top_k_grid=top_k_grid,
                train_pool=train_pool,
            )
            oof[test] = prediction
            if audit is not None:
                audit.append({"arm": arm_name, **receipt})
        if not np.isfinite(oof).all():
            raise HeadFitError("out-of-fold prediction coverage differs")
        predictions[seed] = oof
    return predictions


def block_bootstrap(
    observed: np.ndarray,
    arms: Mapping[str, np.ndarray],
    folds: np.ndarray,
    blocks: np.ndarray,
    *,
    resamples: int,
    seed: int,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    by_fold_blocks = {}
    for fold in range(5):
        rows = np.flatnonzero(folds == fold)
        unique = sorted(set(blocks[rows].tolist()))
        by_fold_blocks[fold] = (unique, {block: rows[blocks[rows] == block] for block in unique})
    rng = np.random.default_rng(seed)
    macro = {name: np.full(resamples, np.nan) for name in arms}
    pooled = {name: np.full(resamples, np.nan) for name in arms}
    names = list(arms)
    for iteration in range(resamples):
        per_fold_indices = []
        for fold in range(5):
            unique, lookup = by_fold_blocks[fold]
            sampled = rng.integers(0, len(unique), size=len(unique))
            per_fold_indices.append(np.concatenate([lookup[unique[index]] for index in sampled]))
        all_indices = np.concatenate(per_fold_indices)
        fold_observed = [observed[indices] for indices in per_fold_indices]
        pooled_observed = observed[all_indices]
        for name in names:
            values = arms[name]
            macro[name][iteration] = fisher_macro(
                [fast_spearman(fold_observed[fold], values[per_fold_indices[fold]]) for fold in range(5)]
            )
            pooled[name][iteration] = fast_spearman(pooled_observed, values[all_indices])
    return macro, pooled


def interval(samples: np.ndarray) -> tuple[float, float, float]:
    valid = samples[np.isfinite(samples)]
    if valid.size < int(0.95 * samples.size):
        return float("nan"), float("nan"), float("nan")
    low, high = (float(value) for value in np.quantile(valid, (0.025, 0.975)))
    return low, high, float(np.std(valid, ddof=1))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resamples", type=int, default=10000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260914)
    parser.add_argument("--skip-curve", action="store_true")
    arguments = parser.parse_args()
    out = arguments.output
    (out / "tables").mkdir(parents=True, exist_ok=True)

    labels = pd.read_csv(arguments.inputs / "currin_lead_labels.tsv.gz", sep="\t")
    labels = labels.sort_values("row_index").reset_index(drop=True)
    y = labels["beta_alt"].to_numpy(np.float64)
    folds = labels["heldout_fold"].to_numpy(np.int64)
    blocks = labels["block_1mb"].to_numpy(dtype=str)
    ambiguous = labels["strand_ambiguous"].to_numpy(bool)
    n = len(labels)

    raw: dict[str, np.ndarray] = {}
    for model in ("hyenadna", "caduceus"):
        path = arguments.raw / f"{model}_allele_embeddings.npz"
        with np.load(path, allow_pickle=False) as data:
            order = tuple(data["allele_order"].astype(str).tolist())
            if order != ("REF", "ALT", "REF_RC", "ALT_RC") or not np.array_equal(
                data["row_index"], labels["row_index"].to_numpy(np.int64)
            ):
                raise HeadFitError(f"{model} embedding alignment differs")
            raw[model] = data["embeddings"]
        if raw[model].shape[:2] != (n, 4) or not np.isfinite(raw[model]).all():
            raise HeadFitError(f"{model} embedding geometry differs")

    allele_features = allele_identity_features(
        labels["ref"].to_numpy(dtype=str), labels["alt"].to_numpy(dtype=str)
    )

    arms_config = [
        ("hyenadna", "delta_ridge", DELTA_SLICE, DELTA_TOP_K),
        ("hyenadna", "full_ridge", FULL_SLICE, FULL_TOP_K),
        ("caduceus", "delta_ridge", DELTA_SLICE, DELTA_TOP_K),
        ("caduceus", "full_ridge", FULL_SLICE, FULL_TOP_K),
    ]

    def campaign(pool: np.ndarray, tag: str) -> tuple[dict[str, dict[int, np.ndarray]], list[dict]]:
        audit: list[dict[str, object]] = []
        per_seed: dict[str, dict[int, np.ndarray]] = {}
        for model in ("hyenadna", "caduceus"):
            inner_by_fold, outer_by_fold = {}, {}
            for held_fold in range(5):
                inner_by_fold[held_fold], outer_by_fold[held_fold] = build_fold_features(
                    raw[model], folds, held_fold, pool
                )
            for model_id, head, feature_slice, top_k in arms_config:
                if model_id != model:
                    continue
                name = f"{model}/{head}"
                print(f"[{tag}] fitting {name}", flush=True)
                per_seed[name] = run_arm(
                    inner_features_by_fold=inner_by_fold,
                    outer_features_by_fold=outer_by_fold,
                    outcomes=y,
                    folds=folds,
                    blocks=blocks,
                    feature_slice=feature_slice,
                    top_k_grid=top_k,
                    train_pool=pool,
                    audit=audit,
                    arm_name=name,
                )
            del inner_by_fold, outer_by_fold
        identity_by_fold = {fold: allele_features for fold in range(5)}
        print(f"[{tag}] fitting allele_identity_ridge control", flush=True)
        per_seed["available_simple_controls/allele_identity_ridge"] = run_arm(
            inner_features_by_fold=identity_by_fold,
            outer_features_by_fold=identity_by_fold,
            outcomes=y,
            folds=folds,
            blocks=blocks,
            feature_slice=slice(0, 16),
            top_k_grid=ALLELE_TOP_K,
            train_pool=pool,
            audit=audit,
            arm_name="available_simple_controls/allele_identity_ridge",
        )
        return per_seed, audit

    summaries: list[dict[str, object]] = []
    per_seed_rows: list[dict[str, object]] = []
    per_fold_rows: list[dict[str, object]] = []

    for tag, pool in (("primary_all_variants", np.ones(n, bool)), ("no_strand_ambiguous", ~ambiguous)):
        if tag == "no_strand_ambiguous" and int(ambiguous.sum()) == 0:
            continue
        per_seed, audit = campaign(pool, tag)
        evaluated = pool if tag == "no_strand_ambiguous" else np.ones(n, bool)
        rows = np.flatnonzero(evaluated)
        ensemble = {name: np.mean(np.vstack([values[s] for s in SEEDS]), axis=0) for name, values in per_seed.items()}
        control = "available_simple_controls/allele_identity_ridge"
        for name, values in per_seed.items():
            for seed in SEEDS:
                per_seed_rows.append(
                    {
                        "stratum": tag,
                        "arm": name,
                        "seed": seed,
                        "macro_spearman": macro_over_folds(y[rows], values[seed][rows], folds[rows]),
                        "pooled_spearman": fast_spearman(y[rows], values[seed][rows]),
                    }
                )
            for fold in range(5):
                keep = rows[folds[rows] == fold]
                per_fold_rows.append(
                    {
                        "stratum": tag,
                        "arm": name,
                        "held_out_fold": fold,
                        "variants": int(keep.size),
                        "blocks": int(len(set(blocks[keep].tolist()))),
                        "spearman_seed_ensemble": fast_spearman(y[keep], ensemble[name][keep]),
                    }
                )
        macro_samples, pooled_samples = block_bootstrap(
            y[rows],
            {name: values[rows] for name, values in ensemble.items()},
            folds[rows],
            blocks[rows],
            resamples=arguments.resamples,
            seed=arguments.bootstrap_seed,
        )
        for name in ensemble:
            macro_point = macro_over_folds(y[rows], ensemble[name][rows], folds[rows])
            pooled_point = fast_spearman(y[rows], ensemble[name][rows])
            low, high, se = interval(macro_samples[name])
            plow, phigh, pse = interval(pooled_samples[name])
            gain_samples = macro_samples[name] - macro_samples[control]
            gain_low, gain_high, gain_se = interval(gain_samples)
            pooled_gain_samples = pooled_samples[name] - pooled_samples[control]
            pgain_low, pgain_high, _ = interval(pooled_gain_samples)
            control_macro = macro_over_folds(y[rows], ensemble[control][rows], folds[rows])
            control_pooled = fast_spearman(y[rows], ensemble[control][rows])
            comparable = 0
            positive = 0
            for s in SEEDS:
                arm_value = macro_over_folds(y[rows], per_seed[name][s][rows], folds[rows])
                control_value = macro_over_folds(y[rows], per_seed[control][s][rows], folds[rows])
                if not (np.isfinite(arm_value) and np.isfinite(control_value)):
                    continue
                comparable += 1
                positive += int(arm_value > control_value)
            summaries.append(
                {
                    "stratum": tag,
                    "arm": name,
                    "variants": int(rows.size),
                    "blocks": int(len(set(blocks[rows].tolist()))),
                    "macro_spearman": macro_point,
                    "macro_ci_low": low,
                    "macro_ci_high": high,
                    "macro_bootstrap_se": se,
                    "gain_vs_allele_identity": macro_point - control_macro,
                    "gain_ci_low": gain_low,
                    "gain_ci_high": gain_high,
                    "gain_bootstrap_se": gain_se,
                    "positive_gain_seeds": int(positive),
                    "comparable_seeds": int(comparable),
                    "pooled_spearman": pooled_point,
                    "pooled_ci_low": plow,
                    "pooled_ci_high": phigh,
                    "pooled_gain_vs_allele_identity": pooled_point - control_pooled,
                    "pooled_gain_ci_low": pgain_low,
                    "pooled_gain_ci_high": pgain_high,
                    "bootstrap_resamples": arguments.resamples,
                    "bootstrap_seed": arguments.bootstrap_seed,
                    "bootstrap_unit": "1Mb_block_within_fold",
                }
            )
        if tag == "primary_all_variants":
            frame = pd.DataFrame({"row_index": labels["row_index"], "variant_id": labels["lead_variant_id"],
                                  "chr": labels["chr"], "pos_hg38": labels["pos_hg38"],
                                  "block_1mb": blocks, "heldout_fold": folds, "beta_alt": y,
                                  "strand_ambiguous": ambiguous})
            for name, values in ensemble.items():
                frame[name.replace("/", "__") + "__ensemble"] = values
            for name, values in per_seed.items():
                for seed in SEEDS:
                    frame[f"{name.replace('/', '__')}__seed{seed}"] = values[seed]
            frame.to_csv(out / "tables/oof_predictions.tsv.gz", sep="\t", index=False, compression="gzip")
            pd.DataFrame(audit).to_csv(out / "tables/selection_audit.tsv", sep="\t", index=False)
            np.savez_compressed(
                out / "tables/bootstrap_samples.npz",
                **{name.replace("/", "__"): macro_samples[name] for name in macro_samples},
            )

    pd.DataFrame(summaries).to_csv(out / "tables/head_summary.tsv", sep="\t", index=False)
    pd.DataFrame(per_seed_rows).to_csv(out / "tables/per_seed_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(per_fold_rows).to_csv(out / "tables/per_fold_metrics.tsv", sep="\t", index=False)

    # ---------------------------------------------------------- C3 learning curve
    curve_rows: list[dict[str, object]] = []
    if not arguments.skip_curve:
        full_pool = np.ones(n, bool)
        for model in ("hyenadna", "caduceus"):
            # The full-size arm uses the whole training pool, which does not depend on the seed, so
            # its two projections are fitted once per fold and reused across seeds.
            full_cache: dict[int, tuple[np.ndarray, np.ndarray]] = {}
            for size in CURVE_SIZES:
                for seed in SEEDS:
                    oof = np.full(n, np.nan)
                    achieved: list[int] = []
                    for held_fold in range(5):
                        test = folds == held_fold
                        candidates = np.flatnonzero(~test)
                        rng = np.random.default_rng(np.random.SeedSequence([seed, 70_000 + held_fold, size]))
                        permuted = rng.permutation(candidates)
                        take = len(candidates) if size == 0 else min(size, len(candidates))
                        pool = np.zeros(n, bool)
                        pool[permuted[:take]] = True
                        pool[test] = True  # test rows are never used for fitting; mask only gates training
                        achieved.append(int(take))
                        if size == 0:
                            if held_fold not in full_cache:
                                full_cache[held_fold] = build_fold_features(raw[model], folds, held_fold, pool & ~test)
                            inner, outer = full_cache[held_fold]
                        else:
                            inner, outer = build_fold_features(raw[model], folds, held_fold, pool & ~test)
                        prediction, _ = fit_fold(
                            inner_features=inner[:, DELTA_SLICE],
                            outer_features=outer[:, DELTA_SLICE],
                            outcomes=y,
                            folds=folds,
                            block_ids=blocks,
                            held_fold=held_fold,
                            seed=seed,
                            top_k_grid=DELTA_TOP_K,
                            train_pool=pool,
                        )
                        oof[test] = prediction
                    curve_rows.append(
                        {
                            "model": model,
                            "head": "delta_ridge",
                            "requested_training_size": size if size else -1,
                            "mean_training_variants_per_fold": float(np.mean(achieved)),
                            "seed": seed,
                            "macro_spearman": macro_over_folds(y, oof, folds),
                            "pooled_spearman": fast_spearman(y, oof),
                        }
                    )
                    print(f"[curve] {model} size={size} seed={seed} macro={curve_rows[-1]['macro_spearman']:.4f}", flush=True)
            full_cache.clear()
        pd.DataFrame(curve_rows).to_csv(out / "tables/learning_curve.tsv", sep="\t", index=False)

    metrics = {
        "schema_version": "agp-c2-endogenous-head-metrics-v1",
        "variants": int(n),
        "blocks": int(len(set(blocks.tolist()))),
        "strand_ambiguous": int(ambiguous.sum()),
        "seeds": list(SEEDS),
        "bootstrap_resamples": arguments.resamples,
        "bootstrap_seed": arguments.bootstrap_seed,
        "bootstrap_unit": "1Mb_block_resampled_within_fold_paired_across_arms",
        "primary_metric": "signed Spearman within held-out fold, Fisher z-macro over 5 folds, seed-ensemble prediction",
        "secondary_metric": "pooled out-of-fold signed Spearman",
        "reporter_reference_gain": 0.082007807505673547,
        "curve_sizes": list(CURVE_SIZES),
    }
    (out / "metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    print(json.dumps(metrics, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
