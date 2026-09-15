#!/usr/bin/env python3
"""E2 step 2: does reporter supervision improve endogenous prediction?

Two matched arms on frozen HyenaDNA REF/ALT embeddings of the same 4,096-bp windows:

* arm A `endogenous_only`  - delta-ridge head trained on the Currin caQTL label alone.
* arm B `plus_reporter_oof` - identical, plus one force-included column: the out-of-fold prediction
  of a reporter-trained delta-ridge head fitted on the GSE281364 MPRA allele effect.

Leakage rule.  Held-out folds are whole chromosome groups (ChromBPNet `hepatocyte_5fold_v2`), so no
1-Mb block straddles a fold boundary.  For evaluation fold `f` the endogenous head sees only rows
with fold != f, and the reporter head that supplies the extra column for fold `f` is fitted only on
MPRA elements with fold != f.  The extra column for arm B's own training rows is produced by nested
out-of-fold prediction: for inner fold `g != f` the reporter head is fitted on elements in folds not
in {f, g}.  No row ever receives a reporter prediction from a head that saw it.

Head recipe is the shipped reporter recipe (`fit_gse281364_dna_language_seeded_heads.py` and the
helpers it imports), reused here through `c2_03_fit_heads.py`: fold-safe whitened SVD projection on
training rows only, head features [REF, ALT, ALT-REF, |ALT-REF|] with reverse-complement forms
averaged into each allele, delta slice [512:768], absolute-training-Pearson top-k selection, ridge
alpha by inner-validation RMSE against a constant-training-mean floor, 1-Mb block bootstrap of the
training rows per seed.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge

sys.path.insert(0, str(Path(__file__).resolve().parent))
from c2_03_fit_heads import (  # noqa: E402  (module lives beside this file)
    ALPHAS,
    DELTA_TOP_K,
    allele_identity_features,
    apply_projection,
    bootstrap_indices,
    fast_spearman,
    fisher_macro,
    fit_projection,
    head_features,
    interval,
)

SEEDS = (1103, 2909, 4721, 6673, 8111)
DELTA_SLICE = slice(512, 768)
ALLELE_TOP_K = (16,)
INNER_SALT = 10_000
OUTER_SALT = 20_000
REPORTER_INNER_SALT = 30_000
REPORTER_OUTER_SALT = 40_000
STRATA = (("500", 500.0), ("1000", 1000.0), ("5000", 5000.0), ("10000", 10000.0),
          ("100000", 100000.0), ("unlimited", np.inf))
PRIMARY_STRATUM = "5000"
C1_FAMILY_P = {"c1_endpoint1_contrast": 1.0e-4, "c1_endpoint2_contrast": 0.0102}


class TransferError(RuntimeError):
    """Raised when a fold, leakage or feature contract differs."""


# ---------------------------------------------------------------- head machinery


def rank_columns(features: np.ndarray, outcomes: np.ndarray, rows: np.ndarray,
                 candidate: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Training mean, scale and |Pearson|-ordered candidate columns, fitted on `rows` only."""
    training = features[rows]
    mean = training.mean(axis=0)
    scale = np.maximum(training.std(axis=0), 1.0e-8)
    usable = candidate[training[:, candidate].std(axis=0) > 1.0e-8]
    if usable.size == 0:
        return mean, scale, usable
    normalized = (training[:, usable] - mean[usable]) / scale[usable]
    centered = outcomes[rows] - outcomes[rows].mean()
    denominator = np.sqrt(np.sum(normalized**2, axis=0) * np.sum(centered**2))
    correlation = np.divide(
        np.abs(normalized.T @ centered), denominator,
        out=np.zeros(usable.size, dtype=np.float64), where=denominator > 0,
    )
    order = np.lexsort((usable, -correlation))
    return mean, scale, usable[order]


def fit_head(
    *,
    inner_features: np.ndarray,
    outer_features: np.ndarray,
    outcomes: np.ndarray,
    inner_train: np.ndarray,
    inner_valid: np.ndarray,
    outer_train: np.ndarray,
    apply_rows: np.ndarray,
    candidate: np.ndarray,
    forced: np.ndarray,
    top_k_grid: Sequence[int],
    locked: tuple[int, float] | None = None,
    constant_floor: bool = True,
) -> tuple[np.ndarray, dict[str, object]]:
    """Select (k, alpha) on the inner split, refit on the outer split, predict `apply_rows`.

    `constant_floor=True` is the shipped reporter recipe: the constant training mean competes in the
    inner selection and wins when no feature set beats the intercept, which makes the fold's
    prediction constant and its Spearman undefined.  `constant_floor=False` removes that candidate,
    identically in every arm, so the head always returns a fitted ridge.
    """
    mean, scale, rank = rank_columns(inner_features, outcomes, inner_train, candidate)
    inner_valid_y = outcomes[inner_valid]
    inner_train_y = outcomes[inner_train]
    constant_rmse = float(np.sqrt(np.mean((inner_valid_y - inner_train_y.mean()) ** 2)))
    if locked is not None:
        best_k, best_alpha = locked
        best_rmse = float("nan")
    else:
        choices: list[tuple[float, int, float]] = []
        for k in sorted({min(int(v), int(rank.size)) for v in top_k_grid if int(v) > 0}):
            if k == 0:
                continue
            selected = np.concatenate([forced, rank[:k]])
            x_train = (inner_features[np.ix_(inner_train, selected)] - mean[selected]) / scale[selected]
            x_valid = (inner_features[np.ix_(inner_valid, selected)] - mean[selected]) / scale[selected]
            for alpha in ALPHAS:
                model = Ridge(alpha=float(alpha), fit_intercept=True, solver="lsqr", tol=1.0e-7)
                model.fit(x_train, inner_train_y)
                prediction = model.predict(x_valid)
                choices.append(
                    (float(np.sqrt(np.mean((inner_valid_y - prediction) ** 2))), k, float(alpha))
                )
        floor = [(constant_rmse, 0, math.inf)] if (constant_floor or not choices) else []
        best_rmse, best_k, best_alpha = min(
            [*choices, *floor], key=lambda v: (v[0], v[1], v[2])
        )
    outer_mean, outer_scale, outer_rank = rank_columns(
        outer_features, outcomes, outer_train, candidate
    )
    if best_k == 0 and forced.size == 0:
        prediction = np.full(apply_rows.size, float(outcomes[outer_train].mean()))
        used_alpha: float | str = "constant_training_mean"
        used_size = 0
    else:
        selected = np.concatenate([forced, outer_rank[: min(int(best_k), int(outer_rank.size))]])
        x_train = (outer_features[np.ix_(outer_train, selected)] - outer_mean[selected]) / outer_scale[selected]
        x_apply = (outer_features[np.ix_(apply_rows, selected)] - outer_mean[selected]) / outer_scale[selected]
        alpha = 1.0e6 if not np.isfinite(best_alpha) else float(best_alpha)
        model = Ridge(alpha=alpha, fit_intercept=True, solver="lsqr", tol=1.0e-7)
        model.fit(x_train, outcomes[outer_train])
        prediction = np.asarray(model.predict(x_apply), dtype=np.float64)
        used_alpha = float(best_alpha) if np.isfinite(best_alpha) else "inf_floored_to_1e6"
        used_size = int(selected.size)
    receipt = {
        "selected_k": int(best_k),
        "selected_alpha": used_alpha,
        "constant_floor_offered": bool(constant_floor),
        "selected_feature_count": used_size,
        "forced_columns": int(forced.size),
        "inner_validation_rmse": best_rmse,
        "constant_inner_rmse": constant_rmse,
        "inner_train_rows": int(inner_train.size),
        "inner_valid_rows": int(inner_valid.size),
        "outer_train_rows": int(np.unique(outer_train).size),
        "apply_rows": int(apply_rows.size),
    }
    return prediction, (int(best_k), float(best_alpha)), receipt


class ProjectionCache:
    """Whitened SVD projections keyed by the exact training row set."""

    def __init__(self, embeddings: np.ndarray) -> None:
        self.embeddings = embeddings
        self.store: dict[bytes, np.ndarray] = {}
        self.fits = 0

    def features(self, train_mask: np.ndarray) -> np.ndarray:
        # The shipped projection keeps 256 components unconditionally, so it requires at least 256
        # training allele rows (4 per variant).  Every training pool here is far above that; the
        # guard makes the requirement explicit instead of failing inside the SVD.
        if 4 * int(train_mask.sum()) < 256:
            raise TransferError("projection training pool has fewer than 256 allele rows")
        key = np.packbits(train_mask).tobytes()
        cached = self.store.get(key)
        if cached is None:
            self.fits += 1
            cached = head_features(
                apply_projection(self.embeddings, fit_projection(self.embeddings, train_mask))
            ).astype(np.float64)[:, DELTA_SLICE]
            self.store[key] = cached
        return cached


# ---------------------------------------------------------------- bootstrap


def block_bootstrap(
    observed: np.ndarray,
    arms: Mapping[str, np.ndarray],
    folds: np.ndarray,
    blocks: np.ndarray,
    *,
    resamples: int,
    seed: int,
) -> dict[str, np.ndarray]:
    present = sorted(set(folds.tolist()))
    layout = {}
    for fold in present:
        rows = np.flatnonzero(folds == fold)
        unique = sorted(set(blocks[rows].tolist()))
        layout[fold] = (unique, [rows[blocks[rows] == block] for block in unique])
    rng = np.random.default_rng(seed)
    names = list(arms)
    macro = {name: np.full(resamples, np.nan) for name in names}
    for iteration in range(resamples):
        per_fold = []
        for fold in present:
            unique, groups = layout[fold]
            sampled = rng.integers(0, len(unique), size=len(unique))
            per_fold.append(np.concatenate([groups[index] for index in sampled]))
        fold_observed = [observed[indices] for indices in per_fold]
        for name in names:
            values = arms[name]
            macro[name][iteration] = fisher_macro(
                [fast_spearman(fold_observed[i], values[per_fold[i]]) for i in range(len(present))]
            )
    return macro


def bootstrap_p(samples: np.ndarray) -> float:
    valid = samples[np.isfinite(samples)]
    if valid.size == 0:
        return float("nan")
    below = float(np.mean(valid <= 0.0))
    above = float(np.mean(valid >= 0.0))
    return float(min(1.0, max(2.0 * min(below, above), 1.0 / (valid.size + 1))))


def macro_over_folds(observed: np.ndarray, predicted: np.ndarray, folds: np.ndarray) -> float:
    present = sorted(set(folds.tolist()))
    return fisher_macro(
        [fast_spearman(observed[folds == f], predicted[folds == f]) for f in present]
    )


# ---------------------------------------------------------------- campaign


def build_reporter_oof(
    *,
    cache: ProjectionCache,
    reporter_y: np.ndarray,
    reporter_pool: np.ndarray,
    folds: np.ndarray,
    blocks: np.ndarray,
    seed: int,
    audit: list[dict[str, object]],
) -> tuple[dict[int, np.ndarray], np.ndarray]:
    """Nested out-of-fold reporter predictions, one full vector per evaluation fold."""
    n = folds.size
    per_eval_fold: dict[int, np.ndarray] = {}
    self_eval = np.full(n, np.nan)
    for held in range(5):
        column = np.full(n, np.nan)
        for target in range(5):
            if target == held:
                pool = reporter_pool & (folds != held)
                valid_fold = (held + 1) % 5
                apply_rows = np.flatnonzero(folds == held)
                tag = "outer"
            else:
                pool = reporter_pool & (folds != held) & (folds != target)
                valid_fold = next(
                    f for f in [(target + step) % 5 for step in range(1, 5)]
                    if f not in {held, target}
                )
                apply_rows = np.flatnonzero(folds == target)
                tag = "nested"
            inner_train_mask = pool & (folds != valid_fold)
            inner_valid_mask = pool & (folds == valid_fold)
            if inner_valid_mask.sum() < 20 or inner_train_mask.sum() < 50:
                raise TransferError("reporter inner split is too small")
            inner_features = cache.features(inner_train_mask)
            outer_features = cache.features(pool)
            inner_train = bootstrap_indices(
                inner_train_mask, blocks, seed=seed, salt=REPORTER_INNER_SALT + 5 * held + target
            )
            outer_train = bootstrap_indices(
                pool, blocks, seed=seed, salt=REPORTER_OUTER_SALT + 5 * held + target
            )
            candidate = np.arange(inner_features.shape[1])
            prediction, _, receipt = fit_head(
                inner_features=inner_features,
                outer_features=outer_features,
                outcomes=reporter_y,
                inner_train=inner_train,
                inner_valid=np.flatnonzero(inner_valid_mask),
                outer_train=outer_train,
                apply_rows=apply_rows,
                candidate=candidate,
                forced=np.empty(0, dtype=int),
                top_k_grid=DELTA_TOP_K,
            )
            column[apply_rows] = prediction
            if tag == "outer":
                self_eval[apply_rows] = prediction
            audit.append({
                "component": "reporter_head", "seed": seed, "evaluation_fold": held,
                "predicted_fold": target, "role": tag, "validation_fold": valid_fold,
                "reporter_training_elements": int(pool.sum()),
                "reporter_training_blocks": int(len(set(blocks[pool].tolist()))),
                "reporter_training_shares_block_with_evaluation_fold": int(
                    len(set(blocks[pool].tolist()) & set(blocks[folds == held].tolist()))
                ),
                **receipt,
            })
        if np.isnan(column).any():
            raise TransferError("reporter out-of-fold column has gaps")
        per_eval_fold[held] = column
    return per_eval_fold, self_eval


def run_campaign(
    *,
    cache: ProjectionCache,
    identity: np.ndarray,
    endo_y: np.ndarray,
    endo_pool: np.ndarray,
    reporter_oof: Mapping[int, dict[int, np.ndarray]],
    folds: np.ndarray,
    blocks: np.ndarray,
    audit: list[dict[str, object]],
    tag: str,
) -> dict[str, dict[int, np.ndarray]]:
    n = folds.size
    arms = ("endogenous_only", "plus_reporter_oof", "plus_reporter_oof_free_selection",
            "plus_reporter_oof_locked_hparams", "allele_identity_ridge")
    per_seed: dict[str, dict[int, np.ndarray]] = {arm: {} for arm in arms}
    for seed in SEEDS:
        predictions = {arm: np.full(n, np.nan) for arm in arms}
        for held in range(5):
            inner_fold = (held + 1) % 5
            outer_mask = endo_pool & (folds != held)
            inner_mask = outer_mask & (folds != inner_fold)
            inner_valid_mask = endo_pool & (folds == inner_fold)
            evaluation = np.flatnonzero(endo_pool & (folds == held))
            if set(blocks[outer_mask].tolist()) & set(blocks[evaluation].tolist()):
                raise TransferError("held-out block overlaps endogenous training")
            inner_delta = cache.features(inner_mask)
            outer_delta = cache.features(outer_mask)
            inner_train = bootstrap_indices(inner_mask, blocks, seed=seed, salt=INNER_SALT + held)
            outer_train = bootstrap_indices(outer_mask, blocks, seed=seed, salt=OUTER_SALT + held)
            inner_valid = np.flatnonzero(inner_valid_mask)
            width = inner_delta.shape[1]
            reporter_column = reporter_oof[seed][held]

            base, locked, receipt = fit_head(
                inner_features=inner_delta, outer_features=outer_delta, outcomes=endo_y,
                inner_train=inner_train, inner_valid=inner_valid, outer_train=outer_train,
                apply_rows=evaluation, candidate=np.arange(width),
                forced=np.empty(0, dtype=int), top_k_grid=DELTA_TOP_K,
            )
            predictions["endogenous_only"][evaluation] = base
            audit.append({"component": "endogenous_head", "arm": "endogenous_only", "stratum": tag,
                          "seed": seed, "evaluation_fold": held, "validation_fold": inner_fold,
                          "endogenous_training_variants": int(outer_mask.sum()),
                          "endogenous_training_blocks": int(len(set(blocks[outer_mask].tolist()))),
                          "evaluation_variants": int(evaluation.size),
                          "evaluation_blocks": int(len(set(blocks[evaluation].tolist()))),
                          "training_blocks_shared_with_evaluation": 0, **receipt})

            inner_plus = np.column_stack([inner_delta, reporter_column])
            outer_plus = np.column_stack([outer_delta, reporter_column])
            forced = np.asarray([width], dtype=int)
            for arm, forced_columns, candidate, lock in (
                ("plus_reporter_oof", forced, np.arange(width), None),
                ("plus_reporter_oof_free_selection", np.empty(0, int), np.arange(width + 1), None),
                ("plus_reporter_oof_locked_hparams", forced, np.arange(width), locked),
            ):
                value, _, extra = fit_head(
                    inner_features=inner_plus, outer_features=outer_plus, outcomes=endo_y,
                    inner_train=inner_train, inner_valid=inner_valid, outer_train=outer_train,
                    apply_rows=evaluation, candidate=candidate, forced=forced_columns,
                    top_k_grid=DELTA_TOP_K, locked=lock,
                )
                predictions[arm][evaluation] = value
                audit.append({"component": "endogenous_head", "arm": arm, "stratum": tag,
                              "seed": seed, "evaluation_fold": held, "validation_fold": inner_fold,
                              "endogenous_training_variants": int(outer_mask.sum()),
                              "endogenous_training_blocks": int(len(set(blocks[outer_mask].tolist()))),
                              "evaluation_variants": int(evaluation.size),
                              "evaluation_blocks": int(len(set(blocks[evaluation].tolist()))),
                              "training_blocks_shared_with_evaluation": 0, **extra})

            value, _, extra = fit_head(
                inner_features=identity, outer_features=identity, outcomes=endo_y,
                inner_train=inner_train, inner_valid=inner_valid, outer_train=outer_train,
                apply_rows=evaluation, candidate=np.arange(identity.shape[1]),
                forced=np.empty(0, dtype=int), top_k_grid=ALLELE_TOP_K,
            )
            predictions["allele_identity_ridge"][evaluation] = value
            audit.append({"component": "endogenous_head", "arm": "allele_identity_ridge",
                          "stratum": tag, "seed": seed, "evaluation_fold": held,
                          "validation_fold": inner_fold,
                          "endogenous_training_variants": int(outer_mask.sum()),
                          "endogenous_training_blocks": int(len(set(blocks[outer_mask].tolist()))),
                          "evaluation_variants": int(evaluation.size),
                          "evaluation_blocks": int(len(set(blocks[evaluation].tolist()))),
                          "training_blocks_shared_with_evaluation": 0, **extra})
        for arm in arms:
            per_seed[arm][seed] = predictions[arm]
        print(f"[{tag}] seed {seed} done", flush=True)
    return per_seed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resamples", type=int, default=10000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260914)
    arguments = parser.parse_args()
    out = arguments.output
    (out / "tables").mkdir(parents=True, exist_ok=True)

    units = pd.read_csv(arguments.inputs / "e2_units.tsv.gz", sep="\t")
    units = units.sort_values("row_index").reset_index(drop=True)
    n = len(units)
    folds = units.chrom_fold.to_numpy(np.int64)
    blocks = units.block_1mb.to_numpy(dtype=str)

    with np.load(arguments.raw / "hyenadna_allele_embeddings.npz", allow_pickle=False) as data:
        if tuple(data["allele_order"].astype(str).tolist()) != ("REF", "ALT", "REF_RC", "ALT_RC"):
            raise TransferError("embedding allele order differs")
        if not np.array_equal(data["row_index"], units.row_index.to_numpy(np.int64)):
            raise TransferError("embedding row alignment differs")
        embeddings = data["embeddings"]
    if embeddings.shape != (n, 4, 256) or not np.isfinite(embeddings).all():
        raise TransferError("embedding geometry differs")

    cache = ProjectionCache(embeddings)
    identity = allele_identity_features(
        units.mpra_ref.to_numpy(dtype=str), units.mpra_alt.to_numpy(dtype=str)
    )

    endo_pool_full = units.has_endogenous.to_numpy(bool)
    reporter_variants = {
        "hepg2_mean": "d_hepg2_mean",
        "hepg2_control": "d_HepG2_control",
        "hepg2_paoa": "d_HepG2_PAOA",
        "all_context_mean": "d_all_context_mean",
    }

    summary_rows: list[dict[str, object]] = []
    per_seed_rows: list[dict[str, object]] = []
    per_fold_rows: list[dict[str, object]] = []
    reporter_rows: list[dict[str, object]] = []
    audit: list[dict[str, object]] = []
    fold_rows: list[dict[str, object]] = []

    # ------------------------------------------------ reporter heads, one per label variant
    reporter_oof_by_label: dict[str, dict[int, dict[int, np.ndarray]]] = {}
    for label_name, column in reporter_variants.items():
        values = units[column].to_numpy(np.float64)
        pool = np.isfinite(values)
        filled = np.where(pool, values, 0.0)
        per_seed_oof: dict[int, dict[int, np.ndarray]] = {}
        self_eval_by_seed: dict[int, np.ndarray] = {}
        for seed in SEEDS:
            clock = time.time()
            per_eval, self_eval = build_reporter_oof(
                cache=cache, reporter_y=filled, reporter_pool=pool, folds=folds, blocks=blocks,
                seed=seed, audit=audit,
            )
            per_seed_oof[seed] = per_eval
            self_eval_by_seed[seed] = self_eval
            print(f"[reporter/{label_name}] seed {seed} in {time.time() - clock:.1f}s", flush=True)
        reporter_oof_by_label[label_name] = per_seed_oof
        ensemble = np.mean(np.vstack([self_eval_by_seed[s] for s in SEEDS]), axis=0)
        rows = np.flatnonzero(pool)
        reporter_rows.append({
            "reporter_label": label_name,
            "elements_with_label": int(pool.sum()),
            "macro_spearman_vs_reporter_label": macro_over_folds(
                values[rows], ensemble[rows], folds[rows]
            ),
            "pooled_spearman_vs_reporter_label": fast_spearman(values[rows], ensemble[rows]),
            **{
                f"seed_{s}_macro": macro_over_folds(
                    values[rows], self_eval_by_seed[s][rows], folds[rows]
                )
                for s in SEEDS
            },
        })
        print(f"[reporter/{label_name}] self-evaluation macro "
              f"{reporter_rows[-1]['macro_spearman_vs_reporter_label']:.4f}", flush=True)

    pd.DataFrame(reporter_rows).to_csv(
        out / "tables" / "reporter_head_self_evaluation.tsv", sep="\t", index=False
    )

    # ------------------------------------------------ leakage evidence
    for label_name in reporter_variants:
        pool = np.isfinite(units[reporter_variants[label_name]].to_numpy(np.float64))
        for held in range(5):
            evaluation = endo_pool_full & (folds == held)
            endo_train = endo_pool_full & (folds != held)
            rep_train = pool & (folds != held)
            fold_rows.append({
                "reporter_label": label_name,
                "evaluation_fold": held,
                "evaluation_chromosomes": ",".join(
                    sorted(set(units.contig[evaluation].tolist()))
                ),
                "evaluation_variants": int(evaluation.sum()),
                "evaluation_blocks_1mb": int(units.block_1mb[evaluation].nunique()),
                "endogenous_training_variants": int(endo_train.sum()),
                "endogenous_training_blocks_1mb": int(units.block_1mb[endo_train].nunique()),
                "reporter_training_elements": int(rep_train.sum()),
                "reporter_training_blocks_1mb": int(units.block_1mb[rep_train].nunique()),
                "elements_in_both_evaluation_and_endogenous_training": int(
                    (evaluation & endo_train).sum()
                ),
                "elements_in_both_evaluation_and_reporter_training": int(
                    (evaluation & rep_train).sum()
                ),
                "blocks_1mb_shared_evaluation_and_endogenous_training": int(
                    len(set(blocks[evaluation]) & set(blocks[endo_train]))
                ),
                "blocks_1mb_shared_evaluation_and_reporter_training": int(
                    len(set(blocks[evaluation]) & set(blocks[rep_train]))
                ),
                "blocks_lr239_shared_evaluation_and_reporter_training": int(
                    len(set(units.block_lr239[evaluation]) & set(units.block_lr239[rep_train]))
                ),
            })
    pd.DataFrame(fold_rows).to_csv(out / "tables" / "fold_membership.tsv", sep="\t", index=False)

    # ------------------------------------------------ campaigns
    campaigns: list[tuple[str, str, str, np.ndarray, np.ndarray]] = []
    beta_nearest = units.beta_nearest.to_numpy(np.float64)
    beta_minp = units.beta_minp.to_numpy(np.float64)
    campaigns.append(("primary", "nearest_peak", "hepg2_mean", beta_nearest, endo_pool_full))
    for label_name in ("hepg2_control", "hepg2_paoa", "all_context_mean"):
        campaigns.append((f"reporter_label_{label_name}", "nearest_peak", label_name,
                          beta_nearest, endo_pool_full))
    campaigns.append(("outcome_selected_minp_peak", "minp_peak", "hepg2_mean",
                      beta_minp, np.isfinite(beta_minp)))
    absolute = units.abs_distance_nearest.to_numpy(np.float64)
    within_primary = endo_pool_full & (absolute <= 5000.0)
    campaigns.append(("train_within_primary_stratum", "nearest_peak", "hepg2_mean",
                      beta_nearest, within_primary))

    for campaign, endo_label, reporter_label, raw_y, pool in campaigns:
        pool = pool & np.isfinite(raw_y)
        y = np.where(pool, raw_y, 0.0)
        per_seed = run_campaign(
            cache=cache, identity=identity, endo_y=y, endo_pool=pool,
            reporter_oof=reporter_oof_by_label[reporter_label], folds=folds, blocks=blocks,
            audit=audit, tag=campaign,
        )
        # Every pooled row is predicted by every seed, so a plain mean is defined on the pool and
        # NaN elsewhere; nanmean would only hide a coverage gap.
        ensemble = {
            arm: np.mean(np.vstack([values[s] for s in SEEDS]), axis=0)
            for arm, values in per_seed.items()
        }
        reporter_alone = np.full(n, np.nan)
        for held in range(5):
            rows = np.flatnonzero(pool & (folds == held))
            reporter_alone[rows] = np.mean(
                np.vstack([reporter_oof_by_label[reporter_label][s][held][rows] for s in SEEDS]),
                axis=0,
            )
        ensemble["reporter_oof_alone"] = reporter_alone

        strata = STRATA if campaign == "primary" else ((PRIMARY_STRATUM, 5000.0),)
        for stratum_name, threshold in strata:
            rows = np.flatnonzero(pool & (absolute <= threshold))
            if rows.size == 0:
                continue
            observed = raw_y[rows]
            arms = {name: values[rows] for name, values in ensemble.items()}
            if not all(np.isfinite(value).all() for value in arms.values()):
                raise TransferError("arm prediction coverage differs")
            clock = time.time()
            samples = block_bootstrap(
                observed, arms, folds[rows], blocks[rows],
                resamples=arguments.resamples, seed=arguments.bootstrap_seed,
            )
            print(f"[{campaign}/{stratum_name}] bootstrap in {time.time() - clock:.1f}s", flush=True)
            reference = samples["endogenous_only"]
            for name, values in arms.items():
                point = macro_over_folds(observed, values, folds[rows])
                low, high, standard_error = interval(samples[name])
                gain = samples[name] - reference
                gain_low, gain_high, gain_se = interval(gain)
                seed_gain = []
                for seed in SEEDS:
                    if name == "reporter_oof_alone":
                        continue
                    arm_value = macro_over_folds(observed, per_seed[name][seed][rows], folds[rows])
                    base_value = macro_over_folds(
                        observed, per_seed["endogenous_only"][seed][rows], folds[rows]
                    )
                    per_seed_rows.append({
                        "campaign": campaign, "stratum_bp": stratum_name, "arm": name,
                        "seed": seed, "macro_spearman": arm_value,
                        "gain_vs_endogenous_only": arm_value - base_value,
                    })
                    seed_gain.append(arm_value - base_value)
                summary_rows.append({
                    "campaign": campaign,
                    "endogenous_label": endo_label,
                    "reporter_label": reporter_label,
                    "stratum_max_abs_distance_bp": stratum_name,
                    "arm": name,
                    "evaluation_variants": int(rows.size),
                    "evaluation_blocks_1mb": int(len(set(blocks[rows].tolist()))),
                    "training_variants": int(pool.sum()),
                    "macro_spearman": point,
                    "macro_ci_low": low,
                    "macro_ci_high": high,
                    "macro_bootstrap_se": standard_error,
                    "gain_vs_endogenous_only": point
                    - macro_over_folds(observed, arms["endogenous_only"], folds[rows]),
                    "gain_ci_low": gain_low,
                    "gain_ci_high": gain_high,
                    "gain_bootstrap_se": gain_se,
                    "gain_bootstrap_p_two_sided": bootstrap_p(gain)
                    if name != "endogenous_only" else float("nan"),
                    "gain_ci_excludes_zero": bool(
                        np.isfinite(gain_low) and np.isfinite(gain_high)
                        and (gain_low > 0 or gain_high < 0)
                    ),
                    "positive_gain_seeds": int(sum(1 for v in seed_gain if v > 0)),
                    "comparable_seeds": int(sum(1 for v in seed_gain if np.isfinite(v))),
                })
                for fold in sorted(set(folds[rows].tolist())):
                    keep = rows[folds[rows] == fold]
                    per_fold_rows.append({
                        "campaign": campaign, "stratum_bp": stratum_name, "arm": name,
                        "evaluation_fold": fold, "variants": int(keep.size),
                        "blocks_1mb": int(len(set(blocks[keep].tolist()))),
                        "spearman": fast_spearman(raw_y[keep], ensemble[name][keep]),
                    })
            if campaign == "primary" and stratum_name == PRIMARY_STRATUM:
                np.savez_compressed(
                    out / "tables" / "primary_bootstrap_samples.npz",
                    **{name: samples[name] for name in samples},
                )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out / "tables" / "transfer_summary.tsv", sep="\t", index=False)
    pd.DataFrame(per_seed_rows).to_csv(out / "tables" / "per_seed_macro.tsv", sep="\t", index=False)
    pd.DataFrame(per_fold_rows).to_csv(out / "tables" / "per_fold_spearman.tsv", sep="\t", index=False)
    pd.DataFrame(audit).to_csv(out / "tables" / "fit_audit.tsv.gz", sep="\t", index=False,
                               compression="gzip")

    primary = summary[
        (summary.campaign == "primary")
        & (summary.stratum_max_abs_distance_bp == PRIMARY_STRATUM)
        & (summary.arm == "plus_reporter_oof")
    ]
    if len(primary) != 1:
        raise TransferError("primary contrast row is not unique")
    row = primary.iloc[0]
    family = dict(C1_FAMILY_P)
    family["e2_transfer_contrast"] = float(row.gain_bootstrap_p_two_sided)
    order = sorted(family.items(), key=lambda item: item[1])
    holm = []
    previous_rejected = True
    for index, (name, value) in enumerate(order):
        threshold = 0.05 / (len(order) - index)
        rejected = bool(previous_rejected and value <= threshold)
        previous_rejected = rejected
        holm.append({"member": name, "raw_p": value, "holm_threshold": threshold,
                     "rejected_at_fwer_0.05": rejected, "rank": index + 1})
    pd.DataFrame(holm).to_csv(out / "tables" / "holm_family.tsv", sep="\t", index=False)

    reporter_frame = pd.DataFrame(reporter_rows).set_index("reporter_label")
    alone = summary[
        (summary.campaign == "primary")
        & (summary.stratum_max_abs_distance_bp == PRIMARY_STRATUM)
        & (summary.arm == "reporter_oof_alone")
    ].iloc[0]
    predictions = [
        {
            "prediction": "E2.1",
            "statement": "gain of arm B over arm A is below 0.02 and its 95% interval includes 0",
            "observed": (
                f"gain {row.gain_vs_endogenous_only:.4f} "
                f"[{row.gain_ci_low:.4f}, {row.gain_ci_high:.4f}]"
            ),
            "met": bool(
                row.gain_vs_endogenous_only < 0.02
                and row.gain_ci_low <= 0.0 <= row.gain_ci_high
            ),
        },
        {
            "prediction": "E2.2",
            "statement": "the reporter head is not dead: out-of-fold macro Spearman on the "
                         "reporter label it was trained on is at least 0.10",
            "observed": f"{reporter_frame.loc['hepg2_mean', 'macro_spearman_vs_reporter_label']:.4f}",
            "met": bool(
                reporter_frame.loc["hepg2_mean", "macro_spearman_vs_reporter_label"] >= 0.10
            ),
        },
        {
            "prediction": "E2.3",
            "statement": "reporter_oof_alone has |signed Spearman| below 0.05 against the caQTL "
                         "label in the primary stratum",
            "observed": f"{alone.macro_spearman:.4f} "
                        f"[{alone.macro_ci_low:.4f}, {alone.macro_ci_high:.4f}]",
            "met": bool(abs(alone.macro_spearman) < 0.05),
        },
    ]
    pd.DataFrame(predictions).to_csv(out / "tables" / "predictions.tsv", sep="\t", index=False)

    receipt = {
        "schema_version": "agp-e2-transfer-v1",
        "seeds": list(SEEDS),
        "bootstrap_seed": arguments.bootstrap_seed,
        "bootstrap_resamples": arguments.resamples,
        "primary_stratum_max_abs_distance_bp": int(PRIMARY_STRATUM),
        "projection_fits": cache.fits,
        "primary_gain": float(row.gain_vs_endogenous_only),
        "primary_gain_ci": [float(row.gain_ci_low), float(row.gain_ci_high)],
        "primary_gain_bootstrap_p": float(row.gain_bootstrap_p_two_sided),
        "holm": holm,
        "predictions": predictions,
    }
    (out / "tables" / "run_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
