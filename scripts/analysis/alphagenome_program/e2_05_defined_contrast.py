#!/usr/bin/env python3
"""E2 step 5: the transfer contrast under a head that is defined in every fold.

Why this exists.  Under the shipped reporter head recipe the inner validation offers the constant
training mean as a candidate, and on this substrate the endogenous-only arm selects it in 3 of the 5
held-out folds at every seed (`tables/fit_audit.tsv.gz`).  A constant prediction has no rank, so the
fold's Spearman is undefined and the reporter's Fisher macro rule, which is undefined if any
component is undefined, makes arm A - and therefore the paired contrast - undefined.  That is a real
property of the endogenous head on this label, and it is reported as such; it is not a reason to
leave the prespecified contrast unanswered.

This script therefore runs both head modes over the same folds, blocks, seeds and bootstrap:

* `shipped_floor`      - the constant training mean competes (the published recipe).
* `no_constant_floor`  - that one candidate is removed, identically in both arms, so every fold
                         returns a fitted ridge and the macro is defined.

Removing the floor is applied to arm A and arm B alike, so the two arms still differ in exactly one
component: the presence of the reporter out-of-fold column.  Three estimators are reported for each
mode: the Fisher macro over all five folds, the Fisher macro over the folds where BOTH arms are
defined, and the pooled Spearman over all held-out rows.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from c2_03_fit_heads import allele_identity_features, fast_spearman, fisher_macro, interval  # noqa: E402
from e2_02_fit_transfer import (  # noqa: E402
    ALLELE_TOP_K,
    DELTA_TOP_K,
    PRIMARY_STRATUM,
    SEEDS,
    STRATA,
    C1_FAMILY_P,
    ProjectionCache,
    TransferError,
    bootstrap_indices,
    bootstrap_p,
    build_reporter_oof,
    fit_head,
    macro_over_folds,
)

MODES = ("shipped_floor", "no_constant_floor")


def macro_over(observed, predicted, folds, keep_folds) -> float:
    return fisher_macro(
        [fast_spearman(observed[folds == f], predicted[folds == f]) for f in keep_folds]
    )


def paired_bootstrap(
    observed: np.ndarray,
    arms: Mapping[str, np.ndarray],
    folds: np.ndarray,
    blocks: np.ndarray,
    keep_folds: list[int],
    *,
    resamples: int,
    seed: int,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Block bootstrap giving macro over all folds, macro over `keep_folds`, and pooled."""
    present = sorted(set(folds.tolist()))
    layout = {}
    for fold in present:
        rows = np.flatnonzero(folds == fold)
        unique = sorted(set(blocks[rows].tolist()))
        layout[fold] = (unique, [rows[blocks[rows] == block] for block in unique])
    rng = np.random.default_rng(seed)
    names = list(arms)
    macro_all = {name: np.full(resamples, np.nan) for name in names}
    macro_keep = {name: np.full(resamples, np.nan) for name in names}
    pooled = {name: np.full(resamples, np.nan) for name in names}
    for iteration in range(resamples):
        per_fold = {}
        for fold in present:
            unique, groups = layout[fold]
            sampled = rng.integers(0, len(unique), size=len(unique))
            per_fold[fold] = np.concatenate([groups[index] for index in sampled])
        stacked = np.concatenate([per_fold[fold] for fold in present])
        pooled_observed = observed[stacked]
        fold_observed = {fold: observed[per_fold[fold]] for fold in present}
        for name in names:
            values = arms[name]
            per = {
                fold: fast_spearman(fold_observed[fold], values[per_fold[fold]])
                for fold in present
            }
            macro_all[name][iteration] = fisher_macro([per[f] for f in present])
            macro_keep[name][iteration] = fisher_macro([per[f] for f in keep_folds])
            pooled[name][iteration] = fast_spearman(pooled_observed, values[stacked])
    return macro_all, macro_keep, pooled


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

    units = pd.read_csv(arguments.inputs / "e2_units.tsv.gz", sep="\t").sort_values("row_index")
    n = len(units)
    folds = units.chrom_fold.to_numpy(np.int64)
    blocks = units.block_1mb.to_numpy(dtype=str)
    absolute = units.abs_distance_nearest.to_numpy(np.float64)
    beta = units.beta_nearest.to_numpy(np.float64)
    pool = units.has_endogenous.to_numpy(bool) & np.isfinite(beta)
    endo_y = np.where(pool, beta, 0.0)

    with np.load(arguments.raw / "hyenadna_allele_embeddings.npz", allow_pickle=False) as data:
        if not np.array_equal(data["row_index"], units.row_index.to_numpy(np.int64)):
            raise TransferError("embedding row alignment differs")
        embeddings = data["embeddings"]
    cache = ProjectionCache(embeddings)
    identity = allele_identity_features(
        units.mpra_ref.to_numpy(dtype=str), units.mpra_alt.to_numpy(dtype=str)
    )

    reporter_values = units.d_hepg2_mean.to_numpy(np.float64)
    reporter_pool = np.isfinite(reporter_values)
    audit: list[dict[str, object]] = []
    reporter_oof: dict[int, dict[int, np.ndarray]] = {}
    for seed in SEEDS:
        reporter_oof[seed], _ = build_reporter_oof(
            cache=cache, reporter_y=np.where(reporter_pool, reporter_values, 0.0),
            reporter_pool=reporter_pool, folds=folds, blocks=blocks, seed=seed, audit=audit,
        )
        print(f"[reporter] seed {seed} done", flush=True)

    predictions: dict[tuple[str, str], dict[int, np.ndarray]] = {}
    for mode in MODES:
        arms = ("endogenous_only", "plus_reporter_oof", "allele_identity_ridge")
        for arm in arms:
            predictions[(mode, arm)] = {}
        for seed in SEEDS:
            store = {arm: np.full(n, np.nan) for arm in arms}
            for held in range(5):
                inner_fold = (held + 1) % 5
                outer_mask = pool & (folds != held)
                inner_mask = outer_mask & (folds != inner_fold)
                inner_valid = np.flatnonzero(pool & (folds == inner_fold))
                evaluation = np.flatnonzero(pool & (folds == held))
                inner_delta = cache.features(inner_mask)
                outer_delta = cache.features(outer_mask)
                width = inner_delta.shape[1]
                inner_train = bootstrap_indices(inner_mask, blocks, seed=seed, salt=10_000 + held)
                outer_train = bootstrap_indices(outer_mask, blocks, seed=seed, salt=20_000 + held)
                column = reporter_oof[seed][held]
                base, _, receipt = fit_head(
                    inner_features=inner_delta, outer_features=outer_delta, outcomes=endo_y,
                    inner_train=inner_train, inner_valid=inner_valid, outer_train=outer_train,
                    apply_rows=evaluation, candidate=np.arange(width),
                    forced=np.empty(0, dtype=int), top_k_grid=DELTA_TOP_K,
                    constant_floor=(mode == "shipped_floor"),
                )
                store["endogenous_only"][evaluation] = base
                audit.append({"component": "defined_contrast", "mode": mode,
                              "arm": "endogenous_only", "seed": seed,
                              "evaluation_fold": held, **receipt})
                value, _, receipt = fit_head(
                    inner_features=np.column_stack([inner_delta, column]),
                    outer_features=np.column_stack([outer_delta, column]),
                    outcomes=endo_y, inner_train=inner_train, inner_valid=inner_valid,
                    outer_train=outer_train, apply_rows=evaluation,
                    candidate=np.arange(width), forced=np.asarray([width], dtype=int),
                    top_k_grid=DELTA_TOP_K, constant_floor=(mode == "shipped_floor"),
                )
                store["plus_reporter_oof"][evaluation] = value
                audit.append({"component": "defined_contrast", "mode": mode,
                              "arm": "plus_reporter_oof", "seed": seed,
                              "evaluation_fold": held, **receipt})
                value, _, receipt = fit_head(
                    inner_features=identity, outer_features=identity, outcomes=endo_y,
                    inner_train=inner_train, inner_valid=inner_valid, outer_train=outer_train,
                    apply_rows=evaluation, candidate=np.arange(identity.shape[1]),
                    forced=np.empty(0, dtype=int), top_k_grid=ALLELE_TOP_K,
                    constant_floor=(mode == "shipped_floor"),
                )
                store["allele_identity_ridge"][evaluation] = value
                audit.append({"component": "defined_contrast", "mode": mode,
                              "arm": "allele_identity_ridge", "seed": seed,
                              "evaluation_fold": held, **receipt})
            for arm in arms:
                predictions[(mode, arm)][seed] = store[arm]
            print(f"[{mode}] seed {seed} done", flush=True)

    reporter_alone = np.full(n, np.nan)
    for held in range(5):
        rows = np.flatnonzero(pool & (folds == held))
        reporter_alone[rows] = np.mean(
            np.vstack([reporter_oof[s][held][rows] for s in SEEDS]), axis=0
        )

    summary: list[dict[str, object]] = []
    per_fold_rows: list[dict[str, object]] = []
    per_seed_rows: list[dict[str, object]] = []
    saved = {}
    for mode in MODES:
        ensemble = {
            arm: np.mean(np.vstack([predictions[(mode, arm)][s] for s in SEEDS]), axis=0)
            for arm in ("endogenous_only", "plus_reporter_oof", "allele_identity_ridge")
        }
        ensemble["reporter_oof_alone"] = reporter_alone
        for stratum, threshold in STRATA:
            rows = np.flatnonzero(pool & (absolute <= threshold))
            observed = beta[rows]
            arms = {name: values[rows] for name, values in ensemble.items()}
            present = sorted(set(folds[rows].tolist()))
            defined = {
                name: [
                    fold for fold in present
                    if np.isfinite(fast_spearman(observed[folds[rows] == fold],
                                                 values[folds[rows] == fold]))
                ]
                for name, values in arms.items()
            }
            keep_folds = [
                fold for fold in present
                if fold in defined["endogenous_only"] and fold in defined["plus_reporter_oof"]
            ]
            macro_all, macro_keep, pooled = paired_bootstrap(
                observed, arms, folds[rows], blocks[rows], keep_folds,
                resamples=arguments.resamples, seed=arguments.bootstrap_seed,
            )
            if mode == "no_constant_floor" and stratum == PRIMARY_STRATUM:
                saved = {f"{k}__{name}": v[name]
                         for k, v in (("macro_all", macro_all), ("macro_keep", macro_keep),
                                      ("pooled", pooled))
                         for name in v}
            base = "endogenous_only"
            for name, values in arms.items():
                row: dict[str, object] = {
                    "mode": mode,
                    "stratum_max_abs_distance_bp": stratum,
                    "arm": name,
                    "evaluation_variants": int(rows.size),
                    "evaluation_blocks_1mb": int(len(set(blocks[rows].tolist()))),
                    "folds_with_defined_spearman": ",".join(str(f) for f in defined[name]),
                    "folds_used_by_macro_keep": ",".join(str(f) for f in keep_folds),
                }
                for tag, samples, point in (
                    ("macro_all_folds", macro_all,
                     macro_over(observed, values, folds[rows], present)),
                    ("macro_defined_folds", macro_keep,
                     macro_over(observed, values, folds[rows], keep_folds)),
                    ("pooled", pooled, fast_spearman(observed, values)),
                ):
                    low, high, standard_error = interval(samples[name])
                    gain_samples = samples[name] - samples[base]
                    gain_low, gain_high, _ = interval(gain_samples)
                    row[f"{tag}_spearman"] = point
                    row[f"{tag}_ci_low"] = low
                    row[f"{tag}_ci_high"] = high
                    row[f"{tag}_bootstrap_se"] = standard_error
                    row[f"{tag}_gain"] = point - macro_over(
                        observed, arms[base], folds[rows],
                        present if tag == "macro_all_folds" else keep_folds
                    ) if tag != "pooled" else point - fast_spearman(observed, arms[base])
                    row[f"{tag}_gain_ci_low"] = gain_low
                    row[f"{tag}_gain_ci_high"] = gain_high
                    row[f"{tag}_gain_p"] = (
                        float("nan") if name == base else bootstrap_p(gain_samples)
                    )
                summary.append(row)
                for fold in present:
                    keep = rows[folds[rows] == fold]
                    per_fold_rows.append({
                        "mode": mode, "stratum_max_abs_distance_bp": stratum, "arm": name,
                        "evaluation_fold": fold, "variants": int(keep.size),
                        "blocks_1mb": int(len(set(blocks[keep].tolist()))),
                        "spearman": fast_spearman(beta[keep], ensemble[name][keep]),
                    })
                if name in ("endogenous_only", "plus_reporter_oof", "allele_identity_ridge"):
                    for seed in SEEDS:
                        values_seed = predictions[(mode, name)][seed][rows]
                        per_seed_rows.append({
                            "mode": mode, "stratum_max_abs_distance_bp": stratum, "arm": name,
                            "seed": seed,
                            "macro_all_folds": macro_over(observed, values_seed, folds[rows], present),
                            "macro_defined_folds": macro_over(observed, values_seed, folds[rows], keep_folds),
                            "pooled": fast_spearman(observed, values_seed),
                        })
            print(f"[{mode}/{stratum}] done", flush=True)

    frame = pd.DataFrame(summary)
    frame.to_csv(out / "tables" / "defined_contrast_summary.tsv", sep="\t", index=False)
    pd.DataFrame(per_fold_rows).to_csv(
        out / "tables" / "defined_contrast_per_fold.tsv", sep="\t", index=False
    )
    pd.DataFrame(per_seed_rows).to_csv(
        out / "tables" / "defined_contrast_per_seed.tsv", sep="\t", index=False
    )
    pd.DataFrame(audit).to_csv(
        out / "tables" / "defined_contrast_audit.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    if saved:
        np.savez_compressed(out / "tables" / "defined_contrast_bootstrap_samples.npz", **saved)

    primary = frame[
        (frame["mode"] == "no_constant_floor")
        & (frame.stratum_max_abs_distance_bp == PRIMARY_STRATUM)
        & (frame.arm == "plus_reporter_oof")
    ]
    if len(primary) != 1:
        raise TransferError("primary row is not unique")
    row = primary.iloc[0]
    family = dict(C1_FAMILY_P)
    family["e2_transfer_contrast"] = float(row.macro_all_folds_gain_p)
    order = sorted(family.items(), key=lambda item: item[1])
    holm = []
    previous = True
    for index, (name, value) in enumerate(order):
        threshold = 0.05 / (len(order) - index)
        rejected = bool(previous and value <= threshold)
        previous = rejected
        holm.append({"member": name, "raw_p": value, "holm_threshold": threshold,
                     "rank": index + 1, "rejected_at_fwer_0.05": rejected})
    pd.DataFrame(holm).to_csv(out / "tables" / "holm_family_final.tsv", sep="\t", index=False)

    receipt = {
        "schema_version": "agp-e2-defined-contrast-v1",
        "primary_mode": "no_constant_floor",
        "primary_stratum_max_abs_distance_bp": int(PRIMARY_STRATUM),
        "primary_gain_macro_all_folds": float(row.macro_all_folds_gain),
        "primary_gain_ci": [float(row.macro_all_folds_gain_ci_low),
                            float(row.macro_all_folds_gain_ci_high)],
        "primary_gain_bootstrap_p": float(row.macro_all_folds_gain_p),
        "primary_arm_b_macro": float(row.macro_all_folds_spearman),
        "e2_1_met": bool(
            row.macro_all_folds_gain < 0.02
            and row.macro_all_folds_gain_ci_low <= 0.0 <= row.macro_all_folds_gain_ci_high
        ),
        "holm": holm,
        "seeds": list(SEEDS),
        "bootstrap_seed": arguments.bootstrap_seed,
        "bootstrap_resamples": arguments.resamples,
    }
    (out / "tables" / "defined_contrast_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
