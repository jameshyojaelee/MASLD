#!/usr/bin/env python3
"""Fit a fixed regularized combination from fold-excluding Currin predictions.

Fold 0 is excluded before outcome parsing. Four development folds are already
inspected; this analysis cannot supply a new independent confirmation.
"""
import argparse
import hashlib
import json
import os
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model"
NATIVE_LONG = BASE / "native/native_1048576.tsv"
NATIVE_SHORT = BASE / "native_lengths_21775189/native_2048.tsv"
SEEDS = range(20260917, 20260922)
FOLDS = range(1, 5)
RECIPE = "adapter_r16_last5_symmetric"
FIXED = "frozen_2048_symmetric_shared_mlp64"
PENALTY = 0.01  # fixed before these two-fold fits; no scored-fold tuning
SEED = 20260922


def native(path):
    # Scan only fold membership first. The native 1-Mb source also contains the
    # spent fold, whose outcome fields must not enter the analysis dataframe.
    membership = pd.read_csv(path, sep="\t", usecols=["heldout_fold"])
    if not membership.heldout_fold.isin(range(5)).all():
        raise ValueError("Unexpected source fold membership")
    skip = np.flatnonzero(~membership.heldout_fold.isin(FOLDS).to_numpy()) + 1
    frame = pd.read_csv(path, sep="\t", usecols=["key", "chr", "ref", "alt", "block_1mb",
                                                 "heldout_fold", "beta_alt", "local_atac_liver"],
                        skiprows=skip)
    if frame.key.duplicated().any() or not frame.heldout_fold.isin(FOLDS).all():
        raise ValueError("Duplicate key or forbidden fold in native source")
    if not np.isfinite(frame[["beta_alt", "local_atac_liver"]].to_numpy(float)).all():
        raise ValueError("Nonfinite native outcome or prediction")
    return frame.set_index("key")


def ensemble(folder, recipe, expected_fold, expected_training):
    frames = []
    for seed in SEEDS:
        name = f"{recipe}_seed{seed}" if recipe.startswith("outer") else f"{recipe}__seed{seed}"
        path = folder/name
        split = json.loads((path/"split.json").read_text())
        if split["validation_fold"] != expected_fold or sorted(split["training_folds"]) != expected_training:
            raise ValueError(f"Training exposure mismatch: {path}")
        if split["seed"] != seed:
            raise ValueError(f"Seed identity differs from directory: {path}")
        frame = pd.read_csv(path/"validation_predictions.tsv", sep="\t")
        if frame.variant_id.duplicated().any() or set(frame.validation_fold.unique()) != {expected_fold}:
            raise ValueError(f"Variant identity or validation fold mismatch: {path}")
        frame["key"] = frame.variant_id.str.removeprefix("chr")
        frames.append(frame.set_index("key")[["predicted_beta", "observed_beta"]].rename(
            columns={"predicted_beta": f"seed{seed}", "observed_beta": f"observed{seed}"}))
    if any(set(frame.index) != set(frames[0].index) for frame in frames[1:]):
        raise ValueError("Seed prediction coverage differs; an inner join would hide missing variants")
    joined = pd.concat([frame.reindex(frames[0].index) for frame in frames], axis=1)
    if joined.empty or joined.index.has_duplicates:
        raise ValueError("Missing common seed predictions")
    if not np.isfinite(joined.to_numpy(float)).all():
        raise ValueError("Nonfinite seed outcome or prediction")
    for seed in SEEDS:
        if not np.allclose(joined[f"observed{seed}"], joined[f"observed{SEEDS.start}"], atol=1e-6, rtol=1e-6):
            raise ValueError("Seed labels differ")
    return pd.DataFrame({"adapter": joined[[f"seed{s}" for s in SEEDS]].mean(1),
                         "beta_alt": joined[f"observed{SEEDS.start}"]})


def calibrated_native(frame, train_folds, target, column="local_atac_liver"):
    train = frame.loc[frame.heldout_fold.isin(train_folds)]
    x = train[column].to_numpy(float)
    y = train.beta_alt.to_numpy(float)
    slope = np.dot(x, y)/np.dot(x, x)
    if not np.isfinite(slope):
        raise ValueError("Native calibration undefined")
    return target[column].to_numpy(float)*slope, float(slope)


def stack_fit(native_score, adapter, y):
    matrix = np.column_stack([native_score, adapter])
    scale = np.sqrt(np.mean(matrix**2, axis=0))
    if np.any(scale < 1e-10):
        raise ValueError("One input predictor has zero scale")
    z = matrix/scale
    gram = z.T@z
    penalty = PENALTY*np.linalg.eigvalsh(gram).max()
    weight = np.linalg.solve(gram + penalty*np.eye(2), z.T@y)/scale
    return weight


def paired_group_bootstrap(frame, simple, complex_name, group, draws=3000):
    loss_simple = (frame.beta_alt.to_numpy()-frame[simple].to_numpy())**2
    loss_complex = (frame.beta_alt.to_numpy()-frame[complex_name].to_numpy())**2
    diff = loss_simple-loss_complex
    groups = [part.index.to_numpy() for _, part in frame.reset_index(drop=True).groupby(group, sort=True)]
    rng = np.random.default_rng(SEED + 1)
    bs = np.empty(draws)
    for j in range(draws):
        idx = np.concatenate([groups[k] for k in rng.integers(0, len(groups), len(groups))])
        bs[j] = diff[idx].mean()
    obs = float(diff.mean())
    centered = bs-bs.mean()
    p = (1+np.sum(np.abs(centered) >= abs(obs)))/(draws+1)
    return {"MSE_reduction": obs, "CI95": np.quantile(bs, [0.025, 0.975]).tolist(),
            "p_nominal_centred_bootstrap": float(p), "resampling_unit": group,
            "n_groups": len(groups), "draws": draws, "seed": SEED+1}


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute-node allocation required")
    if args.out.exists():
        raise FileExistsError(args.out)
    if json.loads((args.inputs/"completion.json").read_text())["fits"] != 60:
        raise ValueError("Two-fold nested input fits incomplete")
    args.out.mkdir(parents=True)
    long, short = native(NATIVE_LONG), native(NATIVE_SHORT)
    rows, calibration, fit_weights, coverage, meta_coverage = [], [], [], [], []
    for outer in FOLDS:
        meta_parts = []
        for meta in FOLDS:
            if outer == meta:
                continue
            expected_training = sorted(set(FOLDS)-{outer, meta})
            ens = ensemble(args.inputs, f"outer{outer}_meta{meta}", meta, expected_training)
            keys = ens.index.intersection(long.index)
            keys = keys[long.loc[keys, "heldout_fold"].eq(meta)]
            native_keys = long.index[long.heldout_fold.eq(meta)]
            meta_coverage.append({"outer": outer, "meta": meta,
                                  "adapter_rows": len(ens), "native1m_rows": len(native_keys),
                                  "shared_rows": len(keys),
                                  "adapter_rows_without_native1m": len(ens.index.difference(keys)),
                                  "native1m_rows_without_adapter": len(native_keys.difference(keys))})
            n1 = long.loc[keys]
            pred, slope = calibrated_native(long, expected_training, n1)
            part = ens.loc[keys].copy()
            if not np.allclose(part.beta_alt, n1.beta_alt, atol=1e-6, rtol=1e-6):
                raise ValueError("Meta labels disagree with native source")
            part["native1m"] = pred
            part["fold"] = meta
            meta_parts.append(part)
            calibration.append({"outer": outer, "meta": meta, "fit_folds": expected_training,
                                "native1m_slope": slope, "calibration_n": int(long.heldout_fold.isin(expected_training).sum())})
        meta_frame = pd.concat(meta_parts)
        weight = stack_fit(meta_frame.native1m, meta_frame.adapter, meta_frame.beta_alt)
        outer_training = sorted(set(FOLDS)-{outer})
        outer_dir = BASE/f"nested_f{outer}_{21783074+outer}"
        ens = ensemble(outer_dir, RECIPE, outer, outer_training)
        fixed = ensemble(outer_dir, FIXED, outer, outer_training).rename(columns={"adapter": "fixed_head"})
        if set(ens.index) != set(fixed.index):
            raise ValueError("Matched 2-kb adapter and fixed-head coverage differs")
        keys = ens.index.intersection(fixed.index).intersection(long.index).intersection(short.index)
        keys = keys[long.loc[keys, "heldout_fold"].eq(outer)]
        if not len(keys):
            raise ValueError("No common outer-fold variants")
        l, s = long.loc[keys], short.loc[keys]
        identity = ["chr", "ref", "alt", "heldout_fold", "block_1mb"]
        if not np.array_equal(l[identity].to_numpy(), s[identity].to_numpy()):
            raise ValueError("Native variant, fold or genomic grouping differs by context length")
        if not np.allclose(l.beta_alt, s.beta_alt, atol=1e-8, rtol=1e-8):
            raise ValueError("Native source effects differ by input length")
        n1, slope1 = calibrated_native(long, outer_training, l)
        n2, slope2 = calibrated_native(short, outer_training, s)
        part = ens.loc[keys].copy()
        if not np.allclose(part.beta_alt, l.beta_alt, atol=1e-6, rtol=1e-6):
            raise ValueError("Adapter labels differ from source effect")
        if not np.allclose(fixed.loc[keys, "beta_alt"], l.beta_alt, atol=1e-6, rtol=1e-6):
            raise ValueError("Fixed-head labels differ from source effect")
        part["fixed_head"] = fixed.loc[keys, "fixed_head"]
        part["native1m"] = n1
        part["native2k"] = n2
        part["combined"] = np.column_stack([n1, part.adapter.to_numpy()])@weight
        part["chr"] = l.chr.to_numpy()
        part["historical_bin"] = l.block_1mb.to_numpy()
        part["fold"] = outer
        part["ref"] = l.ref.to_numpy()
        part["alt"] = l.alt.to_numpy()
        rows.append(part.reset_index(names="key"))
        fit_weights.append({"held_fold": outer, "meta_folds": [m for m in FOLDS if m != outer],
                            "native1m_weight": float(weight[0]), "adapter_weight": float(weight[1]),
                            "penalty_fraction": PENALTY, "native1m_slope": slope1, "native2k_slope": slope2})
        coverage.append({"fold": outer, "shared_rows": len(part),
                         "adapter_rows": len(ens), "fixed_head_rows": len(fixed),
                         "adapter_rows_outside_shared_population": len(ens.index.difference(keys)),
                         "source_native1m_rows": int(long.heldout_fold.eq(outer).sum()),
                         "source_native2k_rows": int(short.heldout_fold.eq(outer).sum()),
                         "native1m_rows_outside_shared_population": int(long.heldout_fold.eq(outer).sum())-len(part),
                         "native2k_rows_outside_shared_population": int(short.heldout_fold.eq(outer).sum())-len(part)})
    all_rows = pd.concat(rows, ignore_index=True)
    if all_rows.key.duplicated().any() or all_rows.fold.eq(0).any():
        raise ValueError("Duplicate variant or fold 0 in evaluation")
    model_names = ["native2k", "fixed_head", "adapter", "native1m", "combined"]
    result = {"evidence_state": "inspected Currin development folds 1-4; no external source",
              "n_variants": len(all_rows), "n_chromosomes": int(all_rows.chr.nunique()),
              "n_historical_bins": int(all_rows.historical_bin.nunique()),
              "endpoint": "signed Currin FastQTL ALT dosage beta",
              "uncertainty_scope": "paired genomic resampling conditional on the saved prediction tables and measured source effects; excludes refitting, recipe selection, fitted-seed variation and runtime numerical variation",
              "fold0_rows_evaluated": 0, "fold0_labels_used": 0,
              "native_source_loading": "metadata-only fold scan; fold 0 skipped before outcome parsing",
              "fixed_penalty_fraction": PENALTY,
              "coverage": coverage, "meta_coverage": meta_coverage,
              "weights": fit_weights, "native_meta_calibration": calibration,
              "MSE": {name: float(np.mean((all_rows.beta_alt-all_rows[name])**2)) for name in model_names},
              "contrasts": {}, "environment": {"python": sys.version, "numpy": np.__version__,
                    "pandas": pd.__version__, "platform": platform.platform(),
                    "slurm_job_id": os.environ["SLURM_JOB_ID"]},
              "inputs_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in (args.inputs/"completion.json", NATIVE_LONG, NATIVE_SHORT)}}
    pairs = [("native1m", "combined"), ("adapter", "combined")]
    for simple, complex_name in pairs:
        result["contrasts"][f"{complex_name}_vs_{simple}"] = {
            "chromosome": paired_group_bootstrap(all_rows, simple, complex_name, "chr"),
            "historical_bin_sensitivity": paired_group_bootstrap(all_rows, simple, complex_name, "historical_bin")}
    p = [result["contrasts"][f"combined_vs_{s}"]["chromosome"]["p_nominal_centred_bootstrap"] for s,_ in pairs]
    # Two predeclared comparisons, Holm family-wise correction.
    order = np.argsort(p)
    q = [1.0, 1.0]
    q[int(order[0])] = min(1, 2*p[int(order[0])])
    q[int(order[1])] = max(q[int(order[0])], p[int(order[1])])
    for (simple, _), corrected in zip(pairs, q):
        result["contrasts"][f"combined_vs_{simple}"]["p_holm_primary_family"] = corrected
    best_simple = min(result["MSE"]["native1m"], result["MSE"]["adapter"])
    result["combined_relative_MSE_gain_over_best_simple"] = (best_simple-result["MSE"]["combined"])/best_simple
    result["development_complexity_margin_met"] = result["combined_relative_MSE_gain_over_best_simple"] >= .05
    all_rows.to_csv(args.out/"matched_development_predictions.tsv.gz", sep="\t", index=False)
    (args.out/"results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in ("n_variants", "n_chromosomes", "MSE", "contrasts",
                                             "combined_relative_MSE_gain_over_best_simple",
                                             "development_complexity_margin_met")}, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())
