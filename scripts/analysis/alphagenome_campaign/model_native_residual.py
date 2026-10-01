#!/usr/bin/env python3
"""Training-local native-score, pooled-sequence and combined caQTL comparisons.

All regressions have zero intercept and use training RMS without centering.
Only folds 2--4 select recipes; fold 1 labels are loaded after states are saved.
Fold 0 is excluded. The chromosome bootstrap is conditional on the fitted
models and deposited association estimates, not biological-donor uncertainty.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import sys
import time

import numpy as np
import pandas as pd
import scipy
from scipy.linalg import eigh
from scipy.stats import rankdata

PROJ = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJ/"scripts/analysis/alphagenome_program"))
import i1_common as C  # numpy-only import; no JAX/model initialization

BASE = PROJ / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model"
LAMBDAS = (0.001, 0.01, 0.1, 1., 10., 100.)
LENGTHS = (2048, 16384)
POOLS = ("variant", "symmetric", "target")
FAMILIES = ("sequence", "residual", "joint")
TRAIN = (2, 3, 4)
RMS_MIN = 1e-8
SEED = 1103
IDENTITY = ["key", "chr", "pos_hg38", "ref", "alt", "peak_id",
            "peak_start_hg38", "peak_stop_hg38", "heldout_fold", "block_1mb"]
NATIVE_COLUMNS = ["local_atac_liver", "local_dnase_liver"]
ATAC_TRACKS = ["UBERON:0001114 ATAC-seq", "UBERON:0001115 ATAC-seq", "UBERON:0002107 ATAC-seq"]
DNASE_TRACKS = ["UBERON:0001114 DNase-seq", "UBERON:0001115 DNase-seq"]


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite(value, name):
    if not np.isfinite(value).all():
        raise ValueError(f"Nonfinite {name}")


def stats(x, z, y):
    """Sufficient statistics; callers pass only the permitted training rows."""
    x = np.asarray(x, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    for value, name in ((x, "features"), (z, "native scores"), (y, "training targets")):
        finite(value, name)
    if len(y) < 2 or z.shape != (len(y), 2) or x.shape[0] != len(y):
        raise ValueError("Invalid training dimensions")
    return {"xx": x.T @ x, "xz": x.T @ z, "xy": x.T @ y,
            "zz": z.T @ z, "zy": z.T @ y, "n": len(y)}


def add_stats(items):
    return {key: sum(item[key] for item in items) for key in items[0]}


def scale_from_gram(gram, n):
    diagonal = np.diag(gram)
    if (diagonal < -1e-12).any():
        raise ValueError("Negative feature squared norm")
    rms = np.sqrt(np.maximum(diagonal, 0) / n)
    admitted = rms > RMS_MIN
    scale = np.where(admitted, rms, 1.)
    return scale, admitted


def state(b, a, xs, xa, zs, za, metadata):
    b = np.asarray(b, dtype=float).copy()
    a = np.asarray(a, dtype=float).copy()
    b[~xa] = 0
    a[~za] = 0
    result = {"feature_coef_scaled": b, "native_coef_scaled": a,
              "feature_scale": xs.copy(), "native_scale": zs.copy(),
              "feature_admitted": xa.copy(), "native_admitted": za.copy(),
              "feature_coef_raw": b / xs, "native_coef_raw": a / zs,
              "metadata": dict(metadata, intercept=0., centering="none",
                               rms_threshold=RMS_MIN)}
    for key, value in result.items():
        if key != "metadata":
            finite(value, key)
    return result


def predict(model, x, z):
    z = np.asarray(z, dtype=np.float64)
    finite(z, "inference native scores")
    if z.shape[1] != 2:
        raise ValueError("Expected native ATAC,DNase columns")
    result = z @ model["native_coef_raw"]
    if len(model["feature_coef_raw"]):
        x = np.asarray(x)
        finite(x, "inference feature contrasts")
        if x.shape != (len(z), len(model["feature_coef_raw"])):
            raise ValueError("Inference feature dimensions differ from saved recipe")
        result = result + x @ model["feature_coef_raw"]
    finite(result, "predictions")
    return result


def save_state(path, model):
    arrays = {key: value for key, value in model.items() if key != "metadata"}
    arrays["metadata_json"] = np.asarray(json.dumps(model["metadata"], sort_keys=True))
    np.savez_compressed(path, **arrays)


def load_state(path):
    with np.load(path, allow_pickle=False) as archive:
        result = {key: archive[key] for key in archive.files if key != "metadata_json"}
        result["metadata"] = json.loads(str(archive["metadata_json"]))
    for key in ("feature", "native"):
        if not np.array_equal(result[key + "_coef_raw"],
                              result[key + "_coef_scaled"] / result[key + "_scale"]):
            raise ValueError("Saved preprocessing/coefficients disagree")
    return result


def native_candidates(z, y, metadata=None):
    """OLS and ridge use the exact training rows, including their RMS scales."""
    n = len(y)
    g = np.asarray(z, dtype=float).T @ z
    zy = np.asarray(z, dtype=float).T @ y
    zs, za = scale_from_gram(g, n)
    d = g / (n * np.outer(zs, zs))
    rhs = zy / (n * zs)
    d[~za, :] = 0
    d[:, ~za] = 0
    rhs[~za] = 0
    metadata = metadata or {}
    fitted = {}
    for name, dims, penalty in [("native_atac", [0], 0.),
                                ("native_dnase", [1], 0.),
                                ("native_ad_ols", [0, 1], 0.)] + [
            (f"native_ad_ridge:{lam:g}", [0, 1], lam) for lam in (0., *LAMBDAS)]:
        dims = np.asarray([i for i in dims if za[i]], dtype=int)
        a = np.zeros(2)
        if len(dims):
            small = d[np.ix_(dims, dims)] + penalty * np.eye(len(dims))
            a[dims] = np.linalg.pinv(small, rcond=1e-12) @ rhs[dims]
        fitted[name] = state(np.empty(0), a, np.empty(0), np.empty(0, bool), zs, za,
                             dict(metadata, family=name, penalty=penalty, train_rows=n))
    return fitted


def native_cv(z, y, folds):
    oof, receipts, saved = {}, [], []
    for held in TRAIN:
        train, test = np.isin(folds, [f for f in TRAIN if f != held]), folds == held
        candidates = native_candidates(z[train], y[train], {"inner_held_fold": held})
        for name, model in candidates.items():
            if name not in oof:
                oof[name] = np.full(len(y), np.nan)
            pred = predict(model, None, z[test])
            oof[name][test] = pred
            receipts.append({"candidate": name, "inner_held_fold": held,
                             "n": int(test.sum()), "SSE": float(np.sum((pred-y[test])**2))})
            saved.append((f"native_fold{held}_{name.replace(':', '_')}", model))
    mask = np.isin(folds, TRAIN)
    scores = {name: float(np.mean((value[mask]-y[mask])**2)) for name, value in oof.items()}
    selected = min(scores, key=scores.get)  # registered order resolves exact ties
    ridge = min((name for name in scores if name.startswith("native_ad_ridge:")), key=scores.get)
    return {"selected": selected, "ridge": ridge, "scores": scores,
            "oof": oof, "receipts": receipts, "states": saved}


def feature_path(sufficient, metadata=None):
    """Exact ridge penalty paths, sharing one feature Gram eigendecomposition.

    Joint ridge solves the 2x2 Schur complement for the native channels. This
    is algebraically identical to primal ridge with alpha=n*lambda, including
    its native-channel penalty. Residual ridge leaves native OLS unchanged.
    """
    s = sufficient
    n = s["n"]
    xs, xa = scale_from_gram(s["xx"], n)
    zs, za = scale_from_gram(s["zz"], n)
    idx = np.flatnonzero(xa)
    if not len(idx):
        raise ValueError("No nonzero training feature RMS")
    k = s["xx"][np.ix_(idx, idx)] / (n * np.outer(xs[idx], xs[idx]))
    if not np.allclose(k, k.T, rtol=1e-12, atol=1e-12):
        raise ValueError("Asymmetric feature Gram")
    eigen, u = eigh((k+k.T)/2, check_finite=True, driver="evd")
    tolerance = 1e-10 * max(1., float(np.max(np.abs(eigen))))
    minimum = float(eigen.min())
    if minimum < -tolerance:
        raise ValueError(f"Feature Gram not PSD: {minimum}, tolerance={tolerance}")
    clipped = int((eigen < 0).sum())
    eigen = np.maximum(eigen, 0)
    b = s["xz"][idx] / (n * xs[idx, None] * zs[None, :])
    gy = s["xy"][idx] / (n * xs[idx])
    d = s["zz"] / (n * np.outer(zs, zs))
    gz = s["zy"] / (n * zs)
    b[:, ~za] = 0
    d[~za, :] = 0
    d[:, ~za] = 0
    gz[~za] = 0
    a0 = np.linalg.pinv(d, rcond=1e-12) @ gz
    transformed = u.T @ np.column_stack([gy, b])
    result = {}
    for penalty in LAMBDAS:
        solved = u @ (transformed / (eigen[:, None] + penalty))
        avy, avb = solved[:, 0], solved[:, 1:]
        schur = d + penalty*np.eye(2) - b.T @ avb
        schur = (schur+schur.T)/2
        schur_eigen = np.linalg.eigvalsh(schur)
        # The Schur complement of full ridge is at least lambda*I. A much
        # smaller eigenvalue indicates numerical cancellation or bad inputs.
        if schur_eigen.min() < penalty * (1-1e-5):
            raise ValueError(f"Joint Schur positivity guard failed: {schur_eigen}")
        aj = np.linalg.solve(schur, gz - b.T @ avy)
        for family, active, a in (("sequence", avy, np.zeros(2)),
                                  ("residual", avy-avb@a0, a0),
                                  ("joint", avy-avb@aj, aj)):
            full = np.zeros(len(xs))
            full[idx] = active
            meta = dict(metadata or {}, family=family, penalty=penalty, train_rows=n,
                        admitted_feature_count=len(idx), eigen_minimum=minimum,
                        eigen_negative_clipped=clipped, eigen_tolerance=tolerance,
                        schur_minimum=float(schur_eigen.min()))
            result[(family, penalty)] = state(full, a, xs, xa, zs, za, meta)
    return result


def feature_cv(x, z, y, folds, length, pooling, out=None, after_first=None):
    started = time.monotonic()
    by_fold = {fold: stats(x[folds == fold], z[folds == fold], y[folds == fold]) for fold in TRAIN}
    stat_seconds = time.monotonic()-started
    records, oof = [], {}
    for held in TRAIN:
        tick = time.monotonic()
        train_stats = add_stats([by_fold[f] for f in TRAIN if f != held])
        path = feature_path(train_stats, {"length": length, "pooling": pooling,
                                         "inner_held_fold": held})
        test = folds == held
        matrix = np.column_stack([m["feature_coef_raw"] for m in path.values()])
        native = np.column_stack([m["native_coef_raw"] for m in path.values()])
        predicted = x[test] @ matrix + z[test] @ native
        for j, ((family, penalty), model) in enumerate(path.items()):
            oof.setdefault((family, penalty), np.full(len(y), np.nan))[test] = predicted[:, j]
            records.append({"length": length, "pooling": pooling, "family": family,
                            "penalty": penalty, "inner_held_fold": held, "n": int(test.sum()),
                            "SSE": float(np.sum((predicted[:, j]-y[test])**2)),
                            "admitted_features": int(model["feature_admitted"].sum())})
            if out:
                save_state(out/f"inner_{length}_{pooling}_{held}_{family}_{penalty:g}.npz", model)
        elapsed = time.monotonic()-tick
        if after_first and held == TRAIN[0]:
            after_first(stat_seconds, elapsed, len(path), x.shape[1])
    mask = np.isin(folds, TRAIN)
    candidates = [{"length": length, "pooling": pooling, "family": family,
                   "penalty": penalty, "inner_MSE": float(np.mean((pred[mask]-y[mask])**2))}
                  for (family, penalty), pred in oof.items()]
    return candidates, records, add_stats(list(by_fold.values()))


def choose_recipes(candidates):
    def choose(rows):
        return min(rows, key=lambda r: (r["inner_MSE"], -r["penalty"], r["length"], POOLS.index(r["pooling"])))
    selected = {}
    for family in FAMILIES:
        for length in LENGTHS:
            selected[f"{family}_{length}"] = choose([r for r in candidates if r["family"] == family and r["length"] == length])
        selected[f"{family}_selected"] = choose([r for r in candidates if r["family"] == family])
    return selected


def contrast_pairs():
    pairs = [(f"{family}_{length}", "native_selected") for family in FAMILIES for length in LENGTHS]
    pairs += [(f"{family}_selected", "native_selected") for family in FAMILIES]
    pairs += [("residual_selected", "native_atac"), ("joint_selected", "native_ad_ridge"),
              ("joint_2048", "residual_2048"), ("joint_16384", "residual_16384")]
    pairs += [(f"{family}_16384", f"{family}_2048") for family in FAMILIES]
    pairs += [("native_ad_ridge", "native_atac"), ("residual_selected", "sequence_selected")]
    if len(pairs) != 18 or len(set(pairs)) != 18:
        raise AssertionError("Exploratory contrast family changed")
    return pairs


def rho(y, prediction):
    a, b = rankdata(y), rankdata(prediction)
    a -= a.mean()
    b -= b.mean()
    denom = np.sqrt(a@a * (b@b))
    return float(a@b/denom) if denom > 0 else np.nan


def bootstrap(y, predictions, blocks, repeats, seed):
    """Use identical cluster draws for every model and every endpoint."""
    groups = [np.flatnonzero(blocks == b) for b in np.unique(blocks)]
    if len(groups) < 2:
        raise ValueError("Fewer than two uncertainty units")
    rng = np.random.default_rng(seed)
    mse = np.empty((repeats, predictions.shape[1]))
    ranks = np.empty_like(mse)
    for repeat in range(repeats):
        indices = np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))])
        yy, pp = y[indices], predictions[indices]
        mse[repeat] = np.mean((pp-yy[:, None])**2, axis=0)
        for j in range(pp.shape[1]):
            ranks[repeat, j] = rho(yy, pp[:, j])
    return mse, ranks


def bh(values):
    p = np.asarray(values)
    order = np.argsort(p)
    adjusted = np.empty(len(p))
    adjusted[order] = np.minimum(1., np.minimum.accumulate((p[order]*len(p)/np.arange(1, len(p)+1))[::-1])[::-1])
    return adjusted


def summarize(frame, names, margin, out, repeats):
    y = frame.beta_alt.to_numpy(float)
    p = frame[names].to_numpy(float)
    metrics = []
    for j, name in enumerate(names):
        mse = float(np.mean((p[:, j]-y)**2))
        metrics.append({"model": name, "n": len(y), "chromosomes": frame.chr.nunique(),
                        "historical_1Mb_bins": frame.block_1mb.nunique(), "MSE": mse,
                        "RMSE_beta_units": np.sqrt(mse), "MAE": float(np.mean(abs(p[:, j]-y))),
                        "signed_spearman": rho(y, p[:, j])})
    pd.DataFrame(metrics).to_csv(out/"performance.tsv", sep="\t", index=False)
    index = {name: j for j, name in enumerate(names)}
    for block_column, label in (("chr", "chromosome_primary"), ("block_1mb", "historical_bin_sensitivity")):
        mse, ranks = bootstrap(y, p, frame[block_column].to_numpy(str), repeats, SEED)
        np.savez_compressed(out/f"bootstrap_{label}.npz", names=np.asarray(names), MSE=mse, signed_spearman=ranks)
        rows = []
        for arm, comparator in contrast_pairs():
            ai, ci = index[arm], index[comparator]
            identical = np.array_equal(p[:, ai], p[:, ci])
            for endpoint, draws, point in (
                ("MSE_gain", mse[:, ci]-mse[:, ai], metrics[ci]["MSE"]-metrics[ai]["MSE"]),
                ("signed_spearman_gain", ranks[:, ai]-ranks[:, ci], metrics[ai]["signed_spearman"]-metrics[ci]["signed_spearman"])):
                defined = bool(np.isfinite(draws).all() and np.isfinite(point))
                lo, hi = np.quantile(draws, [0.025, 0.975]) if defined else (np.nan, np.nan)
                pvalue = min(1., 2*min((np.sum(draws <= 0)+1)/(repeats+1),
                                     (np.sum(draws >= 0)+1)/(repeats+1))) if defined else 1.
                disposition = "degenerate_identical_predictions" if identical else "uncertain"
                if endpoint == "MSE_gain" and not identical and defined:
                    disposition = "CI_above_operational_margin" if lo > margin else (
                        "conditional_CI_below_operational_margin" if hi < margin else "operational_margin_unresolved")
                rows.append({"arm": arm, "comparator": comparator, "endpoint": endpoint,
                             "gain": point, "CI_low": lo, "CI_high": hi, "nominal_p": pvalue,
                             "all_bootstrap_draws_defined": defined, "identical_predictions": identical,
                             "margin_beta_squared": margin if endpoint == "MSE_gain" else np.nan,
                             "relative_MSE_gain": point/metrics[ci]["MSE"] if endpoint == "MSE_gain" else np.nan,
                             "precision_disposition": disposition, "planned_family_n": 18})
        result = pd.DataFrame(rows)
        for endpoint in result.endpoint.unique():
            mask = result.endpoint == endpoint
            result.loc[mask, "BH_q_complete18"] = bh(result.loc[mask, "nominal_p"])
        result.to_csv(out/f"contrasts_{label}.tsv", sep="\t", index=False)


def read_labels(path, manifest, admitted_folds):
    """Skip all unrequested label strings before numeric conversion or use."""
    wanted = manifest.loc[manifest.heldout_fold.isin(admitted_folds)].set_index("key")
    result = pd.Series(np.nan, index=manifest.key)
    seen = set()
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row["key"]
            if key not in wanted.index:
                continue
            if key in seen:
                raise ValueError("Duplicate source label")
            seen.add(key)
            target = wanted.loc[key]
            for field in IDENTITY:
                if str(row[field]) != str(target[field] if field != "key" else key):
                    raise ValueError(f"Label identity mismatch: {key}, {field}")
            result.loc[key] = float(row["beta_alt"])
    if seen != set(wanted.index):
        raise ValueError("Source label identities missing")
    finite(result.loc[wanted.index].to_numpy(), "admitted source beta")
    return result.to_numpy()


def load_inputs(args):
    manifests = {length: pd.read_csv(args.features/str(length)/"manifest.tsv", sep="\t") for length in LENGTHS}
    manifest = manifests[2048]
    if manifest.key.duplicated().any() or not manifests[16384][IDENTITY].equals(manifest[IDENTITY]):
        raise ValueError("Feature identities differ or repeat")
    coverage = np.ones(len(manifest), bool)
    census, metadata = [], {}
    for length in LENGTHS:
        directory = args.features/str(length)
        meta = json.loads((directory/"complete.json").read_text())
        metadata[length] = meta
        if (meta["allele_order"] != ["REF", "ALT"] or not meta["complete"] or
            meta["card_tag"] != "l40s" or meta["outcomes_read"] or
            meta["manifest_sha256"] != sha256(directory/"manifest.tsv")):
            raise ValueError("Feature extraction identity/protocol mismatch")
        with np.load(directory/"coverage.npz", allow_pickle=False) as archive:
            for pool in POOLS:
                allowed = archive[pool]
                if allowed.dtype != bool or allowed.shape != (len(manifest),):
                    raise ValueError("Invalid coverage vector")
                coverage &= allowed
                census.append({"length": length, "pooling": pool, "total_rows": len(manifest),
                               "covered_rows": int(allowed.sum())})
    if metadata[2048]["checkpoint"] != metadata[16384]["checkpoint"]:
        raise ValueError("Feature checkpoints differ")
    run_path = args.native.with_name("native_1048576_run.json")
    native_meta = json.loads(run_path.read_text())
    reuse_path = PROJ/"GWAS/finemapping/results/alphagenome_program/ad1-gate-20260915T114615Z/raw/gate_1048576_run.json"
    reused_meta = json.loads(reuse_path.read_text())
    for meta in (native_meta, reused_meta):
        if (meta["checkpoint"] != metadata[2048]["checkpoint"] or meta["checkpoint"] != str(C.CHECKPOINT) or
            meta["fasta"] != C.FASTA_PATH or meta["atac_liver_track_names"] != ATAC_TRACKS or
            meta["dnase_liver_track_names"] != DNASE_TRACKS or meta["length_bp"] != 1048576):
            raise ValueError("Native checkpoint/reference/track correspondence differs from cached feature producer")
    metadata["native_protocol"] = {"current_run_sha256": sha256(run_path), "reused_run_sha256": sha256(reuse_path),
        "checkpoint": native_meta["checkpoint"], "fasta": native_meta["fasta"],
        "ATAC_tracks": ATAC_TRACKS, "DNase_tracks": DNASE_TRACKS,
        "feature_producer_sha256": sha256(Path(__file__).with_name("model_frozen.py")),
        "common_helper_sha256": sha256(C.__file__),
        "fasta_binding": "Cached extraction calls i1_common.FASTA_PATH; both native run records name that exact reference",
        "limitation": "Original feature receipt stores checkpoint path and producer reference route, not original content hashes"}
    # Never read the raw-score file's label column or any old fitted prediction.
    track_columns = ["local_atac__"+name.split()[0] for name in ATAC_TRACKS]+["local_dnase__"+name.split()[0] for name in DNASE_TRACKS]
    fields = ["key", "chr", "pos_hg38", "ref", "alt", "block_1mb", "heldout_fold",
              "length_bp", "mask_width_bp", "aggregation", "card_tag", *NATIVE_COLUMNS, *track_columns]
    native = pd.read_csv(args.native, sep="\t", usecols=fields)
    if native.key.duplicated().any():
        raise ValueError("Native variant identity repeats")
    if not ((native.length_bp == 1048576) & (native.mask_width_bp == 501) &
            (native.aggregation == "DIFF_LOG2_SUM") & (native.card_tag == "l40s")).all():
        raise ValueError("Native input/readout protocol changed")
    if not (np.allclose(native.local_atac_liver, native[track_columns[:3]].mean(axis=1), atol=1e-10, rtol=1e-10) and
            np.allclose(native.local_dnase_liver, native[track_columns[3:]].mean(axis=1), atol=1e-10, rtol=1e-10)):
        raise ValueError("Pooled native scores differ from registered liver track means")
    native = native.set_index("key")
    indexed = native.reindex(manifest.key)
    has_native = np.isfinite(indexed[NATIVE_COLUMNS].to_numpy(float)).all(axis=1)
    common = coverage & has_native
    chosen = common & manifest.heldout_fold.isin([1, *TRAIN]).to_numpy()
    frame = manifest.loc[chosen].copy().reset_index(drop=True)
    for field in ("chr", "pos_hg38", "ref", "alt", "block_1mb", "heldout_fold"):
        if not np.array_equal(frame[field].to_numpy(), native.loc[frame.key, field].to_numpy()):
            raise ValueError(f"Native/feature identity mismatch: {field}")
    if set(frame.heldout_fold) != {1, 2, 3, 4}:
        raise ValueError("Missing a required training/development fold")
    for field in ("chr", "block_1mb"):
        if frame.groupby(field).heldout_fold.nunique().max() != 1:
            raise ValueError("A chromosome/window group crosses folds")
    census.append({"stage": "six_feature_common", "covered_rows": int(coverage.sum())})
    census.append({"stage": "six_features_and_native_common_all_folds", "covered_rows": int(common.sum())})
    census += [{"stage": f"admitted_fold_{f}", "covered_rows": int((frame.heldout_fold == f).sum())} for f in (1, *TRAIN)]
    return frame, native.loc[frame.key, NATIVE_COLUMNS].to_numpy(float), np.flatnonzero(chosen), census, metadata


def load_feature(directory, indices):
    array = np.load(directory, mmap_mode="r", allow_pickle=False)
    if array.ndim != 3 or array.shape[1:] != (2, 3072) or array.dtype != np.float32:
        raise ValueError("Unexpected pooled feature shape/dtype")
    delta = np.asarray(array[indices, 1], np.float32) - np.asarray(array[indices, 0], np.float32)
    finite(delta, "pooled ALT-REF contrasts")
    return delta


def fit(args):
    started = time.monotonic()
    args.out.mkdir(parents=True, exist_ok=False)
    states_dir = args.out/"states"
    states_dir.mkdir()
    design = json.loads(args.design.read_text())
    if (tuple(design["ridge"]["lambda"]) != LAMBDAS or
        tuple(design["representation_candidates"]["length_bp"]) != LENGTHS or
        tuple(design["representation_candidates"]["pooling"]) != POOLS or
        tuple(design["folds"]["training"]) != TRAIN or
        design["folds"]["development_evaluation"] != 1 or
        design["exploratory_contrasts"]["family_size"] != len(contrast_pairs())):
        raise ValueError("Registered design differs from executable recipes")
    write_json(args.out/"registered_design.json", design)
    write_json(args.out/"environment.json", {"python": sys.version, "numpy": np.__version__,
        "scipy": scipy.__version__, "pandas": pd.__version__, "platform": platform.platform(),
        "job_id": os.environ.get("SLURM_JOB_ID"), "seed": SEED,
        "source_sha256": sha256(__file__), "design_sha256": sha256(args.design)})
    frame, z, indices, census, source_meta = load_inputs(args)
    folds = frame.heldout_fold.to_numpy(int)
    train, valid = np.isin(folds, TRAIN), folds == 1
    y = read_labels(args.labels, frame, TRAIN)  # development labels remain NaN
    pd.DataFrame(census).to_csv(args.out/"coverage.tsv", sep="\t", index=False)
    frame.to_csv(args.out/"input_population.tsv.gz", sep="\t", index=False)
    native = native_cv(z, y, folds)
    for name, model in native["states"]:
        save_state(states_dir/(name+".npz"), model)
    write_json(args.out/"native_inner_selection.json", {k: native[k] for k in ("selected", "ridge", "scores")})
    pd.DataFrame(native["receipts"]).to_csv(args.out/"native_inner_errors.tsv", sep="\t", index=False)
    # Benchmark the complete statistical workload using synthetic outcomes,
    # never the held development targets, before scientific fits continue.
    rng = np.random.default_rng(SEED)
    pilot_y = rng.normal(size=int(valid.sum()))
    pilot_p = rng.normal(size=(int(valid.sum()), 15))
    tick = time.monotonic()
    for col in ("chr", "block_1mb"):
        bootstrap(pilot_y, pilot_p, frame.loc[valid, col].to_numpy(str), 20, SEED)
    bootstrap_seconds = (time.monotonic()-tick)*args.resamples/20
    features, sufficient, candidates, inner_rows, loads = {}, {}, [], [], []

    def admit(stat_seconds, path_seconds, path_count, dimension):
        # Eigh is the expensive operation. There are 18 inner paths and at
        # most six outer paths; feature loading/statistics occur six times.
        projected = 2*(6*(loads[0]+stat_seconds) + 24*path_seconds + bootstrap_seconds + 90)
        # Six retained feature arrays and full training Grams, temporary
        # per-fold Grams, double copies and eigensolver workspaces, plus 2 GiB.
        memory = (6*len(frame)*dimension*4 + 9*dimension**2*8 +
                  4*len(frame)*dimension*8 + 2*1024**3)
        receipt = {"first_feature_load_seconds": loads[0], "three_fold_Gram_seconds": stat_seconds,
            "first_inner_full_path_seconds": path_seconds, "path_models": path_count,
            "bootstrap_both_schemes_projected_seconds": bootstrap_seconds,
            "projected_total_seconds_factor2": projected,
            "projected_allocated_core_hours": projected*args.cpus/3600,
            "projected_peak_GB": memory/1024**3, "wall_budget_seconds": args.max_seconds,
            "admitted": projected < args.max_seconds and memory < 28*1024**3,
            "development_labels_loaded": False}
        write_json(args.out/"throughput_admission.json", receipt)
        print(json.dumps(receipt), flush=True)
        if not receipt["admitted"]:
            raise RuntimeError("Measured full-comparison projection exceeds requested allocation")

    for length in LENGTHS:
        for pool in POOLS:
            if time.monotonic()-started > args.max_seconds:
                raise RuntimeError("Runtime guard; full comparison unfinished")
            tick = time.monotonic()
            x = load_feature(args.features/str(length)/(pool+".npy"), indices)
            loads.append(time.monotonic()-tick)
            features[length, pool] = x
            cc, rr, ss = feature_cv(x, z, y, folds, length, pool, states_dir,
                                   admit if len(features) == 1 else None)
            candidates.extend(cc)
            inner_rows.extend(rr)
            sufficient[length, pool] = ss
            pd.DataFrame(inner_rows).to_csv(args.out/"sequence_inner_errors.tsv", sep="\t", index=False)
            print(json.dumps({"completed_representation": [length, pool],
                              "elapsed_seconds": time.monotonic()-started}), flush=True)
    selected = choose_recipes(candidates)
    margin = 0.05*native["scores"][native["selected"]]
    write_json(args.out/"selection_before_development_labels.json", {
        "selected_recipes": selected, "native_selected": native["selected"], "native_ridge": native["ridge"],
        "margin_beta_squared": margin, "margin_definition": "5% selected native inner OOF MSE",
        "development_labels_loaded": False, "training_folds": list(TRAIN), "evaluation_fold": 1,
        "candidate_count": len(candidates), "native_candidates": len(native["scores"]),
        "limits": design["limits"], "input_geometry": design["representation_candidates"]["source_geometry"]})
    pd.DataFrame(candidates).to_csv(args.out/"all_candidate_inner_MSE.tsv", sep="\t", index=False)
    outer_native = native_candidates(z[train], y[train], {"training_folds": list(TRAIN)})
    zero = outer_native["native_ad_ols"].copy()
    zero = dict(zero, native_coef_raw=np.zeros(2), native_coef_scaled=np.zeros(2),
                metadata=dict(zero["metadata"], family="zero"))
    fitted = {"zero": zero, "native_atac": outer_native["native_atac"],
              "native_dnase": outer_native["native_dnase"], "native_ad_ols": outer_native["native_ad_ols"],
              "native_ad_ridge": outer_native[native["ridge"]], "native_selected": outer_native[native["selected"]]}
    outer_paths = {}
    for name, recipe in selected.items():
        key = recipe["length"], recipe["pooling"]
        if key not in outer_paths:
            if time.monotonic()-started > args.max_seconds:
                raise RuntimeError("Runtime guard before selected outer refit")
            outer_paths[key] = feature_path(sufficient[key], {"length": key[0], "pooling": key[1], "training_folds": list(TRAIN)})
        fitted[name] = outer_paths[key][recipe["family"], recipe["penalty"]]
    predicted, reload_checks = {}, []
    for name, model in fitted.items():
        model["metadata"].update(native_length=1048576, native_readout_bp=501,
                                 native_score_order=NATIVE_COLUMNS, checkpoint=source_meta[2048]["checkpoint"],
                                 native_aggregation="DIFF_LOG2_SUM", fasta=C.FASTA_PATH,
                                 allele_order=["REF", "ALT"], beta_units="source FastQTL ALT dosage slope")
        path = states_dir/f"final_{name}.npz"
        save_state(path, model)
        x = features[model["metadata"]["length"], model["metadata"]["pooling"]][valid] if name in selected else None
        pred = predict(model, x, z[valid])
        restored = predict(load_state(path), x, z[valid])
        if not np.array_equal(pred, restored):
            raise ValueError("Final saved-state reload changed predictions")
        if not np.array_equal(-pred, predict(model, -x if x is not None else None, -z[valid])):
            raise ValueError("Final model allele-swap invariant failed")
        zeros = predict(model, np.zeros_like(x) if x is not None else None, np.zeros_like(z[valid]))
        if np.count_nonzero(zeros):
            raise ValueError("Final model zero-contrast invariant failed")
        predicted[name] = pred
        reload_checks.append({"model": name, "reload_max_abs": 0., "swap_max_abs": 0., "zero_max_abs": 0.})
    write_json(args.out/"actual_model_invariants.json", reload_checks)
    # First numeric access to fold-1 outcomes occurs only here, after every
    # choice, complete fitted state, margin and held-input prediction is fixed.
    development = read_labels(args.labels, frame, [1])[valid]
    result = frame.loc[valid].copy()
    result["beta_alt"] = development
    for name, pred in predicted.items():
        result[name] = pred
    result.to_csv(args.out/"development_predictions.tsv.gz", sep="\t", index=False)
    summarize(result, list(predicted), margin, args.out, args.resamples)
    write_json(args.out/"completion.json", {"complete": True, "models": len(predicted),
        "planned_contrasts_per_endpoint_per_scheme": len(contrast_pairs()), "bootstrap_resamples": args.resamples,
        "training_rows": int(train.sum()), "development_rows": int(valid.sum()),
        "training_chromosomes": int(frame.loc[train, "chr"].nunique()),
        "development_chromosomes": int(frame.loc[valid, "chr"].nunique()), "fold0_rows_loaded_as_labels": 0,
        "elapsed_seconds": time.monotonic()-started, "allocated_core_hours": args.cpus*(time.monotonic()-started)/3600,
        "peak_host_GB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,
        "source_native_sha256": sha256(args.native), "source_labels_sha256": sha256(args.labels),
        "source_metadata": source_meta, "conditional_uncertainty": design["exploratory_contrasts"]["uncertainty"],
        "evidence": "nested training selection, single previously used development fold; not protected evaluation",
        "biological_n": "Not established by variant, chromosome or historical-bin counts"})


def inference(args):
    model = load_state(args.model)
    with np.load(args.inputs, allow_pickle=False) as source:
        if (not np.array_equal(source["native_score_order"], model["metadata"]["native_score_order"]) or
            int(source["native_length"]) != model["metadata"]["native_length"] or
            int(source["native_readout_bp"]) != model["metadata"]["native_readout_bp"] or
            str(source["native_aggregation"]) != model["metadata"]["native_aggregation"] or
            str(source["checkpoint"]) != model["metadata"]["checkpoint"] or
            str(source["fasta"]) != model["metadata"]["fasta"]):
            raise ValueError("Inference native channel order, geometry, checkpoint or reference differs")
        native = source["native_atac_dnase"]
        if len(model["feature_coef_raw"]):
            if not np.array_equal(source["allele_order"], ["REF", "ALT"]):
                raise ValueError("Inference allele order mismatch")
            if int(source["length"]) != model["metadata"]["length"] or str(source["pooling"]) != model["metadata"]["pooling"]:
                raise ValueError("Inference pooling/length differs from fitted recipe")
            x = np.asarray(source["alt"], np.float32)-np.asarray(source["ref"], np.float32)
        else:
            x = None
        predictions = predict(model, x, native)
        keys = source["key"]
    if args.out.exists():
        raise FileExistsError(args.out)
    pd.DataFrame({"key": keys, "predicted_beta_alt": predictions}).to_csv(args.out, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    training = sub.add_parser("fit")
    training.add_argument("--features", type=Path, default=BASE/"frozen")
    training.add_argument("--native", type=Path, default=BASE/"native/native_1048576.tsv")
    training.add_argument("--labels", type=Path, default=BASE/"native/currin_variants.tsv")
    training.add_argument("--design", type=Path, default=Path(__file__).with_name("model_native_residual_design.json"))
    training.add_argument("--out", type=Path, required=True)
    training.add_argument("--cpus", type=int, default=4)
    training.add_argument("--max-seconds", type=float, default=6900)
    training.add_argument("--resamples", type=int, default=2000)
    prediction = sub.add_parser("predict")
    prediction.add_argument("--model", type=Path, required=True)
    prediction.add_argument("--inputs", type=Path, required=True)
    prediction.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "predict":
        inference(args)
        return
    if args.out.exists():
        raise FileExistsError(args.out)
    try:
        fit(args)
    except Exception as exc:
        if args.out.exists() and not (args.out/"completion.json").exists():
            write_json(args.out/"failure.json", {"complete": False, "error": repr(exc),
                       "disposition": "Retain failure; no incomplete-arm performance interpretation"})
        raise


if __name__ == "__main__":
    main()
