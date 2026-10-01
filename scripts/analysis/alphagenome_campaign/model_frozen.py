#!/usr/bin/env python3
"""Cache three allele-pooling views and fit 12 matched development recipes.

All 12 recipes are evaluated on exactly the same eligible rows. Coverage for
every length/pooling is separate. No representation fitting uses held data.
One seed and fixed hyperparameters are development comparisons only; these
outputs do not replace nested five-seed finalist evaluation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

PROJ = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJ/"scripts/analysis/alphagenome_program"))
import i1_common as C

POOLINGS = ("variant", "symmetric", "target")


def extract(args):
    import jax
    import jax.numpy as jnp
    import pysam
    from alphagenome_research.model import one_hot_encoder
    from i1_extract_mpra_embeddings import load_trunk
    import model_scalar as M
    out = args.out/str(args.length)
    out.mkdir(parents=True, exist_ok=False)
    manifest = pd.read_csv(args.manifest, sep="\t")
    if manifest.key.duplicated().any():
        raise ValueError("Manifest keys are not unique")
    if any(c in manifest for c in ("beta_alt", "p_val", "p_nominal")):
        raise ValueError("Extraction manifest must omit outcomes")
    if args.limit:
        manifest = manifest.iloc[:args.limit].copy()
    manifest.to_csv(out/"manifest.tsv", sep="\t", index=False)
    n = len(manifest)
    arrays = {p: np.lib.format.open_memmap(out/(p+".npy"), mode="w+", dtype="float32", shape=(n, 2, 3072)) for p in POOLINGS}
    for array in arrays.values():
        array[:] = np.nan
    coverage = {p: np.zeros(n, dtype=bool) for p in POOLINGS}
    fasta = pysam.FastaFile(C.FASTA_PATH)
    encoder = one_hot_encoder.DNAOneHotEncoder()
    card = C.gpu_card()
    if card.get("card_tag") != "l40s":
        raise ValueError("L40S required for comparable frozen representations")
    trunk, base, state, device, load_seconds = load_trunk()
    organism = jnp.zeros((1,), dtype=jnp.int32)
    elapsed, rejected = [], []
    for i, row in enumerate(manifest.itertuples(index=False)):
        lo, hi = C.window_bounds(row.pos_hg38, args.length)
        check = C.validate_window(fasta, row.chr, lo, hi, row.pos_hg38, row.ref)
        if not check["acgt_ok"] or not check["reference_match"]:
            rejected.append({"key": row.key, "pooling": "all", "reason": "reference_or_alphabet"})
            continue
        weights = {}
        for pooling in POOLINGS:
            try:
                weights[pooling] = M.pool_weights(args.length, pooling, target_start=row.peak_start_hg38,
                                                  target_end=row.peak_stop_hg38, window_start=lo)
            except ValueError as exc:
                rejected.append({"key": row.key, "pooling": pooling, "reason": str(exc)})
        sequences = C.extract_ref_alt(fasta, row.chr, lo, hi, row.pos_hg38, row.ref, row.alt)
        tick = time.monotonic()
        for allele, seq in enumerate(sequences):
            x = jnp.asarray(encoder.encode(seq)[None], dtype=jnp.float32)
            e = np.asarray(trunk(base, state, x, organism).get_sequence_embeddings(128), dtype=np.float32)[0]
            if e.shape != (args.length//128, 3072) or not np.isfinite(e).all():
                raise ValueError("Unexpected or nonfinite trunk representation")
            for pooling, w in weights.items():
                arrays[pooling][i, allele] = np.sum(e*w[:, None], axis=0, dtype=np.float32)
        for pooling in weights:
            coverage[pooling][i] = True
        elapsed.append(time.monotonic()-tick)
        if (i+1) % 100 == 0 or i == 0:
            for array in arrays.values():
                array.flush()
            np.savez(out/"coverage.npz", **coverage)
            print(json.dumps({"length": args.length, "completed": i+1, "total": n,
                  "seconds_per_variant": float(np.median(elapsed[1:] or elapsed))}), flush=True)
    for array in arrays.values():
        array.flush()
    np.savez(out/"coverage.npz", **coverage)
    pd.DataFrame(rejected, columns=["key", "pooling", "reason"]).to_csv(out/"excluded.tsv", sep="\t", index=False)
    receipt = {"length": args.length, "rows": n, "coverage": {p:int(v.sum()) for p,v in coverage.items()},
        "mean_warm_seconds_per_variant": float(np.mean(elapsed[1:])) if len(elapsed)>1 else None,
        "load_seconds": load_seconds, "allele_order": ["REF", "ALT"], "poolings": list(POOLINGS),
        "representation": "local_AlphaGenome_embeddings_128bp", "checkpoint": str(C.CHECKPOINT),
        "manifest_sha256": hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
        "measured_target_policy": "full_interval_coverage_required_overlap_weighted_bins",
        "outcomes_read": False, "complete": True, **card}
    (out/"complete.json").write_text(json.dumps(receipt, indent=2)+"\n")


def fit_shared(ref_train, alt_train, y_train, ref_test, alt_test, recipe):
    import jax
    import jax.numpy as jnp
    import optax
    import model_scalar as M
    t, _ = M.initialize({}, mode="frozen", hidden=64, width=ref_train.shape[1], seed=recipe["seed"])
    optimizer = optax.chain(optax.clip_by_global_norm(1), optax.adamw(3e-4, weight_decay=1e-4))
    state = optimizer.init(t)
    @jax.jit
    def update(weights, opt_state, ref, alt, y):
        loss, gradients = jax.value_and_grad(lambda w:M.effect_loss(M.effect_from_features(w, ref, alt), y))(weights)
        changes, opt_state = optimizer.update(gradients, opt_state, weights)
        return optax.apply_updates(weights, changes), opt_state, loss
    rng = np.random.default_rng(recipe["seed"])
    for _ in range(recipe["steps"]):
        idx = rng.choice(len(y_train), size=min(64, len(y_train)), replace=False)
        t, state, value = update(t, state, jnp.asarray(ref_train[idx]), jnp.asarray(alt_train[idx]), jnp.asarray(y_train[idx]))
        if not np.isfinite(float(value)):
            raise ValueError("Nonfinite frozen MLP loss")
    # Predict in bounded batches to avoid a large transient device allocation.
    return np.concatenate([np.asarray(M.effect_from_features(t, ref_test[i:i+128], alt_test[i:i+128]))
                           for i in range(0, len(ref_test), 128)])


def fit(args):
    from sklearn.linear_model import Ridge
    out = args.out/"comparisons"
    out.mkdir(parents=True, exist_ok=False)
    labels = pd.read_csv(args.labels, sep="\t").set_index("lead_variant_id")
    manifests = {l: pd.read_csv(args.out/str(l)/"manifest.tsv", sep="\t") for l in (2048, 16384)}
    if not np.array_equal(manifests[2048].key, manifests[16384].key):
        raise ValueError("Length manifests differ in identities/order")
    matched = np.ones(len(manifests[2048]), dtype=bool)
    coverage_rows = []
    for length in (2048, 16384):
        meta = json.loads((args.out/str(length)/"complete.json").read_text())
        if not meta["complete"] or meta["card_tag"] != "l40s":
            raise ValueError("Extraction incomplete or hardware differs")
        with np.load(args.out/str(length)/"coverage.npz", allow_pickle=False) as cover:
            for pooling in POOLINGS:
                matched &= cover[pooling]
                coverage_rows.append({"length": length, "pooling": pooling, "eligible_rows": int(cover[pooling].sum()), "total_rows": len(matched)})
    manifest = manifests[2048].loc[matched].copy()
    y = labels.loc["chr"+manifest.key, "beta_alt"].to_numpy(dtype=np.float32)
    folds = manifest.heldout_fold.to_numpy()
    if not len(y) or len(set(folds)) != 5:
        raise ValueError("Matched population does not support all five chromosome folds")
    pd.DataFrame(coverage_rows).to_csv(out/"length_specific_coverage.tsv", sep="\t", index=False)
    recipes = [json.loads(p.read_text()) for p in sorted(args.configs.glob("frozen_*.json"))]
    if len(recipes) != 12:
        raise ValueError("Expected 12 frozen recipes")
    predictions = manifest[["key", "block_1mb", "heldout_fold"]].copy()
    predictions["beta_alt"] = y
    metrics, failures = [], []
    for recipe in recipes:
        tick = time.monotonic()
        arrays = np.load(args.out/str(recipe["length"])/(recipe["pooling"]+".npy"), mmap_mode="r")
        features = np.asarray(arrays[matched], dtype=np.float32)
        pred = np.full(len(y), np.nan, dtype=np.float32)
        for fold in range(5):
            train, test = folds != fold, folds == fold
            try:
                if recipe["head"] == "ridge":
                    delta = features[:,1]-features[:,0]
                    scale = delta[train].std(axis=0)
                    scale[scale < 1e-6] = 1.0
                    ridge = Ridge(alpha=recipe["alpha"], fit_intercept=False, solver="lsqr", tol=1e-6)
                    ridge.fit(delta[train]/scale, y[train])
                    pred[test] = ridge.predict(delta[test]/scale)
                elif recipe["head"] == "shared_mlp64":
                    pred[test] = fit_shared(features[train,0], features[train,1], y[train],
                                            features[test,0], features[test,1], recipe)
                else:
                    raise ValueError("Unsupported head")
            except Exception as exc:
                failures.append({"recipe": recipe["id"], "fold": fold, "error": repr(exc)})
        predictions[recipe["id"]] = pred
        okay = np.isfinite(pred).all()
        metrics.append({"recipe": recipe["id"], "matched_rows": len(y), "blocks": manifest.block_1mb.nunique(),
            "macro_spearman": C.macro_over_folds(y, pred, folds) if okay else np.nan,
            "RMSE_beta_units": float(np.sqrt(np.mean((pred-y)**2))) if okay else np.nan,
            "seconds": time.monotonic()-tick, "failed_folds": sum(f["recipe"] == recipe["id"] for f in failures),
            "evidence": "fixed_recipe_single_seed_development_cross_validation"})
        predictions.to_csv(out/"predictions.tsv.gz", sep="\t", index=False)
        pd.DataFrame(metrics).to_csv(out/"performance.tsv", sep="\t", index=False)
        pd.DataFrame(failures, columns=["recipe", "fold", "error"]).to_csv(out/"failed_fits.tsv", sep="\t", index=False)
        print(json.dumps(metrics[-1]), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("extract", "fit"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--length", type=int, choices=(2048, 16384))
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--labels", type=Path, default=C.C2_LABELS)
    parser.add_argument("--configs", type=Path)
    args = parser.parse_args()
    (extract if args.mode == "extract" else fit)(args)
