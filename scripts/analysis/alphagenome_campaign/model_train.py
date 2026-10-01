#!/usr/bin/env python3
"""Bounded local AlphaGenome scalar-effect training on development chromosomes.

This is a feasibility comparison by default (one seed, one split, 100 steps).
Chromosome 0 is the held fold, the next fold is validation, and all other folds
train. Default validation does not inspect the held fold's measured outcomes.
The caller must reserve resources before running; this script never submits jobs.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import sys
import time

import jax
import jax.numpy as jnp
import numpy as np
import optax
import pandas as pd
import pysam

import model_scalar as M

PROJ = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJ/"scripts/analysis/alphagenome_program"))
import i1_common as C
from i1_extract_mpra_embeddings import load_trunk


def emit(path, payload):
    path.write_text(json.dumps(payload, indent=2, default=str)+"\n")


def train(args):
    out = args.out
    out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    config = json.loads(args.config.read_text())
    if config.get("head", "shared_mlp64") != "shared_mlp64":
        raise ValueError("Scalar sequence trainer requires shared_mlp64; ridge recipes use model_frozen.py")
    emit(out/"config.json", config)
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute-node allocation required")
    gpu = C.gpu_card()
    if gpu.get("card_tag") != "l40s":
        raise RuntimeError("Fixed L40S hardware required")
    emit(out/"environment.json", {"python": sys.version, "platform": platform.platform(),
        "jax": jax.__version__, "numpy": np.__version__, "optax": optax.__version__,
        "job_id": os.environ["SLURM_JOB_ID"], "checkpoint": str(C.CHECKPOINT), **gpu})
    labels = pd.read_csv(args.labels, sep="\t")
    if labels.groupby("chr").heldout_fold.nunique().max() != 1:
        raise ValueError("Training requires whole-chromosome folds")
    if labels.lead_variant_id.duplicated().any():
        raise ValueError("Input contains repeated variants; explicit locus-weighted labels required")
    held = args.held_fold
    valid = args.validation_fold if args.validation_fold is not None else (held+1)%5
    if valid == held:
        raise ValueError("Validation fold must differ from the held fold")
    if args.seed is not None:
        config["seed"] = args.seed
        emit(out/"config.json", config)
    train_rows = labels.loc[~labels.heldout_fold.isin([held, valid])].copy()
    validation_rows = labels.loc[labels.heldout_fold == valid].copy()
    # Deterministic positional thinning is independent of outcomes; use all
    # training rows by default. No validation-derived early stopping in pilots.
    if args.train_limit:
        train_rows = train_rows.iloc[np.linspace(0, len(train_rows)-1, args.train_limit, dtype=int)].copy()
    validation_n = min(args.validation_limit, len(validation_rows)) if args.validation_limit else len(validation_rows)
    validation_rows = validation_rows.iloc[np.linspace(0, len(validation_rows)-1,
        validation_n, dtype=int)].copy()
    if config["loss"] == "gaussian_effect":
        if not args.se_column or args.se_column not in labels:
            raise ValueError("Trustworthy source SE column is required; not reconstructed from p-values")
        if not np.isfinite(train_rows[args.se_column]).all() or (train_rows[args.se_column] <= 0).any():
            raise ValueError("Nonpositive or nonfinite training standard errors")
    native_slope = 0.0
    if args.native_scores:
        native = pd.read_csv(args.native_scores, sep="\t").set_index("key").local_atac_liver
        for frame in (train_rows, validation_rows):
            frame["native_score"] = frame.lead_variant_id.str.replace("chr", "", regex=False).map(native)
        if train_rows.native_score.isna().any() or validation_rows.native_score.isna().any():
            raise ValueError("Residual model requires matched native-score coverage")
        x, y = train_rows.native_score.to_numpy(), train_rows.beta_alt.to_numpy()
        native_slope = float(np.dot(x, y)/max(np.dot(x, x), 1e-20))
    else:
        train_rows["native_score"] = 0.0
        validation_rows["native_score"] = 0.0
    emit(out/"split.json", {"held_fold": held, "validation_fold": valid,
        "training_folds": sorted(int(f) for f in train_rows.heldout_fold.unique()),
        "held_fold_never_trained_or_validated": held not in set(train_rows.heldout_fold) | {valid},
        "training_n": len(train_rows), "validation_n": len(validation_rows),
        "held_outcomes_evaluated": False, "training_native_slope": native_slope,
        "population": "existing_significance_selected_Currin_caQTL_leads",
        "loss_units": "source_FastQTL_ALT_dosage_beta", "seed": config["seed"]})
    fasta = pysam.FastaFile(C.FASTA_PATH)
    from alphagenome_research.model import one_hot_encoder
    encoder = one_hot_encoder.DNAOneHotEncoder()
    rejects = []
    def batch(frame):
        refs, alts, weights, observed, native, se, keys = [], [], [], [], [], [], []
        for row in frame.itertuples(index=False):
            lo, hi = C.window_bounds(int(row.pos_hg38), config["length"])
            check = C.validate_window(fasta, row.chr, lo, hi, int(row.pos_hg38), row.ref)
            if not check["acgt_ok"] or not check["reference_match"]:
                rejects.append({"variant_id": row.lead_variant_id, "reason": "reference_or_alphabet"})
                continue
            try:
                w = M.pool_weights(config["length"], config["pooling"],
                    target_start=int(row.peak_start_hg38), target_end=int(row.peak_stop_hg38), window_start=lo)
            except ValueError as exc:
                rejects.append({"variant_id": row.lead_variant_id, "reason": str(exc)})
                continue
            ref, alt = C.extract_ref_alt(fasta, row.chr, lo, hi, int(row.pos_hg38), row.ref, row.alt)
            refs.append(encoder.encode(ref)); alts.append(encoder.encode(alt)); weights.append(w)
            observed.append(float(row.beta_alt)); native.append(float(row.native_score)); keys.append(row.lead_variant_id)
            se.append(float(getattr(row, args.se_column)) if args.se_column else 1.0)
        if not refs:
            return None
        return tuple(jnp.asarray(x, dtype=jnp.float32) for x in (refs, alts, weights, observed, native, se)), keys
    trunk, base, state, device, restore_s = load_trunk()
    t, paths = M.initialize(base, mode=config["mode"], rank=config.get("rank", 4),
        last_blocks=config.get("last_blocks", 3), seed=config["seed"], hidden=config["hidden"])
    if config["loss"] == "gaussian_effect":
        t["log_residual_variance"] = jnp.asarray(-2.0, jnp.float32)
    emit(out/"trainable_paths.json", [{"module": m, "parameter": n, "shape": list(base[m][n].shape)} for m, n in paths])
    # Hashes are only for unchanged frozen parameter/state checks; not labels.
    base_hash, state_hash = M.tree_hash(base), M.tree_hash(state)
    labels_tree = jax.tree_util.tree_map(lambda _: "head", t)
    labels_tree["adaptation"] = jax.tree_util.tree_map(lambda _: "backbone", t["adaptation"])
    optimizer = optax.chain(optax.clip_by_global_norm(1.0), optax.multi_transform({
        "head": optax.adamw(3e-4, weight_decay=1e-4),
        "backbone": optax.adamw(config.get("backbone_lr", 3e-4), weight_decay=1e-4)}, labels_tree))
    opt_state = optimizer.init(t)
    def predict(weights, backbone, stats, ref, alt, pw, native):
        return M.sequence_effect(trunk, backbone, stats, weights, paths, config["mode"], ref, alt, pw, native, native_slope)
    predict_jit = jax.jit(predict)
    def loss_fn(weights, backbone, stats, ref, alt, pw, y, native, se):
        p = predict(weights, backbone, stats, ref, alt, pw, native)
        return M.effect_loss(p, y, loss=config["loss"], se=se,
            log_residual_variance=weights.get("log_residual_variance"))
    @jax.jit
    def step(weights, optimizer_state, backbone, stats, data):
        value, gradients = jax.value_and_grad(loss_fn)(weights, backbone, stats, *data)
        updates, next_state = optimizer.update(gradients, optimizer_state, weights)
        norms = {"all": optax.global_norm(gradients), "head": optax.global_norm(gradients["head"]),
                 "adaptation": optax.global_norm(gradients["adaptation"])}
        return optax.apply_updates(weights, updates), next_state, value, norms
    first = batch(train_rows.iloc[:args.microbatch])
    if first is None:
        raise ValueError("First input examples ineligible")
    data, _ = first
    ref, alt, pw, y, native, se = data
    original_effect = predict_jit(t, base, state, ref, alt, pw, native)
    same = predict_jit(t, base, state, ref, ref, pw, jnp.zeros_like(native))
    reverse = predict_jit(t, base, state, alt, ref, pw, -native)
    np.testing.assert_allclose(same, 0, atol=0, rtol=0)
    np.testing.assert_allclose(reverse, -original_effect, atol=1e-6, rtol=1e-6)
    zero_params = M.adapted_params(base, t, paths, config["mode"])
    if config["mode"] in ("adapter", "frozen", "partial"):
        assert M.tree_hash(zero_params) == base_hash
    emit(out/"pretraining_checks.json", {"identical_zero": True, "swapped_sign": True,
        "zero_adapter_or_initial_partial_exact": True, "restore_seconds": restore_s})
    del zero_params
    rng = np.random.default_rng(config["seed"])
    history = []
    wall_limit = args.max_hours*3600
    training_start = time.monotonic()
    for iteration in range(args.steps):
        if time.monotonic()-start > wall_limit-600:
            break
        iteration_start = time.monotonic()
        indices = rng.choice(len(train_rows), size=args.microbatch, replace=False)
        current = batch(train_rows.iloc[indices])
        if current is None:
            continue
        tick = time.monotonic()
        t, opt_state, value, norms = step(t, opt_state, base, state, current[0])
        jax.block_until_ready(value)
        record = {"step": iteration+1, "loss": float(value), "seconds": time.monotonic()-tick,
                  "iteration_wall_seconds": time.monotonic()-iteration_start,
                  **{f"gradient_norm_{k}": float(v) for k,v in norms.items()}}
        history.append(record)
        if not np.isfinite(record["loss"]) or not np.isfinite(record["gradient_norm_all"]):
            raise ValueError("Nonfinite loss or gradients")
        if iteration == 0 and (record["gradient_norm_head"] <= 0 or
                (config["mode"] != "frozen" and record["gradient_norm_adaptation"] <= 0)):
            raise ValueError("Intended head/backbone gradients are zero")
        if (iteration+1) % 10 == 0 or iteration == 0:
            print(json.dumps(record), flush=True)
            pd.DataFrame(history).to_csv(out/"training_history.tsv", sep="\t", index=False)
    if not history:
        raise RuntimeError("No training steps completed")
    training_wall_seconds = time.monotonic()-training_start
    if M.tree_hash(base) != base_hash or M.tree_hash(state) != state_hash:
        raise RuntimeError("Frozen backbone or running state changed")
    M.save_trainable(out/"weights.npz", t, {"config": config, "paths": paths, "checkpoint": str(C.CHECKPOINT),
        "terms": "AlphaGenome_noncommercial", "training_native_slope": native_slope})
    restored = M.load_trainable(out/"weights.npz")
    final_p = predict_jit(t, base, state, ref, alt, pw, native)
    reloaded_p = predict_jit(restored, base, state, ref, alt, pw, native)
    np.testing.assert_array_equal(final_p, reloaded_p)
    # A separate held genomic sequence panel assesses retained native function
    # by the original ATAC/DNase heads. It is a drift diagnostic, not measured
    # native-assay accuracy (those labels are unavailable here).
    from alphagenome.models import dna_model as api
    from alphagenome_research.model import dna_model as local_model
    from alphagenome_research.model.metadata import metadata
    md = {o:metadata.load(o) for o in (api.Organism.HOMO_SAPIENS, api.Organism.MUS_MUSCULUS)}
    _, full_apply, _, _, _ = local_model.create_model(md)
    from alphagenome.models import dna_output
    native_indices = {name: np.flatnonzero(~md[api.Organism.HOMO_SAPIENS].padding[kind])
        for name, kind in (("atac", dna_output.OutputType.ATAC), ("dnase", dna_output.OutputType.DNASE))}
    @jax.jit
    def native_function(params, stats, x):
        pred = full_apply(params, stats, x, jnp.zeros(x.shape[0], dtype=jnp.int32))
        center = config["length"]//2
        return {name: jnp.log2(1+jnp.sum(pred[name]["predictions_1bp"][:,center-250:center+251, native_indices[name]].astype(jnp.float32), axis=1))
                for name in ("atac", "dnase")}
    native_panel = batch(validation_rows.iloc[:1])
    retention = {"scope": "separate_validation_region_predicted_native_ATAC_DNase_track_drift", "measured_native_accuracy": "not_available"}
    if native_panel is not None:
        x = native_panel[0][0]
        p0 = native_function(base, state, x)
        p1 = native_function(M.adapted_params(base, t, paths, config["mode"]), state, x)
        for name in p0:
            delta = np.asarray(p1[name], dtype=np.float32)-np.asarray(p0[name], dtype=np.float32)
            if not np.isfinite(delta).all():
                raise ValueError(f"Nonfinite native function drift on real {name} tracks")
            retention[name] = {"rms_log2_track_change": float(np.sqrt(np.mean(delta**2))),
                               "max_abs_log2_track_change": float(np.max(np.abs(delta))), "tracks": delta.size}
    emit(out/"native_function_drift.json", retention)
    predictions = []
    validation_start = time.monotonic()
    validation_batch_seconds = []
    for offset in range(0, len(validation_rows), args.microbatch):
        tick = time.monotonic()
        b = batch(validation_rows.iloc[offset:offset+args.microbatch])
        if b is None:
            continue
        d, keys = b
        actual_n = len(keys)
        # Reuse the trained batch shape. Padding is discarded after prediction;
        # frozen statistics make a repeated padding example independent of its
        # neighbors and it contributes to neither labels nor metrics.
        if actual_n < args.microbatch:
            padded = tuple(jnp.concatenate([v,jnp.repeat(v[-1:],args.microbatch-actual_n,axis=0)],axis=0) for v in d)
        else:
            padded = d
        pred = predict_jit(t, base, state, padded[0], padded[1], padded[2], padded[4])
        pred = np.asarray(pred)
        for i, key in enumerate(keys):
            predictions.append({"variant_id": key, "predicted_beta": float(pred[i]),
                                "observed_beta": float(d[3][i]), "validation_fold": valid})
        validation_batch_seconds.append(time.monotonic()-tick)
    pd.DataFrame(predictions).to_csv(out/"validation_predictions.tsv", sep="\t", index=False)
    pd.DataFrame(rejects, columns=["variant_id", "reason"]).drop_duplicates().to_csv(out/"excluded.tsv", sep="\t", index=False)
    validation_wall_seconds = time.monotonic()-validation_start
    warm = [h["seconds"] for h in history[1:]]
    emit(out/"feasibility.json", {"completed_steps": len(history), "requested_steps": args.steps,
        "microbatch": args.microbatch, "first_step_seconds": history[0]["seconds"],
        "warm_step_seconds_median": float(np.median(warm)) if warm else None,
        "training_wall_seconds_including_batching_and_progress_writes": training_wall_seconds,
        "validation_wall_seconds": validation_wall_seconds, "validation_examples": len(predictions),
        "validation_first_batch_seconds": validation_batch_seconds[0] if validation_batch_seconds else None,
        "validation_warm_seconds_per_example": float(np.mean(validation_batch_seconds[1:]))/args.microbatch if len(validation_batch_seconds)>1 else None,
        "iteration_wall_seconds_median": float(np.median([h["iteration_wall_seconds"] for h in history[1:]])) if warm else None,
        "estimated_training_epoch_seconds": float(np.median(warm))*len(train_rows)/args.microbatch if warm else None,
        "elapsed_seconds_including_restore_checks_evaluation_checkpoint": time.monotonic()-start,
        "trainable_parameter_count": sum(v.size for v in jax.tree_util.tree_leaves(t)),
        "frozen_parameters_unchanged": True, "running_statistics_unchanged": True,
        "save_reload_agreement": True, "memory": C.memory_stats(device),
        "interpretation": "single_split_single_seed_feasibility_only", "formal_finalists": "five_seeds_nested_selection_not_run"})


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--labels", type=Path, default=C.C2_LABELS)
    p.add_argument("--native-scores", type=Path)
    p.add_argument("--se-column")
    p.add_argument("--steps", type=int, default=100)
    p.add_argument("--microbatch", type=int, default=1)
    p.add_argument("--train-limit", type=int, default=0)
    p.add_argument("--validation-limit", type=int, default=32)
    p.add_argument("--held-fold", type=int, default=0)
    p.add_argument("--validation-fold", type=int, default=None,
                   help="Validation fold; default (held+1)%%5. Set explicitly to rotate validation "
                        "while keeping one fold closed to both training and validation.")
    p.add_argument("--seed", type=int, default=None,
                   help="Override the config seed, for repeated fits of the same recipe.")
    p.add_argument("--max-hours", type=float, default=3.8)
    args = p.parse_args()
    try:
        train(args)
    except Exception as exc:
        if args.out.exists():
            emit(args.out/"failure.json", {"type": type(exc).__name__, "message": str(exc),
                                          "status": "failed_not_a_negative_scientific_result"})
        raise
