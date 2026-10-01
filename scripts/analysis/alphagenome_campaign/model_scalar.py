#!/usr/bin/env python3
"""Shared per-allele scalar effect model for measured molecular effects.

The pretrained AlphaGenome backbone remains a separate read-only artifact.
LoRA changes q_layer/w and v_layer/w in the final 3 or 5 of 9 attention
blocks. Frozen RMSBatchNorm statistics come from inference-mode trunk_apply.
All subtraction, pooling, loss and trainable head arithmetic are float32.
Source-native association beta is an association target, not an isolated
causal allele effect. A Gaussian effect likelihood requires supplied SEs.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re

import jax
import jax.numpy as jnp
import numpy as np


def qv_paths(params, last_blocks):
    if last_blocks not in (3, 5):
        raise ValueError("Initial comparison uses exactly last three/five blocks")
    found = []
    for module, leaves in params.items():
        match = re.search(r"/mha_block(?:_(\d+))?/(q_layer|v_layer)$", module)
        if match:
            block = int(match.group(1) or 0)
            if set(leaves) != {"w"} or leaves["w"].ndim != 2:
                raise ValueError(f"Unexpected projection structure: {module}")
            found.append((block, module, "w"))
    if len(found) != 18 or set(b for b, _, _ in found) != set(range(9)):
        raise ValueError(f"Expected 9 pairs of installed q_layer/v_layer paths; found {found}")
    return [(m, n) for b, m, n in sorted(found) if b >= 9-last_blocks]


def initialize(params, *, mode="adapter", rank=4, last_blocks=3, hidden=64, seed=1103, width=3072):
    keys = jax.random.split(jax.random.PRNGKey(seed), 100)
    if mode not in ("adapter", "partial", "frozen"):
        raise ValueError(mode)
    head = {"w1": jax.random.normal(keys[0], (width, hidden), dtype=jnp.float32) / np.sqrt(width),
            "b1": jnp.zeros(hidden, dtype=jnp.float32),
            "w2": jax.random.normal(keys[1], (hidden, 1), dtype=jnp.float32) / np.sqrt(hidden)}
    trainable = {"head": head, "adaptation": {}}
    paths = [] if mode == "frozen" else qv_paths(params, last_blocks)
    for i, (module, leaf) in enumerate(paths):
        w = params[module][leaf]
        if mode == "adapter":
            if rank not in (4, 16):
                raise ValueError("Initial ranks are 4/16")
            trainable["adaptation"][str(i)] = {
                "a": jax.random.normal(keys[i+2], (w.shape[0], rank), dtype=jnp.float32) * 0.01,
                "b": jnp.zeros((rank, w.shape[1]), dtype=jnp.float32)}
        else:
            trainable["adaptation"][str(i)] = {"w": jnp.asarray(w, dtype=jnp.float32)}
    # An unconstrained residual variance belongs to the scalar-effect loss,
    # not to model inference; initialized separately only for Gaussian runs.
    return trainable, paths


def adapted_params(base, trainable, paths, mode):
    """Construct a view with replacements; never mutate the pretrained tree."""
    merged = dict(base)
    for i, (module, leaf) in enumerate(paths):
        update = trainable["adaptation"][str(i)]
        if mode == "adapter":
            # alpha/rank=1. A random, B zero: exact pretrained initial function.
            value = jnp.asarray(base[module][leaf], jnp.float32) + update["a"] @ update["b"]
        elif mode == "partial":
            value = update["w"]
        else:
            raise ValueError(mode)
        merged[module] = dict(base[module], **{leaf: value})
    return merged


def pool_weights(length, pooling, *, target_start=None, target_end=None, window_start=0):
    """128-bp bins; symmetric uses exact overlap of [variant-192,variant+192).

    Measured-target pooling requires the whole measured interval in the window;
    truncation is a coverage failure. Fractional bin overlap prevents rounding a
    target boundary to a different genomic interval.
    """
    if length % 128:
        raise ValueError("Input length must tile 128-bp bins")
    centre = length//2
    if pooling == "variant":
        lo, hi = (centre//128)*128, (centre//128+1)*128
    elif pooling == "symmetric":
        lo, hi = centre-192, centre+192
    elif pooling == "target":
        if target_start is None or target_end is None:
            raise ValueError("Measured-target coordinates missing")
        lo, hi = target_start-window_start, target_end-window_start
        if not 0 <= lo < hi <= length:
            raise ValueError("measured_target_not_fully_observed")
    else:
        raise ValueError(pooling)
    starts = np.arange(length//128)*128
    overlap = np.maximum(0, np.minimum(starts+128, hi)-np.maximum(starts, lo)).astype(np.float32)
    return overlap/overlap.sum()


def shared_head(head, embedding):
    x = jnp.asarray(embedding, jnp.float32)
    # Per-allele layer normalization uses no held-example fitted statistics.
    x = (x-jnp.mean(x, axis=-1, keepdims=True)) * jax.lax.rsqrt(jnp.var(x, axis=-1, keepdims=True)+1e-5)
    return (jax.nn.gelu(x @ head["w1"] + head["b1"]) @ head["w2"])[..., 0]


def effect_from_features(trainable, ref, alt, native_score=None, native_slope=0.0):
    effect = shared_head(trainable["head"], alt).astype(jnp.float32) - shared_head(trainable["head"], ref).astype(jnp.float32)
    if native_score is not None:
        effect = effect + jnp.asarray(native_slope, jnp.float32)*jnp.asarray(native_score, jnp.float32)
    return effect


def sequence_effect(trunk, base, state, trainable, paths, mode, ref, alt, weights, native_score=None, native_slope=0.0):
    params = adapted_params(base, trainable, paths, mode)
    n = ref.shape[0]
    # Identical execution for the two alleles, running statistics frozen.
    both = jnp.concatenate([ref, alt], axis=0)
    embeddings = trunk(params, state, both, jnp.zeros(2*n, dtype=jnp.int32)).get_sequence_embeddings(128).astype(jnp.float32)
    w = jnp.concatenate([weights, weights], axis=0)
    pooled = jnp.sum(embeddings*w[..., None], axis=1, dtype=jnp.float32)
    return effect_from_features(trainable, pooled[:n], pooled[n:], native_score, native_slope)


def effect_loss(prediction, observed, *, loss="mse", se=None, log_residual_variance=None):
    residual = jnp.asarray(prediction, jnp.float32)-jnp.asarray(observed, jnp.float32)
    if loss == "mse":
        return jnp.mean(jnp.square(residual), dtype=jnp.float32)
    if loss == "gaussian_effect":
        if se is None or log_residual_variance is None:
            raise ValueError("Gaussian effect likelihood requires trustworthy source SE and residual variance")
        variance = jnp.square(jnp.asarray(se, jnp.float32)) + jax.nn.softplus(log_residual_variance) + 1e-8
        return 0.5*jnp.mean(jnp.square(residual)/variance+jnp.log(variance), dtype=jnp.float32)
    raise ValueError(loss)


def save_trainable(path, tree, metadata):
    """Portable small derivative weights only, never raw participant data."""
    arrays = {}
    def collect(node, prefix):
        for key, value in node.items():
            name = f"{prefix}/{key}" if prefix else key
            if isinstance(value, dict):
                collect(value, name)
            else:
                arrays[name] = np.asarray(value)
    collect(tree, "")
    np.savez_compressed(path, **arrays)
    Path(str(path)+".json").write_text(json.dumps(metadata, indent=2, default=str)+"\n")


def load_trainable(path):
    tree = {}
    with np.load(path, allow_pickle=False) as archive:
        for key in archive.files:
            node = tree
            parts = key.split("/")
            for part in parts[:-1]:
                node = node.setdefault(part, {})
            node[parts[-1]] = jnp.asarray(archive[key])
    tree.setdefault("adaptation", {})
    return tree


def tree_hash(tree):
    digest = hashlib.sha256()
    for value in jax.tree_util.tree_leaves(tree):
        x = np.asarray(value)
        digest.update(str(x.shape).encode())
        digest.update(x.tobytes())
    return digest.hexdigest()
