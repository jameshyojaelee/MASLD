#!/usr/bin/env python3
"""Small compute-node scientific invariants; full-backbone checks run in GPU pilot."""
import argparse
import json
from pathlib import Path
import tempfile

import jax
import jax.numpy as jnp
import numpy as np
import optax

import model_scalar as M


def run():
    base = {f"alpha_genome/transformer_tower/mha_block{('_'+str(b)) if b else ''}/{q}_layer":
            {"w": jnp.eye(8, dtype=jnp.float32)} for b in range(9) for q in ("q", "v")}
    base["unrelated"] = {"w": jnp.ones((2, 2))}
    original = M.tree_hash(base)
    state = {"running_variance": jnp.ones(8)}
    state_hash = M.tree_hash(state)
    t, paths = M.initialize(base, width=8, hidden=8)
    assert len(paths) == 6
    zero = M.adapted_params(base, t, paths, "adapter")
    assert M.tree_hash(zero) == original
    ref = jax.random.normal(jax.random.PRNGKey(42), (5, 8))
    alt = ref + jax.random.normal(jax.random.PRNGKey(43), (5, 8))*0.2
    same = M.effect_from_features(t, ref, ref)
    forward = M.effect_from_features(t, ref, alt)
    reverse = M.effect_from_features(t, alt, ref)
    np.testing.assert_array_equal(np.asarray(same), 0)
    np.testing.assert_array_equal(np.asarray(forward), -np.asarray(reverse))
    native = jnp.linspace(-0.2, 0.2, 5)
    np.testing.assert_allclose(M.effect_from_features(t, ref, alt, native, 2),
                               -M.effect_from_features(t, alt, ref, -native, 2), atol=0)
    def objective(weights):
        p = M.adapted_params(base, weights, paths, "adapter")
        transformed_ref, transformed_alt = ref, alt
        for module, leaf in paths:
            transformed_ref = transformed_ref @ p[module][leaf]
            transformed_alt = transformed_alt @ p[module][leaf]
        return M.effect_loss(M.effect_from_features(weights, transformed_ref, transformed_alt), jnp.ones(5))
    grads = jax.grad(objective)(t)
    assert sum(float(jnp.linalg.norm(g["b"])) for g in grads["adaptation"].values()) > 0
    assert float(jnp.linalg.norm(grads["head"]["w1"])) > 0
    optimizer = optax.chain(optax.clip_by_global_norm(1), optax.adamw(3e-4, weight_decay=1e-4))
    updates, _ = optimizer.update(grads, optimizer.init(t), t)
    changed = optax.apply_updates(t, updates)
    assert M.tree_hash(base) == original and M.tree_hash(state) == state_hash
    assert M.tree_hash(changed) != M.tree_hash(t)
    with tempfile.TemporaryDirectory() as temp:
        path = Path(temp)/"weights.npz"
        M.save_trainable(path, changed, {"test": True})
        restored = M.load_trainable(path)
        np.testing.assert_array_equal(M.effect_from_features(changed, ref, alt), M.effect_from_features(restored, ref, alt))
    weights = M.pool_weights(2048, "symmetric")
    np.testing.assert_allclose(weights.sum(), 1)
    assert np.count_nonzero(weights) == 4
    try:
        M.pool_weights(2048, "target", target_start=-1, target_end=5)
    except ValueError:
        pass
    else:
        raise AssertionError("Target truncation must not silently change measurement")
    return {"identical_alleles_zero": True, "swapped_alleles_reverse": True,
            "native_residual_antisymmetric": True, "zero_adapter_exact": True,
            "nonzero_head_and_adapter_B_gradient": True, "frozen_base_and_state_unchanged": True,
            "save_reload_agreement": True, "target_coverage_refusal": True,
            "scope": "synthetic_small_tensors_actual_JAX_Optax; full_backbone_invariants_separate"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result))
