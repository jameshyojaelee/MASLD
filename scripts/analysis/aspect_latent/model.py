"""Noise-aware multi-aspect latent model (PRESPEC_B section 4).

Latent aspects z_d = W^T x_d + eps_d (eps ~ N(0, I_K)). Each histology label is
an ordinal-probit readout of z with cutpoints shared across cohorts and a
per-cohort shift and scale. NAS is the exact sum of steatosis (0-3), ballooning
(0-2) and lobular inflammation (0-3); a participant with only NAS contributes the
probability summed over the 48 component combinations. Fibrosis may be
interval-censored (GSE213621 records F0F1 / F2 / F3F4).

Label arrays use -1 for missing. Fibrosis intervals are given as (lo, hi).
"""
import itertools

import jax.numpy as jnp
import numpy as np
import numpyro
import numpyro.distributions as dist
from jax.scipy.stats import norm

# Loading order fixes the rotation: fibrosis anchors axis 1, steatosis axis 2.
# Likelihoods and the reported coefficient vectors W @ lambda_j do not depend on it.
LABELS = ("fibrosis", "steatosis", "inflammation", "ballooning")
N_LEVELS = {"steatosis": 4, "ballooning": 3, "inflammation": 4, "fibrosis": 5}

# 48 component combinations -> NAS value (0..8), as a one-hot matrix (48, 9)
_COMBOS = np.array(list(itertools.product(range(4), range(3), range(4))))
NAS_OF_COMBO = _COMBOS.sum(axis=1)
NAS_ONEHOT = np.eye(9)[NAS_OF_COMBO]


def ordinal_probs(eta, cut, shift, scale):
    """P(y = c) for c = 0..C-1 under an ordinal probit. eta (N,), cut (C-1,)."""
    z = (cut[None, :] - shift[:, None] - eta[:, None]) / scale[:, None]
    cdf = norm.cdf(z)
    upper = jnp.concatenate([cdf, jnp.ones((eta.shape[0], 1))], axis=1)
    lower = jnp.concatenate([jnp.zeros((eta.shape[0], 1)), cdf], axis=1)
    return jnp.clip(upper - lower, 1e-12, 1.0)


def nas_probs(p_s, p_b, p_i):
    """P(NAS = n), n = 0..8, from component probabilities (N,4), (N,3), (N,4)."""
    joint = p_s[:, :, None, None] * p_b[:, None, :, None] * p_i[:, None, None, :]
    return joint.reshape(joint.shape[0], 48) @ jnp.asarray(NAS_ONEHOT)


def lower_triangular_loadings(n_labels, k):
    """Label-by-K loading matrix, lower triangular with positive diagonal."""
    rows = []
    for j in range(n_labels):
        row = []
        for c in range(k):
            if c > j:
                row.append(0.0)
            elif c == j:
                row.append(numpyro.sample(f"lam_{j}_{c}", dist.HalfNormal(1.0)))
            else:
                row.append(numpyro.sample(f"lam_{j}_{c}", dist.Normal(0.0, 1.0)))
        rows.append(jnp.stack([jnp.asarray(v) for v in row]))
    return jnp.stack(rows)


def aspect_model(x, cohort, n_cohorts, labels, nas, fib_lo, fib_hi, k):
    """x (N,P) within-cohort scaled features; cohort (N,) ints; labels dict of (N,)
    arrays for steatosis/ballooning/inflammation (-1 missing); nas (N,) (-1 missing,
    used only where all three components are missing); fibrosis as [fib_lo, fib_hi]
    (both -1 if missing)."""
    n, p = x.shape
    w = numpyro.sample("W", dist.Normal(0.0, 0.5).expand([p, k]).to_event(2))
    with numpyro.plate("participants", n):
        eps = numpyro.sample("eps", dist.Normal(0.0, 1.0).expand([k]).to_event(1))
    z = x @ w + eps
    lam = lower_triangular_loadings(len(LABELS), k)

    probs = {}
    for j, name in enumerate(LABELS):
        c = N_LEVELS[name]
        gaps = numpyro.sample(f"cut_gaps_{name}", dist.HalfNormal(1.0).expand([c - 2]).to_event(1))
        first = numpyro.sample(f"cut0_{name}", dist.Normal(0.0, 2.0))
        cut = jnp.concatenate([jnp.array([first]), first + jnp.cumsum(gaps)])
        shift = numpyro.sample(f"shift_{name}", dist.Normal(0.0, 1.0).expand([n_cohorts]).to_event(1))
        scale = numpyro.sample(f"scale_{name}", dist.LogNormal(0.0, 0.25).expand([n_cohorts]).to_event(1))
        probs[name] = ordinal_probs(z @ lam[j], cut, shift[cohort], scale[cohort])

    ll = jnp.zeros(n)
    for name in ("steatosis", "ballooning", "inflammation"):
        y = labels[name]
        ll += jnp.where(y >= 0, jnp.log(jnp.take_along_axis(probs[name], jnp.clip(y, 0)[:, None], 1)[:, 0]), 0.0)
    comp_missing = (labels["steatosis"] < 0) & (labels["ballooning"] < 0) & (labels["inflammation"] < 0)
    p_nas = nas_probs(probs["steatosis"], probs["ballooning"], probs["inflammation"])
    use_nas = comp_missing & (nas >= 0)
    ll += jnp.where(use_nas, jnp.log(jnp.clip(jnp.take_along_axis(p_nas, jnp.clip(nas, 0)[:, None], 1)[:, 0], 1e-12)), 0.0)
    cum = jnp.cumsum(probs["fibrosis"], axis=1)
    lo, hi = jnp.clip(fib_lo, 0), jnp.clip(fib_hi, 0)
    p_int = jnp.take_along_axis(cum, hi[:, None], 1)[:, 0] - jnp.where(
        lo > 0, jnp.take_along_axis(cum, (lo - 1)[:, None], 1)[:, 0], 0.0)
    ll += jnp.where(fib_lo >= 0, jnp.log(jnp.clip(p_int, 1e-12)), 0.0)
    numpyro.factor("loglik", ll.sum())
    numpyro.deterministic("z", z)
