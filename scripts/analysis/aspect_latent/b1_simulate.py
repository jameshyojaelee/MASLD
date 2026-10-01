#!/usr/bin/env python3
"""B1: design-only checks for Model B, run before PRESPEC_B is hashed.

1. Unit test: the NAS-sum likelihood equals the brute-force marginal over the
   48 component combinations.
2. Recovery: data simulated from a known K = 2 truth at the development design
   (cohort sizes and label types as in PRESPEC_B) are refit with NUTS; the
   rotation-invariant coefficient vectors W @ lambda_j must match the truth
   (cosine >= 0.9 for fibrosis and steatosis).
3. MDE in the sealed GSE281797 design (n = 91, observed label marginals): power
   of the component-specific partial Spearman test when the component carries
   its own signal of a given size beyond the shared axis. The predictor is the
   true component signal plus the development-fit noise, so this MDE is a lower
   bound (the real model predicts less well). Reads no RNA and no label-RNA link.
"""
import argparse
import itertools
import json
from pathlib import Path

import jax

jax.config.update("jax_enable_x64", True)  # float64: ordinal CDF differences underflow in float32
import jax.numpy as jnp  # noqa: E402
import numpy as np
from numpyro.infer import MCMC, NUTS
from scipy import stats

from model import LABELS, N_LEVELS, NAS_OF_COMBO, aspect_model, nas_probs, ordinal_probs

SEED = 20260923
# development design: (cohort, n, label type) per PRESPEC_B section 2
DESIGN = [("GSE130970", 72, "components"), ("GSE267145", 99, "components"),
          ("GSE135251", 206, "nas"), ("GSE162694", 112, "nas"), ("GSE193066", 106, "nas"),
          ("GSE213621", 299, "fib_interval"), ("GSE240729", 67, "fib_only")]


def test_nas_sum(rng):
    n = 7
    p = [rng.dirichlet(np.ones(c), size=n) for c in (4, 3, 4)]
    fast = np.asarray(nas_probs(*[jnp.asarray(a) for a in p]))
    brute = np.zeros((n, 9))
    for a, b, c in itertools.product(range(4), range(3), range(4)):
        brute[:, a + b + c] += p[0][:, a] * p[1][:, b] * p[2][:, c]
    assert np.allclose(fast, brute, atol=1e-12), "NAS-sum likelihood disagrees with brute force"
    assert len(NAS_OF_COMBO) == 48
    return float(np.abs(fast - brute).max())


def draw_ordinal(eta, cut, rng):
    u = eta + rng.standard_normal(eta.shape)
    return np.searchsorted(cut, u)


def simulate(rng, p=117, k=2):
    w = rng.normal(0, 0.15, size=(p, k))
    lam = {"fibrosis": np.array([1.0, 0.0]), "steatosis": np.array([0.4, 0.9]),
           "inflammation": np.array([0.6, 0.5]), "ballooning": np.array([0.7, 0.3])}
    cuts = {name: np.linspace(-1.0, 1.5, N_LEVELS[name] - 1) for name in LABELS}
    xs, cohort, lab = [], [], {n: [] for n in ("steatosis", "ballooning", "inflammation")}
    nas, flo, fhi = [], [], []
    for ci, (_, n, kind) in enumerate(DESIGN):
        x = rng.standard_normal((n, p))
        z = x @ w + rng.standard_normal((n, k))
        y = {name: draw_ordinal(z @ lam[name], cuts[name], rng) for name in LABELS}
        comp = kind == "components"
        for name in lab:
            lab[name].append(y[name] if comp else -np.ones(n, int))
        nas.append(y["steatosis"] + y["ballooning"] + y["inflammation"] if kind == "nas" else -np.ones(n, int))
        if kind == "fib_interval":
            lo = np.select([y["fibrosis"] <= 1, y["fibrosis"] == 2], [0, 2], 3)
            hi = np.select([y["fibrosis"] <= 1, y["fibrosis"] == 2], [1, 2], 4)
        else:
            lo = hi = y["fibrosis"]
        flo.append(lo); fhi.append(hi); xs.append(x); cohort.append(np.full(n, ci))
    return dict(x=np.vstack(xs), cohort=np.concatenate(cohort),
                labels={k2: np.concatenate(v) for k2, v in lab.items()},
                nas=np.concatenate(nas), fib_lo=np.concatenate(flo), fib_hi=np.concatenate(fhi)), w, lam


def recovery(rng, chains, warmup, draws):
    d, w, lam = simulate(rng)
    mcmc = MCMC(NUTS(aspect_model), num_warmup=warmup, num_samples=draws, num_chains=chains,
                chain_method="sequential", progress_bar=False)
    mcmc.run(jax.random.PRNGKey(SEED), x=jnp.asarray(d["x"]), cohort=jnp.asarray(d["cohort"]),
             n_cohorts=len(DESIGN), labels={k: jnp.asarray(v) for k, v in d["labels"].items()},
             nas=jnp.asarray(d["nas"]), fib_lo=jnp.asarray(d["fib_lo"]), fib_hi=jnp.asarray(d["fib_hi"]), k=2)
    s = mcmc.get_samples()
    w_hat = np.asarray(s["W"]).mean(0)
    out = {}
    for j, name in enumerate(LABELS):
        lam_hat = np.array([np.asarray(s[f"lam_{j}_{c}"]).mean() if f"lam_{j}_{c}" in s else 0.0 for c in range(2)])
        a, b = w_hat @ lam_hat, w @ lam[name]
        out[name] = float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))
    div = int(np.asarray(mcmc.get_extra_fields()["diverging"]).sum())
    return out, div


def partial_spearman(pred, y, others):
    """Spearman correlation of pred and y after regressing both ranks on the ranks of others."""
    design = np.column_stack([np.ones(len(y)), stats.rankdata(others, axis=0)])
    def resid(v):
        r = stats.rankdata(v)
        return r - design @ np.linalg.lstsq(design, r, rcond=None)[0]
    return float(np.corrcoef(resid(pred), resid(y))[0, 1])


def label_marginals(path):
    """GSE281797 label marginals (labels only; no RNA is read)."""
    import pandas as pd
    d = pd.read_csv(path, sep="\t")
    d = d[d["sex_concordant"].astype(str) == "True"]
    out = {}
    for col, name, levels in [("steatosis", "steatosis", 4), ("inflammation", "inflammation", 4),
                              ("ballooning", "ballooning", 3), ("fibrosis", "fibrosis", 5)]:
        v = d[col].dropna().astype(int).clip(0, levels - 1)
        out[name] = (np.bincount(v, minlength=levels) + 0.5) / (len(v) + 0.5 * levels)
    return out, int(len(d))


def mde(rng, marg, n, sims=2000, effects=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5)):
    """Power of the partial Spearman test per component at the sealed design.
    Each label = shared axis + optional own signal + readout noise; the predictor is
    the true signal plus noise, so power is an upper bound and the MDE a lower bound."""
    cuts = {k: stats.norm.ppf(np.cumsum(v)[:-1]) for k, v in marg.items()}
    power = {}
    for comp in marg:
        power[comp] = {}
        null_q95 = None
        for e in effects:
            vals = []
            for _ in range(sims):
                shared, own = rng.standard_normal(n), rng.standard_normal(n)
                lat = {k: 0.7 * shared + 0.3 * rng.standard_normal(n) for k in marg}
                lat[comp] = lat[comp] + e * own
                y = {k: np.searchsorted(cuts[k], lat[k] + rng.standard_normal(n)) for k in marg}
                pred = 0.7 * shared + e * own + 0.5 * rng.standard_normal(n)
                others = np.column_stack([y[k] for k in marg if k != comp]).astype(float)
                vals.append(partial_spearman(pred, y[comp], others))
            vals = np.array(vals)
            if e == 0.0:
                null_q95 = float(np.quantile(vals, 0.95))
            power[comp][str(e)] = float((vals > null_q95).mean())
        power[comp]["null_q95"] = null_q95
        passing = [e for e in effects if e > 0 and power[comp][str(e)] >= 0.8]
        power[comp]["mde_80"] = min(passing) if passing else None
    return power


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--chains", type=int, default=4)
    ap.add_argument("--warmup", type=int, default=500)
    ap.add_argument("--draws", type=int, default=500)
    ap.add_argument("--sims", type=int, default=2000)
    ap.add_argument("--sealed-labels", default="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/GSE281797/metadata/sample_phenotypes_resolved.tsv")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    res = {"nas_sum_max_abs_diff": test_nas_sum(rng)}
    cos, div = recovery(rng, a.chains, a.warmup, a.draws)
    res["recovery_cosine"] = cos
    res["recovery_divergences"] = div
    res["recovery_pass"] = bool(cos["fibrosis"] >= 0.9 and cos["steatosis"] >= 0.9)
    marg, n_sealed = label_marginals(a.sealed_labels)
    res["gse281797_n"] = n_sealed
    res["gse281797_marginals"] = {k: v.round(4).tolist() for k, v in marg.items()}
    res["gse281797_power"] = mde(rng, marg, n_sealed, sims=a.sims)
    (out / "b1_results.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
