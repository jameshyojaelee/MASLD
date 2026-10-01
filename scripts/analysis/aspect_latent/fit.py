"""Model B fitting utilities (PRESPEC_B sections 4-5, Amendment 2).

SVI fits of aspect_model, held-out log predictive density for a whole cohort
(latent noise integrated per participant, the cohort's label shift/scale integrated
jointly over their priors), and simulation of labels from a fitted model for the
parametric-bootstrap null. The same code computes observed and null statistics.
"""
import jax

jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from jax.scipy.stats import norm  # noqa: E402
from numpyro.infer import SVI, Trace_ELBO  # noqa: E402
from numpyro.infer.autoguide import AutoNormal  # noqa: E402
from numpyro.optim import Adam  # noqa: E402

from model import LABELS, N_LEVELS, NAS_ONEHOT, aspect_model  # noqa: E402

SEED = 20260923
COMPONENTS = ("steatosis", "ballooning", "inflammation")
LABEL_COLS = ["steatosis", "ballooning", "inflammation", "nas", "fib_lo", "fib_hi"]


def load_inputs(path):
    d = pd.read_csv(path, sep="\t")
    prog = [c for c in d.columns if c.startswith("hotspot_")]
    return d, prog


def arrays(d, prog, cohorts):
    """Model arrays for the rows of d; cohort index follows the order in `cohorts`."""
    x = d[prog].to_numpy(float)
    x = np.where(np.isfinite(x), x, 0.0)  # program untestable in a cohort -> cohort mean
    lab = {c: d[c].fillna(-1).astype(int).to_numpy() for c in COMPONENTS}
    return dict(x=jnp.asarray(x), cohort=jnp.asarray(d["cohort"].map({c: i for i, c in enumerate(cohorts)}).to_numpy()),
                n_cohorts=len(cohorts), labels={k: jnp.asarray(v) for k, v in lab.items()},
                nas=jnp.asarray(d["nas"].fillna(-1).astype(int).to_numpy()),
                fib_lo=jnp.asarray(d["fib_lo"].fillna(-1).astype(int).to_numpy()),
                fib_hi=jnp.asarray(d["fib_hi"].fillna(-1).astype(int).to_numpy()))


def fit_svi(a, k, steps=4000, seed=SEED, draws=200):
    guide = AutoNormal(aspect_model)
    svi = SVI(aspect_model, guide, Adam(0.01), Trace_ELBO())
    res = svi.run(jax.random.PRNGKey(seed), steps, progress_bar=False, k=k, **a)
    post = guide.sample_posterior(jax.random.PRNGKey(seed + 1), res.params, sample_shape=(draws,), k=k, **a)
    return post, float(res.losses[-1])


def _cut(post, name, s):
    first = post[f"cut0_{name}"][s]
    return jnp.concatenate([jnp.array([first]), first + jnp.cumsum(post[f"cut_gaps_{name}"][s])])


def _lam(post, k):
    rows = []
    for j in range(len(LABELS)):
        rows.append(jnp.stack([post[f"lam_{j}_{c}"] if f"lam_{j}_{c}" in post else jnp.zeros_like(post["lam_0_0"])
                               for c in range(k)], axis=-1))
    return jnp.stack(rows, axis=-2)  # draws x labels x k


def _probs(eta, cut, shift, scale):
    z = (cut - shift - eta[..., None]) / scale
    cdf = norm.cdf(z)
    ones = jnp.ones(eta.shape + (1,))
    return jnp.clip(jnp.concatenate([cdf, ones], -1) - jnp.concatenate([0 * ones, cdf], -1), 1e-12, 1.0)


def cohort_logpred(post, k, a, rng, n_eps=32, n_shift=10, draws=None):
    """log p(labels of one held-out cohort | x, fitted global parameters)."""
    x = a["x"]
    n = x.shape[0]
    lam = _lam(post, k)
    S = post["W"].shape[0] if draws is None else draws
    vals = []
    for s in range(S):
        z_mean = x @ post["W"][s]  # n x k
        eps = rng.standard_normal((n_eps, n, k))
        z = z_mean[None] + eps  # e x n x k
        for _ in range(n_shift):
            ll = jnp.zeros((n_eps, n))
            probs = {}
            for j, name in enumerate(LABELS):
                shift = rng.normal(0.0, 1.0)
                scale = float(np.exp(rng.normal(0.0, 0.25)))
                probs[name] = _probs(z @ lam[s, j], _cut(post, name, s), shift, scale)  # e x n x C
            for name in COMPONENTS:
                y = a["labels"][name]
                p = jnp.take_along_axis(probs[name], jnp.clip(y, 0)[None, :, None], -1)[..., 0]
                ll += jnp.where(y >= 0, jnp.log(p), 0.0)
            miss = (a["labels"]["steatosis"] < 0) & (a["labels"]["ballooning"] < 0) & (a["labels"]["inflammation"] < 0)
            joint = probs["steatosis"][..., :, None, None] * probs["ballooning"][..., None, :, None] * probs["inflammation"][..., None, None, :]
            pn = joint.reshape(joint.shape[:-3] + (48,)) @ jnp.asarray(NAS_ONEHOT)
            use = miss & (a["nas"] >= 0)
            ll += jnp.where(use, jnp.log(jnp.clip(jnp.take_along_axis(pn, jnp.clip(a["nas"], 0)[None, :, None], -1)[..., 0], 1e-12)), 0.0)
            cum = jnp.cumsum(probs["fibrosis"], -1)
            lo, hi = jnp.clip(a["fib_lo"], 0), jnp.clip(a["fib_hi"], 0)
            p_hi = jnp.take_along_axis(cum, hi[None, :, None], -1)[..., 0]
            p_lo = jnp.where(lo > 0, jnp.take_along_axis(cum, jnp.clip(lo - 1, 0)[None, :, None], -1)[..., 0], 0.0)
            ll += jnp.where(a["fib_lo"] >= 0, jnp.log(jnp.clip(p_hi - p_lo, 1e-12)), 0.0)
            per_participant = jax.scipy.special.logsumexp(ll, axis=0) - np.log(n_eps)  # integrate eps
            vals.append(float(per_participant.sum()))
    vals = np.array(vals)
    return float(np.logaddexp.reduce(vals) - np.log(len(vals)))


def simulate_labels(post_mean, k, a, rng):
    """New labels from a fitted model with the same missingness pattern as `a`."""
    x = np.asarray(a["x"])
    n = x.shape[0]
    lam = np.asarray(_lam({kk: jnp.asarray(v)[None] for kk, v in post_mean.items()}, k))[0]
    z = x @ np.asarray(post_mean["W"]) + rng.standard_normal((n, k))
    coh = np.asarray(a["cohort"])
    out = {}
    for j, name in enumerate(LABELS):
        cut = np.asarray(_cut({kk: jnp.asarray(v)[None] for kk, v in post_mean.items()}, name, 0))
        shift = np.asarray(post_mean[f"shift_{name}"])[coh]
        scale = np.asarray(post_mean[f"scale_{name}"])[coh]
        u = z @ lam[j] + shift + scale * rng.standard_normal(n)
        out[name] = np.searchsorted(cut, u)
    lab = {c: np.where(np.asarray(a["labels"][c]) >= 0, out[c], -1) for c in COMPONENTS}
    nas = np.where(np.asarray(a["nas"]) >= 0, out["steatosis"] + out["ballooning"] + out["inflammation"], -1)
    lo0, hi0 = np.asarray(a["fib_lo"]), np.asarray(a["fib_hi"])
    interval = lo0 != hi0
    f = out["fibrosis"]
    lo = np.where(lo0 < 0, -1, np.where(interval, np.select([f <= 1, f == 2], [0, 2], 3), f))
    hi = np.where(lo0 < 0, -1, np.where(interval, np.select([f <= 1, f == 2], [1, 2], 4), f))
    b = dict(a)
    b.update(labels={c: jnp.asarray(v) for c, v in lab.items()}, nas=jnp.asarray(nas),
             fib_lo=jnp.asarray(lo), fib_hi=jnp.asarray(hi))
    return b


def subset(a, mask):
    m = np.asarray(mask)
    return dict(x=a["x"][m], cohort=a["cohort"][m], n_cohorts=a["n_cohorts"],
                labels={c: v[m] for c, v in a["labels"].items()}, nas=a["nas"][m],
                fib_lo=a["fib_lo"][m], fib_hi=a["fib_hi"][m])


def posterior_mean(post):
    return {k: np.asarray(v).mean(0) for k, v in post.items() if k not in ("eps", "z")}


def loco_statistic(a, cohorts, k_pair, rng, steps=4000):
    """Mean per-participant gain in held-out log predictive density, K2 over K1,
    leave-one-cohort-out over `cohorts` (indices into a['cohort'])."""
    out = {}
    for c in cohorts:
        test = np.asarray(a["cohort"]) == c
        tr, te = subset(a, ~test), subset(a, test)
        lp = {}
        for k in k_pair:
            post, _ = fit_svi(tr, k, steps=steps)
            lp[k] = cohort_logpred(post, k, te, rng)
        out[int(c)] = {"n": int(test.sum()), "lpd": {str(k): v for k, v in lp.items()},
                       "gain_per_participant": (lp[k_pair[1]] - lp[k_pair[0]]) / int(test.sum())}
    n_tot = sum(v["n"] for v in out.values())
    stat = sum(v["gain_per_participant"] * v["n"] for v in out.values()) / n_tot
    return stat, out
