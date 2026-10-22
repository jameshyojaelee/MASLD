"""Latent-severity model for stage misclassification (shared by CGA1 and CGB3).

L ~ N(mu, s^2) in Kleiner units; a reading is K = clip(round(L + e), 0, 4),
e ~ N(0, sd_e^2); a cohort may report K binned into groups. sd_e is calibrated so
two independent readings of one L have within-person SD sigma, defined as
sd(K1 - K2) / sqrt(2), the definition behind RESULTS.md:1934-1956.
"""
import numpy as np
from scipy import optimize, stats

BINS_GSE213621 = {0: 1, 1: 1, 2: 2, 3: 3, 4: 3}
NMC = 400_000


def group_map(dataset):
    if dataset == "GSE213621":
        return BINS_GSE213621
    return {k: k for k in range(5)}


def read_stage(L, sd_e, rng):
    return np.clip(np.rint(L + rng.normal(0, sd_e, L.shape) if sd_e > 0 else L), 0, 4).astype(int)


def to_group(K, gmap):
    out = np.empty_like(K)
    for k, v in gmap.items():
        out[K == k] = v
    return out


def group_probs(mu, tau, gmap):
    """P(group) when L + e ~ N(mu, tau^2)."""
    edges = np.array([-np.inf, 0.5, 1.5, 2.5, 3.5, np.inf])
    pk = np.diff(stats.norm.cdf((edges - mu) / tau))
    groups = sorted(set(gmap.values()))
    return groups, np.array([pk[[k for k, v in gmap.items() if v == g]].sum() for g in groups])


def fit_mu_tau(group_counts, gmap):
    groups = sorted(set(gmap.values()))
    n = np.array([group_counts.get(g, 0) for g in groups], float)

    def nll(par):
        mu, ltau = par
        _, p = group_probs(mu, np.exp(ltau), gmap)
        return -np.sum(n * np.log(np.clip(p, 1e-12, None)))

    r = optimize.minimize(nll, x0=[1.5, 0.0], method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-9})
    return r.x[0], float(np.exp(r.x[1]))


def within_person_sd(mu, s, sd_e, rng, n=NMC):
    L = rng.normal(mu, s, n)
    k1, k2 = read_stage(L, sd_e, rng), read_stage(L, sd_e, rng)
    return np.std(k1 - k2) / np.sqrt(2)


def fit_latent(group_counts, gmap, sigma, seed):
    """Return mu, s, sd_e with s^2 + sd_e^2 = tau^2 and within-person SD = sigma."""
    mu, tau = fit_mu_tau(group_counts, gmap)
    if sigma == 0:
        return dict(mu=mu, s=tau, sd_e=0.0, tau=tau, sigma_achieved=0.0, capped=False)
    rng = np.random.default_rng(seed)

    def gap(sd_e):
        s = np.sqrt(max(tau ** 2 - sd_e ** 2, 1e-6))
        return within_person_sd(mu, s, sd_e, np.random.default_rng(seed)) - sigma

    hi = tau * 0.999
    if gap(hi) < 0:
        sd_e, capped = hi, True
    else:
        sd_e, capped = optimize.brentq(gap, 1e-3, hi, xtol=1e-3), False
    s = np.sqrt(max(tau ** 2 - sd_e ** 2, 1e-6))
    return dict(mu=mu, s=s, sd_e=sd_e, tau=tau,
                sigma_achieved=within_person_sd(mu, s, sd_e, rng), capped=capped)


def group_moments(par, gmap, seed, n=NMC):
    """E[L | group], Var(L | group), P(group) by Monte Carlo."""
    rng = np.random.default_rng(seed)
    L = rng.normal(par["mu"], par["s"], n)
    g = to_group(read_stage(L, par["sd_e"], rng), gmap)
    out = {}
    for k in sorted(set(gmap.values())):
        m = g == k
        out[k] = (L[m].mean(), L[m].var(), m.mean()) if m.sum() > 10 else (np.nan, np.nan, m.mean())
    return out


def fit_a(obs_means, obs_n, moments):
    """Weighted LS of observed group means of C on model E[L | group]."""
    ks = [k for k in obs_means if obs_n[k] > 0 and np.isfinite(moments[k][0])]
    x = np.array([moments[k][0] for k in ks])
    y = np.array([obs_means[k] for k in ks])
    w = np.array([obs_n[k] for k in ks], float)
    X = np.column_stack([np.ones_like(x), x])
    W = np.diag(w)
    beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ y)
    return float(beta[0]), float(beta[1])
