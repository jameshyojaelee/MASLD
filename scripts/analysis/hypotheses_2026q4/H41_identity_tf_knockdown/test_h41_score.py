"""Unit checks for h41_score.py core statistics. Run: python3 -I test_h41_score.py (from this directory)."""
import os
import sys

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h41_score as h  # noqa: E402


def es_brute(z, mask):
    order = sorted(range(len(z)), key=lambda i: (-z[i], i))
    hit_w = sum(abs(z[i]) for i in range(len(z)) if mask[i])
    n_miss = len(z) - sum(mask)
    run, best = 0.0, 0.0
    for i in order:
        run += abs(z[i]) / hit_w if mask[i] else -1.0 / n_miss
        if abs(run) > abs(best):
            best = run
    return best


def test_es_matches_brute():
    rng = np.random.default_rng(1)
    for _ in range(20):
        z = rng.normal(size=60)
        mask = np.zeros(60, bool)
        mask[rng.choice(60, 8, replace=False)] = True
        assert np.isclose(h.es_matrix(z[None], mask)[0], es_brute(z, mask))


def test_es_extremes():
    z = np.arange(20, 0, -1).astype(float)
    top = np.zeros(20, bool)
    top[:5] = True
    assert np.isclose(h.es_matrix(z[None], top)[0], 1.0)
    assert np.isclose(h.es_matrix(z[None], top[::-1].copy())[0], -1.0)


def test_wtcs_sign_rule():
    assert h.wtcs(np.array([0.6]), np.array([-0.4]))[0] == 0.5
    assert h.wtcs(np.array([0.6]), np.array([0.4]))[0] == 0.0
    assert h.wtcs(np.array([-0.6]), np.array([0.4]))[0] == -0.5


def test_rank_std_spearman():
    rng = np.random.default_rng(2)
    a, b = rng.normal(size=100), rng.normal(size=100)
    assert np.isclose(h.rank_std(a)[0] @ h.rank_std(b)[0], spearmanr(a, b)[0])


def test_modz():
    rng = np.random.default_rng(3)
    z2 = rng.normal(size=(1, 2, 50))
    c, w, _ = h.modz_batch(z2, h.rank_std(z2[0])[None])
    assert np.allclose(w, 0.5) and np.allclose(c[0], z2[0].mean(axis=0))
    base = rng.normal(size=50)
    z3 = np.stack([base + 0.1 * rng.normal(size=50), base + 0.1 * rng.normal(size=50), -base])[None]
    c, w, _ = h.modz_batch(z3, h.rank_std(z3[0])[None])
    assert w[0, 2] < 0.05 and np.isclose(w.sum(), 1.0)


def test_bh():
    p = np.array([0.01, 0.04, 0.03, 0.5])
    q = h.bh(p)
    assert np.allclose(q, [0.04, 0.05333333, 0.05333333, 0.5])


def test_p_upper():
    assert h.p_upper(np.array([0.0, 1.0, 2.0]), 1.0) == 3.0 / 4.0
    assert h.p_upper(np.zeros(9999), 1.0) == 1.0 / 10000.0


def test_planted_signal():
    """A consensus that copies the disease direction gets a high WTCS and rho; a random one does not."""
    rng = np.random.default_rng(4)
    g = 978
    lfc = rng.normal(size=g)
    up = np.zeros(g, bool)
    dn = np.zeros(g, bool)
    o = np.argsort(lfc)
    up[o[-50:]] = True
    dn[o[:50]] = True
    u = h.Universe(np.arange(g), up, dn, lfc)
    s = u.score(lfc + 0.5 * rng.normal(size=g))
    assert s["wtcs"][0] > 0.5 and s["rho"][0] > 0.8
    r = u.score(rng.normal(size=g))
    assert abs(r["rho"][0]) < 0.15


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
