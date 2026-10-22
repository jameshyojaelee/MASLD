"""Unit checks for the H47 scoring core on synthetic data (no LINCS data used).

Usage: python3 -I test_h47_score.py
"""
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h47_score as h  # noqa: E402


def naive_es(z_row, mask):
    order = np.argsort(-z_row, kind="stable")
    n = z_row.size
    k = mask.sum()
    nr = np.abs(z_row[mask]).sum()
    run, best = 0.0, 0.0
    for g in order:
        run += abs(z_row[g]) / nr if mask[g] else -1.0 / (n - k)
        if abs(run) > abs(best):
            best = run
    return best


def main():
    rng = np.random.default_rng(1)
    z = rng.normal(size=(200, 300))
    up = np.zeros(300, bool)
    dn = np.zeros(300, bool)
    up[rng.choice(300, 30, replace=False)] = True
    dn[rng.choice(np.flatnonzero(~up), 25, replace=False)] = True
    ref = h.RankedRef(z)
    es_up = ref.es(up, chunk=64)
    naive = np.array([naive_es(z[i], up) for i in range(200)])
    assert np.allclose(es_up, naive), "vectorized ES differs from loop"

    # self-connectivity: query from signature 0's extremes ranks signature 0 first
    o = np.argsort(-z[0])
    q_up = np.zeros(300, bool)
    q_dn = np.zeros(300, bool)
    q_up[o[:25]] = True
    q_dn[o[-25:]] = True
    w = h.wtcs(ref.es(q_up), ref.es(q_dn))
    assert w.argmax() == 0 and w[0] > 0.9, w[0]
    # its mirror is the strongest reverser
    w_rev = h.wtcs(ref.es(q_dn), ref.es(q_up))
    assert w_rev.argmin() == 0 and w_rev[0] < -0.9

    # WTCS sign rule
    assert h.wtcs(np.array([0.3]), np.array([0.2]))[0] == 0
    assert np.isclose(h.wtcs(np.array([0.3]), np.array([-0.2]))[0], 0.25)

    # NCS: mean of positive NCS within each group is 1, negative -1
    grp = rng.integers(0, 3, size=200)
    wt = h.wtcs(ref.es(up), ref.es(dn))
    n = h.ncs(wt, grp, 3)
    for g in range(3):
        m = grp == g
        if (wt[m] > 0).any():
            assert np.isclose(n[m & (wt > 0)].mean(), 1.0)
        if (wt[m] < 0).any():
            assert np.isclose(n[m & (wt < 0)].mean(), -1.0)

    # AUC equals Mann-Whitney U / (n1 n0)
    s = rng.normal(size=500)
    mem = np.zeros(500, bool)
    mem[:20] = True
    s[:20] += 0.5
    u = stats.mannwhitneyu(s[mem], s[~mem]).statistic
    assert np.isclose(h.auc(s, mem), u / (20 * 480))

    # compound summary: median within line, mean across lines
    r = np.array([1.0, 3.0, 5.0, 10.0, -2.0])
    comp = np.array([0, 0, 0, 0, 1])
    line = np.array([0, 0, 0, 2, 1])
    pooled, per_line = h.compound_scores(r, comp, line, 2)
    assert np.isclose(per_line[0, 0], 3.0) and np.isclose(per_line[0, 2], 10.0) and np.isnan(per_line[0, 1])
    assert np.isclose(pooled[0], 6.5) and np.isclose(pooled[1], -2.0)

    # BH
    assert np.allclose(h.bh(np.array([0.01, 0.04, 0.03])), [0.03, 0.04, 0.04])
    print("all H47 scoring checks passed")


if __name__ == "__main__":
    main()
