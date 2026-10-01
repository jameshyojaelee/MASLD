"""Checks for t5a helpers on simulated data with known answers."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import t5a_allelic_table as t  # noqa: E402

RULE = {"min_dp": 20, "min_minor_reads": 2, "min_minor_frac": 0.02}


def sim_bb(rng, p, rho, n):
    s = (1 - rho) / rho
    return rng.binomial(n, rng.beta(p * s, (1 - p) * s, size=len(n)))


def test_bb_recovers_effect():
    rng = np.random.default_rng(1)
    n = rng.poisson(60, 400) + 20
    k = sim_bb(rng, 0.7, 0.05, n)
    par, _ = t.bb_fit(k.astype(float), n.astype(float))
    est = par[0] / t.LN2
    assert abs(est - np.log2(0.7 / 0.3)) < 0.1, est
    se = t.bb_se_intercept(par, k.astype(float), n.astype(float)) / t.LN2
    assert 0.01 < se < 0.2, se


def test_bb_lrt_null_calibrated():
    rng = np.random.default_rng(2)
    ps = []
    for _ in range(400):
        n = rng.poisson(40, 40) + 20
        k = sim_bb(rng, 0.5, 0.05, n).astype(float)
        _, l1 = t.bb_fit(k, n.astype(float))
        _, l0 = t.bb_fit(k, n.astype(float), fix_b=0.0)
        ps.append(t.stats.chi2.sf(max(0, 2 * (l0 - l1)), 1))
    frac = np.mean(np.array(ps) < 0.05)
    assert 0.02 < frac < 0.09, frac


def test_paule_mandel_recovers_tau2():
    rng = np.random.default_rng(3)
    n = 3000
    v = rng.uniform(0.02, 0.2, n)
    X = np.column_stack([np.ones(n), rng.normal(size=n)])
    y = X @ np.array([0.3, 0.1]) + rng.normal(0, np.sqrt(v + 0.05))
    tau2 = t.paule_mandel(y, v, X)
    assert abs(tau2 - 0.05) < 0.01, tau2
    assert t.paule_mandel(X @ np.array([0.3, 0.1]) + rng.normal(0, np.sqrt(v)), v, X) < 0.01


def test_n_min_monotone():
    assert t.n_min_for(0.5, RULE) == 20
    a, b = t.n_min_for(0.3, RULE), t.n_min_for(0.08, RULE)
    assert a <= b, (a, b)


def test_stage_coding():
    meta = pd.DataFrame({
        "sample_id": ["a", "b", "c", "d", "e", "f", "g"],
        "dataset": ["GSE213621", "GSE213621", "GSE135251", "GSE135251", "GSE130970", "GSE130970", "GSE126848"],
        "condition": ["Fibrosis_F3F4", "Control", "NASH", "Control", "Control", "Control", "NASH"],
        "fibrosis_stage": [3, 0, 2, 0, 0, 0, np.nan]})
    s = t.stage_coding(meta, {"f"}).set_index("run")
    assert s.loc["a", "S"] == 2 and s.loc["a", "C"] == 0
    assert np.isnan(s.loc["b", "S"]) and s.loc["b", "C"] == 1
    assert s.loc["c", "S"] == 1
    assert s.loc["d", "C"] == 1
    assert s.loc["e", "S"] == 0 and s.loc["e", "C"] == 0  # GSE130970 F0 NAFLD is staged
    assert s.loc["f", "C"] == 1                          # source control
    assert np.isnan(s.loc["g", "S"])


def test_het_call():
    ref = np.array([10, 18, 30, 29, 100])
    alt = np.array([10, 2, 1, 1, 1])
    got = t.het_call(ref, alt, RULE)
    assert got.tolist() == [True, True, False, False, False]


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f(); print("ok", name)
