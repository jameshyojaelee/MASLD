"""Freedman-Lane calibration and power for t5b on simulated genes with known answers."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import t5b_stage_lineage as t  # noqa: E402

LIN = ["Hepatocytes", "Fibroblasts", "Cholangiocytes"]


def sim_gene(rng, n=120, stage_eff=0.0, lin_eff=0.0, tau=0.3, fib_noise=0.01):
    coh = rng.choice(["A", "B", "C"], n)
    S = rng.integers(0, 3, n).astype(float)
    hep = np.clip(0.9 - 0.08 * S + rng.normal(0, 0.05, n), 0.3, 0.99)  # composition tracks stage
    fib = np.clip(0.02 + 0.03 * S + rng.normal(0, fib_noise, n), 0.001, 0.5)
    cho = np.clip(0.02 + rng.normal(0, 0.01, n), 0.001, 0.5)
    depth = rng.poisson(60, n) + 20
    mu = 0.5 + stage_eff * S + lin_eff * (hep - hep.mean()) * 5 + rng.normal(0, tau, n)
    p = 1 / (1 + 2.0 ** (-mu))
    a = rng.binomial(depth, p)
    return pd.DataFrame({"cohort": coh, "S": S, "Hepatocytes": hep, "Fibroblasts": fib, "Cholangiocytes": cho,
                         "a_o": a, "r_o": depth - a, "n": depth})


def test_stage_null_calibrated_with_composition_confounding():
    t.N_PERM = 200
    rng = np.random.default_rng(7)
    ps = []
    for _ in range(200):
        g = sim_gene(rng, lin_eff=0.4)  # composition effect, no stage effect
        yw, Xf, Xr, coh, _ = t.whitened(g, LIN, "S")
        ps.append(t.freedman_lane(yw, Xf, Xr, coh, 1, rng)[1])
    frac = np.mean(np.array(ps) < 0.05)
    assert 0.02 <= frac <= 0.10, frac


def test_stage_power():
    t.N_PERM = 200
    rng = np.random.default_rng(8)
    ps, bs = [], []
    for _ in range(50):
        yw, Xf, Xr, coh, _ = t.whitened(sim_gene(rng, stage_eff=0.25, fib_noise=0.06), LIN, "S")
        ps.append(t.freedman_lane(yw, Xf, Xr, coh, 1, rng)[1])
        bs.append(np.linalg.lstsq(Xf, yw, rcond=None)[0][0])
    assert abs(np.mean(bs) - 0.25) < 0.05, np.mean(bs)
    assert np.mean(np.array(ps) < 0.05) > 0.6, np.mean(np.array(ps) < 0.05)


def test_lineage_null_and_power():
    t.N_PERM = 200
    rng = np.random.default_rng(9)
    p0 = [t.freedman_lane(*t.whitened(sim_gene(rng), LIN, "lineage")[:4], 3, rng)[1] for _ in range(150)]
    assert 0.02 <= np.mean(np.array(p0) < 0.05) <= 0.10
    p1 = [t.freedman_lane(*t.whitened(sim_gene(rng, lin_eff=0.6), LIN, "lineage")[:4], 3, rng)[1] for _ in range(40)]
    assert np.mean(np.array(p1) < 0.05) > 0.5


def test_halves_keep_relatives_together():
    ids = [f"IND{i:05d}" for i in range(200)]
    xw = pd.DataFrame({"run": [f"R{i}" for i in range(200)], "individual_id": ids})
    kin = pd.DataFrame({"#IID1": ["R1", "R5"], "IID2": ["R2", "R9"], "relation": ["first_degree", "first_degree"]})
    h = t.halves(ids, kin, xw)
    assert h["IND00001"] == h["IND00002"] and h["IND00005"] == h["IND00009"]
    assert 0.35 < np.mean(list(h.values())) < 0.65


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f(); print("ok", name)
