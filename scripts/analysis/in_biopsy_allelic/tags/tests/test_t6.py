"""L1 lineage mixture recovers lineage effects and beats L0 when composition varies; ties it when it does not."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import t6_predict as t  # noqa: E402

LIN = ["Hepatocytes", "Macrophages", "Endothelial cells", "T cells"]


def sim(rng, n=300, vary=True):
    hep = rng.uniform(0.5, 0.95, n) if vary else np.full(n, 0.75)
    mac = (1 - hep) * 0.6
    theta = pd.DataFrame({"Hepatocytes": hep, "Macrophages": mac, "Endothelial cells": (1 - hep) * 0.3,
                          "T cells": (1 - hep) * 0.1})
    e_g = pd.Series({"Hepatocytes": 10.0, "Macrophages": 100.0, "Endothelial cells": 5.0, "T cells": 1.0})
    comps = ["Hepatocytes", "Macrophages", "Endothelial cells"]
    W = t.lineage_weights(theta, e_g, comps)
    pi = np.array([0.75, 0.25, 0.5])
    depth = rng.poisson(80, n) + 20
    k = rng.binomial(depth, W @ pi)
    return k.astype(float), depth.astype(float), W


def test_weights_sum_to_one():
    rng = np.random.default_rng(0)
    _, _, W = sim(rng)
    assert np.allclose(W.sum(1), 1)


def test_l1_recovers_and_beats_l0():
    rng = np.random.default_rng(1)
    k, n, W = sim(rng)
    fold = np.arange(len(k)) % t.K_FOLD
    g = pd.DataFrame({"a_o": k, "n": n})
    ll = t.crossfit_gene(g, W, fold, np.nan)
    assert np.sum(ll["l1"]) > np.sum(ll["l0"]) + 20, (np.sum(ll["l1"]), np.sum(ll["l0"]))
    p0, _ = t.fit_l0(k, n)
    pi, _, _ = t.fit_l1(k, n, W, p0)
    assert abs(pi[0] - 0.75) < 0.05 and abs(pi[1] - 0.25) < 0.1, pi


def test_l1_no_gain_without_variation():
    rng = np.random.default_rng(2)
    k, n, W = sim(rng, vary=False)
    fold = np.arange(len(k)) % t.K_FOLD
    ll = t.crossfit_gene(pd.DataFrame({"a_o": k, "n": n}), W, fold, np.nan)
    assert abs(np.sum(ll["l1"]) - np.sum(ll["l0"])) < 5, (np.sum(ll["l1"]), np.sum(ll["l0"]))


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f(); print("ok", name)
