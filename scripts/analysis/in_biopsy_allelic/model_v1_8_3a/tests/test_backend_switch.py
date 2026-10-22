"""v1.8.3-gpu: the torch backend against the CPU evaluators on the B-FIX fixture and a synthetic problem.

Device: MODEL_A_TORCH_DEVICE (default cuda when available, else cpu). Tolerances are rounding-level: the two
backends sum in different orders, so values agree to about 1e-13 relative, not bitwise.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    import torch
except ImportError:  # the rnaseq env has no torch
    pytest.skip("torch not installed in this env", allow_module_level=True)
if "MODEL_A_TORCH_DEVICE" not in os.environ:
    os.environ["MODEL_A_TORCH_DEVICE"] = "cuda" if torch.cuda.is_available() else "cpu"
import fixture_io as fio  # noqa: E402
import model_a as ma  # noqa: E402
from test_model_a import FIXTURE, _random_theta, _synthetic_problem  # noqa: E402


def both(fun, *a, **k):
    ma.set_backend("cpu")
    try:
        c = fun(*a, **k)
    finally:
        ma.set_backend("torch")
    try:
        t = fun(*a, **k)
    finally:
        ma.set_backend("cpu")
    return c, t


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    assert np.array_equal(np.isfinite(a), np.isfinite(b))
    return float(np.max(np.abs(a[m] - b[m])) / max(np.max(np.abs(a[m])), 1e-300)) if m.any() else 0.0


@pytest.fixture(scope="module")
def fx():
    prob, _ = fio.load_problem(FIXTURE)
    return prob, prob.data(), fio.load_params(FIXTURE, prob.L)


def test_default_backend_is_cpu():
    assert os.environ.get("MODEL_A_BACKEND", "cpu") == "cpu" or ma.BACKEND == os.environ["MODEL_A_BACKEND"]
    with pytest.raises(ValueError):
        ma.set_backend("gpu")


@pytest.mark.parametrize("name", ["theta1_start_like", "theta2_truth_v0", "theta4_perturbed"])
def test_loglik3_agrees(fx, name):
    prob, D, params = fx
    th = params[name].copy()
    th[prob.L.v] = 0.0
    for order in (0, 1, 2):
        c, t = both(ma.loglik3, th, D, order)
        if order == 0:
            c, t = (c,), (t,)
        assert abs(c[0] - t[0]) <= 1e-12 * abs(c[0])
        for x, y in zip(c[1:], t[1:]):
            assert rel(x, y) <= 1e-10
    c, t = both(ma.boundary_score, th, D)
    assert abs(c - t) <= 1e-10 * max(1.0, abs(c))


@pytest.mark.parametrize("name", ["theta3_truth_vpos", "theta4_perturbed"])
def test_loglik4_agrees(fx, name):
    prob, D, params = fx
    th = params[name]
    (lc, gc, Hc, ic), (lt, gt, Ht, it) = both(ma.loglik4, th, D, 2, detail=True)
    assert abs(lc - lt) <= 1e-12 * abs(lc)
    assert rel(gc, gt) <= 1e-9
    assert rel(Hc, Ht) <= 1e-9
    assert np.array_equal(ic["integrator"], it["integrator"])
    assert np.array_equal(ic["n_modes"], it["n_modes"])
    assert np.allclose(ic["u_hat"], it["u_hat"], rtol=1e-8, atol=1e-10)
    c, t = both(ma.loglik4, th, D, 0, all_reference=True)
    assert abs(c - t) <= 1e-12 * abs(c)


def test_loglik4_v0_branch(fx):
    prob, D, params = fx
    th = params["theta2_truth_v0"].copy()
    th[prob.L.v] = 0.0
    (lc, gc, Hc), (lt, gt, Ht) = both(ma.loglik4, th, D, 2)
    assert abs(lc - lt) <= 1e-12 * abs(lc)
    assert rel(gc, gt) <= 1e-9
    assert np.array_equal(np.isnan(Hc), np.isnan(Ht))
    assert rel(Hc, Ht) <= 1e-9


def test_synthetic_reference_genes():
    """A synthetic problem with large genes and a v at which some genes use the reference integrator, so the
    device fixed-node terms are exercised."""
    prob, th = _synthetic_problem([3, 40, 400, 1500], seed=5)
    D = prob.data()
    L = prob.L
    th = th.copy()
    th[L.v] = 0.05
    (lc, gc, Hc, ic), (lt, gt, Ht, it) = both(ma.loglik4, th, D, 2, detail=True, all_reference=True)
    assert abs(lc - lt) <= 1e-12 * abs(lc)
    assert rel(gc, gt) <= 1e-9
    assert rel(Hc, Ht) <= 1e-9
    th2 = _random_theta(L, np.random.default_rng(3))
    (lc, gc), (lt, gt) = both(ma.loglik4, th2, D, 1)
    assert abs(lc - lt) <= 1e-12 * abs(lc)
    assert rel(gc, gt) <= 1e-9


def test_device_cache_follows_data_identity(fx):
    import copy

    import backend_torch as bt
    prob, D, params = fx
    th = params["theta4_perturbed"].copy()
    th[prob.L.v] = 0.0
    D2 = copy.copy(D)
    D2.omega = np.zeros_like(D.omega)
    ma.set_backend("torch")
    try:
        a = ma.loglik3(th, D, 0)
        b = ma.loglik3(th, D2, 0)
        assert bt.device_model(D) is not bt.device_model(D2)
    finally:
        ma.set_backend("cpu")
    ma.set_backend("cpu")
    assert abs(b - ma.loglik3(th, D2, 0)) <= 1e-12 * abs(b)
    assert abs(a - ma.loglik3(th, D, 0)) <= 1e-12 * abs(a)
