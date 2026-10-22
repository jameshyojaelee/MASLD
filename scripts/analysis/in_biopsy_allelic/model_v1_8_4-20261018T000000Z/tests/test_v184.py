"""v1.8.4 unit tests (SPEC_v184_ADDENDUM.md, H-18): K = 3 selection (item 3), probe start (item 4), seed
entropies (item 7), the section-4 pair without an extra rule-4 restart (item 6.1), l3* and LR (items 6.3-6.4)."""
import itertools
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import drivers as dv  # noqa: E402
import fitting as ft  # noqa: E402
import fixture_io as fio  # noqa: E402

FIXTURE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/"
               "MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z")


class F:
    def __init__(self, conv, ll, v=0.0, ll_ref=None):
        self.converged, self.loglik = conv, ll
        self.loglik_ref = ll if ll_ref is None else ll_ref
        self.theta = np.array([0.3, v])


# ------------------------------------------------------------------ item 3
def test_select_k_reduces_to_v18_table_at_k2():
    vals = [-10.0, -10.0 + 5e-7, -10.0 + 1e-6, -10.0 + 2e-6, -9.0, -np.inf, np.nan]
    for lw, lc, cw, cc in itertools.product(vals, vals, (True, False), (True, False)):
        if (cw and not np.isfinite(lw)) or (cc and not np.isfinite(lc)):
            continue   # item 3 drops converged starts with non-finite l from C; v1.8 select did not say
        w, c = F(cw, lw), F(cc, lc)
        kept, _ = ft.select(w, c, 1e-6)
        k, _ = ft.select_k([w, c], 1e-6)
        old = None if kept is None else (0 if kept is w else 1)
        assert k == old, (lw, lc, cw, cc, k, old)


def test_select_k_ties_and_failures():
    e = 1e-6
    assert ft.select_k([F(True, -10.0), F(True, -10.0 + 5e-7), F(True, -10.0 + 9e-7)], e)[0] == 0
    assert ft.select_k([F(True, -10.0), F(True, -9.0), F(True, -9.0 + 5e-7)], e)[0] == 1
    assert ft.select_k([F(True, -10.0), F(True, -9.5), F(True, -9.0)], e)[0] == 2
    # a difference of exactly eps keeps the higher fit
    assert ft.select_k([F(True, -1.0), F(True, -1.0 + 1.0 / 1024)], 1.0 / 1024)[0] == 1
    # unconverged probe higher by > eps fails; within eps does not
    assert ft.select_k([F(True, -10.0), F(True, -10.0), F(False, -10.0 + 2e-6)], e)[0] is None
    assert ft.select_k([F(True, -10.0), F(True, -10.0), F(False, -10.0 + 5e-7)], e)[0] == 0
    # non-finite unconverged ignored; none converged fails
    assert ft.select_k([F(True, -10.0), F(False, np.nan), F(False, np.inf)], e)[0] == 0
    k, rule = ft.select_k([F(False, -10.0), F(False, -9.0), F(False, -8.0)], e)
    assert k is None and rule == "failed, none converged"
    # kept is the lowest index within eps of the max even when the max is a later index
    assert ft.select_k([F(True, -10.0 + 1e-7), F(True, -20.0), F(True, -10.0 + 8e-7)], e)[0] == 0


# ------------------------------------------------------------------ item 7
def test_entropies():
    assert ft.probe_entropy_gateb_observed(27, 4) == [20261015, 27, 4, 10, 0, 0, 1]
    assert ft.probe_entropy_gateb_refit(27, 4, 3, 11) == [20261015, 27, 4, 10, 3, 11, 1]
    assert ft.probe_entropy_full_crt_observed(5) == [20261015, 99, 5, 10, 0, 0, 1]
    assert ft.probe_entropy_full_crt_refit(5, 999) == [20261015, 99, 5, 10, 1, 999, 1]
    assert ft.probe_entropy_fixed_observed() == [20260923, 13, 0, 0, 1]
    assert ft.probe_entropy_fixed_refit(3, 24) == [20260923, 13, 3, 24, 1]
    assert dv.probe_entropies_fixed(1, 0) == [[20260923, 13, 1, 0, 1]]
    with pytest.raises(ValueError):
        ft.probe_entropy_fixed_refit(0, 1)
    with pytest.raises(ValueError):
        ft.probe_entropy_gateb_refit(1, 1, 4, 0)


# ------------------------------------------------------------------ item 4 (on the fixture)
@pytest.fixture(scope="module")
def fx():
    prob, glob = fio.load_problem(FIXTURE)
    return prob, prob.data(), glob


def test_probe_start_recipe(fx):
    prob, D, _ = fx
    L = prob.L
    e = [20260923, 13, 3, 24, 1]
    x = ft.probe_start(D, e)
    rng = np.random.default_rng(np.random.SeedSequence(e))
    y = ft.start_values(D)
    y[L.alpha] = rng.uniform(-1.5, 1.5, L.G)
    y[L.r] = rng.uniform(-6.0, -2.0, L.G)
    y[L.kappa] = rng.uniform(-0.5, 0.5)
    y[[L.rho_S, L.rho_b, L.rho_d, L.omega_S]] = rng.uniform(-1.0, 1.0, 4)
    y[L.v] = 0.0
    hold = ~D.free3
    y[L.alpha[hold[L.alpha]]] = 0.0
    y[L.r[hold[L.r]]] = ft.START_R
    y[L.delta[hold[L.delta]]] = 0.0
    assert np.array_equal(x, np.clip(y, L.lb, L.ub))
    assert np.array_equal(x, y)                     # the clip does nothing with these ranges
    assert x[L.v] == 0.0 and np.all(x[L.delta] == 0.0)
    # 2G + 5 doubles consumed in that order
    r2 = np.random.default_rng(np.random.SeedSequence(e))
    u = r2.uniform(0, 1, 2 * L.G + 5)
    free_a = D.free3[L.alpha]
    assert np.allclose(x[L.alpha][free_a], (-1.5 + 3.0 * u[:L.G])[free_a])
    assert np.isclose(x[L.omega_S], -1.0 + 2.0 * u[2 * L.G + 4])
    # what fit_section3 builds from it is the same point (holds already applied)
    assert np.array_equal(x, ft.probe_start(D, list(e)))


# ------------------------------------------------------------------ items 6.1, 6.3, 6.4
def test_pair_has_no_rule4_restart(monkeypatch):
    calls = []

    def fake(D, start, v_start, v_scale, fix_v=False, tol=ft.TOL, nodes=20, rule4_restart=True):
        calls.append((v_start, rule4_restart))
        return F(True, -1.0, v=v_start)
    monkeypatch.setattr(ft, "fit_section4", fake)
    monkeypatch.setattr(ft, "eps_tie4", lambda D: 1e-6)
    ft.fit_section4_two_starts(None, np.zeros(2), 0.0025)
    assert calls == [(0.0, False), (0.0025, False)]


class _L:
    v = 1


class _D:
    L = _L()


def test_lr_v184():
    D = _D()
    assert ft.lr_v184(F(True, -5.0, v=0.0, ll_ref=-1.0), -10.0, D) == 0.0
    assert ft.lr_v184(F(True, -5.0, v=1e-8, ll_ref=-1.0), -10.0, D) == 0.0
    assert ft.lr_v184(F(True, -5.0, v=2e-8, ll_ref=-9.0), -10.0, D) == pytest.approx(2.0)
    # kept v > 0 with l4 <= l3*: LR = 0 (item 6.4, conflict C2)
    assert ft.lr_v184(F(True, -5.0, v=0.01, ll_ref=-10.5), -10.0, D) == 0.0


def test_ell3_star(monkeypatch):
    D = _D()
    monkeypatch.setattr(ft.ma, "loglik3", lambda th, D, order=1: -9.9 if th[1] == 0.0 else np.nan)
    m0 = F(True, -9.0, v=5e-9)
    m1 = F(True, -8.0, v=0.01)
    l3s, c = ft.ell3_star(-10.0, [m0, m1], D)
    assert l3s == -9.9 and c == [(0, -9.9), (1, None)]
    l3s, c = ft.ell3_star(-10.0, [F(False, -1.0, v=0.0), m1], D)   # unconverged member contributes nothing
    assert l3s == -10.0 and c == [(0, None), (1, None)]
