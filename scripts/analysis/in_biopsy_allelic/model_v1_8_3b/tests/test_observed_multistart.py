"""v1.8.3 unit tests: the H-11 M1 observed section-3 selection rule (n_obs = 2) and its seed.
v1.8.3a: the seed and only-one-converged tests follow H-12 m1 and m2 (see test_v183a.py for the new rules).
v1.8.3b: the only-one-converged test uses the 1e-6 failure margin (H-14 M1)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import fitting as ft  # noqa: E402
import fixture_io as fio  # noqa: E402

FIXTURE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/"
               "MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z")


class F:
    def __init__(self, conv, ll):
        self.converged, self.loglik = conv, ll


def sel(*fits):
    return ft.select_observed(list(fits))[0]


def test_both_converged_higher_wins():
    assert sel(F(True, -10.0), F(True, -9.0)) == 1
    assert sel(F(True, -9.0), F(True, -10.0)) == 0


def test_ties_within_1e9_go_to_section7_start():
    assert sel(F(True, -10.0), F(True, -10.0 + 5e-10)) == 0
    assert sel(F(True, -10.0), F(True, -10.0 + 9e-10)) == 0
    assert sel(F(True, -10.0), F(True, -10.0 + 3e-9)) == 1
    assert sel(F(True, -10.0 + 5e-10), F(True, -10.0)) == 0


def test_only_one_converged_is_reported_unless_the_other_is_higher():
    # v1.8.3a (H-12 m1), margin per v1.8.3b (H-14 M1): an unconverged start more than 1e-6 above the converged
    # one fails the observed fit; within 1e-6 the converged start is reported
    assert sel(F(False, -11.0), F(True, -10.0)) == 1
    assert sel(F(True, -10.0), F(False, -10.0 + 5e-10)) == 0
    assert sel(F(True, -10.0), F(False, -10.0 + 5e-7)) == 0
    assert sel(F(True, -10.0), F(False, -10.0 + 2e-6)) is None
    assert sel(F(False, -1.0), F(True, -10.0)) is None
    assert sel(F(True, -10.0), F(False, -1.0)) is None


def test_neither_converged_fails():
    idx, rule = ft.select_observed([F(False, -10.0), F(False, -9.0)])
    assert idx is None and "failed" in rule


def test_nonfinite_converged_start_is_not_eligible():
    assert sel(F(True, -np.inf), F(True, -10.0)) == 1
    assert sel(F(True, np.nan), F(False, -10.0)) is None


def test_seed_entropy():
    assert ft.OBS_START_STREAM == 7          # v1.8.3a, H-12 m2
    e = ft.obs_start_seed(20261015, 2, 16)
    assert e == [20261015, 2, 16, 7, 1]
    a = np.random.default_rng(np.random.SeedSequence(e)).normal(size=5)
    b = np.random.default_rng(np.random.SeedSequence(list(e))).normal(size=5)
    assert np.array_equal(a, b)
    c = np.random.default_rng(np.random.SeedSequence(ft.obs_start_seed(20261015, 2, 17))).normal(size=5)
    assert not np.array_equal(a, c)


def test_n_obs_is_fixed_at_two():
    with pytest.raises(ValueError):
        ft.fit_section3_observed(None, ft.FIXTURE_OBS_START_SEED, n_obs=3)


@pytest.fixture(scope="module")
def fx():
    prob, _ = fio.load_problem(FIXTURE)
    return prob, prob.data()


def test_wrapper_on_fixture(fx):
    """Start 0 of the wrapper is fit_section3(D) bit for bit; the jittered start is reproducible and lies
    inside the bounds; the record describes both starts."""
    prob, D = fx
    L = prob.L
    f0 = ft.fit_section3(D)
    y1 = ft.jittered_start(D, ft.FIXTURE_OBS_START_SEED)
    y2 = ft.jittered_start(D, ft.FIXTURE_OBS_START_SEED)
    assert np.array_equal(y1, y2)
    assert np.all(y1 >= L.lb) and np.all(y1 <= L.ub)
    x0 = ft.section3_start(D)
    assert np.array_equal(y1[~D.free3], x0[~D.free3])
    assert not np.array_equal(y1[D.free3], x0[D.free3])
    f = ft.fit_section3_observed(D, ft.FIXTURE_OBS_START_SEED)
    m = f.multistart
    assert m["n_obs"] == 2 and len(m["starts"]) == 2
    assert m["starts"][0]["dl_vs_start0"] == 0.0 and m["starts"][0]["bitwise_equal_start0"]
    k = m["reported_start"]
    if k == 0:
        assert np.array_equal(f.theta, f0.theta) and f.loglik == f0.loglik
    else:
        assert f.converged and m["starts"][1]["dl_vs_start0"] > ft.OBS_TIE
    assert f.converged == (not m["failed"])
