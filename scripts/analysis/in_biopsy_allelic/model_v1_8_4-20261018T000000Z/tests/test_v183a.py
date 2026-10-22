"""v1.8.3a unit tests (v1.8.3b: the m1 failure margin is 1e-6 per H-14 M1): LEAD_RULINGS_W4 H-12 m1 (failure rule for an unconverged start or restart of higher l),
m2 (Gate B stream 7, full-CRT cell_id 99), m3 (B-FIX and real-data second start [20260923, 7, 0]), m4/m6 (timing
field and run header of run_fixture.py)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import fitting as ft  # noqa: E402

HERE = Path(__file__).resolve().parents[1]


class F:
    def __init__(self, conv, ll, log=None):
        self.converged, self.loglik = conv, ll
        self.attempt_log = [] if log is None else log


# ------------------------------------------------------------------ H-12 m1, between starts
def test_m1_unconverged_start_higher_fails():
    idx, rule = ft.select_observed([F(True, -10.0), F(False, -9.0)])
    assert idx is None and "H-12 m1" in rule and "failed" in rule
    idx, rule = ft.select_observed([F(False, -9.0), F(True, -10.0)])
    assert idx is None and "H-12 m1" in rule


def test_m1_failure_margin_is_1e6_not_the_tie():
    # v1.8.3b (H-14 M1; spec v1.8.3 item 3.4 step 3, lines 1759-1772): the failure margin is 1e-6, not OBS_TIE
    assert ft.OBS_FAIL_MARGIN == 1e-6 and ft.OBS_TIE == 1e-9
    assert ft.select_observed([F(True, -10.0), F(False, -10.0 + 0.9e-9)])[0] == 0
    assert ft.select_observed([F(True, -10.0), F(False, -10.0 + 1.1e-9)])[0] == 0
    assert ft.select_observed([F(True, -10.0), F(False, -10.0 + 0.9e-6)])[0] == 0
    assert ft.select_observed([F(True, -10.0), F(False, -10.0 + 1.1e-6)])[0] is None
    assert ft.select_observed([F(True, -10.0), F(False, -np.inf)])[0] == 0
    assert ft.select_observed([F(True, -10.0), F(False, np.nan)])[0] == 0


def test_m1_compares_against_the_reported_start():
    # start 1 reported (higher); start 0 converged lower: no failure
    assert ft.select_observed([F(True, -10.0), F(True, -9.0)])[0] == 1
    # reported start 1 at -9; an unconverged attempt of start 0 at -9.5 is below it: no failure
    assert ft.select_observed([F(True, -10.0, [(-9.5, False), (-10.0, True)]), F(True, -9.0)])[0] == 1


def test_v183b_unconverged_start_5e7_higher_does_not_fail_2e6_does():
    # v1.8.3b (H-14 M1): an unconverged START 5e-7 above the candidate is inside the 1e-6 margin (no failure);
    # 2e-6 above fails. Both orders of the starts.
    assert ft.select_observed([F(True, -10.0), F(False, -10.0 + 5e-7)])[0] == 0
    assert ft.select_observed([F(False, -10.0 + 5e-7), F(True, -10.0)])[0] == 1
    idx, rule = ft.select_observed([F(True, -10.0), F(False, -10.0 + 2e-6)])
    assert idx is None and "H-12 m1" in rule and "1e-06" in rule
    assert ft.select_observed([F(False, -10.0 + 2e-6), F(True, -10.0)])[0] is None


def test_v183b_unconverged_restart_5e7_higher_does_not_fail_2e6_does():
    # v1.8.3b (H-14 M1): the same margin for an unconverged restart (attempt) inside a converged start
    s1 = F(True, -10.5, [(-10.5, True)])
    assert ft.select_observed([F(True, -10.0, [(-10.0 + 5e-7, False), (-10.0, True)]), s1])[0] == 0
    idx, rule = ft.select_observed([F(True, -10.0, [(-10.0 + 2e-6, False), (-10.0, True)]), s1])
    assert idx is None and "restart within start 0" in rule


def test_v183b_tie_between_converged_starts_still_1e9():
    # OBS_TIE (1e-9) still decides between two converged starts; the 1e-6 margin is not used there
    assert ft.select_observed([F(True, -10.0), F(True, -10.0 + 5e-7)])[0] == 1
    assert ft.select_observed([F(True, -10.0), F(True, -10.0 + 0.9e-9)])[0] == 0


# ------------------------------------------------------------------ H-12 m1, restarts within a start
def test_m1_unconverged_restart_within_reported_start_fails():
    s0 = F(True, -10.0, [(-9.0, False), (-10.0, True)])
    s1 = F(True, -10.5, [(-10.5, True)])
    idx, rule = ft.select_observed([s0, s1])
    assert idx is None and "restart within start 0" in rule


def test_m1_unconverged_restart_within_other_start_fails():
    s0 = F(True, -10.0, [(-10.0, True)])
    s1 = F(True, -10.5, [(-9.9, False), (-10.5, True)])
    idx, rule = ft.select_observed([s0, s1])
    assert idx is None and "start 1" in rule


def test_m1_converged_restarts_only_no_failure():
    s0 = F(True, -10.0, [(-10.0, True)])
    s1 = F(True, -10.0, [(-10.0, True)])
    assert ft.select_observed([s0, s1])[0] == 0


def test_run_records_attempt_log(monkeypatch):
    """_run writes (l, converged) for every attempt; an earlier unconverged attempt of higher l is kept in the
    log although the converged attempt is the one returned (the within-fit rule of v1.8.2 is unchanged)."""
    outs = iter([
        (np.zeros(3), -1.0, np.zeros(3), 1.0, "m", True, "ok", False),      # attempt 1: pg 1.0 >= tol
        (np.ones(3), -2.0, np.zeros(3), 1e-9, "m", True, "ok", False),      # attempt 2: converged
    ])
    monkeypatch.setattr(ft, "maximize", lambda *a, **k: next(outs))
    monkeypatch.setattr(ft, "hessian_scale", lambda ev, x0, free, base: base)

    class Ev:
        is4 = False
        last_choice = None
    lb, ub = np.full(3, -5.0), np.full(3, 5.0)
    f = ft._run(Ev(), np.zeros(3), lb, ub, np.ones(3, bool), np.ones(3), ft.TOL)
    assert f.converged and f.loglik == -2.0 and f.attempts == 2
    assert f.attempt_log == [(-1.0, False), (-2.0, True)]
    idx, rule = ft.select_observed([f, F(True, -3.0)])
    assert idx is None and "restart within start 0" in rule


def test_failed_observed_fit_is_marked_not_converged(monkeypatch):
    fits = iter([ft.Fit(theta=np.zeros(2), loglik=-10.0, grad=np.zeros(2), proj_max_grad=0.0, converged=True,
                        attempts=2, seconds=1.0, attempt_log=[(-9.0, False), (-10.0, True)]),
                 ft.Fit(theta=np.ones(2), loglik=-11.0, grad=np.zeros(2), proj_max_grad=0.0, converged=True,
                        attempts=1, seconds=2.0, attempt_log=[(-11.0, True)])])
    monkeypatch.setattr(ft, "fit_section3", lambda D, start=None, tol=ft.TOL: next(fits))
    monkeypatch.setattr(ft, "jittered_start", lambda D, seed: None)
    # v1.8.4: K_obs = 3; the third (probe) start is a converged fit far below
    fits = iter(list(fits) + [ft.Fit(theta=np.ones(2), loglik=-50.0, grad=np.zeros(2), proj_max_grad=0.0,
                                     converged=True, attempts=1, seconds=0.0, attempt_log=[(-50.0, True)])])
    monkeypatch.setattr(ft, "probe_start", lambda D, e: None)
    f = ft.fit_section3_observed(None, ft.FIXTURE_OBS_START_SEED, ft.probe_entropy_fixed_observed(1))
    assert not f.converged and f.multistart["failed"]
    assert "H-12 m1" in f.status
    assert f.multistart["seconds_total"] == 3.0
    assert f.multistart["starts"][0]["attempts_dl_vs_start0"] == [[1.0, False], [0.0, True]]


# ------------------------------------------------------------------ H-12 m2, m3 seeds
def test_m2_gate_b_stream_7():
    assert ft.OBS_START_STREAM == 7
    assert ft.obs_start_seed(20261015, 3, 12) == [20261015, 3, 12, 7, 1]
    assert ft.FULL_CRT_CELL_ID == 99
    assert ft.obs_start_seed_full_crt(20261015, 5) == [20261015, 99, 5, 7, 1]
    a = np.random.default_rng(np.random.SeedSequence([20261015, 3, 12, 7, 1])).normal(size=4)
    b = np.random.default_rng(np.random.SeedSequence(ft.obs_start_seed(20261015, 3, 12))).normal(size=4)
    assert np.array_equal(a, b)


def test_m3_bfix_and_real_data_seed():
    assert ft.FIXTURE_OBS_START_SEED == [20260923, 7, 0]
    assert ft.REAL_DATA_OBS_START_SEED == [20260923, 7, 0]
    assert ft.FIXTURE_OBS_START_SEED[:2] != ft.RESTART_SEED[:2]


def test_old_72bit_stream_is_gone():
    src = (HERE / "fitting.py").read_text()
    assert 'b"obs_start"' not in src
    assert "2054682710969277313652" not in src


# ------------------------------------------------------------------ m4, m6 run_fixture header and timing
def test_run_fixture_header_and_timing_fields():
    src = (HERE / "run_fixture.py").read_text()
    assert "header=run_header()" in src
    assert "section3_observed_fit_both_starts" in src
    sys.path.insert(0, str(HERE))
    import run_fixture as rf
    h = rf.run_header()
    for k in ("host", "cpu_model", "glibc", "env", "backend", "numpy", "scipy"):
        assert k in h
    assert set(h["env"]) >= {"OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"}
