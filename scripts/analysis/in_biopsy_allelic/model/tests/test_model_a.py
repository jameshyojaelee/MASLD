"""Unit tests for the Model A-final likelihood (spec v1 sections 1.3, 2.2, 4 and 7)."""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.special import expit, logit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import betabinom as bbm  # noqa: E402
import fitting as ft  # noqa: E402
import fixture_io as fio  # noqa: E402
import model_a as ma  # noqa: E402
from het_rule import het, truncation_bounds  # noqa: E402

FIXTURE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/"
               "MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z")
TAU_MAX = 0.05
RULES = ([(d, k, F) for d in (10, 20, 30) for k in (1, 2, 3, 4, 5) for F in (1, 2, 5, 10, 15, 20)]
         + [(10, 2, 10), (0, 0, 0), (1, 1, 50), (0, 0, 50), (5, 0, 33)])


@pytest.fixture(scope="module")
def fx():
    prob, glob = fio.load_problem(FIXTURE)
    params = fio.load_params(FIXTURE, prob.L)
    return prob, prob.data(), params


def fd4(f, x, i, h):
    """Fourth-order central difference of f along coordinate i."""
    def at(s):
        y = x.copy()
        y[i] += s * h
        return f(y)
    return (-at(2) + 8 * at(1) - 8 * at(-1) + at(-2)) / (12 * h)


# --- section 1.3 ------------------------------------------------------------------------------

def test_truncation_set_matches_het_enumeration():
    n = np.arange(0, 2001)
    seg, a = bbm.seg_ranges(np.zeros_like(n), n + 1)
    N = n[seg]
    for d, k, F in RULES:
        lo, hi, ne = truncation_bounds(n, d, k, F)
        in_k = ne[N] & (a >= lo[N]) & (a <= hi[N])
        assert np.array_equal(in_k, het(N - a, a, d, k, F)), (d, k, F)
        assert np.array_equal(in_k, het(a, N - a, d, k, F)), (d, k, F)


# --- section 2.2 ------------------------------------------------------------------------------

def test_hand_example():
    lo, hi, ne = truncation_bounds(10, 10, 2, 10)
    assert (int(lo), int(hi), bool(ne)) == (2, 8, True)
    # orientation -1 -> s_t = +1; a = tag-REF reads = 7; logit p = eta + s_t * omega
    p = expit(0.5 + 1 * 0.1)
    assert abs(p - 0.64566) < 1e-5
    f_binom = np.exp(bbm.bb_logpmf(7, 10, p, 0.0))
    assert abs(f_binom - 0.24973) < 1e-4
    # production path: BB at the rho floor (x2 -> -inf gives rho = 1e-6)
    st = ma._RowState(np.array([0.6]), np.array([-40.0]), np.array([10]), 0)
    f_bb = np.exp(st.terms(np.array([0]), np.array([7]), bbm.lchoose(10, 7), 0)[0])[0]
    assert abs(f_bb - 0.24973) < 1e-4


def _random_rows(rng, nrow, rule=(10, 2, 15), orient=None, a=None, nmax=400):
    n = rng.integers(rule[0], nmax, nrow)
    lo, hi, _ = truncation_bounds(n, *rule)
    if a is None:
        a = lo + (rng.random(nrow) * (hi - lo + 1)).astype(int)
    if orient is None:
        orient = rng.choice([-1, 1], nrow)
    p_x = rng.uniform(0, 0.05, nrow)
    p_ref = rng.uniform(0.05, 0.9, nrow) * (1 - p_x)
    p_alt = 1 - p_x - p_ref
    e = rng.uniform(1e-4, 0.01, nrow)
    return dict(n=n, a=a, orient=orient, p_ref=p_ref, p_alt=p_alt, p_x=p_x, e=e, rule=rule)


def _rows(d, phi=0.0439, flip=False):
    a = d["n"] - d["a"] if flip else d["a"]
    o = -d["orient"] if flip else d["orient"]
    return ma.Rows(d["n"], a, *d["rule"], o, d["p_ref"], d["p_alt"], d["p_x"], d["e"], phi)


def test_orientation_invariance_rows():
    rng = np.random.default_rng(11)
    d = _random_rows(rng, 300)
    x1 = rng.uniform(-3, 3, 300)
    x2 = rng.uniform(-6, -1, 300)
    r1 = ma.row_eval(x1, x2, _rows(d), 2)
    r2 = ma.row_eval(-x1, x2, _rows(d, flip=True), 2)
    assert np.allclose(r1.logL, r2.logL, rtol=0, atol=1e-11)
    assert np.allclose(r1.g1, -r2.g1, rtol=1e-9, atol=1e-11)
    assert np.allclose(r1.g2, r2.g2, rtol=1e-9, atol=1e-11)
    assert np.allclose(r1.g11, r2.g11, rtol=1e-9, atol=1e-11)


def test_orientation_invariance_fixture(fx):
    prob, D, params = fx
    rows = fio.read_tsv(FIXTURE / "rows.tsv")
    col = lambda k, f=float: np.array([f(r[k]) for r in rows])  # noqa: E731
    cohorts = {r["cohort"]: r for r in fio.read_tsv(FIXTURE / "cohorts.tsv")}
    rule = np.array([[int(cohorts[r["cohort"]][k]) for k in ("min_dp", "min_minor_reads",
                                                              "min_minor_pct")] for r in rows])
    n, a, o = col("n", int), col("a", int), col("orientation", int)
    Rf = ma.Rows(n, n - a, rule[:, 0], rule[:, 1], rule[:, 2], -o, col("p_ref"), col("p_alt"),
                 col("p_x"), col("e_i"), 0.0439)
    Df = prob.data(rows=Rf)
    L = prob.L
    for name in ("theta2_truth_v0", "theta4_perturbed"):
        th = params[name].copy()
        thf = th.copy()
        thf[L.alpha] = -th[L.alpha]
        thf[L.omega_S] = th[L.omega_S]      # s_t flips with orientation, so omega terms flip too
        assert abs(ma.loglik3(th, D, 0) - ma.loglik3(thf, Df, 0)) < 1e-9
        th[L.v] = thf[L.v] = 0.0036
        assert abs(ma.loglik4(th, D, 0) - ma.loglik4(thf, Df, 0)) < 1e-9


# --- section 7 normalizers --------------------------------------------------------------------

def test_C_R_equals_C_A():
    rng = np.random.default_rng(5)
    for rule in [(10, 1, 10), (10, 2, 15), (10, 1, 5), (20, 3, 20)]:
        n = rng.integers(rule[0], 3000, 200)
        lo, hi, _ = truncation_bounds(n, *rule)
        e = rng.uniform(1e-4, 0.02, 200)
        for phi in (0.0, 0.0439, 0.2):
            cr = bbm.log_mass_over_K(n, lo, hi, e, phi)
            ca = bbm.log_mass_over_K(n, lo, hi, 1 - e, phi)
            assert np.allclose(cr, ca, rtol=1e-12, atol=1e-12)


def test_gamma_family_differences_against_exact_sums():
    # for integer m: lgamma(x+m) - lgamma(x) = sum_j log(x+j), psi difference = sum_j 1/(x+j),
    # trigamma difference = -sum_j 1/(x+j)^2; math.fsum gives correctly rounded sums
    import math
    for x in (0.0021, 0.37, 12.5, 29.99, 30.0, 999.0, 2.5e4, 7.7e5, 1e6):
        for m in (0, 1, 3, 40, 2118, 6821):
            j = range(m)
            lg = math.fsum(math.log(x + k) for k in j)
            dg = math.fsum(1.0 / (x + k) for k in j)
            tg = -math.fsum(1.0 / (x + k) ** 2 for k in j)
            assert abs(bbm.lgamma_diff(x, m) - lg) <= 1e-13 * max(1.0, abs(lg)), (x, m)
            assert abs(bbm.digamma_diff(x, m) - dg) <= 1e-11 * max(1e-3, abs(dg)), (x, m)
            assert abs(bbm.trigamma_diff(x, m) - tg) <= 1e-11 * abs(tg) + 1e-30, (x, m)


def test_recurrence_matches_direct():
    for n, mu, rho in [(10, 0.3, 0.02), (300, 0.64, 0.03), (2500, 0.5, 1e-5), (6821, 0.2, 0.1)]:
        rec = bbm.bb_logpmf_recurrence(n, mu, rho, 0, n)
        direct = bbm.bb_logpmf(np.arange(n + 1), n, mu, rho)
        assert np.allclose(rec, direct, rtol=0, atol=1e-9)


def test_normalizer_tail_and_direct_agree_and_derivatives():
    rng = np.random.default_rng(7)
    d = _random_rows(rng, 400, rule=(10, 2, 15), nmax=3000)
    d["n"][:40] = rng.integers(10, 16, 40)             # small n with extreme p -> T >= 0.5 rows
    lo, hi, _ = truncation_bounds(d["n"], *d["rule"])
    d["a"] = lo
    R = _rows(d)
    x1 = rng.uniform(-3, 3, 400)
    x1[:40] = rng.choice([-6.0, 6.0], 40)
    x2 = rng.uniform(-8, -0.5, 400)
    st = ma._RowState(x1, x2, R.n, 2)
    lc = ma.log_C_H(st, R, 2)
    rho = ma.EPS_RHO + (1 - 2 * ma.EPS_RHO) * expit(x2)
    direct = bbm.log_mass_over_K(R.n, R.lo, R.hi, expit(x1), rho)
    tl = st.terms(R.tail_row, R.tail_a, R.tail_lchoose, 0)[0]
    T = np.exp(bbm.seg_logsumexp(tl, R.tail_row, 400))
    assert (T >= 0.5).sum() >= 5 and (T < 0.5).sum() >= 50   # both branches are exercised
    assert np.allclose(lc[0], direct, rtol=1e-11, atol=1e-11)
    h = 1e-4
    for k, (dx1, dx2) in enumerate([(h, 0), (0, h)]):
        fp = ma.log_C_H(ma._RowState(x1 + dx1, x2 + dx2, R.n, 0), R, 0)[0]
        fm = ma.log_C_H(ma._RowState(x1 - dx1, x2 - dx2, R.n, 0), R, 0)[0]
        assert np.allclose((fp - fm) / (2 * h), lc[1 + k], rtol=1e-5, atol=1e-7)


def test_row_derivatives_fd():
    rng = np.random.default_rng(3)
    d = _random_rows(rng, 200, nmax=2000)
    R = _rows(d)
    x1 = rng.uniform(-4, 4, 200)
    x2 = rng.uniform(-9, -0.5, 200)
    r = ma.row_eval(x1, x2, R, 2)
    h = 1e-3
    f = lambda a, b: ma.row_eval(a, b, R, 1)  # noqa: E731

    def d4(fun, dx1, dx2):
        return (-fun(x1 + 2 * dx1, x2 + 2 * dx2) + 8 * fun(x1 + dx1, x2 + dx2)
                - 8 * fun(x1 - dx1, x2 - dx2) + fun(x1 - 2 * dx1, x2 - 2 * dx2)) / (12 * h)
    ll = lambda a, b: ma.row_eval(a, b, R, 0).logL  # noqa: E731
    # FD rounding: log L carries ~1e-14 * n absolute rounding from lgamma terms of size ~n log n,
    # which the stencil amplifies by ~1.5/h; the first-derivative atol scales with n for that
    at_n = 1e-8 + 1e-10 * R.n
    checks = [(d4(ll, h, 0), r.g1, 1e-6, at_n), (d4(ll, 0, h), r.g2, 1e-6, at_n),
              (d4(lambda a, b: f(a, b).g1, h, 0), r.g11, 1e-5, 1e-7),
              (d4(lambda a, b: f(a, b).g1, 0, h), r.g12, 1e-5, 1e-7),
              (d4(lambda a, b: f(a, b).g2, 0, h), r.g22, 1e-5, 1e-7)]
    for k, (num, ana, rt, at) in enumerate(checks):
        bad = np.abs(num - ana) > at + rt * np.abs(ana)
        assert not bad.any(), (k, num[bad][:3], ana[bad][:3], x1[bad][:3], x2[bad][:3])


# --- gradients and Hessians in section 7 coordinates -----------------------------------------

def _random_theta(L, rng):
    th = np.zeros(L.size)
    th[L.alpha] = rng.uniform(-1.5, 1.5, L.G)
    th[L.r] = rng.uniform(-5, -2, L.G)
    th[L.kappa] = 0.2
    th[L.delta] = rng.normal(0, 0.3, L.p)
    th[L.omega_S], th[L.rho_S], th[L.rho_b], th[L.rho_d] = 0.05, -0.1, 0.2, 0.3
    th[L.v] = 0.003
    return th


def test_gradient_section3_fd(fx):
    prob, D, params = fx
    L = prob.L
    thetas = [params["theta1_start_like"], params["theta2_truth_v0"], params["theta4_perturbed"],
              _random_theta(L, np.random.default_rng(2))]
    for th in thetas:
        ll, g = ma.loglik3(th, D, 1)
        for i in np.flatnonzero(D.free3):
            num = fd4(lambda x: ma.loglik3(x, D, 0), th, i, 1e-3)
            assert abs(num - g[i]) <= 1e-8 + 1e-6 * abs(g[i]), (L.names[i], num, g[i])


def test_gradient_section4_fd(fx):
    prob, D, params = fx
    L = prob.L
    thetas = [params["theta3_truth_vpos"], params["theta4_perturbed"],
              _random_theta(L, np.random.default_rng(4))]
    for th in thetas:
        ll, g = ma.loglik4(th, D, 1)
        for i in np.flatnonzero(D.free4):
            h = 1e-3 * th[L.v] if i == L.v else 1e-3
            num = fd4(lambda x: ma.loglik4(x, D, 0), th, i, h)
            assert abs(num - g[i]) <= 1e-8 + 1e-6 * abs(g[i]), (L.names[i], num, g[i])


def test_kappa_A_gradient_fd():
    prob, _ = fio.load_problem(FIXTURE)
    probA = ma.Problem(prob.ind_ids, prob.cohort, prob.S, prob.zraw, prob.base_names, prob.cluster,
                       prob.rows, prob.row_ind, prob.row_gene, prob.L.gene_names, prob.omega_row,
                       prob.w_fib, prob.w_mye, kappa_A=True)
    D = probA.data()
    L = probA.L
    th = _random_theta(L, np.random.default_rng(8))
    th[L.gamma_fib], th[L.gamma_mye] = 0.3, -0.2
    ll, g, H = ma.loglik3(th, D, 2)
    for i in (L.gamma_fib, L.gamma_mye, L.kappa, L.alpha[2]):
        num = fd4(lambda x: ma.loglik3(x, D, 0), th, i, 1e-3)
        assert abs(num - g[i]) <= 1e-8 + 1e-6 * abs(g[i]), L.names[i]
        numh = np.array([fd4(lambda x: ma.loglik3(x, D, 1)[1][j], th, i, 1e-3)
                         for j in (L.gamma_fib, L.alpha[2])])
        assert np.allclose(numh, H[i, [L.gamma_fib, L.alpha[2]]], rtol=1e-6, atol=1e-6)


def test_hessian_section3_fd(fx):
    prob, D, params = fx
    th = params["theta4_perturbed"]
    ll, g, H = ma.loglik3(th, D, 2)
    idx = np.flatnonzero(D.free3)
    for i in idx:
        col = fd4(lambda x: ma.loglik3(x, D, 1)[1], th, i, 1e-3)
        assert np.allclose(col[idx], H[idx, i], rtol=1e-6, atol=1e-6), prob.L.names[i]


def test_hessian_section4_fd(fx):
    prob, D, params = fx
    L = prob.L
    th = params["theta3_truth_vpos"]
    ll, g, H = ma.loglik4(th, D, 2)
    idx = np.flatnonzero(D.free4)
    for i in idx:
        h = 1e-3 * th[L.v] if i == L.v else 1e-3
        col = fd4(lambda x: ma.loglik4(x, D, 1)[1], th, i, h)
        scale = np.sqrt(np.abs(np.diag(H))[idx] * abs(H[i, i]))
        assert np.all(np.abs(col[idx] - H[idx, i]) <= 1e-5 * scale + 1e-6), L.names[i]


def test_row_scores_sum_to_gradient(fx):
    prob, D, params = fx
    th = params["theta4_perturbed"]
    ll, g = ma.loglik3(th, D, 1)
    S = ma.row_scores(th, D)
    assert np.allclose(np.asarray(S.sum(axis=0)).ravel(), g, rtol=1e-12, atol=1e-10)


# --- section 4 --------------------------------------------------------------------------------

def test_tau_to_zero_reduction(fx):
    prob, D, params = fx
    L = prob.L
    th = params["theta2_truth_v0"].copy()
    l3 = ma.loglik3(th, D, 0)
    score = ma.boundary_score(th, D)
    th[L.v] = 1e-12
    assert abs(ma.loglik4(th, D, 0) - l3) < 1e-9
    v = 1e-8
    th[L.v] = v
    l4, g4 = ma.loglik4(th, D, 1)
    assert abs((l4 - l3) / v - score) <= 1e-4 * max(1.0, abs(score))
    assert abs(g4[L.v] - score) <= 1e-4 * max(1.0, abs(score))
    th[L.v] = 0.0
    l0, g0 = ma.loglik4(th, D, 1)
    assert l0 == l3 and g0[L.v] == score


def _synthetic_problem(sizes, seed=21):
    """Genes with the given row counts; BB counts drawn and kept in K (a known-answer design,
    not the B-FIX generator)."""
    rng = np.random.default_rng(seed)
    n_ind = max(sizes)
    coh = np.array(["A", "B", "C"])[rng.integers(0, 3, n_ind)]
    S = rng.integers(0, 4, n_ind).astype(float)
    zraw = rng.normal(0, 1, (n_ind, 6))
    zraw[:, 0] = rng.normal(-6, 0.5, n_ind)
    ind_ids = [f"I{i:04d}" for i in range(n_ind)]
    rules = {"A": (10, 1, 10), "B": (10, 2, 15), "C": (10, 1, 5)}
    rows = []
    for g, m in enumerate(sizes):
        for i in np.sort(rng.choice(n_ind, m, replace=False)):
            rows.append((i, g))
    rows.sort()
    ri = np.array([r[0] for r in rows])
    rg = np.array([r[1] for r in rows])
    rule = np.array([rules[coh[i]] for i in ri])
    n = np.maximum(10, np.round(np.exp(rng.normal(3.8, 0.9, len(ri))))).astype(int)
    lo, hi, _ = truncation_bounds(n, rule[:, 0], rule[:, 1], rule[:, 2])
    alpha = rng.uniform(-1, 1, len(sizes))
    mu = expit(alpha[rg])
    a = np.empty(len(ri), int)
    for k in range(len(ri)):
        while True:
            p = rng.beta(mu[k] * 32, (1 - mu[k]) * 32)
            x = rng.binomial(n[k], p)
            if lo[k] <= x <= hi[k]:
                a[k] = x
                break
    orient = rng.choice([-1, 1], len(ri))
    p_ref = rng.uniform(0.2, 0.8, len(ri))
    R = ma.Rows(n, a, rule[:, 0], rule[:, 1], rule[:, 2], orient, p_ref, 1 - p_ref,
                np.zeros(len(ri)), np.full(len(ri), 0.002), 0.0439)
    omega = np.array([{"A": 0.2, "B": 0.1, "C": 0.05}[coh[i]] for i in ri])
    prob = ma.Problem(ind_ids, coh, S, zraw, fio.BASE_Z, np.array(ind_ids), R, ri, rg,
                      [f"S{g}" for g in range(len(sizes))], omega)
    th = np.zeros(prob.L.size)
    L = prob.L
    th[L.alpha] = alpha
    th[L.r] = logit(0.03)
    th[L.kappa], th[L.omega_S], th[L.rho_S] = -0.1, 0.02, 0.1
    th[L.delta] = rng.normal(0, 0.05, L.p)
    return prob, th


def test_adaptive_gh_matches_quad():
    sizes = [5, 20, 50, 200, 500]
    prob, th = _synthetic_problem(sizes)
    D = prob.data()
    L = prob.L
    for mult in (0.25, 1.0, 4.0):
        v = (mult * TAU_MAX) ** 2
        th[L.v] = v
        _, det = ma.loglik4(th, D, 0, detail=True)
        for g in range(len(sizes)):
            ref, err = ma.gh_gene_integral_quad(th, D, g, v)
            assert abs(det["ell_g"][g] - ref) <= 1e-8 * abs(ref), (sizes[g], mult, det["ell_g"][g], ref)


def test_adaptive_gh_matches_quad_fixture_genes(fx):
    prob, D, params = fx
    L = prob.L
    th = params["theta4_perturbed"].copy()
    for mult in (0.25, 1.0, 4.0):
        v = (mult * TAU_MAX) ** 2
        th[L.v] = v
        _, det = ma.loglik4(th, D, 0, detail=True)
        for g in range(L.G):
            ref, err = ma.gh_gene_integral_quad(th, D, g, v)
            assert abs(det["ell_g"][g] - ref) <= 1e-8 * abs(ref), (L.gene_names[g], mult)


def test_adaptive_gh_20_vs_40_fixture(fx):
    prob, D, params = fx
    L = prob.L
    for name in ("theta3_truth_vpos", "theta4_perturbed"):
        th = params[name].copy()
        for mult in (0.25, 1.0, 4.0):
            th[L.v] = (mult * TAU_MAX) ** 2
            _, d20 = ma.loglik4(th, D, 0, nodes=20, detail=True)
            _, d40 = ma.loglik4(th, D, 0, nodes=40, detail=True)
            assert np.allclose(d20["ell_g"], d40["ell_g"], rtol=1e-8, atol=0)


# --- design -----------------------------------------------------------------------------------

def test_centring_with_multiplicities_equals_duplication():
    rng = np.random.default_rng(9)
    coh = np.array(["A"] * 6 + ["B"] * 5)
    S = rng.integers(0, 4, 11).astype(float)
    z = rng.normal(0, 1, (11, 2))
    z[1, 1] = z[7, 1] = np.nan
    mult = np.array([0, 1, 2, 1, 3, 1, 1, 0, 2, 1, 1], float)
    spec = ma.design_spec(coh, z, ["z1", "z2"])
    assert [s[2] for s in spec] == ["z1", "z2", "z2_missing"]
    rep = np.repeat(np.arange(11), mult.astype(int))
    St = ma.centre_stage(coh, S, mult)
    St_dup = ma.centre_stage(coh[rep], S[rep], np.ones(len(rep)))
    assert np.allclose(St[rep], St_dup, atol=1e-14)
    Z = ma.centre_z(coh, z, mult, spec)
    Z_dup = ma.centre_z(coh[rep], z[rep], np.ones(len(rep)), spec)
    assert np.allclose(Z[rep], Z_dup, atol=1e-14)


def test_fixture_design_columns(fx):
    prob, D, params = fx
    assert prob.design_columns == ["z_b_log_e", "z_d_dup", "z_ffpe", "z_sex_female", "z_age10",
                                   "z_c", "z_ffpe_missing"]
    # z_age10 is all-missing in C3 and must be 0 there after centring
    c3 = prob.cohort[D.ind] == "C3"
    assert np.all(D.Z[c3, 4] == 0.0)


def test_start_values_fixture(fx):
    prob, D, params = fx
    th = ft.start_values(D)
    L = prob.L
    assert set(np.unique(th[L.alpha])) <= {-1.0, 0.0, 1.0}
    assert np.all(th[L.r] == logit(0.02))
