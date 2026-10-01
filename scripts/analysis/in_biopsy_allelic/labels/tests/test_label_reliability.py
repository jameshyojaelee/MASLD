"""Checks of the label-reliability probability model on simulated data (no project data read).

Run: PYTHONNOUSERSITE=1 python -m pytest -q scripts/analysis/in_biopsy_allelic/labels/tests
"""
import importlib.util
from pathlib import Path

import numpy as np
from scipy.stats import multivariate_normal, norm

spec = importlib.util.spec_from_file_location(
    "label_reliability", Path(__file__).resolve().parents[1] / "label_reliability.py")
lr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lr)

SEED = 20260926


def simulate_pairs(rng, n, mu, sd, rho, delta):
    var_T, var_e = rho * sd ** 2, (1 - rho) * sd ** 2
    T = rng.normal(mu, np.sqrt(var_T), n)
    z1 = T + rng.normal(0, np.sqrt(var_e), n)
    z2 = T + delta + rng.normal(0, np.sqrt(var_e), n)
    return np.digitize(z1, lr.CUTS_Y), np.digitize(z2, lr.CUTS_Y)


def test_bvn_cdf_matches_scipy_including_zero_arguments():
    rng = np.random.default_rng(SEED)
    pts = [(0.0, 0.7, 0.3), (-1.2, 0.0, -0.6), (0.0, 0.0, 0.5), (2.5, -2.5, 0.9)]
    pts += list(zip(rng.uniform(-4, 4, 40), rng.uniform(-4, 4, 40), rng.uniform(-0.95, 0.95, 40)))
    for h, k, r in pts:
        ref = multivariate_normal.cdf([h, k], mean=[0, 0], cov=[[1, r], [r, 1]], abseps=1e-12, releps=1e-12)
        assert abs(lr.bvn_cdf(h, k, r) - ref) < 1e-7, (h, k, r)


def test_cell_probs_sum_to_one_and_have_probit_marginals():
    P = lr.cell_probs(lr.CUTS_Y, 1.3, 1.6, 1.1, 0.55)
    assert abs(P.sum() - 1) < 1e-10
    m1 = np.diff(np.concatenate([[0], norm.cdf((lr.CUTS_Y - 1.3) / 1.1), [1]]))
    m2 = np.diff(np.concatenate([[0], norm.cdf((lr.CUTS_Y - 1.6) / 1.1), [1]]))
    assert np.allclose(P.sum(1), m1, atol=1e-10) and np.allclose(P.sum(0), m2, atol=1e-10)


def test_fit_recovers_planted_parameters():
    rng = np.random.default_rng(SEED)
    truth = dict(mu=1.4, sd=1.05, rho=0.62, delta=0.25)
    y1, y2 = simulate_pairs(rng, 40000, **truth)
    fit = lr.fit_pairs(lr.pair_table(y1, y2))
    mu, sd, rho, delta = lr.unpack(fit["theta"])
    assert fit["converged"] and not fit["at_bound"]
    assert abs(mu - truth["mu"]) < 0.03 and abs(sd - truth["sd"]) < 0.03
    assert abs(rho - truth["rho"]) < 0.03 and abs(delta - truth["delta"]) < 0.03


def test_implied_within_pair_sd_and_lambda_match_monte_carlo():
    rng = np.random.default_rng(SEED + 1)
    mu, sd, rho, delta = 1.2, 1.0, 0.55, 0.1
    y1, y2 = simulate_pairs(rng, 400000, mu, sd, rho, delta)
    row = lr.summarise_fit(np.array([mu, np.log(sd), np.arctanh(rho), delta]),
                           [dict(key="X", set="F", cohort="X", scale="Y5", n=2, mu_Z=mu, sd_Z=sd)])
    assert abs(row["implied_within_pair_sd"] - np.sqrt(((y2 - y1) ** 2).mean() / 2)) < 0.005
    s1b, s2b = (lr.harmonize(v) for v in simulate_pairs(rng, 400000, mu, sd, rho, 0.0))
    assert abs(row["lambda_in_deposit"] - np.corrcoef(s1b, s2b)[0, 1]) < 0.01
    # the same deposit's var_e carried to a target with the same latent SD gives the same lambda
    assert abs(row["lambda_X"] - row["lambda_in_deposit"]) < 1e-10


def test_lambda_is_the_slope_attenuation_factor():
    """Regress an outcome that depends on E[S|T] on one noisy S read; the slope shrinks by lambda."""
    rng = np.random.default_rng(SEED + 2)
    mu, sd, rho = 1.5, 1.1, 0.6
    n = 400000
    T = rng.normal(mu, np.sqrt(rho) * sd, n)
    e_sd = np.sqrt(1 - rho) * sd
    upper, lower = np.append(lr.CUTS_S, np.inf), np.insert(lr.CUTS_S, 0, -np.inf)
    p = norm.cdf((upper - T[:, None]) / e_sd) - norm.cdf((lower - T[:, None]) / e_sd)   # P(S = s | T)
    es_t = p @ lr.S_VALUES
    s = lr.harmonize(np.digitize(T + rng.normal(0, e_sd, n), lr.CUTS_Y))
    outcome = 0.3 * es_t + rng.normal(0, 0.2, n)
    slope = np.cov(outcome, s)[0, 1] / s.var(ddof=1)
    lam, _ = lr.single_occasion_corr(lr.CUTS_S, lr.S_VALUES, mu, sd, rho)
    assert abs(slope / 0.3 - lam) < 0.01


def test_s3_fit_recovers_planted_parameters():
    rng = np.random.default_rng(SEED + 3)
    truth = dict(mu=1.9, sd=1.1, rho=0.6, delta=0.2)
    y1, y2 = simulate_pairs(rng, 60000, **truth)
    fit = lr.fit_pairs(lr.table_for(y1, y2, "Y5", "S3"))
    mu, sd, rho, delta = lr.unpack(fit["theta"])
    assert fit["converged"] and not fit["at_bound"]
    assert abs(mu - truth["mu"]) < 0.05 and abs(sd - truth["sd"]) < 0.05
    assert abs(rho - truth["rho"]) < 0.04 and abs(delta - truth["delta"]) < 0.04


def test_joint_fit_recovers_shared_var_e():
    rng = np.random.default_rng(SEED + 4)
    var_e = 0.2
    arms = [dict(mu=0.9, var_T=0.5, delta=0.1), dict(mu=1.3, var_T=0.9, delta=-0.1)]
    Ns = []
    for a in arms:
        sd = np.sqrt(a["var_T"] + var_e)
        y1, y2 = simulate_pairs(rng, 40000, a["mu"], sd, a["var_T"] / sd ** 2, a["delta"])
        Ns.append(lr.pair_table(y1, y2))
    fit = lr.fit_joint(Ns)
    ps, ve = lr.joint_arm_params(fit["theta"])
    assert fit["converged"] and abs(ve - var_e) < 0.02
    for (mu, sd, rho, d), a in zip(ps, arms):
        assert abs(mu - a["mu"]) < 0.03 and abs(rho * sd ** 2 - a["var_T"]) < 0.04 and abs(d - a["delta"]) < 0.03


def test_three_level_single_read_fit_equals_closed_form():
    counts = np.array([89, 94, 80])
    f = lr.fit_single(counts, lr.CUTS_S)
    mu, sd = lr.exact_three_level_probit(counts)
    assert abs(f["mu_Z"] - mu) < 1e-4 and abs(f["sd_Z"] - sd) < 1e-4
    assert np.allclose(f["fitted_counts"], counts, atol=1e-3)


def test_latent_shape_bounds_bracket_the_normal_lambda():
    """Feed the exact normal-model marginal: the NPMLE reproduces it and the LP bounds contain the normal lambda."""
    mu, sd, var_e = 1.6, 1.2, 0.3
    p = np.diff(np.concatenate([[0], norm.cdf((lr.CUTS_Y - mu) / sd), [1]]))
    lam_normal, _ = lr.single_occasion_corr(lr.CUTS_S, lr.S_VALUES, mu, sd, 1 - var_e / sd ** 2)
    s = lr.latent_shape_sensitivity(p * 1e6, lr.CUTS_Y, var_e)
    assert s["lp_marginal"] == "observed"
    assert s["lp_min_lambda"] - 1e-6 <= lam_normal <= s["lp_max_lambda"] + 1e-6
    assert s["npmle_max_abs_marginal_error"] < 1e-3
    assert s["lp_min_lambda"] - 1e-6 <= s["npmle_lambda"] <= s["lp_max_lambda"] + 1e-6


def test_lambda_conditional_is_the_slope_attenuation_with_a_covariate():
    """Outcome depends on E[S|T] and on W; W correlates with T; label error is independent of W."""
    rng = np.random.default_rng(SEED + 5)
    mu, var_T, var_e, n = 1.5, 0.9, 0.35, 400000
    W = rng.normal(0, 1, n)
    T = mu + np.sqrt(var_T) * (0.6 * W + 0.8 * rng.normal(0, 1, n))
    e_sd = np.sqrt(var_e)
    upper, lower = np.append(lr.CUTS_S, np.inf), np.insert(lr.CUTS_S, 0, -np.inf)
    es_t = (norm.cdf((upper - T[:, None]) / e_sd) - norm.cdf((lower - T[:, None]) / e_sd)) @ lr.S_VALUES
    s = lr.harmonize(np.digitize(T + rng.normal(0, e_sd, n), lr.CUTS_Y))
    y = 0.3 * es_t + 0.5 * W + rng.normal(0, 0.2, n)
    X = np.column_stack([np.ones(n), s, W])
    slope = np.linalg.lstsq(X, y, rcond=None)[0][1]
    lam = np.corrcoef(s, lr.harmonize(np.digitize(T + rng.normal(0, e_sd, n), lr.CUTS_Y)))[0, 1]
    r2 = np.corrcoef(s, W)[0, 1] ** 2
    assert abs(slope / 0.3 - lr.lambda_conditional(lam, r2)) < 0.015


def test_lambda_falls_as_var_e_rises_so_the_carried_spec_is_the_larger_var_e():
    for mu, sd in ((1.2, 1.3), (1.9, 1.02), (1.3, 1.8)):
        lams = [lr.cohort_lambda(v, mu, sd)[0] for v in np.linspace(0.01, 0.9 * sd ** 2, 40)]
        assert np.all(np.diff(lams) < 0), (mu, sd)


def test_het_model_at_b0_is_the_homoscedastic_model():
    mu, sd_T, sd_e, delta = 1.4, 0.9, 0.5, 0.2
    sd = np.sqrt(sd_T ** 2 + sd_e ** 2)
    P_het = lr.het_cell_probs(lr.CUTS_Y, mu, sd_T, np.log(sd_e), 0.0, delta)
    P_hom = lr.cell_probs(lr.CUTS_Y, mu, mu + delta, sd, sd_T ** 2 / sd ** 2)
    assert np.abs(P_het - P_hom).max() < 1e-7
    lam_het = lr.het_cohort_lambda(mu, sd_T, np.log(sd_e), 0.0)[0]
    lam_hom = lr.cohort_lambda(sd_e ** 2, mu, sd)[0]
    assert abs(lam_het - lam_hom) < 1e-7


def test_het_fit_recovers_a_planted_level_dependent_error():
    rng = np.random.default_rng(SEED + 6)
    mu, sd_T, a, b, delta, n = 1.3, 1.0, np.log(0.4), 0.5, 0.1, 60000
    T = rng.normal(mu, sd_T, n)
    z1 = T + np.exp(a + b * (T - lr.HET_CENTRE)) * rng.normal(0, 1, n)
    z2 = T + delta + np.exp(a + b * (T + delta - lr.HET_CENTRE)) * rng.normal(0, 1, n)
    N = lr.pair_table(np.digitize(z1, lr.CUTS_Y), np.digitize(z2, lr.CUTS_Y))
    f = lr.fit_het(N, [lr.het_from_hom(lr.fit_pairs(N, rho_nonneg=True)["theta"])])
    m, lsd, a_, b_, d_ = f["theta"]
    assert f["converged"] and abs(b_ - b) < 0.08 and abs(a_ - a) < 0.05 and abs(np.exp(lsd) - sd_T) < 0.04


def test_r2_design_drops_all_missing_adds_indicators_and_adjusts():
    rng = np.random.default_rng(SEED + 7)
    n = 300
    Z = rng.normal(0, 1, (n, 3))
    Z[:, 1] = np.nan                          # all missing: dropped
    miss = rng.random(n) < 0.3
    Z[miss, 2] = np.nan                       # partly missing: filled, indicator appended
    S = (Z[:, 0] + rng.normal(0, 1, n) > 0).astype(float) + miss
    X, s = lr.design_matrix(Z, S)
    assert X.shape[1] == 3 and np.allclose(X.mean(0), 0) and abs(s.mean()) < 1e-12
    n_, p, raw, adj = lr.r2_of(Z, S)
    ref = np.column_stack([np.ones(n), Z[:, 0], np.where(miss, np.nanmean(Z[:, 2]), Z[:, 2]), miss])
    res = S - ref @ np.linalg.lstsq(ref, S, rcond=None)[0]
    raw_ref = 1 - res @ res / ((S - S.mean()) @ (S - S.mean()))
    assert p == 3 and abs(raw - raw_ref) < 1e-12
    assert abs(adj - (1 - (1 - raw_ref) * (n - 1) / (n - 4))) < 1e-12
