"""Synthetic checks for REC's independent v0 likelihood draft."""

import unittest

import numpy as np
from scipy.special import expit
from scipy.stats import betabinom

from likelihood_v0 import (
    ProportionalModel, Row, mixture_logpmf, row_loglik_score,
    tau_marginal_loglik, truncation_logmass,
)


def make_rows():
    rows = []
    for gene in range(2):
        for individual in range(5):
            n = 35 + 7 * individual + 4 * gene
            rows.append(
                Row(
                    gene=gene,
                    n=n,
                    a=n // 2 + (-1) ** individual * (2 + gene),
                    q=0.24 + 0.12 * gene,
                    e=0.008 + individual * 0.001,
                    phi=0.01,
                    omega=0.03,
                    orientation=(-1) ** gene,
                    stage=float(individual % 3 - 1),
                    z=(float(individual - 2) / 2,),
                    background=float(individual - 2) / 3,
                    duplication=float(gene - 0.5),
                    min_depth=10,
                    min_minor=3,
                    min_fraction=0.08,
                )
            )
    return rows


class LikelihoodChecks(unittest.TestCase):
    def test_orientation_swap_and_boundary_genotype_prior(self):
        row = make_rows()[0]
        for q in (0.0, 0.37, 1.0):
            forward = Row(**{**row.__dict__, "q": q, "orientation": 1})
            reverse = Row(**{**row.__dict__, "q": q, "orientation": -1})
            # Reversing the lead haplotype reverses the oriented count and eta;
            # tag ALT frequency q still describes the same physical tag ALT.
            a = mixture_logpmf(forward, np.array([forward.a]), 0.23, -3.6)[0]
            b = mixture_logpmf(reverse, np.array([reverse.n - reverse.a]), -0.23, -3.6)[0]
            self.assertTrue(np.isfinite(a))
            self.assertAlmostEqual(a, b, places=10)

    def test_explicit_normalizer_matches_all_count_enumeration(self):
        row = make_rows()[0]
        all_counts = np.arange(row.n + 1)
        logmass = mixture_logpmf(row, all_counts, 0.27, -3.6)
        self.assertAlmostEqual(np.exp(logmass).sum(), 1, places=10)
        allowed_mass = np.exp(logmass[row.allowed]).sum()
        self.assertAlmostEqual(np.exp(truncation_logmass(row, 0.27, -3.6)), allowed_mass, places=10)

    def test_thousands_of_reads_and_rare_tag_frequency(self):
        base = make_rows()[0]
        row = Row(
            gene=base.gene,
            n=5000,
            a=2450,
            q=0.0001,
            e=0.00001,
            phi=base.phi,
            omega=base.omega,
            orientation=base.orientation,
            stage=base.stage,
            z=base.z,
            background=base.background,
            duplication=base.duplication,
            min_depth=base.min_depth,
            min_minor=base.min_minor,
            min_fraction=base.min_fraction,
        )
        self.assertTrue(np.all(np.isfinite(row_loglik_score(row, 0.2, -4.0))))

    def test_gradient_matches_finite_difference(self):
        model = ProportionalModel(make_rows())
        p = model.initial()
        p[:2] = [0.35, -0.42]
        p[2:4] = [-3.4, -3.7]
        p[4:] = [0.11, -0.08, 0.03, -0.02, 0.01]
        _, analytic = model.loglik_gradient(p)
        numeric = np.empty_like(p)
        step = 1e-5
        for index in range(len(p)):
            plus = p.copy()
            minus = p.copy()
            plus[index] += step
            minus[index] -= step
            numeric[index] = (
                model.loglik_gradient(plus)[0] - model.loglik_gradient(minus)[0]
            ) / (2 * step)
        np.testing.assert_allclose(analytic, numeric, atol=1e-6, rtol=2e-5)

    def test_tau_zero_reduces_to_proportional_model(self):
        model = ProportionalModel(make_rows())
        p = model.initial()
        p[0] = 0.2
        p[1] = -0.1
        self.assertEqual(tau_marginal_loglik(model, p, 0), model.loglik_gradient(p)[0])

    def test_truncated_mixture_matches_brute_force(self):
        row = make_rows()[0]
        eta = 0.27
        dispersion_linear = -3.6
        ll = row_loglik_score(row, eta, dispersion_linear)[0]
        rho = expit(dispersion_linear)
        mean_h = expit(eta - row.orientation * row.omega)
        t = (1 - rho) / rho
        h = betabinom.pmf(np.arange(row.n + 1), row.n, mean_h * t, (1 - mean_h) * t)
        err_t = (1 - row.phi) / row.phi
        ref_mean = row.e if row.orientation == 1 else 1 - row.e
        r = betabinom.pmf(
            np.arange(row.n + 1), row.n, ref_mean * err_t, (1 - ref_mean) * err_t
        )
        a = betabinom.pmf(
            np.arange(row.n + 1), row.n, (1 - ref_mean) * err_t, ref_mean * err_t
        )
        mix = 2 * row.q * (1 - row.q) * h + (1 - row.q) ** 2 * r + row.q**2 * a
        self.assertAlmostEqual(ll, np.log(mix[row.a] / mix[row.allowed].sum()), places=10)


if __name__ == "__main__":
    unittest.main()
