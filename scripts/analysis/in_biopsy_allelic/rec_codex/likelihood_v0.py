"""Independent, outcome-blind draft of Model A v0's conditional likelihood.

This is REC code. It uses only the written v0 specification and synthetic inputs.
It deliberately does not read the builder implementation or study observations.
The unresolved v1 choices are recorded in the Codex specification critique.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import minimize
from scipy.special import betaln, digamma, expit, gammaln, logsumexp
from scipy.stats import binom


@dataclass(frozen=True)
class Row:
    gene: int
    n: int
    a: int
    q: float
    e: float
    phi: float
    omega: float
    orientation: int
    stage: float
    z: tuple[float, ...]
    background: float
    duplication: float
    min_depth: int
    min_minor: int
    min_fraction: float

    def __post_init__(self) -> None:
        if self.n < 0 or not 0 <= self.a <= self.n:
            raise ValueError("invalid read counts")
        if not 0 <= self.q <= 1 or not 0 <= self.e <= 1:
            raise ValueError("q and e must be probabilities")
        if not 0 <= self.phi < 1 or self.orientation not in (-1, 1):
            raise ValueError("invalid phi or orientation")
        if not 0 <= self.min_fraction <= 0.5:
            raise ValueError("invalid heterozygote fraction threshold")
        if self.n < self.min_depth or self.a < self.lower or self.a > self.n - self.lower:
            raise ValueError("row does not pass the specified het rule")

    @property
    def lower(self) -> int:
        return max(self.min_minor, int(np.ceil(self.min_fraction * self.n)))

    @property
    def allowed(self) -> np.ndarray:
        return np.arange(self.lower, self.n - self.lower + 1, dtype=np.int64)


def _bb_logpmf_score(k: np.ndarray, n: int, mean: float, rho: float):
    """Full BB log PMF and scores in logit(mean), logit(rho) coordinates."""
    k = np.asarray(k, dtype=np.float64)
    if not 0 <= mean <= 1 or not 0 <= rho < 1:
        raise ValueError("invalid beta-binomial parameter")
    if rho == 0 or mean in (0, 1):
        lp = binom.logpmf(k, n, mean)
        if mean in (0, 1):
            return lp, np.zeros_like(k), np.zeros_like(k)
        return lp, k - n * mean, np.zeros_like(k)
    total = (1 - rho) / rho
    left = mean * total
    right = (1 - mean) * total
    lp = (
        gammaln(n + 1)
        - gammaln(k + 1)
        - gammaln(n - k + 1)
        + betaln(k + left, n - k + right)
        - betaln(left, right)
    )
    left_delta = digamma(k + left) - digamma(left)
    right_delta = digamma(n - k + right) - digamma(right)
    d_logit_mean = mean * (1 - mean) * total * (left_delta - right_delta)
    d_logit_rho = -total * (
        mean * left_delta
        + (1 - mean) * right_delta
        - digamma(n + total)
        + digamma(total)
    )
    return lp, d_logit_mean, d_logit_rho


def mixture_logpmf(row: Row, counts: np.ndarray, eta: float, dispersion_linear: float):
    """Untruncated three-class log PMF, including BB combinatorial constants."""
    counts = np.asarray(counts, dtype=np.int64)
    if np.any((counts < 0) | (counts > row.n)):
        raise ValueError("counts outside 0..n")
    mean_het = expit(eta - row.orientation * row.omega)
    rho = expit(dispersion_linear)
    mean_ref_hom = row.e if row.orientation == 1 else 1 - row.e
    component = np.stack(
        (
            _bb_logpmf_score(counts, row.n, mean_het, rho)[0],
            _bb_logpmf_score(counts, row.n, mean_ref_hom, row.phi)[0],
            _bb_logpmf_score(counts, row.n, 1 - mean_ref_hom, row.phi)[0],
        ), axis=0,
    )
    with np.errstate(divide="ignore"):
        log_prior = np.log([2 * row.q * (1 - row.q), (1 - row.q) ** 2, row.q**2])
    return logsumexp(log_prior[:, None] + component, axis=0)


def truncation_logmass(row: Row, eta: float, dispersion_linear: float):
    """Log mixture mass of the inclusive integer het-call interval."""
    return float(logsumexp(mixture_logpmf(row, row.allowed, eta, dispersion_linear)))


def row_loglik_score(row: Row, eta: float, dispersion_linear: float):
    """Log P(a | n, called het) and derivatives in eta and logit(rho)."""
    allowed = row.allowed
    if allowed.size == 0:
        raise ValueError("heterozygote rule allows no counts")
    s = -row.orientation  # reference-bias sign in the oriented allele ratio
    mean_het = expit(eta + row.omega * s)
    rho = expit(dispersion_linear)
    chosen = np.array([row.a], dtype=np.int64)
    het_a, d_mean_a, d_rho_a = _bb_logpmf_score(chosen, row.n, mean_het, rho)
    het_k, d_mean_k, d_rho_k = _bb_logpmf_score(allowed, row.n, mean_het, rho)
    mean_ref_hom = row.e if row.orientation == 1 else 1 - row.e
    mean_alt_hom = 1 - mean_ref_hom
    ref_a = _bb_logpmf_score(chosen, row.n, mean_ref_hom, row.phi)[0]
    alt_a = _bb_logpmf_score(chosen, row.n, mean_alt_hom, row.phi)[0]
    ref_k = _bb_logpmf_score(allowed, row.n, mean_ref_hom, row.phi)[0]
    alt_k = _bb_logpmf_score(allowed, row.n, mean_alt_hom, row.phi)[0]
    with np.errstate(divide="ignore"):
        log_prior = np.log(np.array([2 * row.q * (1 - row.q), (1 - row.q) ** 2, row.q**2]))

    log_num_terms = log_prior + np.array([het_a[0], ref_a[0], alt_a[0]])
    log_num = logsumexp(log_num_terms)
    num_het_weight = np.exp(log_num_terms[0] - log_num)
    log_het_allowed = log_prior[0] + het_k
    log_den = logsumexp(
        np.concatenate((log_het_allowed, log_prior[1] + ref_k, log_prior[2] + alt_k))
    )
    den_het_weights = np.exp(log_het_allowed - log_den)
    d_eta = num_het_weight * d_mean_a[0] - np.dot(den_het_weights, d_mean_k)
    d_dispersion = num_het_weight * d_rho_a[0] - np.dot(den_het_weights, d_rho_k)
    return float(log_num - log_den), float(d_eta), float(d_dispersion)


class ProportionalModel:
    """Section 3 joint likelihood with an analytic gradient.

    Vector order: alpha[g], r[g], kappa, delta[z], rho_stage,
    rho_background, rho_duplication. Gene IDs must be 0..G-1.
    """

    def __init__(self, rows: list[Row]):
        if not rows:
            raise ValueError("at least one row is required")
        self.rows = tuple(rows)
        self.genes = max(row.gene for row in rows) + 1
        self.z_dim = len(rows[0].z)
        if any(row.gene < 0 or len(row.z) != self.z_dim for row in rows):
            raise ValueError("invalid gene ID or nuisance-vector dimension")
        if set(row.gene for row in rows) != set(range(self.genes)):
            raise ValueError("gene IDs must be dense")
        self.size = 2 * self.genes + self.z_dim + 4

    def initial(self) -> np.ndarray:
        p = np.zeros(self.size)
        p[self.genes : 2 * self.genes] = np.log(0.02 / 0.98)
        return p

    def loglik_gradient(self, parameters: np.ndarray):
        p = np.asarray(parameters, dtype=float)
        if p.shape != (self.size,):
            raise ValueError("wrong parameter-vector length")
        g_count, z_dim = self.genes, self.z_dim
        kappa = p[2 * g_count]
        delta = p[2 * g_count + 1 : 2 * g_count + 1 + z_dim]
        rho_stage, rho_background, rho_duplication = p[-3:]
        gradient = np.zeros_like(p)
        loglik = 0.0
        for row in self.rows:
            alpha = p[row.gene]
            factor = 1 + kappa * row.stage + np.dot(delta, row.z)
            eta = alpha * factor
            disp_linear = (
                p[g_count + row.gene]
                + rho_stage * row.stage
                + rho_background * row.background
                + rho_duplication * row.duplication
            )
            ll, d_eta, d_disp = row_loglik_score(row, eta, disp_linear)
            loglik += ll
            gradient[row.gene] += d_eta * factor
            gradient[g_count + row.gene] += d_disp
            gradient[2 * g_count] += d_eta * alpha * row.stage
            gradient[2 * g_count + 1 : 2 * g_count + 1 + z_dim] += d_eta * alpha * np.asarray(row.z)
            gradient[-3] += d_disp * row.stage
            gradient[-2] += d_disp * row.background
            gradient[-1] += d_disp * row.duplication
        return float(loglik), gradient

    def fit(self, initial: np.ndarray | None = None):
        start = self.initial() if initial is None else np.asarray(initial, dtype=float)

        def objective(p):
            ll, grad = self.loglik_gradient(p)
            return -ll, -grad

        return minimize(objective, start, method="L-BFGS-B", jac=True)

    def with_stages(self, stages: np.ndarray):
        if len(stages) != len(self.rows):
            raise ValueError("one stage value per row is required")
        return ProportionalModel(
            [replace(row, stage=float(stage)) for row, stage in zip(self.rows, stages)]
        )


def tau_marginal_loglik(model: ProportionalModel, parameters: np.ndarray, tau: float, nodes=20):
    """Gene-level marginal log likelihood under section 4; no τ fit yet."""
    if tau < 0:
        raise ValueError("tau must be nonnegative")
    if tau == 0:
        return model.loglik_gradient(parameters)[0]
    roots, weights = np.polynomial.hermite.hermgauss(nodes)
    p = np.asarray(parameters, dtype=float)
    g_count, z_dim = model.genes, model.z_dim
    kappa = p[2 * g_count]
    delta = p[2 * g_count + 1 : 2 * g_count + 1 + z_dim]
    rho_stage, rho_background, rho_duplication = p[-3:]
    total = 0.0
    for gene in range(g_count):
        gene_rows = [row for row in model.rows if row.gene == gene]
        node_loglik = np.zeros(nodes)
        for index, root in enumerate(roots):
            u = np.sqrt(2) * tau * root
            for row in gene_rows:
                alpha = p[gene]
                eta = alpha * (1 + kappa * row.stage + np.dot(delta, row.z)) + u * row.stage
                disp_linear = (
                    p[g_count + gene]
                    + rho_stage * row.stage
                    + rho_background * row.background
                    + rho_duplication * row.duplication
                )
                node_loglik[index] += row_loglik_score(row, eta, disp_linear)[0]
        total += logsumexp(np.log(weights) + node_loglik) - 0.5 * np.log(np.pi)
    return float(total)
