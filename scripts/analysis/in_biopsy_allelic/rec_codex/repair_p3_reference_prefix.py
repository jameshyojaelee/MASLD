"""Value-only, full-support BB reference using exact rising-factorial prefixes.

This independent module does not import the candidate likelihood or its BB code.
The probability identity is exact; floating cumulative sums need equivalence checks.
Numerical imports/execution are restricted to an individual compute allocation.
"""
from functools import lru_cache
import os

if not os.environ.get('SLURM_JOB_ID', '').isdigit() or os.environ.get('SLURM_ARRAY_TASK_ID'):
    raise RuntimeError('Individual compute allocation required for numerical reference')

import numpy as np
from scipy.special import expit, gammaln

EPS = 1e-6


@lru_cache(maxsize=16)
def support_constants(n):
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or not 0 <= n <= 100000:
        raise ValueError('Invalid bounded integer support')
    n = int(n)
    k = np.arange(n + 1, dtype=float)
    choose = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    logj = np.log(np.arange(1, n, dtype=float))
    choose.flags.writeable = False
    logj.flags.writeable = False
    return choose, logj


def rising_prefix(logshape, n, logj):
    """Entry k is log((shape)_k), including the exact zero-count value."""
    if not np.isfinite(logshape):
        raise ValueError('Nonfinite log shape')
    prefix = np.empty(n + 1, dtype=float)
    prefix[0] = 0.
    if n:
        terms = np.empty(n, dtype=float)
        terms[0] = logshape
        terms[1:] = np.logaddexp(logshape, logj)
        np.cumsum(terms, out=prefix[1:])
    return prefix


class PreparedBB:
    """Prepare only n and dispersion, which do not vary with the random effect.

    Callers must create a new object whenever the linear dispersion changes,
    including finite-difference perturbations. No parameter-dependent cache.
    """
    def __init__(self, n, linear_r):
        self.choose, self.logj = support_constants(n)
        self.n = int(n)
        if not np.isfinite(linear_r):
            raise ValueError('Nonfinite linear dispersion')
        self.linear_r = float(linear_r)
        rho = EPS + (1 - 2 * EPS) * expit(self.linear_r)
        self.logt = float(np.log((1 - rho) / rho))
        self.denominator = float(rising_prefix(self.logt, self.n, self.logj)[self.n])

    def logpmf(self, eta):
        if not np.isfinite(eta):
            raise ValueError('Nonfinite predictor')
        loga = self.logt - np.logaddexp(0., -eta)
        logb = self.logt - np.logaddexp(0., eta)
        a = rising_prefix(loga, self.n, self.logj)
        b = rising_prefix(logb, self.n, self.logj)
        return self.choose + a + b[::-1] - self.denominator


def logpmf_full(n, eta, linear_r):
    return PreparedBB(n, linear_r).logpmf(eta)
