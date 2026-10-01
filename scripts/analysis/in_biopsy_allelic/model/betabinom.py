"""Beta-binomial pieces for Model A-final (spec v1 sections 2.2 and 7).

BB(n, mu, rho) has shapes alpha = mu(1-rho)/rho and beta = (1-mu)(1-rho)/rho and
log pmf = lchoose(n, a) + lbeta(a + alpha, n - a + beta) - lbeta(alpha, beta), full constants.
The binomial is used only when rho is exactly 0 (the fixed homozygote dispersion phi = 0).
"""
import numpy as np
from scipy.special import digamma, gammaln, xlog1py, xlogy, zeta

# From this shape value up, the differences below use the asymptotic (Stirling) series, so that
# a shape near 1e6 (rho at its 1e-6 floor) does not lose ~1e-8 of absolute precision per row in
# a difference of two large lgamma values. At x = 30 the first omitted series term is < 1e-23.
SERIES_X = 30.0


def _split(x, m):
    x, m = np.broadcast_arrays(np.asarray(x, float), np.asarray(m, float))
    return x, m, x >= SERIES_X


def _lgamma_series(y):
    y2 = y * y
    return (1 / 12 - (1 / 360 - (1 / 1260 - (1 / 1680 - (1 / 1188 - (691 / 360360
            - 1 / (156 * y2)) / y2) / y2) / y2) / y2) / y2) / y


def _digamma_series(y):
    """psi(y) - log(y) + 1/(2y)."""
    y2 = y * y
    return -(1 / 12 - (1 / 120 - (1 / 252 - (1 / 240 - (1 / 132 - (691 / 32760
             - 1 / (12 * y2)) / y2) / y2) / y2) / y2) / y2) / y2


def _trigamma_series(y):
    """psi_1(y) - 1/y - 1/(2y^2) - 1/(6y^3)."""
    y2 = y * y
    return -(1 / 30 - (1 / 42 - (1 / 30 - (5 / 66 - (691 / 2730 - 7 / (6 * y2)) / y2) / y2)
             / y2) / y2) / (y2 * y2 * y)


def lgamma_diff(x, m):
    """log Gamma(x + m) - log Gamma(x) for x > 0, m >= 0."""
    x, m, big = _split(x, m)
    out = np.empty(x.shape)
    out[~big] = gammaln(x[~big] + m[~big]) - gammaln(x[~big])
    xb, mb = x[big], m[big]
    z = xb + mb
    out[big] = ((xb - 0.5) * np.log1p(mb / xb) + mb * np.log(z) - mb
                + (_lgamma_series(z) - _lgamma_series(xb)))
    return out


def digamma_diff(x, m):
    """psi(x + m) - psi(x) for x > 0, m >= 0."""
    x, m, big = _split(x, m)
    out = np.empty(x.shape)
    out[~big] = digamma(x[~big] + m[~big]) - digamma(x[~big])
    xb, mb = x[big], m[big]
    z = xb + mb
    out[big] = (np.log1p(mb / xb) + mb / (2 * xb * z)
                + (_digamma_series(z) - _digamma_series(xb)))
    return out


def trigamma_diff(x, m):
    """psi_1(x + m) - psi_1(x) for x > 0, m >= 0."""
    x, m, big = _split(x, m)
    out = np.empty(x.shape)
    out[~big] = zeta(2.0, x[~big] + m[~big]) - zeta(2.0, x[~big])
    xb, mb = x[big], m[big]
    z = xb + mb
    out[big] = (-mb / (xb * z) - mb * (xb + z) / (2 * xb ** 2 * z ** 2)
                - mb * (xb ** 2 + xb * z + z ** 2) / (6 * xb ** 3 * z ** 3)
                + (_trigamma_series(z) - _trigamma_series(xb)))
    return out


def lchoose(n, a):
    return gammaln(n + 1.0) - gammaln(a + 1.0) - gammaln(n - a + 1.0)


def bb_logpmf(a, n, mu, rho):
    """log BB(a | n, mu, rho); binomial where rho == 0. Arguments broadcast.

    lbeta(a + alpha, n - a + beta) - lbeta(alpha, beta) is evaluated as
    [lgamma(a+alpha) - lgamma(alpha)] + [lgamma(n-a+beta) - lgamma(beta)] - [lgamma(n+t) - lgamma(t)].
    """
    a, n, mu, rho = np.broadcast_arrays(np.asarray(a, float), np.asarray(n, float),
                                        np.asarray(mu, float), np.asarray(rho, float))
    out = np.empty(a.shape)
    binom = rho == 0
    if binom.any():
        ab, nb, mb = a[binom], n[binom], mu[binom]
        out[binom] = lchoose(nb, ab) + xlogy(ab, mb) + xlog1py(nb - ab, -mb)
    bb = ~binom
    if bb.any():
        ab, nb, mb, rb = a[bb], n[bb], mu[bb], rho[bb]
        t = (1.0 - rb) / rb
        al, be = mb * t, (1.0 - mb) * t
        out[bb] = (lchoose(nb, ab) + lgamma_diff(al, ab) + lgamma_diff(be, nb - ab)
                   - lgamma_diff(t, nb))
    return out


def bb_logpmf_recurrence(n, mu, rho, a0, a1):
    """log BB(a) for a = a0..a1 by the section 7 recurrence from a direct value at a0.

    log f(a+1) - log f(a) = log(n-a) + log(a+alpha) - log(a+1) - log(n-a-1+beta).
    Reference implementation for the unit tests (scalar n, mu, rho; rho > 0).
    """
    t = (1.0 - rho) / rho
    al, be = mu * t, (1.0 - mu) * t
    a = np.arange(a0, a1, dtype=float)
    inc = np.log(n - a) + np.log(a + al) - np.log(a + 1.0) - np.log(n - a - 1.0 + be)
    first = bb_logpmf(a0, n, mu, rho)
    return np.concatenate([np.atleast_1d(first), first + np.cumsum(inc)])


def seg_ranges(start, length):
    """Flattened integer ranges: for segment j, start[j], ..., start[j] + length[j] - 1.

    Returns (segment id per element, value per element).
    """
    start = np.asarray(start, dtype=np.int64)
    length = np.asarray(length, dtype=np.int64)
    seg = np.repeat(np.arange(len(start)), length)
    offs = np.arange(length.sum()) - np.repeat(np.cumsum(length) - length, length)
    return seg, np.repeat(start, length) + offs


def seg_logsumexp(x, seg, nseg):
    """Log-sum-exp of x within segments (seg ids need not be sorted). Empty segment -> -inf."""
    m = np.full(nseg, -np.inf)
    np.maximum.at(m, seg, x)
    mf = np.where(np.isfinite(m), m, 0.0)
    s = np.bincount(seg, weights=np.exp(x - mf[seg]), minlength=nseg)
    with np.errstate(divide="ignore"):
        return np.log(s) + mf


def log_mass_over_K(n, lo, hi, mu, rho):
    """log of sum_{a=lo..hi} BB(a | n, mu, rho) per row, by direct log-sum-exp over K.

    Used for the fixed homozygote normalizers C_R and C_A (inputs e_i and phi only).
    Rows with lo > hi get -inf.
    """
    n = np.asarray(n, dtype=np.int64)
    lo = np.asarray(lo, dtype=np.int64)
    hi = np.asarray(hi, dtype=np.int64)
    length = np.maximum(hi - lo + 1, 0)
    seg, a = seg_ranges(lo, length)
    mu = np.broadcast_to(np.asarray(mu, float), n.shape)
    rho = np.broadcast_to(np.asarray(rho, float), n.shape)
    lf = bb_logpmf(a, n[seg], mu[seg], rho[seg])
    return seg_logsumexp(lf, seg, len(n))
