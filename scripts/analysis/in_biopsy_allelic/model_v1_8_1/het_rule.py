"""Integer het predicate and truncation set K_c(n) for Model A-final (spec v1 section 1.3).

`truncation_bounds` is the one function that defines K_c(n). The row builder, the
normalizers and the simulators all call it, so the truncation set is the same everywhere.
`het` is the definitional predicate; a unit test checks that the two agree for n = 0..2,000.
"""
import numpy as np


def het(r, x, d, k, F):
    """het(r, x) under rule (d, k, F): n >= d, min(r, x) >= k and 100*min >= F*n, in integers."""
    r = np.asarray(r, dtype=np.int64)
    x = np.asarray(x, dtype=np.int64)
    n = r + x
    m = np.minimum(r, x)
    return (n >= d) & (m >= k) & (100 * m >= F * n)


def truncation_bounds(n, d, k, F):
    """Return (lo, hi, nonempty) with K_c(n) = {a : lo <= a <= hi} when nonempty.

    lo = max(k, ceil(F*n/100)) by integer ceil division, hi = n - lo. K is empty when n < d
    or lo > hi. Arguments broadcast.
    """
    n = np.asarray(n, dtype=np.int64)
    d = np.asarray(d, dtype=np.int64)
    k = np.asarray(k, dtype=np.int64)
    F = np.asarray(F, dtype=np.int64)
    lo = np.maximum(k, -((-F * n) // 100))
    hi = n - lo
    nonempty = (n >= d) & (lo <= hi)
    return lo, hi, nonempty
