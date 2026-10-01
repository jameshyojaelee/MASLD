"""Chosen donor-loss weights and penalized reduced-rank RNA regression.

Weights express a predictive target preference, not known error precision.
No likelihood-based standard errors follow from this transformation.
"""
import numpy as np


def donor_weights(strata, strength):
    """Equal graph degree among donors with an exact metadata partner."""
    if not np.isfinite(strength) or strength < 0:
        raise ValueError("Contrast strength must be finite and nonnegative")
    n = len(strata)
    if n < 2:
        raise ValueError("At least two training donors required")
    groups = {}
    for i, key in enumerate(strata):
        groups.setdefault(tuple(key), []).append(i)
    laplacian = np.zeros((n, n))
    paired = sum(len(indices) for indices in groups.values() if len(indices) >= 2)
    for indices in groups.values():
        m = len(indices)
        if m >= 2:
            laplacian[np.ix_(indices, indices)] = (m * np.eye(m) - np.ones((m, m))) / (m - 1)
    if paired == 0 and strength > 0:
        raise ValueError("Positive contrast strength has no comparable donors")
    scale = n / paired if paired else 0.0
    return np.eye(n) + strength * scale * laplacian


def fit_predict(x_train, y_train, x_query, weight, penalty_fraction, rank=None):
    """Minimize weighted SSE plus isotropic coefficient penalty at fixed rank.

    RNA means/scales and output intercept use training donors. Kernel scaling
    divides by feature count; the penalty is relative to its largest eigenvalue.
    Rank=None is the exact ridge comparator. Truncation uses the penalized
    objective's right singular vectors, rather than unpenalized fitted SSE.
    """
    x = np.asarray(x_train, float)
    y = np.asarray(y_train, float)
    query = np.asarray(x_query, float)
    w = np.asarray(weight, float)
    if x.ndim != 2 or y.ndim != 2 or query.ndim != 2:
        raise ValueError("Matrices required")
    n, p = x.shape
    if n < 2 or p < 1 or y.shape[0] != n or query.shape[1] != p or w.shape != (n, n):
        raise ValueError("Incompatible donor/feature axes")
    if not all(np.isfinite(a).all() for a in (x, y, query, w)):
        raise ValueError("Nonfinite inputs")
    if not np.allclose(w, w.T, atol=1e-12) or not np.allclose(w @ np.ones(n), np.ones(n), atol=1e-10):
        raise ValueError("Weight must be symmetric and preserve the intercept")
    if not np.isfinite(penalty_fraction) or penalty_fraction <= 0:
        raise ValueError("Strictly positive ridge penalty required")
    if rank is not None and (not isinstance(rank, int) or rank < 1):
        raise ValueError("Positive integer rank required")
    mean = x.mean(0)
    scale = x.std(0)
    scale[scale < 1e-8] = 1.0
    a, q = (x - mean) / scale, (query - mean) / scale
    intercept = y.mean(0)
    # Cholesky returns C with W=C C'; choose A=C' so A'A=W.
    transform = np.linalg.cholesky(w).T
    z = transform @ a
    target = transform @ (y - intercept)
    kernel = z @ z.T / p
    eigenvalues, u = np.linalg.eigh(kernel)
    if eigenvalues.min() < -1e-8 * max(1.0, eigenvalues.max()):
        raise ValueError("Negative kernel curvature")
    eigenvalues = np.maximum(eigenvalues, 0.0)
    penalty = penalty_fraction * max(float(eigenvalues.max()), 1e-12)
    projected = u.T @ target
    alpha = u @ (projected / (eigenvalues + penalty)[:, None])
    prediction = q @ z.T @ alpha / p
    if rank is not None:
        penalized_response = np.sqrt(eigenvalues / (eigenvalues + penalty))[:, None] * projected
        _, singular, vt = np.linalg.svd(penalized_response, full_matrices=False)
        keep = min(rank, int(np.sum(singular > 1e-10)))
        directions = vt[:keep].T
        prediction = (prediction @ directions) @ directions.T
    return prediction + intercept
