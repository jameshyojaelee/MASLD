"""Fits for Model A-final (spec v1 section 7): starting values, L-BFGS-B with the section 3.4
bounds followed by a projected Newton polish, the convergence rule and restarts, and the
individual-cluster sandwich SE for kappa (section 5.1)."""
import time
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize
from scipy.special import logit

import model_a as ma

TOL = 1e-6                  # projected max |gradient| (section 7)
DECREMENT = 1e-12           # Newton decrement g'H^-1 g at which the polish stops
N_RESTARTS = 3
RESTART_SEED = [20260923, 90]
START_ALPHAS = (-1.0, 0.0, 1.0)
START_R = logit(0.02)


@dataclass
class Fit:
    theta: np.ndarray
    loglik: float
    grad: np.ndarray
    proj_max_grad: float
    converged: bool
    attempts: int
    seconds: float
    at_bound: list = field(default_factory=list)
    lbfgs_message: str = ""


def projected_gradient(g, x, lb, ub, free):
    pg = np.where(free, g, 0.0)
    pg[(x <= lb) & (pg < 0)] = 0.0
    pg[(x >= ub) & (pg > 0)] = 0.0
    return pg


def start_values(D, include_omega_c=True):
    """Section 7 start: per gene, the truncated mixture with kappa = delta = omega_S = 0 and
    dispersion slopes 0, r = logit(0.02), best of alpha in {-1, 0, 1} (first on ties)."""
    L = D.L
    th = np.zeros(L.size)
    th[L.r] = START_R
    if not include_omega_c:
        D = _without_omega(D)
    ll = np.stack([ma.gene_loglik3(_with_alpha(th, L, a), D) for a in START_ALPHAS])
    th[L.alpha] = np.asarray(START_ALPHAS)[np.argmax(ll, axis=0)]
    th[L.alpha[~D.gene_has_rows]] = 0.0
    return th


def _with_alpha(th, L, a):
    t = th.copy()
    t[L.alpha] = a
    return t


def _without_omega(D):
    import copy
    D2 = copy.copy(D)
    D2.omega = np.zeros_like(D.omega)
    return D2


def _newton_direction(H, g):
    """Ascent direction for a maximum: solve (-H) d = g with a Jacobi-scaled eigen solve;
    eigenvalues are replaced by their absolute value, floored at 1e-10 of the largest."""
    A = -H
    s = 1.0 / np.sqrt(np.maximum(np.abs(np.diag(A)), 1e-300))
    As = A * s[:, None] * s[None, :]
    lam, V = np.linalg.eigh(0.5 * (As + As.T))
    lam = np.maximum(np.abs(lam), 1e-10 * np.max(np.abs(lam)))
    return s * (V @ ((V.T @ (s * g)) / lam))


def maximize(fun, x0, lb, ub, free, scale, tol=TOL, newton_iter=60):
    """Maximize fun over the free coordinates: L-BFGS-B (scaled coordinates) then a projected
    Newton polish. fun(x, order) returns ll, (ll, g) or (ll, g, H) in the full coordinates."""
    idx = np.flatnonzero(free)
    x = np.clip(np.array(x0, float), lb, ub)

    def neg(y):
        xx = x.copy()
        xx[idx] = y * scale[idx]
        try:
            f, g = fun(xx, 1)
        except (RuntimeError, FloatingPointError, np.linalg.LinAlgError):
            return np.inf, np.zeros(len(idx))
        if not np.isfinite(f) or not np.all(np.isfinite(g[idx])):
            return np.inf, np.zeros(len(idx))
        return -f, -g[idx] * scale[idx]

    res = minimize(neg, x[idx] / scale[idx], jac=True, method="L-BFGS-B",
                   bounds=list(zip(lb[idx] / scale[idx], ub[idx] / scale[idx])),
                   options=dict(maxiter=5000, maxfun=10000, maxcor=30, ftol=1e-15, gtol=1e-10))
    x[idx] = np.clip(res.x * scale[idx], lb[idx], ub[idx])
    msg = str(res.message)
    f, g, H = fun(x, 2)
    for _ in range(newton_iter):
        pg = projected_gradient(g, x, lb, ub, free)
        fr = free & ~(((x <= lb) & (g < 0)) | ((x >= ub) & (g > 0)))
        if not np.all(np.isfinite(H[np.ix_(fr, fr)])):
            break
        d = _newton_direction(H[np.ix_(fr, fr)], g[fr])
        # polish past the convergence rule until the Newton decrement is negligible
        if np.max(np.abs(pg)) < tol and g[fr] @ d < DECREMENT:
            break
        t, accepted = 1.0, False
        while t > 1e-12:
            xn = x.copy()
            xn[fr] = np.clip(x[fr] + t * d, lb[fr], ub[fr])
            try:
                fn = fun(xn, 0)
            except (RuntimeError, FloatingPointError):
                fn = -np.inf
            if np.isfinite(fn) and fn >= f - 1e-13 * max(1.0, abs(f)):
                accepted = True
                break
            t *= 0.5
        if not accepted:
            break
        x = xn
        f, g, H = fun(x, 2)
    pg = projected_gradient(g, x, lb, ub, free)
    return x, f, g, float(np.max(np.abs(pg))), msg


def _jitter(x, lb, ub, free, scale, attempt):
    rng = np.random.default_rng(np.random.SeedSequence(RESTART_SEED + [attempt]))
    y = x + np.where(free, rng.normal(0.0, 0.1, len(x)) * scale, 0.0)
    return np.clip(y, lb, ub)


def hessian_scale(fun, x0, free, base):
    """L-BFGS-B coordinate scale 1/sqrt(|H_ii|) at the start (the z columns differ in scale by
    ~1e3, so unscaled L-BFGS-B needs many iterations); base where H_ii is unusable."""
    try:
        d = np.abs(np.diag(fun(x0, 2)[2]))
    except (RuntimeError, FloatingPointError, np.linalg.LinAlgError):
        return base
    ok = free & np.isfinite(d) & (d > 1e-12)
    return np.where(ok, 1.0 / np.sqrt(np.where(ok, d, 1.0)), base)


def _run(fun, start, lb, ub, free, scale, tol):
    t0 = time.perf_counter()
    scale = hessian_scale(fun, np.clip(start, lb, ub), free, scale)
    x0 = start
    best = None
    for attempt in range(N_RESTARTS + 1):
        if attempt > 0:
            x0 = _jitter(best[0] if best is not None else start, lb, ub, free, scale, attempt)
        try:
            x, f, g, pgmax, msg = maximize(fun, x0, lb, ub, free, scale, tol)
        except (RuntimeError, FloatingPointError, np.linalg.LinAlgError) as err:
            x, f, g, pgmax, msg = x0, -np.inf, np.full(len(x0), np.nan), np.inf, repr(err)
        if best is None or f > best[1]:
            best = (x, f, g, pgmax, msg)
        if pgmax < tol:
            best = (x, f, g, pgmax, msg)
            break
    x, f, g, pgmax, msg = best
    on = [i for i in np.flatnonzero(free) if x[i] <= lb[i] or x[i] >= ub[i]]
    return Fit(theta=x, loglik=float(f), grad=g, proj_max_grad=pgmax, converged=pgmax < tol,
               attempts=attempt + 1, seconds=time.perf_counter() - t0, at_bound=on,
               lbfgs_message=msg)


def fit_section3(D, start=None, tol=TOL):
    """Section-3 MLE (v = 0). start None -> section 7 starting values."""
    L = D.L
    x0 = start_values(D) if start is None else np.array(start, float)
    x0[L.v] = 0.0
    # coordinates without information (genes with no row, all-zero z columns) are held fixed
    x0[L.alpha[~D.free3[L.alpha]]] = 0.0
    x0[L.r[~D.free3[L.r]]] = START_R
    x0[L.delta[~D.free3[L.delta]]] = 0.0
    scale = np.ones(L.size)
    return _run(lambda x, o: ma.loglik3(x, D, o), x0, L.lb, L.ub, D.free3, scale, tol)


def fit_section4(D, start, v_start, v_scale, fix_v=False, tol=TOL, nodes=20):
    """Section-4 fit from `start` with v = v_start. fix_v gives the profile at that v.
    v is optimized in units of v_scale inside L-BFGS-B; convergence is judged in v."""
    L = D.L
    x0 = np.array(start, float)
    x0[L.v] = v_start
    free = D.free4.copy()
    if fix_v:
        free[L.v] = False
    scale = np.ones(L.size)
    scale[L.v] = v_scale
    fun = lambda x, o: ma.loglik4(x, D, o, nodes)  # noqa: E731
    fit = _run(fun, x0, L.lb, L.ub, free, scale, tol)
    if not fix_v and fit.theta[L.v] == 0.0 and not fit.converged:
        # left at v = 0 with a positive boundary score: restart from v = v_scale
        x1 = fit.theta.copy()
        x1[L.v] = v_scale
        fit2 = _run(fun, x1, L.lb, L.ub, free, scale, tol)
        if fit2.loglik > fit.loglik:
            fit2.attempts += fit.attempts
            fit = fit2
    return fit


def fit_section4_two_starts(D, theta3, v_max, tol=TOL, nodes=20):
    """Section-4 fits started at v = 0 and v = tau_max^2 from the section-3 MLE; the better
    one is kept (section 5.3)."""
    f0 = fit_section4(D, theta3, 0.0, v_max, tol=tol, nodes=nodes)
    f1 = fit_section4(D, theta3, v_max, v_max, tol=tol, nodes=nodes)
    return (f0 if f0.loglik >= f1.loglik else f1), f0, f1


def lr_statistic(fit3, best4):
    """LR = 2 (max(l4(v_hat), l3) - l3). When v_hat = 0, l4(v_hat) is the section-3 maximum,
    so LR is 0 exactly (differences between the two fits' end points are rounding)."""
    v = best4.theta[-1]
    if v == 0.0:
        return 0.0
    return 2.0 * (max(best4.loglik, fit3.loglik) - fit3.loglik)


def sandwich_se(theta, D, free=None):
    """SE^2 = sum_c (e_kappa' H^-1 s_c)^2 with H the observed information of the free section-3
    parameters (bound coordinates removed; kappa always kept) and s_c the score summed over
    the rows of cluster c (section 5.1)."""
    L = D.L
    free = D.free3 if free is None else free
    ll, g, Hs = ma.loglik3(theta, D, 2)
    on_bound = ((theta <= L.lb) | (theta >= L.ub)) & free
    on_bound[L.kappa] = False
    idx = np.flatnonzero(free & ~on_bound)
    info = -Hs[np.ix_(idx, idx)]
    ek = (idx == L.kappa).astype(float)
    a = np.linalg.solve(info, ek)
    S_rows = ma.row_scores(theta, D)
    C = sp.csr_matrix((np.ones(D.nr), (D.row_cluster, np.arange(D.nr))),
                      shape=(len(D.cluster_codes), D.nr))
    Sc = (C @ S_rows).toarray()[:, idx]
    contrib = Sc @ a
    return float(np.sqrt(np.sum(contrib ** 2))), float(np.linalg.cond(info))
