"""Fits for Model A-final (spec v1 section 7 with v1.8 section 2): starting values, L-BFGS-B with
the section 3.4 bounds (v in [0, 1]) followed by a projected Newton polish with the v1.8 step
acceptance, the v1.8 convergence rule (projected gradient and Cholesky of the negative
projected Hessian), the final-point integrator check, restarts, the warm/cold start selection,
the individual-cluster sandwich SE for kappa (section 5.1), and (v1.8.3) the two-start observed
section-3 fit of LEAD_RULINGS_W4 H-11 M1 (fit_section3_observed)."""
import time
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize
from scipy.special import logit

import model_a as ma

TOL = 1e-6                  # projected max |gradient| (section 7; v1.8 section 2 item 1)
DECREMENT = 1e-12           # Newton decrement g'H^-1 g at which the polish stops
BOUND_TOL = 1e-8            # v1.8 section 2 item 1: on a bound when |x - bound| <= 1e-8 max(1, |bound|)
N_RESTARTS = 3
RESTART_SEED = [20260923, 90]
START_ALPHAS = (-1.0, 0.0, 1.0)
START_R = logit(0.02)
# L-BFGS-B ftol. v1.8.2: section 3 back at 1e-15 (LEAD_RULINGS_W4 H-9 rejected the v1.8.1 value 1e-10);
# section 4 1e-15 as in v1.8.
FTOL3 = 1e-15
FTOL4 = 1e-15


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
    hessian_ok: bool = False
    status: str = ""              # "ok", or why the fit is not certified
    repolished: bool = False      # v1.8 section 1 item 8
    n_reference_genes: int = 0    # genes on the reference integrator at the final point
    loglik_ref: float = np.nan    # section 4: every gene on the reference integrator
    v_at_upper: bool = False      # v1.8 section 1 item 10
    events: dict = field(default_factory=dict)
    multistart: dict = None       # v1.8.3: set only on the observed section-3 fit (fit_section3_observed)


def on_lower(x, lb):
    return np.isfinite(lb) & (np.abs(x - lb) <= BOUND_TOL * np.maximum(1.0, np.abs(lb)))


def on_upper(x, ub):
    return np.isfinite(ub) & (np.abs(x - ub) <= BOUND_TOL * np.maximum(1.0, np.abs(ub)))


def projected_gradient(g, x, lb, ub, free):
    """v1.8 section 2 item 1: on a lower bound a negative component is set to 0, on an upper
    bound a positive one."""
    pg = np.where(free, g, 0.0)
    pg[on_lower(x, lb) & (pg < 0)] = 0.0
    pg[on_upper(x, ub) & (pg > 0)] = 0.0
    return pg


def hessian_ok(H, x, lb, ub, free):
    """v1.8 section 2 item 1: Cholesky of the negative Hessian over the free coordinates that
    are not on a bound."""
    idx = np.flatnonzero(free & ~on_lower(x, lb) & ~on_upper(x, ub))
    if len(idx) == 0:
        return True
    A = -H[np.ix_(idx, idx)]
    if not np.all(np.isfinite(A)):
        return False
    try:
        np.linalg.cholesky(0.5 * (A + A.T))
        return True
    except np.linalg.LinAlgError:
        return False


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


class Eval3:
    """Section-3 objective (no integrator state)."""
    is4 = False
    last_choice = None

    def __init__(self, D):
        self.D = D

    def new_fit(self):
        pass

    def force(self, choice):
        pass

    def __call__(self, x, order):
        return ma.loglik3(x, self.D, order)

    def final(self, x):
        return ma.loglik3(x, self.D, 2)


class Eval4:
    """Section-4 objective with the v1.8 per-fit state (previous modes, first-evaluation and
    final-point global checks, integrator choice)."""
    is4 = True

    def __init__(self, D, nodes=20):
        self.D, self.nodes = D, nodes
        self.state = ma.Sec4State(D.L.G)
        self.forced = None

    def new_fit(self):
        self.state = ma.Sec4State(self.D.L.G)
        self.forced = None

    def force(self, choice):
        self.forced = None if choice is None else np.asarray(choice).copy()

    @property
    def last_choice(self):
        c = self.state.last_choice
        return None if c is None else c.copy()

    def __call__(self, x, order):
        return ma.loglik4(x, self.D, order, self.nodes, state=self.state, force=self.forced)

    def final(self, x):
        return ma.loglik4(x, self.D, 2, self.nodes, state=self.state, final=True)


def _same_choice(a, b):
    if a is None or b is None:
        return a is None and b is None
    return bool(np.array_equal(a, b))


def _polish(ev, x, f, g, H, lb, ub, free, tol, newton_iter):
    """Projected Newton polish. v1.8.2 step acceptance (v1.8 section 2 item 2 as amended under
    LEAD_RULINGS_W4 H-10): with eps_l = max(1e-9, 1e-13 |l_best|), a trial step is accepted if
    (l_trial >= l_best - eps_l and the step lowers the projected max |gradient|) or
    l_trial > l_best + eps_l; l_best is the highest l among accepted iterates (the polish start
    included). Rejected trials do not change l_best."""
    f_best = f
    for _ in range(newton_iter):
        pg = projected_gradient(g, x, lb, ub, free)
        pgmax = float(np.max(np.abs(pg)))
        fr = free & ~((on_lower(x, lb) & (g < 0)) | (on_upper(x, ub) & (g > 0)))
        if not np.all(np.isfinite(H[np.ix_(fr, fr)])):
            break
        d = _newton_direction(H[np.ix_(fr, fr)], g[fr])
        if pgmax < tol and g[fr] @ d < DECREMENT:
            break
        eps = max(1e-9, 1e-13 * abs(f_best))
        t, accepted = 1.0, False
        while t > 1e-12:
            xn = x.copy()
            xn[fr] = np.clip(x[fr] + t * d, lb[fr], ub[fr])
            try:
                fn, gn, Hn = ev(xn, 2)
            except (RuntimeError, FloatingPointError, np.linalg.LinAlgError):
                fn = -np.inf
            if np.isfinite(fn):
                pgn = float(np.max(np.abs(projected_gradient(gn, xn, lb, ub, free))))
                if (fn >= f_best - eps and pgn < pgmax) or fn > f_best + eps:
                    accepted = True
                    f_best = max(f_best, fn)
                    break
            t *= 0.5
        if not accepted:
            break
        x, f, g, H = xn, fn, gn, Hn
    return x, f, g, H


def maximize(ev, x0, lb, ub, free, scale, tol=TOL, newton_iter=60, ftol=FTOL4):
    """Maximize over the free coordinates: L-BFGS-B (scaled coordinates), a projected Newton
    polish, and the v1.8 final-point check (section 1 item 8, section 2 item 1).
    Returns x, l, gradient, projected max |gradient|, L-BFGS-B message, Hessian condition,
    status and whether a re-polish ran."""
    ev.new_fit()
    idx = np.flatnonzero(free)
    x = np.clip(np.array(x0, float), lb, ub)

    def neg(y):
        xx = x.copy()
        xx[idx] = np.clip(y * scale[idx], lb[idx], ub[idx])
        try:
            f, g = ev(xx, 1)
        except (RuntimeError, FloatingPointError, np.linalg.LinAlgError):
            return np.inf, np.zeros(len(idx))
        if not np.isfinite(f) or not np.all(np.isfinite(g[idx])):
            return np.inf, np.zeros(len(idx))
        return -f, -g[idx] * scale[idx]

    res = minimize(neg, x[idx] / scale[idx], jac=True, method="L-BFGS-B",
                   bounds=list(zip(lb[idx] / scale[idx], ub[idx] / scale[idx])),
                   options=dict(maxiter=5000, maxfun=10000, maxcor=30, ftol=ftol, gtol=1e-10))
    x[idx] = np.clip(res.x * scale[idx], lb[idx], ub[idx])
    msg = str(res.message)
    f, g, H = ev(x, 2)
    x, f, g, H = _polish(ev, x, f, g, H, lb, ub, free, tol, newton_iter)
    status, repolished = "ok", False
    choice_polish = ev.last_choice
    try:
        f, g, H = ev.final(x)
        choice_final = ev.last_choice
        if not _same_choice(choice_polish, choice_final):
            # item 8: re-polish with each gene's integrator held at its final-point choice
            repolished = True
            ev.force(choice_final)
            f, g, H = ev(x, 2)
            x, f, g, H = _polish(ev, x, f, g, H, lb, ub, free, tol, newton_iter)
            ev.force(None)
            f, g, H = ev.final(x)
            if not _same_choice(choice_final, ev.last_choice):
                status = "integrator changed twice at the final point"
    except (RuntimeError, FloatingPointError, np.linalg.LinAlgError) as err:
        status = f"final-point evaluation failed: {err!r}"
    pg = projected_gradient(g, x, lb, ub, free)
    hok = hessian_ok(H, x, lb, ub, free)
    return x, f, g, float(np.max(np.abs(pg))), msg, hok, status, repolished


def _jitter(x, lb, ub, free, scale, attempt):
    rng = np.random.default_rng(np.random.SeedSequence(RESTART_SEED + [attempt]))
    y = x + np.where(free, rng.normal(0.0, 0.1, len(x)) * scale, 0.0)
    return np.clip(y, lb, ub)


def hessian_scale(ev, x0, free, base):
    """L-BFGS-B coordinate scale 1/sqrt(|H_ii|) at the start (the z columns differ in scale by
    ~1e3, so unscaled L-BFGS-B needs many iterations); base where H_ii is unusable."""
    ev.new_fit()
    try:
        d = np.abs(np.diag(ev(x0, 2)[2]))
    except (RuntimeError, FloatingPointError, np.linalg.LinAlgError):
        return base
    ok = free & np.isfinite(d) & (d > 1e-12)
    return np.where(ok, 1.0 / np.sqrt(np.where(ok, d, 1.0)), base)


def _run(ev, start, lb, ub, free, scale, tol, ftol=FTOL4):
    """Fit with up to N_RESTARTS restarts (section 7). A fit converges when the projected max
    |gradient| < tol, the Hessian condition holds and the final-point check passed."""
    t0 = time.perf_counter()
    scale = hessian_scale(ev, np.clip(start, lb, ub), free, scale)
    x0 = start
    best = None
    for attempt in range(N_RESTARTS + 1):
        if attempt > 0:
            x0 = _jitter(best[0] if best is not None else start, lb, ub, free, scale, attempt)
        try:
            out = maximize(ev, x0, lb, ub, free, scale, tol, ftol=ftol)
        except (RuntimeError, FloatingPointError, np.linalg.LinAlgError) as err:
            out = (x0, -np.inf, np.full(len(x0), np.nan), np.inf, repr(err), False,
                   f"fit failed: {err!r}", False)
        x, f, g, pgmax, msg, hok, status, rep = out
        conv = pgmax < tol and hok and status == "ok"
        out = out + (conv, ev.last_choice, dict(ev.state.events) if ev.is4 else {})
        if best is None or f > best[1] or conv:
            best = out
        if conv:
            break
    x, f, g, pgmax, msg, hok, status, rep, conv, choice, events = best
    if status == "ok" and not conv:
        status = "projected gradient >= tol" if pgmax >= tol else "Hessian condition failed"
    on = [i for i in np.flatnonzero(free) if on_lower(x[i:i + 1], lb[i:i + 1])[0]
          or on_upper(x[i:i + 1], ub[i:i + 1])[0]]
    return Fit(theta=x, loglik=float(f), grad=g, proj_max_grad=pgmax, converged=conv,
               attempts=attempt + 1, seconds=time.perf_counter() - t0, at_bound=on,
               lbfgs_message=msg, hessian_ok=hok, status=status, repolished=rep,
               n_reference_genes=0 if choice is None else int(np.sum(choice == ma.REF)),
               events=events)


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
    fit = _run(Eval3(D), x0, L.lb, L.ub, D.free3, scale, tol, ftol=FTOL3)
    fit.loglik_ref = fit.loglik
    return fit


# v1.8.3 (LEAD_RULINGS_W4 H-11 M1): observed section-3 fit from n_obs = 2 starts.
N_OBS = 2
OBS_TIE = 1e-9                # ties in l within 1e-9 go to the section 7 start (start 0)
# The ruling's seed is SeedSequence([20261015, cell, r, "obs_start", 1]). SeedSequence takes non-negative
# integers only, so the stream label "obs_start" is entered as the integer of its ASCII bytes (big-endian).
# This cannot collide with the integer streams 0-6 and 9 of the Gate B scheme.
OBS_START_STREAM = int.from_bytes(b"obs_start", "big")
# B-FIX and synthetic runs have no (cell, r): their jittered start uses this fixed entropy.
FIXTURE_OBS_START_SEED = [20260923, 92, OBS_START_STREAM, 1]


def obs_start_seed(seed0, cell_id, r, k=1):
    """Gate B entropy of observed start k: [seed0, cell, r, OBS_START_STREAM, k] (seed0 = 20261015)."""
    return [int(seed0), int(cell_id), int(r), OBS_START_STREAM, int(k)]


def section3_start(D):
    """The x0 that fit_section3 builds from the section 7 start."""
    L = D.L
    x0 = start_values(D)
    x0[L.v] = 0.0
    x0[L.alpha[~D.free3[L.alpha]]] = 0.0
    x0[L.r[~D.free3[L.r]]] = START_R
    x0[L.delta[~D.free3[L.delta]]] = 0.0
    return x0


def jittered_start(D, seed):
    """The section 7 start jittered as in the section 7 restart recipe: x0 + N(0, 0.1) s_j on the free
    section-3 coordinates, s_j = 1/sqrt(|H_jj|) at the start (hessian_scale), clipped to the bounds; rng
    default_rng(SeedSequence(seed))."""
    L = D.L
    x0 = section3_start(D)
    free = D.free3
    scale = hessian_scale(Eval3(D), np.clip(x0, L.lb, L.ub), free, np.ones(L.size))
    rng = np.random.default_rng(np.random.SeedSequence(list(seed)))
    y = x0 + np.where(free, rng.normal(0.0, 0.1, len(x0)) * scale, 0.0)
    return np.clip(y, L.lb, L.ub)


def select_observed(fits, tie=OBS_TIE):
    """H-11 M1 selection among the observed starts (index 0 = section 7 start).
    Returns (index of the reported start, rule). Index None means the observed fit failed (no start
    converged). Among converged starts the highest l is reported; a start within `tie` of that maximum
    is a tie, and ties go to the lowest index (so to the section 7 start). Unconverged starts are ignored
    when any start converged."""
    conv = [i for i, f in enumerate(fits) if f.converged and np.isfinite(f.loglik)]
    if not conv:
        return None, "no start converged: observed fit failed"
    if len(conv) == 1:
        return conv[0], f"only start {conv[0]} converged"
    lmax = max(fits[i].loglik for i in conv)
    best = min(i for i in conv if fits[i].loglik >= lmax - tie)
    if best == 0:
        rule = "all converged, start 0 highest or tied within %g" % tie
    else:
        rule = f"all converged, start {best} higher by > {tie:g}"
    return best, rule


def fit_section3_observed(D, seed, n_obs=N_OBS, tol=TOL):
    """v1.8.3 observed section-3 fit (H-11 M1): start 0 is the section 7 start (fit_section3(D)); start 1
    is the jittered section 7 start with entropy `seed` (Gate B: obs_start_seed(20261015, cell, r)). Each
    start is a full fit with its own section 7 restarts. The reported fit is chosen by select_observed;
    if no start converged, the higher-l start is returned marked not converged. The returned Fit is the
    chosen start's Fit object, unchanged, with `multistart` set to a per-start record."""
    if n_obs != 2:
        raise ValueError("H-11 fixes n_obs = 2")
    fits = [fit_section3(D, tol=tol), fit_section3(D, start=jittered_start(D, seed), tol=tol)]
    idx, rule = select_observed(fits)
    if idx is None:
        idx = int(np.argmax([f.loglik if np.isfinite(f.loglik) else -np.inf for f in fits]))
        fits[idx].status = f"observed multistart: {rule}"
    l0 = fits[0].loglik
    fits[idx].multistart = dict(
        n_obs=n_obs, seed=[int(s) for s in seed], reported_start=int(idx), rule=rule,
        failed=not fits[idx].converged,
        starts=[dict(start=k, converged=bool(f.converged), attempts=int(f.attempts), status=str(f.status),
                     seconds=float(f.seconds), dl_vs_start0=float(f.loglik - l0),
                     bitwise_equal_start0=bool(np.array_equal(f.theta, fits[0].theta)))
                for k, f in enumerate(fits)],
        seconds_total=float(sum(f.seconds for f in fits)))
    return fits[idx]


def reference_loglik(theta, D):
    """Section-4 l with the reference integrator for every gene (v1.8 section 2 item 3);
    the exact section-3 value at v = 0."""
    if theta[D.L.v] == 0.0:
        return ma.loglik3(theta, D, 0)
    return ma.loglik4(theta, D, 0, all_reference=True)


def fit_section4(D, start, v_start, v_scale, fix_v=False, tol=TOL, nodes=20):
    """Section-4 fit from `start` with v = v_start. fix_v gives the profile at that v.
    v is optimized in units of v_scale inside L-BFGS-B; convergence is judged in v.
    The end point is also evaluated with the reference integrator for every gene."""
    L = D.L
    x0 = np.array(start, float)
    x0[L.v] = v_start
    free = D.free4.copy()
    if fix_v:
        free[L.v] = False
    scale = np.ones(L.size)
    scale[L.v] = v_scale
    ev = Eval4(D, nodes)
    fit = _run(ev, x0, L.lb, L.ub, free, scale, tol)
    if not fix_v and on_lower(fit.theta[L.v:L.v + 1], L.lb[L.v:L.v + 1])[0] and not fit.converged:
        # left at v = 0 with a positive boundary score: restart from v = v_scale
        x1 = fit.theta.copy()
        x1[L.v] = v_scale
        fit2 = _run(ev, x1, L.lb, L.ub, free, scale, tol)
        if fit2.converged or fit2.loglik > fit.loglik:
            fit2.attempts += fit.attempts
            fit = fit2
    try:
        fit.loglik_ref = float(reference_loglik(fit.theta, D))
    except (RuntimeError, FloatingPointError) as err:
        fit.loglik_ref = np.nan
        if fit.converged:
            fit.converged = False
            fit.status = f"reference evaluation at the end point failed: {err!r}"
    fit.v_at_upper = bool(on_upper(fit.theta[L.v:L.v + 1], L.ub[L.v:L.v + 1])[0])
    return fit


def eps_tie4(D):
    """v1.8 section 2 item 3: max(1e-6, 1e-8 G_fit) for section-4 comparisons."""
    return max(1e-6, 1e-8 * D.n_genes_with_rows)


def select(warm, cold, eps, use_ref=False):
    """v1.8 section 2 item 3 start selection. Returns (kept fit or None for a failed refit,
    rule). Section-4 comparisons use the reference-integrator l of both end points."""
    def val(f):
        x = f.loglik_ref if use_ref else f.loglik
        return x if np.isfinite(x) else -np.inf
    lw, lc = val(warm), val(cold)
    if warm.converged and cold.converged:
        if abs(lw - lc) >= eps:
            return (warm, "both converged, warm higher") if lw > lc else \
                (cold, "both converged, cold higher")
        return warm, "both converged, tie: warm kept"
    if warm.converged or cold.converged:
        conv, other, name = (warm, cold, "warm") if warm.converged else (cold, warm, "cold")
        if val(conv) >= val(other) - eps:
            return conv, f"only {name} converged, kept"
        return None, f"only {name} converged, unconverged fit higher by > eps: failed refit"
    return None, "neither converged: failed refit"


def fit_section4_two_starts(D, theta3, v_max, tol=TOL, nodes=20):
    """Section-4 fits started at v = 0 and v = tau_max^2 from theta3 (section 5.3). The kept
    one follows the v1.8 section 2 item 3 rule with the v = 0 start in the warm role. If that
    rule fails, the higher of the two is returned marked not converged."""
    f0 = fit_section4(D, theta3, 0.0, v_max, tol=tol, nodes=nodes)
    f1 = fit_section4(D, theta3, v_max, v_max, tol=tol, nodes=nodes)
    kept, rule = select(f0, f1, eps_tie4(D), use_ref=True)
    if kept is None:
        import copy
        kept = copy.copy(f0 if np.nan_to_num(f0.loglik_ref, nan=-np.inf)
                         >= np.nan_to_num(f1.loglik_ref, nan=-np.inf) else f1)
        kept.converged = False
        kept.status = "two v starts: " + rule
    return kept, f0, f1, rule


def lr_statistic(fit3, best4):
    """LR = 2 (max(l4(v_hat), l3) - l3) with l4 the reference-integrator value (v1.8 section 2
    item 3). When v_hat = 0, l4(v_hat) is the section-3 maximum, so LR is 0 exactly."""
    v = best4.theta[best4_v_index(best4)]
    if v == 0.0:
        return 0.0
    return 2.0 * (max(best4.loglik_ref, fit3.loglik) - fit3.loglik)


def lr_statistic_gh(fit3, best4):
    """The same LR with the GH20 l of the kept fit (reported beside the primary)."""
    v = best4.theta[best4_v_index(best4)]
    if v == 0.0:
        return 0.0
    return 2.0 * (max(best4.loglik, fit3.loglik) - fit3.loglik)


def best4_v_index(fit):
    return len(fit.theta) - 1


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
