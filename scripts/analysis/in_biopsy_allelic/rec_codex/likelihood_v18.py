"""REC Model A under spec v1.8 sections 1-2. Synthetic fixtures only until Gate B.

The likelihood, row arithmetic, predictors, bounds (except the new v <= 1),
starts and sandwich come unchanged from likelihood_v1.py. This module adds:

  v1.8 §1  bracketed mode search for adaptive GH20, the global grid check,
           the composite Gauss-Legendre reference integrator, the GH20 vs
           reference comparison, the final-point integrator re-polish rule
           and the upper bound v <= 1;
  v1.8 §2  the convergence criterion (projected score < 1e-6 and a Cholesky
           factor of the negative projected Hessian), Newton-polish step
           acceptance against l_best - eps_l, and the start-selection table.

Implementation choices not fixed by v1.8 (each is logged in model.stats):
  * l_best is the highest l among ACCEPTED polish points (the start counts).
  * A Newton step in the local search is kept when h falls by no more than
    1e-10*max(1,|h|) (rounding); golden steps use the strict comparison.
  * Two maxima count as one mode when h at their midpoint is not below the
    lower one by more than 1e-10*max(1,|h|) (no valley between them).
  * Local search stops after 400 iterations; a search that has not met the
    stop rule is treated like trigger (ii) (a global check runs).
  * When no local mode has a defined sigma (h'' >= 0 everywhere found), the
    reference uses sigma_min = U_g / 4096.
  * The reference gradient/Hessian sum drops nodes whose normalised weight
    is below 1e-18 (below double rounding of the weighted sum).
  * Double precision throughout (v1.8 §1.9 allows this if §4 tests pass).
  * Grid scans and reference-integrator values use a node-vectorised copy of
    the likelihood_v1 row value (row_values_vec); GH20 values, all gradients
    and Hessians use likelihood_v1.conditional_terms unchanged. The two agree
    to rounding (checked in the accuracy job).
"""
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import betaln, expit, gammaln, logsumexp

import likelihood_v1 as L
from likelihood_v1 import conditional_terms

GH_X, GH_W = np.polynomial.hermite.hermgauss(20)
GL_X, GL_W = np.polynomial.legendre.leggauss(20)
GOLD = (3 - math.sqrt(5)) / 2
V_UPPER = 1.0
SCORE_TOL = 1e-6
PANEL_CAP = 8192
GRID_CAP = 8193
SECOND_MODE_GAP = 25.0
MAX_LOCAL_ITER = 400
NODE_WEIGHT_FLOOR = 1e-18
ROUND_REL = 1e-10


class IntegrationFailure(ValueError):
    """A v1.8 evaluation failure (reference did not converge at 8,192 panels,
    or a gene changed integrator twice at the final point). Counts as a
    failed fit, like any ValueError in the v1 drivers."""


def on_bound(x, bound):
    return bound is not None and abs(x - bound) <= 1e-8 * max(1., abs(bound))


def projection(p, g, bounds):
    """v1.8 §2.1. Returns projected gradient, Hessian-check mask, step mask."""
    proj = np.array(g, float, copy=True)
    hess_free = np.ones(len(p), bool)
    step_free = np.ones(len(p), bool)
    for j, (lo, hi) in enumerate(bounds):
        if lo is not None and hi is not None and lo == hi:
            proj[j] = 0.; hess_free[j] = False; step_free[j] = False
            continue
        if on_bound(p[j], lo):
            hess_free[j] = False
            if g[j] < 0:
                proj[j] = 0.; step_free[j] = False
        if on_bound(p[j], hi):
            hess_free[j] = False
            if g[j] > 0:
                proj[j] = 0.; step_free[j] = False
    return proj, hess_free, step_free


def clip(p, bounds):
    out = np.array(p, float, copy=True)
    for j, (lo, hi) in enumerate(bounds):
        out[j] = min(max(out[j], -np.inf if lo is None else lo), np.inf if hi is None else hi)
    return out


def hessian_condition(h, mask):
    if not mask.any():
        return True
    sub = -np.asarray(h)[np.ix_(mask, mask)]
    if not np.isfinite(sub).all():
        return False
    try:
        np.linalg.cholesky((sub + sub.T) / 2)
        return True
    except np.linalg.LinAlgError:
        return False


def polish_directions(h, g, mask, p, bounds):
    """Newton direction (eigen-floored) then a bounded quadratic direction.
    Directions are proposals only; acceptance uses the exact objective."""
    h = np.asarray(h)
    usable = mask & np.isfinite(np.diag(h))
    free = np.flatnonzero(usable)
    if not len(free):
        return
    curv = -h[np.ix_(free, free)]
    if not np.isfinite(curv).all():
        return
    curv = (curv + curv.T) / 2
    vals, vecs = np.linalg.eigh(curv)
    floor = max(1e-9, float(np.max(np.abs(vals))) * 1e-10)
    curv = (vecs * np.maximum(vals, floor)) @ vecs.T
    score = np.asarray(g)[free]
    step = np.zeros(len(p))
    step[free] = np.linalg.solve(curv, score)
    yield 'newton', step
    limits = [(None if bounds[j][0] is None else bounds[j][0] - p[j],
               None if bounds[j][1] is None else bounds[j][1] - p[j]) for j in free]
    quad = minimize(lambda s: (.5 * s @ curv @ s - score @ s, curv @ s - score), np.zeros(len(free)),
                    jac=True, method='SLSQP', bounds=limits, options={'ftol': 1e-14, 'maxiter': 1000})
    step = np.zeros(len(p))
    step[free] = quad.x
    yield 'bounded_quadratic', step


def newton_polish(fun, p, bounds, maxiter=60, stats=None):
    """v1.8 §2.2. fun(p, hessian) -> (l, g, H or None)."""
    p = clip(p, bounds)
    ll, g, h = fun(p, True)
    best = ll
    history = []
    for iteration in range(maxiter):
        proj, _, step_free = projection(p, g, bounds)
        pmax = float(np.max(np.abs(proj), initial=0.))
        history.append({'iteration': iteration, 'loglik': ll, 'max_projected_gradient': pmax})
        if pmax < SCORE_TOL:
            break
        eps = max(1e-9, 1e-13 * abs(best))
        accepted = False
        for kind, step in polish_directions(h, g, step_free, p, bounds):
            for damping in 2. ** -np.arange(31):
                trial = clip(p + damping * step, bounds)
                if np.array_equal(trial, p):
                    continue
                lt, gt, _ = fun(trial, False)
                if not (np.isfinite(lt) and np.isfinite(gt).all()):
                    continue
                pt = float(np.max(np.abs(projection(trial, gt, bounds)[0]), initial=0.))
                if lt >= best - eps and pt < pmax:
                    p, accepted = trial, True
                    best = max(best, lt)
                    if stats is not None:
                        stats['polish_steps_' + kind] = stats.get('polish_steps_' + kind, 0) + 1
                    break
            if accepted:
                break
        if not accepted:
            if stats is not None:
                stats['polish_stalls'] = stats.get('polish_stalls', 0) + 1
            break
        ll, g, h = fun(p, True)
    return p, history


def select_start(fits, eps, key='loglik'):
    """v1.8 §2.3 table. Returns (kept name or None, rule)."""
    w, c = fits['warm'], fits['cold']
    lw, lc = w[key], c[key]
    if w['converged'] and c['converged']:
        if abs(lw - lc) >= eps:
            return ('warm' if lw > lc else 'cold'), 'both converged, differ >= eps: higher'
        return 'warm', 'both converged, differ < eps: warm'
    if w['converged'] or c['converged']:
        name, other = ('warm', 'cold') if w['converged'] else ('cold', 'warm')
        if fits[name][key] >= fits[other][key] - eps:
            return name, 'one converged, not more than eps below the other: converged'
        return None, 'one converged but more than eps below the unconverged start: failed refit'
    return None, 'neither start converged: failed refit'


def _log1p_sum(m, x):
    """sum_{i<m} log1p(i/x) for each element of x (likelihood_v1 rising-factorial pieces)."""
    if m <= 0:
        return np.zeros(len(x))
    i = np.arange(m)
    out = np.empty(len(x))
    chunk = max(1, 4_000_000 // m)
    for s in range(0, len(x), chunk):
        out[s:s + chunk] = np.log1p(i[None, :] / x[s:s + chunk, None]).sum(axis=1)
    return out


def bb_logmass_vec(k, n, aa, bb):
    """likelihood_v1.bb_logmass for a vector of (aa, bb); rows = nodes, columns = k."""
    k = np.atleast_1d(np.asarray(k, dtype=float))
    out = np.empty((len(aa), len(k)))
    if not len(k):
        return out
    t = aa + bb
    big = t > 1e4
    breaks = np.r_[0, np.flatnonzero(np.diff(k) != 1) + 1, len(k)]
    for start, stop in zip(breaks[:-1], breaks[1:]):
        first = k[start]
        choose = gammaln(n + 1) - gammaln(first + 1) - gammaln(n - first + 1)
        col = np.empty(len(aa))
        if big.any():
            left, right = int(first), int(n - first)
            a_, b_, t_ = aa[big], bb[big], t[big]
            col[big] = (choose + left * np.log(a_ / t_) + right * np.log(b_ / t_)
                        + _log1p_sum(left, a_) + _log1p_sum(right, b_) - _log1p_sum(int(n), t_))
        if (~big).any():
            a_, b_ = aa[~big], bb[~big]
            col[~big] = choose + betaln(first + a_, n - first + b_) - betaln(a_, b_)
        out[:, start] = col
        counts = k[start:stop - 1]
        if len(counts):
            inc = (np.log(n - counts)[None, :] + np.log(counts[None, :] + aa[:, None])
                   - np.log(counts + 1)[None, :] - np.log(n - counts - 1 + bb[:, None]))
            out[:, start + 1:stop] = col[:, None] + np.cumsum(inc, axis=1)
    return out


def row_values_vec(row, eta, linear_r):
    """likelihood_v1.conditional_terms(...)[0] for a vector of eta (value only)."""
    eta = np.asarray(eta, float) - row.orientation * row.omega
    m, one_m = expit(eta), expit(-eta)
    v = expit(linear_r)
    rho = L.EPS + (1 - 2 * L.EPS) * v
    t = (1 - rho) / rho
    aa, bb = m * t, one_m * t
    with np.errstate(divide='ignore'):
        logprior = np.log(row.prior)
    lo = max(row.minor, (row.percent * row.n + 99) // 100)
    if row.n < row.depth or not lo <= row.a <= row.n - lo:
        raise ValueError('row fails integer het predicate')
    het_a = bb_logmass_vec([row.a], row.n, aa, bb)[:, 0]
    hom_a = (L.fixed_bb(np.array([row.a]), row.n, row.e, row.phi)[0],
             L.fixed_bb(np.array([row.a]), row.n, 1 - row.e, row.phi)[0])
    mix = logsumexp(np.stack((het_a + logprior[0], np.full_like(het_a, hom_a[0] + logprior[1]),
                              np.full_like(het_a, hom_a[1] + logprior[2]))), axis=0)
    tails = np.r_[np.arange(lo), np.arange(row.n - lo + 1, row.n + 1)]
    if len(tails):
        logtail = logsumexp(bb_logmass_vec(tails, row.n, aa, bb), axis=1)
        logc = np.empty(len(eta))
        small = logtail < np.log(.5)
        logc[small] = np.log1p(-np.exp(logtail[small]))
        if (~small).any():
            logc[~small] = logsumexp(bb_logmass_vec(np.arange(lo, row.n - lo + 1), row.n,
                                                    aa[~small], bb[~small]), axis=1)
    else:
        logc = np.zeros(len(eta))
    hom = L.hom_log_inclusion(row.n, lo, row.e, row.phi)
    den = logsumexp(np.stack((logprior[0] + logc, np.full_like(logc, logprior[1] + hom),
                              np.full_like(logc, logprior[2] + hom))), axis=0)
    return mix - den


class ModelV18(L.Model):
    def __init__(self, rows, gene_count=None):
        super().__init__(list(rows), gene_count)
        self.gene_rows = {g: [r for r in self.rows if r.gene == g] for g in range(self.G)}
        self.genes_with_rows = sum(1 for g in range(self.G) if self.gene_rows[g])
        self.stats = {}
        self.events = []
        self._reset_fit_state()

    @classmethod
    def from_model(cls, model):
        return cls(model.rows, model.G)

    # ---------------------------------------------------------------- state
    def _reset_fit_state(self):
        self._prev_mode = {}
        self._seen = set()
        self._held = None
        self.last_choices = {}

    def _count(self, key, k=1):
        self.stats[key] = self.stats.get(key, 0) + k

    def _event(self, record):
        self._count('events_total')
        if len(self.events) < 2000:
            self.events.append(record)

    # ------------------------------------------------------------ gene h(u)
    def _gene_pre(self, p, g):
        pre = []
        for row in self.gene_rows[g]:
            eta, r, _, _ = self.predictors(row, p, 0., need_hessian=False)
            pre.append((row, eta, r))
        return pre

    @staticmethod
    def _h_vec(pre, v, u):
        u = np.asarray(u, float)
        total = -u * u / (2 * v)
        for row, eta, r in pre:
            total = total + row_values_vec(row, eta + u * row.stage, r)
        return total

    @staticmethod
    def _h(pre, v, u, derivs=False):
        if not derivs:
            total = -u * u / (2 * v)
            for row, eta, r in pre:
                total += conditional_terms(row, eta + u * row.stage, r, False)[0]
            return total
        total, d1, d2 = -u * u / (2 * v), -u / v, -1. / v
        for row, eta, r in pre:
            val, s, hh = conditional_terms(row, eta + u * row.stage, r, True)
            total += val; d1 += s[0] * row.stage; d2 += hh[0, 0] * row.stage ** 2
        return total, d1, d2

    def _local(self, pre, v, a, c, b, ha, hc, hb):
        """v1.8 §1.3 bracketed Newton / golden-section search."""
        nonconcave = False
        last_du = None
        converged = False
        for _ in range(MAX_LOCAL_ITER):
            hc, d1, d2 = self._h(pre, v, c, True)
            if d2 >= 0:
                nonconcave = True
            if last_du is not None and last_du < 1e-10 and abs(d1) < 1e-8 * max(1., abs(hc)):
                converged = True
                break
            x = None
            newton = False
            if d2 < 0:
                xn = c - d1 / d2
                if a < xn < b and xn != c:
                    x, newton = xn, True
            if x is None:
                x = c + GOLD * (b - c) if (b - c) > (c - a) else c - GOLD * (c - a)
            if not a < x < b or x == c:
                # Bracket below floating resolution.
                last_du = 0.
                if abs(d1) < 1e-8 * max(1., abs(hc)):
                    converged = True
                break
            hx = self._h(pre, v, x)
            last_du = abs(x - c)
            # A Newton step is kept when h falls by no more than rounding (ROUND_REL);
            # otherwise a rounding-high point can block the stationarity stop rule.
            if hx >= hc or (newton and hx >= hc - ROUND_REL * max(1., abs(hc))):
                if x > c:
                    a, ha = c, hc
                else:
                    b, hb = c, hc
                c, hc = x, hx
            else:
                if x > c:
                    b, hb = x, hx
                else:
                    a, ha = x, hx
        hc, d1, d2 = self._h(pre, v, c, True)
        if d2 >= 0:
            nonconcave = True
        return {'u': c, 'h': hc, 'd1': d1, 'd2': d2, 'nonconcave': nonconcave, 'converged': converged}

    @staticmethod
    def _best(modes):
        top = max(m['h'] for m in modes)
        tied = [m for m in modes if m['h'] >= top - 1e-10]
        return min(tied, key=lambda m: (abs(m['u']), m['u']))

    def _distinct(self, modes, pre, v):
        """Two maxima are one mode if they are within 1e-8*max(1,|u|) or if h at
        their midpoint is not below the lower of the two by more than rounding."""
        out = []
        for m in sorted(modes, key=lambda m: (not m['converged'], -m['h'])):
            same = False
            for o in out:
                if abs(m['u'] - o['u']) <= 1e-8 * max(1., abs(o['u'])):
                    same = True
                    break
                mid = self._h(pre, v, .5 * (m['u'] + o['u']))
                low = min(m['h'], o['h'])
                if mid >= low - ROUND_REL * max(1., abs(low)):
                    same = True
                    self._count('modes_merged_by_midpoint')
                    break
            if not same:
                out.append(m)
        return out

    def _sigma(self, mode):
        return (-mode['d2']) ** -.5 if mode['d2'] < 0 else None

    def _mode_search(self, g, pre, v, final):
        h0 = self._h(pre, v, 0.)
        U = math.sqrt(2 * v * (-h0)) if h0 < 0 else math.sqrt(2 * v)
        hl, hr = self._h(pre, v, -U), self._h(pre, v, U)
        results = [self._local(pre, v, -U, 0., U, hl, h0, hr)]
        prev = self._prev_mode.get(g)
        if prev is not None and -U < prev < U and prev != 0.:
            hp = self._h(pre, v, prev)
            if hp >= max(hl, hr):
                results.append(self._local(pre, v, -U, prev, U, hl, hp, hr))
        triggers = []
        if len(results) == 2 and abs(results[0]['u'] - results[1]['u']) > 1e-8 * max(1., abs(results[0]['u'])):
            triggers.append('i')
        if any(r['nonconcave'] for r in results):
            triggers.append('ii')
        best = self._best(results)
        if best['d2'] > -1e-12:
            triggers.append('iii')
        if g not in self._seen:
            triggers.append('iv')
        if final:
            triggers.append('v')
        if not all(r['converged'] for r in results):
            triggers.append('local_search_unconverged')
        modes = list(results)
        grid = None
        force_reference = False
        if triggers:
            sigma = self._sigma(best)
            n = None if sigma is None else max(129, math.ceil(4 * U / sigma) + 1)
            if n is None or n > GRID_CAP:
                force_reference = True
                grid = {'points': None, 'reason': 'sigma undefined' if n is None else 'grid cap'}
                self._count('grid_cap_or_sigma_undefined')
            else:
                us = np.linspace(-U, U, n)
                hs = self._h_vec(pre, v, us)
                peaks = [j for j in range(1, n - 1) if hs[j] > hs[j - 1] and hs[j] > hs[j + 1]]
                for j in peaks:
                    modes.append(self._local(pre, v, float(us[j - 1]), float(us[j]), float(us[j + 1]),
                                             float(hs[j - 1]), float(hs[j]), float(hs[j + 1])))
                grid = {'points': n, 'peaks': len(peaks)}
                self._count('grid_scans')
        distinct = self._distinct(modes, pre, v)
        best = self._best(distinct)
        others = [m for m in distinct if m is not best]
        multimodal = any(m['h'] >= best['h'] - SECOND_MODE_GAP for m in others)
        flat = best['d2'] > -1e-12
        return {'U': U, 'best': best, 'modes': distinct, 'triggers': triggers, 'grid': grid,
                'multimodal': multimodal, 'flat': flat,
                'use_reference': force_reference or multimodal or flat}

    # ------------------------------------------------------------ integrals
    def _node_integral(self, p, g, pre, v, u, logw, need_gradient, with_hessian, floor=0., hv=None):
        if hv is None:
            hv = np.array([self._h(pre, v, float(x)) for x in u])
        terms = logw + hv
        logsum = logsumexp(terms)
        value = float(logsum - .5 * math.log(2 * math.pi * v))
        if not need_gradient:
            return value, None, None
        prob = np.exp(terms - logsum)
        keep = np.flatnonzero(prob > floor)
        rs = self.gene_rows[g]
        grad = np.zeros(self.size + 1)
        hess = np.zeros((self.size + 1, self.size + 1)) if with_hessian else None
        second = np.zeros((self.size + 1, self.size + 1)) if with_hessian else None
        for i in keep:
            x = float(u[i])
            _, gg, hh = self.evaluate(p, rs, x, need_hessian=with_hessian)
            s = np.r_[gg, (x * x - v) / (2 * v * v)]
            grad += prob[i] * s
            if with_hessian:
                full = np.zeros((self.size + 1, self.size + 1))
                full[:-1, :-1] = hh
                full[-1, -1] = 1 / (2 * v * v) - x * x / (v ** 3)
                hess += prob[i] * full
                second += prob[i] * np.outer(s, s)
        if with_hessian:
            hess += second - np.outer(grad, grad)
        return value, grad, hess

    def _gh_nodes(self, mode, nodes=20):
        x, w = (GH_X, GH_W) if nodes == 20 else np.polynomial.hermite.hermgauss(nodes)
        sigma = self._sigma(mode)
        u = mode['u'] + math.sqrt(2) * sigma * x
        logw = np.log(w) + x * x + math.log(math.sqrt(2) * sigma)
        return u, logw

    def _reference_nodes(self, pre, v, search):
        modes = search['modes']
        best = search['best']
        sigmas = [s for s in (self._sigma(m) for m in modes) if s is not None]
        sigma_min = min(sigmas) if sigmas else search['U'] / 4096
        if not sigmas:
            self._count('reference_sigma_fallback')
        sigma_best = self._sigma(best) or sigma_min
        W = math.sqrt(2 * v * (-best['h'] + 37 + max(0., math.log(math.sqrt(v) / sigma_best))))
        edges = sorted({-W, W} | {m['u'] for m in modes if -W < m['u'] < W})
        lengths = np.diff(edges)
        counts = np.maximum(1, np.ceil(lengths / sigma_min)).astype(int)
        previous = None
        while True:
            if counts.sum() > PANEL_CAP:
                raise IntegrationFailure('reference integrator: %d panels without 1e-10 agreement' % counts.sum())
            us, lw = [], []
            for lo, length, k in zip(edges[:-1], lengths, counts):
                width = length / k
                centres = lo + width * (np.arange(k) + .5)
                us.append((centres[:, None] + .5 * width * GL_X[None, :]).ravel())
                lw.append(np.tile(np.log(GL_W * .5 * width), k))
            u, logw = np.concatenate(us), np.concatenate(lw)
            hv = self._h_vec(pre, v, u)
            value = float(logsumexp(logw + hv) - .5 * math.log(2 * math.pi * v))
            if previous is not None and abs(value - previous) <= 1e-10:
                return value, u, logw, int(counts.sum()), W, hv
            previous = value
            counts = counts * 2

    def _gene_integral(self, p, g, v, need_gradient, with_hessian, final, force_reference, nodes=20):
        pre = self._gene_pre(p, g)
        search = self._mode_search(g, pre, v, final)
        self._seen.add(g)
        self._prev_mode[g] = search['best']['u']
        choice = 'REF' if search['use_reference'] else 'GH20'
        record = {'gene': g, 'v': v, 'triggers': search['triggers'], 'modes': len(search['modes']),
                  'multimodal': search['multimodal'], 'flat': search['flat']}
        ref = None
        held = self._held is not None and not force_reference
        if search['triggers'] and choice == 'GH20' and not force_reference and not held:
            gh_value = self._node_integral(p, g, pre, v, *self._gh_nodes(search['best'], nodes), False, False)[0]
            ref = self._reference_nodes(pre, v, search)
            record['gh20_minus_reference'] = gh_value - ref[0]
            self._count('reference_comparisons')
            if abs(gh_value - ref[0]) > 1e-8:
                choice = 'REF'
                self._count('item7_switch_to_reference')
                self._event({'event': 'item7_switch', **record})
        if force_reference:
            choice = 'REF'
        if held:
            choice = self._held.get(g, choice)
        if choice == 'GH20' and self._sigma(search['best']) is None:
            choice = 'REF'
            self._count('gh20_without_sigma_to_reference')
        self.last_choices[g] = choice
        self._count('integrator_' + choice)
        if choice == 'GH20':
            u, logw = self._gh_nodes(search['best'], nodes)
            value, grad, hess = self._node_integral(p, g, pre, v, u, logw, need_gradient, with_hessian)
        else:
            if ref is None:
                ref = self._reference_nodes(pre, v, search)
            value, u, logw, panels, W, hv = ref
            record['reference_panels'] = panels
            if need_gradient:
                _, grad, hess = self._node_integral(p, g, pre, v, u, logw, True, with_hessian,
                                                    floor=NODE_WEIGHT_FLOOR, hv=hv)
            else:
                grad = hess = None
        if search['multimodal'] or search['flat'] or choice == 'REF':
            self._event({'event': 'gene_integrator', 'choice': choice, **record})
        return value, grad, hess, choice

    def marginal(self, p, v, nodes=20, with_hessian=False, gene_terms=None, final=False,
                 force_reference=False, need_gradient=True):
        if v < 0:
            raise ValueError('negative variance')
        if v > V_UPPER * (1 + 1e-12):
            raise ValueError('variance above the v1.8 upper bound')
        if v == 0:
            out = super().marginal(p, 0., nodes, with_hessian, gene_terms)
            self.last_choices = {g: 'EXACT' for g in range(self.G)}
            return out
        total = 0.
        grad = np.zeros(self.size + 1)
        hess = np.zeros((self.size + 1, self.size + 1)) if with_hessian else None
        choices = {}
        for g in range(self.G):
            if not self.gene_rows[g]:
                if gene_terms is not None:
                    gene_terms[g] = 0.
                continue
            value, gg, hh, choice = self._gene_integral(p, g, v, need_gradient, with_hessian, final,
                                                        force_reference, nodes)
            choices[g] = choice
            total += value
            if gene_terms is not None:
                gene_terms[g] = value
            if need_gradient:
                grad += gg
                if with_hessian:
                    hess += hh
        self.last_choices = choices
        if not need_gradient:
            return float(total)
        if with_hessian:
            return float(total), grad, hess
        return float(total), grad

    def reference_loglik(self, p, v):
        """Every gene by the reference integrator (v1.8 §2.3 comparisons, P2 profile values)."""
        if v == 0:
            return float(super().marginal(p, 0.)[0])
        saved = (dict(self._prev_mode), set(self._seen), self._held, dict(self.last_choices))
        try:
            return self.marginal(p, v, final=True, force_reference=True, need_gradient=False)
        finally:
            self._prev_mode, self._seen, self._held, self.last_choices = saved

    # ----------------------------------------------------- section 3 fits
    def _section3_bounds(self, start):
        informative = self.informative()
        return [b if informative[j] else (start[j], start[j]) for j, b in enumerate(self.bounds())]

    def _certify3(self, p, bounds):
        ll, g, h = self.evaluate(p)
        proj, hess_free, _ = projection(p, g, bounds)
        pmax = float(np.max(np.abs(proj), initial=0.))
        pd = hessian_condition(h, hess_free)
        return {'parameters': p, 'loglik': ll, 'gradient': g, 'hessian': h,
                'max_projected_gradient': pmax, 'hessian_negative_definite': pd,
                'converged': bool(pmax < SCORE_TOL and pd),
                'on_bound': [j for j, (lo, hi) in enumerate(bounds) if lo != hi and (on_bound(p[j], lo) or on_bound(p[j], hi))]}

    def fit_single(self, start):
        start = np.asarray(start, float)
        bounds = self._section3_bounds(start)

        def objective(p):
            ll, g, _ = self.evaluate(p, need_hessian=False)
            return -ll, -g
        candidates = []
        for shift in (0., -.1, .1):
            p = start.copy()
            p[2 * self.G] = np.clip(p[2 * self.G] + shift, -.9, .9)
            candidates.append(minimize(objective, p, jac=True, method='L-BFGS-B', bounds=bounds,
                                       options={'ftol': 1e-13, 'gtol': 1e-7, 'maxiter': 1000}))
        tried = []
        for k, cand in enumerate(sorted(candidates, key=lambda f: f.fun)):
            p, history = newton_polish(lambda q, hess: self.evaluate(q, need_hessian=hess), cand.x, bounds,
                                       stats=self.stats)
            result = self._certify3(p, bounds)
            result.update(restart_index=k, polish_iterations=len(history))
            tried.append(result)
            if result['converged']:
                break
            self._count('section3_restart')
        best = next((t for t in tried if t['converged']), max(tried, key=lambda t: t['loglik']))
        return {**best, 'restarts_tried': len(tried)}

    def fit(self, start=None):
        cold = self.cold_start()
        if start is None:
            fit = self.fit_single(cold)
            return {**fit, 'kept_start': 'cold', 'start_logliks': {'cold': fit['loglik']},
                    'start_rule': 'single cold start (observed fit)'}
        warm = np.asarray(start, float)
        informative = self.informative()
        cold[~informative] = warm[~informative]
        fits = {'warm': self.fit_single(warm), 'cold': self.fit_single(cold)}
        kept, rule = select_start(fits, 1e-6)
        if kept is None:
            name = max(fits, key=lambda k: fits[k]['loglik'])
            out = {**fits[name], 'converged': False}
        else:
            out = dict(fits[kept])
        return {**out, 'kept_start': kept, 'start_rule': rule,
                'start_logliks': {k: f['loglik'] for k, f in fits.items()},
                'start_converged': {k: f['converged'] for k, f in fits.items()}}

    # ----------------------------------------------------- section 4 fits
    def _certify4(self, p, bounds, final):
        out = self.marginal(p[:-1], p[-1], with_hessian=True, final=final)
        ll, g, h = out
        proj, hess_free, _ = projection(p, g, bounds)
        pmax = float(np.max(np.abs(proj), initial=0.))
        pd = hessian_condition(h, hess_free)
        return {'parameters': p[:-1], 'v': float(p[-1]), 'loglik': ll, 'gradient': g,
                'max_projected_gradient': pmax, 'hessian_negative_definite': pd,
                'converged': bool(pmax < SCORE_TOL and pd), 'integrators': dict(self.last_choices),
                'v_on_upper_bound': bool(abs(p[-1] - V_UPPER) <= 1e-8)}

    def fit_marginal_single(self, start, fixed_v):
        self._reset_fit_state()
        start = np.asarray(start, float)
        informative = self.informative()
        bounds = [b if informative[j] else (start[j], start[j]) for j, b in enumerate(self.bounds())]
        bounds += [(0., V_UPPER) if fixed_v is None else (fixed_v, fixed_v)]
        start = clip(start, bounds)

        def objective(p):
            ll, g = self.marginal(p[:-1], p[-1])
            return -ll, -g

        def fun(q, hess):
            if hess:
                return self.marginal(q[:-1], q[-1], with_hessian=True)
            ll, g = self.marginal(q[:-1], q[-1])
            return ll, g, None
        candidates = []
        for shift in (0., -.1, .1):
            p = start.copy()
            p[2 * self.G] = np.clip(p[2 * self.G] + shift, -.9, .9)
            candidates.append(minimize(objective, p, jac=True, method='L-BFGS-B', bounds=bounds,
                                       options={'ftol': 1e-13, 'gtol': 1e-7, 'maxiter': 500}))
        tried = []
        for k, cand in enumerate(sorted(candidates, key=lambda f: f.fun)):
            self._held = None
            p, history = newton_polish(fun, cand.x, bounds, stats=self.stats)
            if p[-1] > 0:
                fun(p, True)                      # last polish-step integrators at p
            last = dict(self.last_choices)
            result = self._certify4(p, bounds, final=True)
            note = None
            if p[-1] > 0 and result['integrators'] != last:
                # v1.8 §1.8: re-polish with integrators held at the final-point choice.
                self._count('item8_repolish')
                held = dict(result['integrators'])
                self._held = held
                p, history2 = newton_polish(fun, p, bounds, stats=self.stats)
                history += history2
                self._held = None
                check = self._certify4(p, bounds, final=True)
                if check['integrators'] != held:
                    self._count('item8_failed')
                    check['converged'] = False
                    note = 'gene changed integrator again at the re-polished final point'
                else:
                    note = 're-polished with held integrators'
                result = check
            result.update(restart_index=k, polish_iterations=len(history), item8=note)
            tried.append(result)
            if result['converged']:
                break
            self._count('section4_restart')
        best = next((t for t in tried if t['converged']), max(tried, key=lambda t: t['loglik']))
        return {**best, 'restarts_tried': len(tried)}

    def fit_marginal(self, tau_max, start=None, fixed_v=None, null_fit=None, profile_start=None,
                     section4_from_null=False):
        null = self.fit(start) if null_fit is None else null_fit
        p0 = np.asarray(null['parameters'], float)
        if fixed_v is not None:
            warm = np.r_[p0 if profile_start is None else profile_start, fixed_v]
            cold = np.r_[p0, fixed_v]
        else:
            warm = np.r_[p0, 0.]
            cold = np.r_[p0 if section4_from_null else self.cold_start(), tau_max ** 2]
        cold[:-1][~self.informative()] = p0[~self.informative()]
        fits = {'warm': self.fit_marginal_single(warm, fixed_v),
                'cold': self.fit_marginal_single(cold, fixed_v)}
        for f in fits.values():
            f['loglik_reference'] = self.reference_loglik(f['parameters'], f['v'])
        eps = max(1e-6, 1e-8 * self.genes_with_rows)
        kept, rule = select_start(fits, eps, key='loglik_reference')
        if kept is None:
            name = max(fits, key=lambda k: fits[k]['loglik_reference'])
            fit = {**fits[name], 'converged': False}
        else:
            fit = dict(fits[kept])
        label = kept
        if fixed_v is None and fit['v'] == 0. and fit['gradient'][-1] > SCORE_TOL:
            restart = self.fit_marginal_single(np.r_[fit['parameters'], tau_max ** 2], None)
            restart['loglik_reference'] = self.reference_loglik(restart['parameters'], restart['v'])
            if restart['loglik_reference'] > fit['loglik_reference']:
                fit, label = restart, 'boundary_restart'
        if fixed_v is None and fit['v'] == 0. and null['loglik'] >= fit['loglik'] - eps:
            ll, g = self.marginal(p0, 0.)
            fit = {**fit, 'parameters': p0, 'loglik': ll, 'loglik_reference': ll, 'gradient': g,
                   'converged': bool(null['converged'] and g[-1] <= SCORE_TOL and
                                     (kept is not None or fit['converged']))}
        return {**fit, 'null': null, 'kept_start': label, 'start_rule': rule, 'eps_tie': eps,
                'start_logliks': {k: f['loglik'] for k, f in fits.items()},
                'start_logliks_reference': {k: f['loglik_reference'] for k, f in fits.items()},
                'start_converged': {k: f['converged'] for k, f in fits.items()},
                'start_v': {k: f['v'] for k, f in fits.items()}}
