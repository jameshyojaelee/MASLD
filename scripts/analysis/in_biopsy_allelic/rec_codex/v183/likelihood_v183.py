"""REC Model A under spec v1.8.2 + v1.8.3 (2026-10-16). Synthetic fixtures only until Gate B.

Built from the spec text alone (spec sha256 285c254c...883a9b, lines 1604-1891). No B-MODEL code
or output was read. Subclass of likelihood_v18.ModelV18; the likelihood, the section-4 integral
and its integrators (v1.8 §1), the sandwich, bounds, §7 start and the v1.8 §2 item 3 table are
inherited unchanged. This module replaces:

  v1.8.2 item 1   step acceptance in every Newton polish (section 3, section 4, item-8 re-polish);
  v1.8.3 item 2   attempts within one fit: attempt 0 from the fit's start, up to 3 restarts, each
                  from a jitter of the highest-l attempt so far; stop at the first converged one;
  v1.8.3 item 3   the observed section-3 fit from two starts, selection rule 3.4 (1e-9 tie,
                  1e-6 "unconverged higher means failed"), record 3.5;
  v1.8.3 item 4   failed observed fit (selection failure, or sandwich SE not finite and positive);
  v1.8.3 item 7   evaluation failures reject a polish trial or make an attempt not converged;
                  they never stop the task. MemoryError and other system errors are not caught.

Implementation choices left open by the spec (v1.7 rule 8), recorded here:
  * Polish direction: (-H) d = g on the step coordinates after symmetric Jacobi scaling, each
    eigenvalue replaced by its absolute value and floored at 1e-10 times the largest. One
    direction per iteration (no second fallback direction). Iteration cap 60.
  * L-BFGS-B: section 3 ftol 1e-15, gtol 1e-7, maxiter 1000; section 4 ftol 1e-13, gtol 1e-7,
    maxiter 500 (the REC v1.8 section-4 settings, kept for cost). An evaluation failure inside
    L-BFGS-B returns f = 1e100 with a zero gradient, so the line search backs off.
  * Restart jitter (item 2.4): the item 3.2 form. s_j = 1/sqrt|H_jj| at the fit's own start
    (section-3 Hessian for section-3 fits, section-4 Hessian for section-4 fits; s_j = 1 where
    undefined), e ~ N(0, 0.1^2) from default_rng(SeedSequence([20260923, 90, a])) for restart a,
    drawn with one entry per theta coordinate including v, added on the free coordinates of the
    highest-l attempt so far, then clipped.
  * Section-4 integrator state (_prev_mode, _seen) is reset once per fit, as in REC v1.8; the
    Hessian used for restart scales is evaluated with that state saved and restored.
"""
import math
import time

import numpy as np
from scipy.optimize import minimize

import likelihood_v18 as V18
from likelihood_v18 import ModelV18, clip, hessian_condition, projection, select_start

SCORE_TOL = V18.SCORE_TOL            # 1e-6, v1.8 §2 item 1
MAX_RESTARTS = 3                     # §7; attempts 0..3
POLISH_CAP = 60
T_MIN = 1e-12                        # item 1.4: trials while t > 1e-12
DECREMENT_TOL = 1e-12                # item 1.5
OBS_TIE = 1e-9                       # item 3.4 step 2
OBS_FAIL_MARGIN = 1e-6               # item 3.4 step 3 (section-3 eps_tie)
EPS_TIE3 = 1e-6
JITTER_SD = 0.1
RESTART_ENTROPY = (20260923, 90)     # + (a,)
EVAL_ERRORS = (ValueError, FloatingPointError, np.linalg.LinAlgError, ZeroDivisionError, OverflowError)


def free_mask(bounds):
    return np.array([not (lo is not None and hi is not None and lo == hi) for lo, hi in bounds])


def newton_direction(h, g, step_free):
    """Ascent direction from the projected Newton system (Jacobi-scaled, |eigenvalue| floored)."""
    h = np.asarray(h, float)
    free = np.flatnonzero(step_free & np.isfinite(np.diag(h)))
    if not len(free):
        return None
    a = -h[np.ix_(free, free)]
    if not np.isfinite(a).all():
        return None
    a = (a + a.T) / 2
    dg = np.sqrt(np.abs(np.diag(a)))
    dg[~(dg > 0)] = 1.
    scaled = a / dg[:, None] / dg[None, :]
    vals, vecs = np.linalg.eigh(scaled)
    top = float(np.max(np.abs(vals)))
    if not (top > 0 and np.isfinite(top)):
        return None
    vals = np.maximum(np.abs(vals), 1e-10 * top)
    rhs = np.asarray(g, float)[free] / dg
    y = vecs @ ((vecs.T @ rhs) / vals)
    d = np.zeros(len(g))
    d[free] = y / dg
    return d


def newton_polish(fun, p, bounds, maxiter=POLISH_CAP, stats=None):
    """v1.8.2 item 1. fun(p, hessian) -> (l, g, H or None); may raise an EVAL_ERRORS error.
    Returns (end point, l at end point or nan, accepted steps, stop reason). An error evaluating the
    polish start propagates to the caller (attempt not converged)."""
    p = clip(p, bounds)
    ll, g, h = fun(p, True)
    best = ll                                   # 1.1: the polish start counts as accepted
    stop = 'iteration_cap'
    steps = 0
    for _ in range(maxiter):
        proj, _, step_free = projection(p, g, bounds)
        pg = float(np.max(np.abs(proj), initial=0.))
        d = newton_direction(h, g, step_free) if h is not None else None
        if d is None:
            stop = 'no_direction'
            break
        decrement = float(np.asarray(g)[step_free] @ d[step_free])
        if pg < SCORE_TOL and decrement < DECREMENT_TOL:
            stop = 'stationary'
            break
        eps = max(1e-9, 1e-13 * abs(best))       # 1.2, fixed for this iteration
        t = 1.
        accepted = None
        while t > T_MIN:
            trial = clip(p + t * d, bounds)
            t /= 2
            if np.array_equal(trial, p):
                continue
            try:
                lt, gt, _ = fun(trial, False)
            except EVAL_ERRORS:
                if stats is not None:
                    stats['polish_trial_eval_failure'] = stats.get('polish_trial_eval_failure', 0) + 1
                continue
            if not (np.isfinite(lt) and np.isfinite(gt).all()):
                continue
            pt = float(np.max(np.abs(projection(trial, gt, bounds)[0]), initial=0.))
            clause_a = lt >= best - eps and pt < pg
            clause_b = lt > best + eps
            if clause_a or clause_b:
                accepted = (trial, lt, 'a' if clause_a else 'b')
                break
        if accepted is None:
            stop = 'no_trial_accepted'
            break
        p, lt, clause = accepted
        steps += 1
        best = max(best, lt)
        if stats is not None:
            key = 'polish_accept_clause_' + clause
            stats[key] = stats.get(key, 0) + 1
        try:
            ll, g, h = fun(p, True)
        except EVAL_ERRORS:
            stop = 'evaluation_error_at_iterate'
            ll, h = float('nan'), None
            break
    if stats is not None:
        stats['polish_stop_' + stop] = stats.get('polish_stop_' + stop, 0) + 1
    return p, ll, steps, stop


class ModelV183(ModelV18):

    # ------------------------------------------------------------ helpers
    def _lbfgs(self, objective, x, bounds, ftol, maxiter):
        def safe(q):
            try:
                f, g = objective(q)
                if np.isfinite(f) and np.isfinite(g).all():
                    return f, g
            except EVAL_ERRORS:
                pass
            self._count('lbfgs_eval_failure')
            return 1e100, np.zeros(len(q))
        res = minimize(safe, x, jac=True, method='L-BFGS-B', bounds=bounds,
                       options={'ftol': ftol, 'gtol': 1e-7, 'maxiter': maxiter})
        return np.asarray(res.x, float)

    @staticmethod
    def _scales(hdiag, free):
        s = np.ones(len(free))
        if hdiag is None:
            return s
        hdiag = np.asarray(hdiag, float)
        ok = free & np.isfinite(hdiag) & (np.abs(hdiag) > 1e-12)
        s[ok] = 1. / np.sqrt(np.abs(hdiag[ok]))
        return s

    def _normal(self, entropy):
        """One N(0, 0.1^2) entry per theta coordinate in §7 order, v included (size + 1)."""
        rng = np.random.default_rng(np.random.SeedSequence(list(entropy)))
        return rng.normal(0.0, JITTER_SD, self.size + 1)

    def _jitter(self, base, scales, free, entropy, bounds):
        e = self._normal(entropy)[:len(base)]
        out = np.array(base, float, copy=True)
        out[free] = out[free] + scales[free] * e[free]
        return clip(out, bounds)

    @staticmethod
    def _summary(att):
        keys = ('loglik', 'converged', 'max_projected_gradient', 'hessian_negative_definite',
                'polish_iterations', 'polish_stop', 'evaluation_error', 'item8', 'v')
        out = {k: att.get(k) for k in keys if k in att}
        out['parameters'] = np.asarray(att['parameters'], float).tolist()
        return out

    def _run_attempts(self, start, bounds, attempt, hdiag_at_start):
        start = clip(np.asarray(start, float), bounds)
        free = free_mask(bounds)
        attempts = []
        scales = None
        x = start
        for a in range(MAX_RESTARTS + 1):
            if a > 0:
                finite = [t for t in attempts if np.isfinite(t['loglik'])]
                base = (np.asarray(max(finite, key=lambda t: t['loglik'])['point'], float)
                        if finite else start)
                if scales is None:
                    scales = self._scales(hdiag_at_start(start), free)
                x = self._jitter(base, scales, free, RESTART_ENTROPY + (a,), bounds)
                self._count('restart')
            res = attempt(x, bounds)
            res['attempt'] = a
            attempts.append(res)
            if res['converged']:
                break
        conv = [t for t in attempts if t['converged']]
        if conv:
            best = conv[0]
        else:
            finite = [t for t in attempts if np.isfinite(t['loglik'])]
            best = max(finite, key=lambda t: t['loglik']) if finite else attempts[-1]
        out = {k: v for k, v in best.items() if k != 'point'}
        out['restarts_tried'] = len(attempts)
        out['attempts'] = [self._summary(t) for t in attempts]
        out['unconverged_attempt_logliks'] = [t['loglik'] for t in attempts if not t['converged']]
        return out

    # ----------------------------------------------------- section 3 fits
    def _attempt3(self, x, bounds):
        def objective(q):
            ll, g, _ = self.evaluate(q, need_hessian=False)
            return -ll, -g
        p = np.array(x, float, copy=True)
        stop, iters = None, 0
        try:
            p = self._lbfgs(objective, p, bounds, 1e-15, 1000)
            p, _, iters, stop = newton_polish(lambda q, hess: self.evaluate(q, need_hessian=hess), p, bounds,
                                              stats=self.stats)
            res = self._certify3(p, bounds)
            if not np.isfinite(res['loglik']):
                res['converged'] = False
        except EVAL_ERRORS as exc:
            self._count('section3_attempt_eval_error')
            res = {'parameters': p, 'loglik': float('-inf'), 'converged': False,
                   'max_projected_gradient': float('nan'), 'hessian_negative_definite': False,
                   'evaluation_error': repr(exc)}
        res.update(point=np.asarray(res['parameters'], float), polish_iterations=iters, polish_stop=stop)
        return res

    def _hdiag3(self, start):
        try:
            return np.diag(self.evaluate(start)[2])
        except EVAL_ERRORS:
            return None

    def fit_single(self, start):
        start = np.asarray(start, float)
        bounds = self._section3_bounds(start)
        return self._run_attempts(start, bounds, self._attempt3, self._hdiag3)

    def fit_observed(self, entropy):
        """v1.8.3 item 3: the observed section-3 fit from two starts."""
        c0 = time.process_time()
        x0 = clip(self.cold_start(), self.bounds())
        bounds = self._section3_bounds(x0)
        free = free_mask(bounds)
        scales = self._scales(self._hdiag3(x0), free)
        x1 = self._jitter(x0, scales, free, entropy, bounds)
        fits = [self.fit_single(x0), self.fit_single(x1)]
        l0, l1 = fits[0]['loglik'], fits[1]['loglik']
        conv = [k for k in (0, 1) if fits[k]['converged'] and np.isfinite(fits[k]['loglik'])]
        failed, case, cand = False, None, None
        if not conv:
            failed, case, rule = True, 'no converged start', 'no converged start: observed fit fails'
        elif len(conv) == 2:
            cand = 1 if l1 > l0 + OBS_TIE else 0
            rule = 'both converged: start 1 above start 0 + 1e-9' if cand else \
                'both converged: start 0 (higher or within 1e-9)'
        else:
            cand = conv[0]
            rule = 'only start %d converged' % cand
        U = [l for f in fits for l in f['unconverged_attempt_logliks'] if np.isfinite(l)]
        u_gap = None
        if cand is not None:
            lc = fits[cand]['loglik']
            u_gap = (max(U) - lc) if U else None
            if U and max(U) > lc + OBS_FAIL_MARGIN:
                failed, case = True, 'unconverged end point more than 1e-6 above the candidate'
                rule += '; then failed by step 3'
        pick = cand if cand is not None else max((0, 1), key=lambda k: fits[k]['loglik']
                                                 if np.isfinite(fits[k]['loglik']) else -np.inf)
        p0 = np.asarray(fits[0]['parameters'], float)
        record = []
        for k, f in enumerate(fits):
            record.append({'start': k, 'converged': bool(f['converged']), 'attempts': f['restarts_tried'],
                           'loglik': f['loglik'], 'loglik_minus_start0': f['loglik'] - l0,
                           'theta_bitwise_equal_start0': bool(np.array_equal(np.asarray(f['parameters'], float), p0)),
                           'attempt_records': f['attempts']})
        out = dict(fits[pick])
        out.update(converged=not failed, observed_fit_failed=failed, failure_case=case,
                   reported_start=pick, selection_rule=rule,
                   max_unconverged_minus_candidate=u_gap,
                   starts=record, start1_point=x1.tolist(), start0_point=x0.tolist(),
                   start1_entropy=list(entropy), cpu_seconds_both_starts=time.process_time() - c0,
                   kept_start='observed start %d' % pick, start_rule=rule)
        return out

    # ----------------------------------------------------- section 4 fits
    def _section4_funs(self):
        def objective(q):
            ll, g = self.marginal(q[:-1], q[-1])
            return -ll, -g

        def fun(q, hess):
            if hess:
                return self.marginal(q[:-1], q[-1], with_hessian=True)
            ll, g = self.marginal(q[:-1], q[-1])
            return ll, g, None
        return objective, fun

    def _attempt4(self, x, bounds):
        objective, fun = self._section4_funs()
        self._held = None
        p = np.array(x, float, copy=True)
        stop, iters, note = None, 0, None
        try:
            p = self._lbfgs(objective, p, bounds, 1e-13, 500)
            p, _, iters, stop = newton_polish(fun, p, bounds, stats=self.stats)
            if p[-1] > 0:
                fun(p, True)                         # integrators of the last polish iterate at p
            last = dict(self.last_choices)
            res = self._certify4(p, bounds, final=True)
            if p[-1] > 0 and res['integrators'] != last:
                self._count('item8_repolish')
                held = dict(res['integrators'])
                self._held = held
                p, _, it2, stop2 = newton_polish(fun, p, bounds, stats=self.stats)
                iters += it2
                stop = stop + '; re-polish ' + stop2
                self._held = None
                res = self._certify4(p, bounds, final=True)
                if res['integrators'] != held:
                    self._count('item8_failed')
                    res['converged'] = False
                    note = 'gene changed integrator again at the re-polished final point'
                else:
                    note = 're-polished with held integrators'
            if not np.isfinite(res['loglik']):
                res['converged'] = False
        except EVAL_ERRORS as exc:
            self._held = None
            self._count('section4_attempt_eval_error')
            res = {'parameters': p[:-1], 'v': float(p[-1]), 'loglik': float('-inf'), 'converged': False,
                   'gradient': np.full(len(p), np.nan), 'max_projected_gradient': float('nan'),
                   'hessian_negative_definite': False, 'integrators': dict(self.last_choices),
                   'v_on_upper_bound': bool(abs(p[-1] - V18.V_UPPER) <= 1e-8), 'evaluation_error': repr(exc)}
        res.update(point=np.r_[np.asarray(res['parameters'], float), res['v']], polish_iterations=iters,
                   polish_stop=stop, item8=note)
        return res

    def _hdiag4(self, start):
        saved = (dict(self._prev_mode), set(self._seen), self._held, dict(self.last_choices))
        try:
            return np.diag(self.marginal(start[:-1], start[-1], with_hessian=True)[2])
        except EVAL_ERRORS:
            return None
        finally:
            self._prev_mode, self._seen, self._held, self.last_choices = saved

    def fit_marginal_single(self, start, fixed_v):
        self._reset_fit_state()
        start = np.asarray(start, float)
        informative = self.informative()
        bounds = [b if informative[j] else (start[j], start[j]) for j, b in enumerate(self.bounds())]
        bounds += [(0., V18.V_UPPER) if fixed_v is None else (fixed_v, fixed_v)]
        return self._run_attempts(start, bounds, self._attempt4, self._hdiag4)

    def _reference_or_fail(self, fit):
        """Reference-integrator l at a section-4 end point; an evaluation failure gives -inf and
        makes that fit not converged (its l cannot be compared)."""
        if not np.isfinite(fit['loglik']):
            return float('-inf')
        try:
            return self.reference_loglik(np.asarray(fit['parameters'], float), fit['v'])
        except EVAL_ERRORS:
            self._count('reference_loglik_eval_error')
            fit['converged'] = False
            return float('-inf')

    def fit_marginal(self, tau_max, start=None, fixed_v=None, null_fit=None, profile_start=None,
                     section4_from_null=False):
        """REC v1.8 fit_marginal with evaluation failures contained (v1.8.3 item 7)."""
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
            f['loglik_reference'] = self._reference_or_fail(f)
        eps = max(1e-6, 1e-8 * self.genes_with_rows)
        kept, rule = select_start(fits, eps, key='loglik_reference')
        if kept is None:
            name = max(fits, key=lambda k: fits[k]['loglik_reference'])
            fit = {**fits[name], 'converged': False}
        else:
            fit = dict(fits[kept])
        label = kept
        if fixed_v is None and fit['v'] == 0. and np.isfinite(fit['gradient'][-1]) and \
                fit['gradient'][-1] > SCORE_TOL:
            restart = self.fit_marginal_single(np.r_[fit['parameters'], tau_max ** 2], None)
            restart['loglik_reference'] = self._reference_or_fail(restart)
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
                'start_v': {k: f['v'] for k, f in fits.items()},
                'start_attempts': {k: f['attempts'] for k, f in fits.items()}}


def observed_failure(fit, se):
    """v1.8.3 item 4 cases 1 and 2 for the observed section-3 fit."""
    if fit.get('observed_fit_failed'):
        return 'case 1: ' + str(fit.get('failure_case'))
    if not (se is not None and math.isfinite(se) and se > 0):
        return 'case 2: sandwich SE not finite and positive'
    return None
