"""REC Model A under spec v1.8.4 (addendum sha256 eaa1a3b8...0978, 2026-10-18). Synthetic fixtures only until Gate B.

Built from the addendum text and the in-force spec (sha256 285c254c...883a9b). No B-MODEL code, output or
README was read. Subclass of v183/likelihood_v183.ModelV183 (unchanged file); the likelihood, the section-4
integral and its integrators, the polish rule, the attempts within one fit (attempt 0 + up to 3 restarts,
REC restart jitter [20260923, 90, a]), the sandwich, bounds and the §7 start are inherited unchanged.

This module adds:
  item 2   K = 3 starts for every section-3 refit: warm (reported observed MLE), cold (§7 start on the
           refit's data), probe start j = 1; holds alpha_g = 0, r_g = logit(0.02) for genes with no row,
           delta_m = 0 for all-zero z columns, applied to every start before the first attempt.
  item 3   selection among the K starts (eps_tie,3 = 1e-6; ties to the lower index; an unconverged start
           more than 1e-6 above the kept one fails the refit).
  item 4   probe start: §7 start, then 2G + 5 uniforms in the fixed order, holds, clip.
  item 5   the observed section-3 fit from three starts (start 2 = probe on stream 10, t = 0); selection
           5.2 (1e-9 tie among converged starts, 1e-6 unconverged-higher margin over every attempt).
  item 6   the section-4 pair for P3 (member 0 at v = 0, member 1 at v = tau_max^2, both from the kept
           section-3 theta), the K = 2 table with eps_tie,4, l3*, and LR by item 6.4. No third fit.
  item 10  records per start and per member, with full end points for the §11 cross-evaluation.
The P2 tau profile fits at tau0 keep v1.7 rule 5 (ModelV183.fit_marginal with fixed_v, unchanged).
"""
import time

import numpy as np

from likelihood_v18 import clip
from likelihood_v183 import EVAL_ERRORS, ModelV183

EPS_TIE3 = 1e-6          # item 1, section-3 refit selection and failure
OBS_TIE = 1e-9           # item 1, eps_obs
OBS_FAIL_MARGIN = 1e-6   # item 1, observed unconverged-higher margin
V_ZERO = 1e-8            # item 1, v lower-bound test
MISSED_V0 = 1e-6         # item 1, "section-3 refit missed a v = 0 point"
LOGIT02 = float(np.log(.02 / .98))
START_NAMES = ('warm', 'cold', 'probe1')


def finite(x):
    return x is not None and bool(np.isfinite(x))


def select_k(fits, eps=EPS_TIE3):
    """Item 3. fits: list of dicts with 'loglik', 'converged'. Returns (kept index or None, rule)."""
    C = [s for s, f in enumerate(fits) if f['converged'] and finite(f['loglik'])]
    if not C:
        return None, 'failed, none converged'
    lmax = max(fits[s]['loglik'] for s in C)
    kept = next(s for s in C if lmax - fits[s]['loglik'] < eps)
    lk = fits[kept]['loglik']
    for s, f in enumerate(fits):
        if s not in C and finite(f['loglik']) and f['loglik'] - lk > eps:
            return None, 'failed, unconverged higher (start %d)' % s
    if len(C) == 1:
        rule = 'only converged start'
    elif lk < lmax:
        rule = 'tie to lower index'
    elif any(s < kept for s in C):
        rule = 'higher'
    else:
        rule = 'higher or tie to lower index (kept is the lowest converged index and the maximum)'
    return kept, rule


def select_observed(fits):
    """Item 5.2. Returns (candidate or None, failed, case, rule, max unconverged minus candidate)."""
    C = [s for s, f in enumerate(fits) if f['converged'] and finite(f['loglik'])]
    if not C:
        return None, True, 'no converged start', 'no converged start: observed fit fails', None
    lmax = max(fits[s]['loglik'] for s in C)
    cand = next(s for s in C if lmax - fits[s]['loglik'] <= OBS_TIE)
    rule = 'candidate start %d (lowest converged index within 1e-9 of the maximum %s)' % (
        cand, 'over starts ' + ','.join(map(str, C)))
    U = [l for f in fits for l in f['unconverged_attempt_logliks'] if finite(l)]
    lc = fits[cand]['loglik']
    gap = (max(U) - lc) if U else None
    if U and max(U) - lc > OBS_FAIL_MARGIN:
        return cand, True, 'unconverged end point more than 1e-6 above the candidate', rule + '; failed by step 3', gap
    return cand, False, None, rule, gap


class ModelV184(ModelV183):

    def cold_start(self):
        if getattr(self, '_cold_cache', None) is None:
            self._cold_cache = np.asarray(super().cold_start(), float)
        return self._cold_cache.copy()

    # ------------------------------------------------------------ item 2 holds, item 4 probe
    def apply_holds(self, x):
        x = np.array(x, float, copy=True)
        present = {r.gene for r in self.rows}
        for g in range(self.G):
            if g not in present:
                x[g] = 0.
                x[self.G + g] = LOGIT02
        for j in range(self.Z):
            if all(r.z[j] == 0 for r in self.rows):
                x[2 * self.G + 1 + j] = 0.
        return x

    def probe_start(self, entropy):
        G, k = self.G, 2 * self.G
        x = self.cold_start()
        rng = np.random.default_rng(np.random.SeedSequence(list(entropy)))
        x[:G] = rng.uniform(-1.5, 1.5, G)
        x[G:2 * G] = rng.uniform(-6.0, -2.0, G)
        x[k] = rng.uniform(-0.5, 0.5)
        rho_s, rho_b, rho_d, omega_s = rng.uniform(-1.0, 1.0, 4)
        x[self.size - 3], x[self.size - 2], x[self.size - 1] = rho_s, rho_b, rho_d
        x[k + 1 + self.Z] = omega_s
        return clip(self.apply_holds(x), self.bounds())

    def _timed_fit3(self, x):
        t0 = time.monotonic()
        f = self.fit_single(x)
        f['seconds'] = time.monotonic() - t0
        return f

    @staticmethod
    def _start_record(s, name, x, f, l_ref):
        return {'index': s, 'name': name, 'start_point': np.asarray(x, float).tolist(),
                'end_point': np.asarray(f['parameters'], float).tolist(), 'loglik': f['loglik'],
                'loglik_minus_warm': f['loglik'] - l_ref if finite(l_ref) and finite(f['loglik']) else None,
                'converged': bool(f['converged']), 'attempts': f['restarts_tried'],
                'unconverged_attempt_logliks': f['unconverged_attempt_logliks'],
                'seconds': f.get('seconds'), 'attempt_records': f['attempts']}

    # ------------------------------------------------------------ items 2-3 section-3 refit
    def fit_refit(self, warm_theta, probe_entropy):
        bounds = self.bounds()
        xs = [clip(self.apply_holds(np.asarray(warm_theta, float)), bounds),
              clip(self.apply_holds(self.cold_start()), bounds),
              self.probe_start(probe_entropy)]
        fits = [self._timed_fit3(x) for x in xs]
        kept, rule = select_k(fits, EPS_TIE3)
        record = [self._start_record(s, START_NAMES[s], xs[s], f, fits[0]['loglik']) for s, f in enumerate(fits)]
        if kept is None:
            fin = [s for s, f in enumerate(fits) if finite(f['loglik'])]
            pick = max(fin, key=lambda s: fits[s]['loglik']) if fin else 0
            out = {**fits[pick], 'converged': False}
        else:
            out = dict(fits[kept])
            out['converged'] = True
        out.update(kept_index=kept, kept_start=None if kept is None else START_NAMES[kept], start_rule=rule,
                   refit_failed=kept is None, starts=record, probe_entropy=list(probe_entropy),
                   start_logliks={START_NAMES[s]: f['loglik'] for s, f in enumerate(fits)},
                   start_converged={START_NAMES[s]: bool(f['converged']) for s, f in enumerate(fits)})
        return out

    # ------------------------------------------------------------ item 5 observed fit
    def fit_observed3(self, start1_entropy, probe_entropy):
        c0 = time.monotonic()
        x0 = clip(self.apply_holds(clip(self.cold_start(), self.bounds())), self.bounds())
        b0 = self._section3_bounds(x0)
        free = np.array([not (lo is not None and hi is not None and lo == hi) for lo, hi in b0])
        scales = self._scales(self._hdiag3(x0), free)
        x1 = self._jitter(x0, scales, free, start1_entropy, b0)          # v1.8.3 item 3.2, unchanged
        x2 = self.probe_start(probe_entropy)
        xs = [x0, x1, x2]
        fits = [self._timed_fit3(x) for x in xs]
        cand, failed, case, rule, gap = select_observed(fits)
        if cand is not None:
            pick = cand
        else:
            fin = [s for s, f in enumerate(fits) if finite(f['loglik'])]
            pick = max(fin, key=lambda s: fits[s]['loglik']) if fin else 0
        p0 = np.asarray(fits[0]['parameters'], float)
        record = []
        for s, f in enumerate(fits):
            r = self._start_record(s, 'observed start %d' % s, xs[s], f, fits[0]['loglik'])
            r['loglik_minus_start0'] = r.pop('loglik_minus_warm')
            r['theta_bitwise_equal_start0'] = bool(np.array_equal(np.asarray(f['parameters'], float), p0))
            record.append(r)
        out = dict(fits[pick])
        out.update(converged=not failed, observed_fit_failed=failed, failure_case=case, reported_start=pick,
                   selection_rule=rule, max_unconverged_minus_candidate=gap, starts=record,
                   start1_entropy=list(start1_entropy), probe_entropy=list(probe_entropy),
                   seconds_all_starts=time.monotonic() - c0, kept_start='observed start %d' % pick,
                   start_rule=rule)
        return out

    # ------------------------------------------------------------ item 6 section-4 pair
    def section4_pair(self, theta3, l3, v_start):
        theta3 = np.asarray(theta3, float)
        members = []
        for m, v0 in enumerate((0., float(v_start))):
            t0 = time.monotonic()
            err = None
            try:
                f = self.fit_marginal_single(np.r_[theta3, v0], None)
                f['loglik_reference'] = self._reference_or_fail(f)
            except EVAL_ERRORS as exc:
                err = repr(exc)
                f = {'parameters': theta3, 'v': v0, 'loglik': float('-inf'), 'loglik_reference': float('-inf'),
                     'converged': False, 'restarts_tried': 0, 'attempts': [], 'unconverged_attempt_logliks': []}
            f['seconds'] = time.monotonic() - t0
            f['member_error'] = err
            members.append(f)
        eps4 = max(1e-6, 1e-8 * self.genes_with_rows)
        # item 6.2: the K = 2 table with eps_tie,4 on reference-integrator l, differences computed as doubles
        kept, rule = select_k([{'loglik': f['loglik_reference'], 'converged': f['converged']} for f in members], eps4)
        l3_star = l3
        member_l3 = []
        for f in members:
            val = None
            if f['converged'] and f['v'] <= V_ZERO:
                try:
                    val = self.reference_loglik(np.asarray(f['parameters'], float), 0.)
                except EVAL_ERRORS:
                    val = None
                if finite(val):
                    l3_star = max(l3_star, val)
            member_l3.append(val)
        if kept is None:
            lr, l4, v_kept = float('nan'), None, None
        else:
            l4 = members[kept]['loglik_reference']
            v_kept = members[kept]['v']
            lr = 0. if v_kept <= V_ZERO else 2 * (max(l4, l3_star) - l3_star)
        rec = []
        for m, f in enumerate(members):
            rec.append({'member': m, 'v_start': (0., float(v_start))[m],
                        'end_point': np.asarray(f['parameters'], float).tolist(), 'v': f['v'],
                        'loglik': f['loglik'], 'loglik_reference': f['loglik_reference'],
                        'loglik_section3_at_v0': member_l3[m], 'converged': bool(f['converged']),
                        'attempts': f.get('restarts_tried'), 'seconds': f['seconds'],
                        'error': f['member_error'], 'attempt_records': f.get('attempts'),
                        'unconverged_attempt_logliks': f.get('unconverged_attempt_logliks')})
        return {'failed': kept is None, 'kept_member': kept, 'rule': rule, 'eps_tie4': eps4,
                'v_start': float(v_start), 'l3': l3, 'l3_star': l3_star, 'l3_star_minus_l3': l3_star - l3,
                'missed_v0_point': bool(l3_star - l3 > MISSED_V0), 'l4_reference': l4, 'v_kept': v_kept,
                'LR': lr, 'members': rec}
