"""One synthetic observed fit plus cold-start invariance; no P1 draws or test."""
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
import scipy

import likelihood_v1 as likelihood
from synthetic_p1_preascertainment import generate


EXPECTED_LIKELIHOOD = '41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6'
EXPECTED_GENERATOR = 'afba32e07526e6785878bef2021d055efffb871ab5a6636be6844151312760a8'


def serial(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Compute allocation required.')
    here = Path(__file__).resolve().parent
    hashes = {name:hashlib.sha256((here/name).read_bytes()).hexdigest()
              for name in ('likelihood_v1.py', 'synthetic_p1_preascertainment.py',
                           'design_v1.py', 'profile_synthetic_p1_fit.py')}
    if hashes['likelihood_v1.py'] != EXPECTED_LIKELIHOOD or hashes['synthetic_p1_preascertainment.py'] != EXPECTED_GENERATOR:
        raise RuntimeError('Reviewed numerical reference or generator changed.')
    out = Path(os.environ['CODEX_REC_OUTPUT']) / ('synthetic_p1_fit_profile_'+os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)

    def save(name, value):
        (out/(name+'.json')).write_text(json.dumps(value, default=serial, indent=2)+'\n')
        print(name+' saved', flush=True)

    started = time.monotonic()
    sample = generate('nuisance_on')
    design = sample.design
    if len(design.participants) != 180 or len(design.rows) != 498:
        raise RuntimeError('Verified fixed synthetic population changed.')
    model = design.model()
    metrics = {}
    phase = ['cold_invariance_original']

    def record(label, elapsed, extra=None):
        entry = metrics.setdefault(phase[0]+'/'+label, {'calls':0, 'seconds':0.})
        entry['calls'] += 1
        entry['seconds'] += elapsed
        if extra:
            for key, value in extra.items():
                entry[key] = entry.get(key, 0)+value

    def instrument(instance):
        for name in ('cold_start', 'fit_single', 'sandwich'):
            original = getattr(instance, name)
            def timed(*args, _original=original, _name=name, **kwargs):
                tick = time.monotonic()
                try:
                    return _original(*args, **kwargs)
                finally:
                    record(_name, time.monotonic()-tick)
            setattr(instance, name, timed)
        original_evaluate = instance.evaluate
        def evaluate(*args, **kwargs):
            with_hessian = kwargs.get('need_hessian', args[3] if len(args) > 3 else True)
            label = 'evaluate_full_hessian' if with_hessian else 'evaluate_score_only'
            tick = time.monotonic()
            try:
                return original_evaluate(*args, **kwargs)
            finally:
                record(label, time.monotonic()-tick)
        instance.evaluate = evaluate

    original_minimize = likelihood.minimize
    def minimize(*args, **kwargs):
        tick = time.monotonic()
        result = original_minimize(*args, **kwargs)
        record('scipy_minimize', time.monotonic()-tick,
               {key:int(getattr(result, key, 0)) for key in ('nit', 'nfev', 'njev')})
        return result
    likelihood.minimize = minimize
    try:
        instrument(model)
        tick = time.monotonic()
        cold = model.cold_start()
        cold_original_seconds = time.monotonic()-tick
        # Fixed within-cohort rotation preserves counts and means. It is an
        # arithmetic invariant check, not a production stage draw or a CRT.
        changed_stage = design.stage.copy()
        for cohort in sorted(set(design.cohort)):
            mask = design.cohort == cohort
            changed_stage[mask] = np.roll(changed_stage[mask], 1)
        if np.array_equal(changed_stage, design.stage):
            raise RuntimeError('Fixed stage rotation did not change stages.')
        changed = design.with_stage(changed_stage).model()
        instrument(changed)
        phase[0] = 'cold_invariance_rotated'
        tick = time.monotonic()
        rotated_cold = changed.cold_start()
        cold_rotated_seconds = time.monotonic()-tick
        if not np.array_equal(cold, rotated_cold):
            save('cold_invariance_failed', {'original':cold, 'rotated':rotated_cold,
                                           'max_absolute_difference':np.max(np.abs(cold-rotated_cold))})
            raise RuntimeError('Cold vectors differ after changing only stages.')
        save('cold_invariance', {'array_equal':True, 'original':cold, 'rotated':rotated_cold,
                                'original_seconds':cold_original_seconds,
                                'rotated_seconds':cold_rotated_seconds,
                                'stage_change':'one-place rotation within each complete cohort',
                                'production_stage_law_used':False, 'fit_performed':False})
        phase[0] = 'observed_full_warm_cold_fit'
        # Generic numerical initialization, never the generating parameters.
        # Explicit start invokes the prescribed warm AND cold fitting routes.
        tick = time.monotonic()
        fit = model.fit(model.initial())
        fit_seconds = time.monotonic()-tick
        save('observed_fit', {'fit':fit, 'elapsed_seconds':fit_seconds,
                              'initialization':'model.initial(); generating truth not supplied',
                              'hashes':hashes})
        if fit['converged']:
            phase[0] = 'observed_sandwich'
            tick = time.monotonic()
            se = model.sandwich(fit['parameters'], design.clusters)
            sandwich_seconds = time.monotonic()-tick
            save('observed_sandwich', {'se':se, 'finite_positive':bool(np.isfinite(se) and se > 0),
                                      'elapsed_seconds':sandwich_seconds})
        else:
            se = None
            sandwich_seconds = None
        save('summary', {'status':'bounded_observed_fit_profile_complete', 'hashes':hashes,
                         'participants':180, 'genes':8, 'selected_rows':498,
                         'regime':'nuisance_on', 'fit_converged':bool(fit['converged']),
                         'cold_stage_invariance':True, 'fit_seconds':fit_seconds,
                         'sandwich_seconds':sandwich_seconds, 'sandwich_se':se,
                         'elapsed_seconds':time.monotonic()-started, 'timings':metrics,
                         'timing_rule':'nested inclusive timings; do not sum overlapping categories',
                         'python':sys.version, 'numpy':np.__version__, 'scipy':scipy.__version__,
                         'stage_refits_performed':0, 'P1_size_or_power_evaluated':False,
                         'real_outcomes_read':False,
                         'limits':'One fixed synthetic observed fit and deterministic cold-vector check. '
                                  'No caching, likelihood, optimizer or test-law changes; no calibration claim.'})
    except Exception as error:
        save('failure', {'exception_type':type(error).__name__, 'message':str(error),
                         'elapsed_seconds':time.monotonic()-started, 'timings':metrics,
                         'hashes':hashes, 'P1_size_or_power_evaluated':False})
        raise
    finally:
        likelihood.minimize = original_minimize


if __name__ == '__main__':
    main()
