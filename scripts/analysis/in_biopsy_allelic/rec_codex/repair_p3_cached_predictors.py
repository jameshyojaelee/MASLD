#!/usr/bin/env python3
"""Instance-scoped predictor reuse during one fixed-p synthetic integration."""
import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import time


@contextmanager
def cached_predictors(model, p):
    """Temporarily reuse fixed-p predictors, restoring the instance on failure.

    The caller must limit this context to one integration at the supplied p.
    Cached Jacobian/Hessian arrays are read-only by contract: callers must not
    mutate returned objects. No BB, selection, integration, or fit rules change.
    """
    if not os.environ.get('SLURM_JOB_ID', '').isdigit():
        raise RuntimeError('Compute allocation required')
    import numpy as np
    fixed = np.asarray(p).copy()
    if fixed.ndim != 1 or not np.isfinite(fixed).all():
        raise ValueError('Finite parameter vector required')
    original = model.predictors
    existed = 'predictors' in model.__dict__
    previous = model.__dict__.get('predictors')
    cache = {}

    def predictor(row, parameters, u=0, need_hessian=True):
        current = np.asarray(parameters)
        if current.shape != fixed.shape or not np.array_equal(current, fixed):
            raise ValueError('Predictor cache cannot span changed parameters')
        key = (id(row),bool(need_hessian))
        if key not in cache:
            eta0, dispersion, jacobian, eta_hessian = original(
                row, fixed, 0, need_hessian=need_hessian)
            # Retain the row object so object identity cannot be recycled.
            cache[key] = (row,eta0,dispersion,jacobian,eta_hessian)
        saved_row,eta0,dispersion,jacobian,eta_hessian = cache[key]
        if saved_row is not row:
            raise ValueError('Row identity changed inside predictor cache')
        return eta0+u*row.stage,dispersion,jacobian,eta_hessian

    model.predictors = predictor
    try:
        yield model
    finally:
        if existed:
            model.predictors = previous
        else:
            del model.predictors


def compare_cached_predictors(model, p, rows):
    """Compute-only comparison against unmodified bound methods at fixed points."""
    if not os.environ.get('SLURM_JOB_ID', '').isdigit():
        raise RuntimeError('Compute allocation required')
    import numpy as np
    rows = tuple(rows)
    if not rows:
        raise ValueError('Nonempty synthetic rows required')
    original_predictor = model.predictors
    original_evaluate = model.evaluate
    existed = 'predictors' in model.__dict__
    previous = model.__dict__.get('predictors')

    def restored():
        return (('predictors' in model.__dict__) == existed and
                (not existed or model.__dict__['predictors'] is previous) and
                model.predictors == original_predictor)

    def differences(reference, candidate):
        errors, equal = [], []
        for a,b in zip(reference,candidate):
            if a is None or b is None:
                errors.append(None)
                equal.append(a is None and b is None)
            else:
                aa,bb = np.asarray(a),np.asarray(b)
                if aa.shape != bb.shape:
                    raise ValueError('Compared predictor/evaluator shapes differ')
                errors.append(float(np.max(np.abs(aa-bb))))
                equal.append(bool(np.array_equal(aa,bb)))
        return dict(component_max_absolute_differences=errors,
                    component_exact_equal=equal,all_components_exact_equal=all(equal))

    predictor_cases, evaluate_cases = [], []
    for need_hessian in (False,True):
        points = (-.5,0.,.5)
        # Save independent uncached outputs before installing any override.
        reference_predictors = [[original_predictor(row,p,u,need_hessian=need_hessian)
                                 for row in rows] for u in points]
        began = time.perf_counter()
        reference_values = [original_evaluate(p,rows=rows,u=u,need_hessian=need_hessian)
                            for u in points]
        reference_seconds = time.perf_counter()-began
        with cached_predictors(model,p):
            for u,reference_rows in zip(points,reference_predictors):
                for ordinal,(row,reference) in enumerate(zip(rows,reference_rows)):
                    actual = model.predictors(row,p,u,need_hessian=need_hessian)
                    predictor_cases.append(dict(row_ordinal=ordinal,u=u,need_hessian=need_hessian,
                                                **differences(reference,actual)))
            began = time.perf_counter()
            cached_values = [model.evaluate(p,rows=rows,u=u,need_hessian=need_hessian)
                             for u in points]
            cached_seconds = time.perf_counter()-began
        if not restored():
            raise ValueError('Predictor override was not restored after comparison')
        for u,reference,actual in zip(points,reference_values,cached_values):
            evaluate_cases.append(dict(u=u,need_hessian=need_hessian,
                                       **differences(reference,actual)))
        evaluate_cases[-1].update(three_point_uncached_seconds=reference_seconds,
                                 three_point_cached_seconds=cached_seconds,
                                 timing_excludes_cache_population=True)

    changed = np.asarray(p).copy()
    changed[0] += 1.
    if np.array_equal(changed,np.asarray(p)):
        raise ValueError('Changed-parameter refusal probe did not change the vector')
    rejected = False
    with cached_predictors(model,p):
        try:
            model.predictors(rows[0],changed,0.)
        except ValueError as error:
            rejected = str(error) == 'Predictor cache cannot span changed parameters'
    changed_restored = restored()

    class IntentionalProbe(Exception):
        pass

    try:
        with cached_predictors(model,p):
            model.predictors(rows[0],p,0.)
            raise IntentionalProbe('Restoration probe')
    except IntentionalProbe:
        pass
    exception_restored = restored()
    return dict(predictor_comparisons=predictor_cases,evaluate_comparisons=evaluate_cases,
                changed_parameters_rejected=rejected,
                restored_after_changed_parameter_probe=changed_restored,
                restored_after_intentional_exception=exception_restored,
                previous_instance_override_preserved=not existed or exception_restored,
                supplied_rows=len(rows),fixed_u_points=[-.5,0.,.5],
                no_model_fits=True,no_integration=True,no_calibration=True,
                single_run_timings_not_end_to_end_speed_claim=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true', required=True)
    parser.parse_args()
    print(json.dumps(dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        numerical_execution=False,data_reads=False,instance_only=True,
        cache_scope='one fixed-p integration',parameter_derivatives_retained=True,
        no_probability_or_quadrature_changes=True)))


if __name__ == '__main__':
    main()
