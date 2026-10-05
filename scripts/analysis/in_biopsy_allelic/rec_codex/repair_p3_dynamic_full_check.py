#!/usr/bin/env python3
"""Full synthetic zero-law scalar comparison, no fits or integration.

The immutable reconstruction supplies the same saved counts/parameters. Seven
declared parameter points exercise each shared dispersion slope and combined
predictor changes, then return exactly to the saved point. All 173 rows remain.
"""
import argparse
import ast
import hashlib
import importlib.util
import os
from pathlib import Path
import sys
import time
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
SOURCE = REC/'repair_p3_fixed_point_failed0gene2_22039132/executed_sources/repair_p3_fixed_point_check.py'
SOURCE_SHA = 'd3f529ad20666c58dc8a84d01785d3f62c42f18a430a3a572a055c180b9b4a9e'
KERNEL_SHA = '461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'
DYNAMIC = Path(__file__).with_name('repair_p3_dynamic_conditional.py')
DYNAMIC_SHA = '57ff8c707f17f5cea2a1d7fd8fd76e2e5b57d666cf0512f17e01d8fc66790cea'
V2 = Path(__file__).with_name('repair_p3_prepared_conditional_v2.py')
V2_SHA = '289e93b651411042093a2b987717d87cdb5ed156ab7c736e1ac5cf820670fd4b'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def compare_cached_predictors(model, p, unused_gene_rows):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Compute allocation required')
    import numpy as np
    from dataclasses import replace
    local_source = Path(__file__).with_name(DYNAMIC.name)
    require(sha(local_source)==DYNAMIC_SHA and sha(local_source.with_name(V2.name))==V2_SHA,
            'Archived dynamic/v2 source differs')
    spec = importlib.util.spec_from_file_location('_repair_p3_full_private_dynamic',local_source)
    dynamic = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(dynamic)
    rows = tuple(model.rows)
    require(len(rows)==173 and model.G==8 and model.Z==7 and model.size==28,
            'Original small synthetic shape differs')
    seed = np.asarray(p,dtype=float).copy()
    require(seed.shape==(28,) and np.isfinite(seed).all(),'Finite saved full parameter vector required')
    counts_sha = hashlib.sha256(np.asarray([r.a for r in rows],dtype='<i8').tobytes()).hexdigest()
    row_signatures = tuple(tuple(vars(r).items()) for r in rows)
    methods = {name:getattr(type(model),name) for name in dynamic.INHERITED}
    setup_started = time.perf_counter()
    context = dynamic.DynamicConditional(model,seed,rows)
    setup_seconds = time.perf_counter()-setup_started
    original = context.original.conditional_terms
    original_bb = context.original.bb_terms
    bounds = model.bounds()
    def add(q,index,amount):
        lo,hi = bounds[index]
        q[index] = np.clip(q[index]+amount,-np.inf if lo is None else lo,
                           np.inf if hi is None else hi)
    points = [('saved_seed',seed.copy())]
    for index,label in ((-3,'rho_stage'),(-2,'rho_background'),(-1,'rho_duplication')):
        q = seed.copy()
        add(q,index,.1)
        points.append((label,q))
    predictor = seed.copy()
    for j in range(model.G): add(predictor,j,.02*(j+1))
    add(predictor,2*model.G,.05)
    for j in range(model.Z): add(predictor,2*model.G+1+j,.01*(j+1))
    points.append(('alpha_kappa_all_delta',predictor))
    combined = predictor.copy()
    for index,amount in ((-3,.1),(-2,-.1),(-1,.1)): add(combined,index,amount)
    points.extend((('predictors_and_all_rho',combined),('exact_return_saved_seed',seed.copy())))
    def retained(value):
        if isinstance(value,np.ndarray): return retained(value.tolist())
        if isinstance(value,(tuple,list)): return [retained(v) for v in value]
        if isinstance(value,(float,np.floating)) and not np.isfinite(value): return str(value)
        return value
    records = []
    for need_hessian in (False,True):
        for label,q in points:
            ref,actual,errors,error = None,None,None,None
            ref_seconds,dynamic_seconds = None,None
            ll_within,g_within,exact = False,False,None
            try:
                require(context.original.conditional_terms is original,'Reference conditional not restored')
                started = time.perf_counter()
                ref = model.evaluate(q,need_hessian=need_hessian)
                ref_seconds = time.perf_counter()-started
                previous = dynamic.install_private(context)
                try:
                    started = time.perf_counter()
                    actual = model.evaluate(q,need_hessian=need_hessian)
                    dynamic_seconds = time.perf_counter()-started
                finally:
                    context.original.conditional_terms = previous
                require(all(np.isfinite(v).all() for v in (*ref,*actual) if v is not None),
                        'Nonfinite zero-law full comparison')
                errors = [None if a is None else float(np.max(np.abs(np.asarray(a)-b)))
                          for a,b in zip(ref,actual)]
                exact = [a is b if a is None else bool(np.array_equal(a,b)) for a,b in zip(ref,actual)]
                ll_within = abs(ref[0]-actual[0])<=1e-6
                g_within = bool(np.allclose(ref[1],actual[1],atol=1e-8,rtol=1e-6))
            except Exception:
                error = traceback.format_exc()
            records.append(dict(case=label,need_hessian=need_hessian,parameters=q.copy(),
                original=retained(ref),dynamic=retained(actual),max_absolute_differences=errors,
                exact_equal=exact,exception=error,LL_within_existing_allowance=ll_within,
                gradient_within_existing_allowance=g_within,original_seconds=ref_seconds,
                warm_dynamic_seconds=dynamic_seconds,
                original_finite_masks=None if ref is None else
                    [None if v is None else np.isfinite(v) for v in ref],
                dynamic_finite_masks=None if actual is None else
                    [None if v is None else np.isfinite(v) for v in actual]))
    row = rows[0]
    eta,r = model.predictors(row,seed,0.,need_hessian=False)[:2]
    guards = []
    def refusal(label,test_row,test_eta,test_r):
        rejected,error = False,None
        try: context.terms(test_row,test_eta,test_r)
        except ValueError as exc: rejected,error = True,str(exc)
        except Exception: error = traceback.format_exc()
        guards.append(dict(case=label,rejected=rejected,error=error))
    refusal('foreign_equal_row',replace(row),eta,r)
    omega = row.omega
    try:
        object.__setattr__(row,'omega',omega+.125)
        refusal('same_row_mutation',row,eta,r)
    finally:
        object.__setattr__(row,'omega',omega)
    for value in (np.nan,np.inf,-np.inf):
        refusal('nonfinite_eta_'+str(value),row,value,r)
        refusal('nonfinite_r_'+str(value),row,eta,value)
    refusal('nonscalar_eta',row,np.array([eta]),r)
    refusal('nonscalar_r',row,eta,np.array([r]))
    valid_after,error = False,None
    try:
        answer = context.terms(row,eta,r)
        valid_after = all(np.isfinite(v).all() for v in answer if v is not None)
    except Exception:
        error = traceback.format_exc()
    preserved = (context.original.conditional_terms is original and context.original.bb_terms is original_bb
                 and all(getattr(type(model),name) is fn for name,fn in methods.items()))
    final_counts_sha = hashlib.sha256(np.asarray([r.a for r in rows],dtype='<i8').tobytes()).hexdigest()
    return dict(full_zero_law_comparisons=records,comparison_records_completed=True,
        comparisons=len(records),parameter_points=len(points),hessian_modes=[False,True],
        exception_count=sum(r['exception'] is not None for r in records),
        all_LL_within_existing_allowance=all(r['LL_within_existing_allowance'] for r in records),
        all_gradients_within_existing_allowance=all(r['gradient_within_existing_allowance'] for r in records),
        existing_fixed_vector_LL_allowance=1e-6,existing_fixed_vector_gradient_atol=1e-8,
        existing_fixed_vector_gradient_rtol=1e-6,hessian_differences_descriptive_no_new_threshold=True,
        setup_seconds=setup_seconds,single_run_not_end_to_end_fit_speed=True,
        registered_rows=173,genes=8,score_coordinates=28,hessian_shape=[28,28],
        count_vector_sha256=counts_sha,final_count_vector_sha256=final_counts_sha,
        counts_unchanged=counts_sha==final_counts_sha,
        row_signatures_unchanged=row_signatures==tuple(tuple(vars(r).items()) for r in rows),
        seed_parameters_unchanged=bool(np.array_equal(seed,p)),
        final_point_exact_saved_seed=bool(np.array_equal(points[-1][1],seed)),
        original_callables_and_inherited_method_identities_unchanged=preserved,
        guard_checks=guards,all_required_refusals_observed=all(r['rejected'] for r in guards),
        valid_context_after_refusal=valid_after,valid_context_after_refusal_exception=error,
        memory_estimate=context.memory_estimate,dispersion_cache_ledger=context.ledger(),
        observed_n_bounds=dict(min=min(r.n for r in rows),max=max(r.n for r in rows)),
        parameter_changes='Each rho separately+.1; all alpha+.02*(gene+1), kappa+.05, delta+.01*(coordinate+1); combined rho(+.1,-.1,+.1); clip to original bounds; exact return',
        no_new_counts=True,no_draw_extension=True,no_integration=True,no_model_fit=True,no_calibration=True,
        completion_means='All comparison records retained; arithmetic comparison only, not model validation')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',choices=('failed0gene2','positive18'),required=True)
    parser.add_argument('--expected-checker-sha',required=True)
    parser.add_argument('--inspect',action='store_true')
    args = parser.parse_args()
    own = Path(__file__).resolve()
    guards = {SOURCE:SOURCE_SHA,own:args.expected_checker_sha,DYNAMIC:DYNAMIC_SHA,V2:V2_SHA}
    for path,digest in guards.items():
        require(sha(path)==digest,'Reviewed full-check dependency differs: '+str(path))
        ast.parse(path.read_text())
    source = SOURCE.read_text()
    labels = {'repair_p3_fixed_point_':'repair_p3_dynamic_full_',
              "save('predictor_cache_comparison', report)":"save('dynamic_full_zero_law_comparison', report)",
              'predictor_comparison_completed=True':'dynamic_full_zero_law_comparison_completed=True'}
    for old,new in labels.items():
        require(source.count(old)==1,'Private label occurrence differs: '+old)
        source = source.replace(old,new)
    namespace = {'__file__':str(SOURCE),'__name__':'private_dynamic_full_reconstruction'}
    exec(compile(ast.parse(source),str(SOURCE),'exec'),namespace)
    namespace['PREDICTOR_CACHE'] = own
    namespace['GUARDS'].update(guards)
    sys.argv = [str(SOURCE),'--case',args.case,'--predictor-points-sha',args.expected_checker_sha,
                '--expected-helper-sha',SOURCE_SHA,'--expected-kernel-sha',KERNEL_SHA]
    if args.inspect: sys.argv.append('--inspect')
    namespace['main']()


if __name__=='__main__':
    main()
