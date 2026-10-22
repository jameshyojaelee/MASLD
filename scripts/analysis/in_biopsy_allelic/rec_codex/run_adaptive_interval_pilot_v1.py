#!/usr/bin/env python3
"""Separate synthetic adaptive-interval pilots with original parallel null starts."""
import argparse
import ast
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import sys
import time
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE = ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
NULL = REC/'p3_parallel_20261003T130323Z/baseline/observed_null.json'
FULL = REC/'adaptive_weighted_reference_full8_22085200'
FULL_REVIEW = REC/'adaptive_weighted_root_review_20261007T111822Z/review.json'
POLISH_REVIEW = REC/'original_first_warm_polish_root_review_20261007T093723Z/review.json'
PARALLEL_REVIEW = REC/'parallel_null_shifts_root_review_20261007T104754Z/review.json'
COUNTREF = {0: REC/'standard_hermite_pilot_draw0_22084966/draw_counts.json',
            11: REC/'standard_hermite_pilot_draw11_22084967/draw_counts.json',
            18: REC/'standard_hermite_parallel_pilot_draw18_22085062/draw_counts.json'}
SOURCES = ('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
           'check_bfix_tests_v13.py','repair_p3_local_moments.py','repair_p3_pathwise_moments_v1.py',
           'repair_p3_longdouble_conditional_v1.py','repair_p3_adaptive_intervals_weighted_v1.py',
           'repair_p3_adaptive_peak_scan_v1.py','repair_p3_adaptive_interval_model_v1.py','repair_p3_bb_gamma_difference_v3.py',
           'repair_p3_score_integral_polish.py','repair_p3_gamma_order_polish_model_v2.py',
           'repair_p3_parallel_null_shifts_v1.py')


def require(ok, message):
    if not ok: raise RuntimeError(message)


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--draw', type=int, choices=(0,11,18), required=True)
    ap.add_argument('--pins', type=Path, required=True)
    ap.add_argument('--expected-pins-sha', required=True)
    ap.add_argument('--expected-runner-sha', required=True)
    ap.add_argument('--inspect', action='store_true')
    a = ap.parse_args()
    require(sha(Path(__file__)) == a.expected_runner_sha and sha(a.pins) == a.expected_pins_sha, 'Runner/pins differ')
    pins = json.loads(a.pins.read_text())
    for path, digest in pins.items():
        p = Path(path)
        require(p.resolve().is_relative_to(ROOT.resolve()) and sha(p) == digest, 'Input differs')
        if p.suffix == '.py': ast.parse(p.read_text())
    required = [*(HERE/n for n in SOURCES), NULL, FULL_REVIEW, POLISH_REVIEW, PARALLEL_REVIEW, COUNTREF[a.draw],
        FULL/'protocol.json', FULL/'summary.json',
        *(FIXTURE/n for n in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv','draws/p3_uniforms.npy'))]
    require(all(str(p) in pins for p in required), 'Missing source/input pin')
    require(pins[str(HERE/'check_bfix_tests_v13.py')] == 'a4ca55a52d28dfdbfeb0fc2611c2dcb3ec412987c8f7618f59002cbf72b08c2d', 'Original P3 source differs')
    read = lambda p: json.loads(p.read_text())
    require(read(FULL_REVIEW)['all_four_fixed_points_passed'] and read(POLISH_REVIEW)['passed'] and read(PARALLEL_REVIEW)['passed'], 'Numerical prerequisites missing')
    source = (HERE/'check_bfix_tests_v13.py').read_text()
    loops = [n for n in ast.parse(source).body if isinstance(n,ast.For) and isinstance(n.target,ast.Tuple)
             and [getattr(x,'id',None) for x in n.target.elts] == ['b','draw']]
    require(len(loops) == 1, 'Original draw loop differs')
    body = loops[0].body
    function = ast.FunctionDef(name='original_draw',args=ast.arguments(posonlyargs=[],args=[ast.arg(arg='b'),ast.arg(arg='draw')],
        kwonlyargs=[],kw_defaults=[],defaults=[]),body=body,decorator_list=[])
    compiled_tree = ast.fix_missing_locations(ast.Module(body=[function],type_ignores=[]))
    require(ast.dump(compiled_tree.body[0].body[0],include_attributes=False) == ast.dump(body[0],include_attributes=False)
        and ast.dump(ast.Module(body=compiled_tree.body[0].body,type_ignores=[]),include_attributes=False)
        == ast.dump(ast.Module(body=body,type_ignores=[]),include_attributes=False), 'P3 body changed')
    if a.inspect:
        print(json.dumps(dict(draw=a.draw, original_draw_body_unchanged=True, original_first_null_polish=True, no_fifty_draw_run=True)))
        return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID')
        and os.environ.get('SLURM_CPUS_PER_TASK') == '8', 'Individual eight-CPU compute required')
    require(all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')), 'Threads differ')
    out = REC/f"adaptive_interval_pilot_draw{a.draw}_{os.environ['SLURM_JOB_ID']}"
    out.mkdir(exist_ok=False); archive = out/'executed_sources'; archive.mkdir()
    for p in (Path(__file__), *(HERE/n for n in SOURCES)):
        q=archive/p.name; q.write_bytes(p.read_bytes())
        require(sha(q) == (a.expected_runner_sha if p == Path(__file__) else pins[str(p)]), 'Archive differs')
    import numpy as np
    import scipy
    sys.path.insert(0,str(archive))
    modules = {n[:-3]:importlib.import_module(n[:-3]) for n in SOURCES if n != 'check_bfix_tests_v13.py'}
    require(all(Path(m.__file__).resolve().parent == archive.resolve() for m in modules.values()), 'Import escaped archive')
    original, loader = modules['likelihood_v1'], modules['bfix_v1']
    require(loader.FIXTURE.resolve() == FIXTURE.resolve(), 'Fixture loader differs')
    design, settings = loader.load(); template = design.model(); null = read(NULL)
    require(null['converged'] and settings['B'] == 50 and len(template.rows) == 173
        and (template.G,template.Z,template.size) == (8,7,28), 'Original synthetic design differs')
    started=time.perf_counter()
    def safe(x):
        if isinstance(x,np.ndarray): return safe(x.tolist())
        if isinstance(x,np.generic): return safe(x.item())
        if isinstance(x,float) and not math.isfinite(x): return {'nonfinite_float':x.hex()}
        if isinstance(x,dict): return {k:safe(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)): return [safe(v) for v in x]
        return x
    def save(name,value):
        with (out/(name+'.json')).open('x') as f:
            json.dump(safe(value),f,indent=2,allow_nan=False); f.write('\n'); f.flush(); os.fsync(f.fileno())
    def event(fields):
        with (out/'progress.jsonl').open('a') as f:
            f.write(json.dumps(safe(dict(elapsed_seconds=time.perf_counter()-started,**fields)),allow_nan=False)+'\n'); f.flush()
    save('protocol',dict(draw=a.draw,B=50,source_and_input_sha256=pins,node=os.environ.get('SLURMD_NODENAME'),
        environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__),
        original_draw_body_unchanged=True,original_generating_null_parameters=True,original_starts_and_tie=True,
        original_convergence_cutoff=1e-6,original_constrained_solver_before_refinement=True,
        parallel_original_null_shifts=True,ordered_null_shift_results=True,
        departure_from_original_adaptive20_approximation=True,separate_bounded_pilot=True,
        no_existing_draw_reclassification=True,no_backend_fallback=True,no_calibration=True,no_human_data=True))
    subject=None
    try:
        wrapper=modules['repair_p3_adaptive_interval_model_v1']
        Base=wrapper.model_class(original,modules,modules['repair_p3_adaptive_intervals_weighted_v1'],
            modules['repair_p3_adaptive_peak_scan_v1'],event_sink=event)
        Base=modules['repair_p3_parallel_null_shifts_v1'].model_class(Base,original)
        preflight=Base(template.rows,gene_count=8); fp=read(FULL/'protocol.json'); fs=read(FULL/'summary.json'); q=np.asarray(fp['parameters'])
        with_h=preflight.marginal(q,fp['variance'],with_hessian=True)
        without_h=preflight.marginal(q,fp['variance'])
        preflight_differences=[abs(with_h[0]-fs['logmass']),float(np.max(np.abs(with_h[1]-fs['gradient']))),
            float(np.max(np.abs(with_h[2]-fs['hessian'])))]
        same_score=with_h[0] == without_h[0] and np.array_equal(with_h[1],without_h[1])
        original_boundary=original.Model(template.rows,gene_count=8).marginal(q,0.)
        new_boundary=preflight.marginal(q,0.)
        boundary_exact=all(np.asarray(x).tobytes() == np.asarray(y).tobytes() for x,y in zip(original_boundary,new_boundary))
        save('preflight',dict(full_reference_differences=preflight_differences,score_and_Hessian_modes_same_LL_gradient=same_score,
                             original_zero_boundary_exact=boundary_exact))
        require(all(x <= tol for x,tol in zip(preflight_differences,(1e-12,1e-12,1e-9))) and same_score and boundary_exact,'Adaptive interval preflight failed')
        Polished=modules['repair_p3_gamma_order_polish_model_v2'].model_class(Base,original,
            modules['repair_p3_bb_gamma_difference_v3'],modules['repair_p3_score_integral_polish'],
            expected_polish_sha=pins[str(HERE/'repair_p3_score_integral_polish.py')],event_sink=event)
        uniforms=np.load(FIXTURE/'draws/p3_uniforms.npy',mmap_mode='r',allow_pickle=False)
        require(uniforms.shape == (50,173,2), 'Original uniform shape differs')
        null_p=np.asarray(null['parameters']); expected=[]
        for row,(u1,u2) in zip(template.rows,uniforms[a.draw]):
            eta,r=template.predictors(row,null_p)[:2]
            expected.append(original.conditional_draw(row,eta,r,u1,u2).a)
        previous=read(COUNTREF[a.draw])
        require(previous['counts'] == expected
            and previous['full_row_keys'] == [[r.individual,settings['gene_order'][r.gene]] for r in template.rows]
            and hashlib.sha256(np.asarray(expected,dtype='<i8').tobytes()).hexdigest() == previous['count_vector_sha256'],
            'Original supplied draw counts or order differ')
        previous=dict(full_row_keys=[[r.individual,settings['gene_order'][r.gene]] for r in template.rows])
        save('draw_counts',dict(counts=expected,full_row_keys=previous['full_row_keys'],
            count_vector_sha256=hashlib.sha256(np.asarray(expected,dtype='<i8').tobytes()).hexdigest(),exact_previous_counts=True))
        caught=[]
        from types import SimpleNamespace
        def factory(rows):
            nonlocal subject
            require([r.a for r in rows] == expected, 'Recipe counts differ')
            subject=Polished(rows,gene_count=8)
            original_null_fit=subject.fit
            def recorded_null_fit(*args,**kw):
                result=original_null_fit(*args,**kw)
                save('null_fit_snapshot',result)
                return result
            subject.fit=recorded_null_fit
            def fit(*args,**kw):
                event(dict(event='fit_started',draw=a.draw))
                try:return subject.fit_marginal(*args,**kw)
                except Exception as exc:
                    known=isinstance(exc,wrapper.IntervalRefusal) or isinstance(exc,(FloatingPointError,np.linalg.LinAlgError))
                    caught.append(dict(error_type=type(exc).__name__,message=str(exc),known_numerical_refusal=known,
                        integration_failure=getattr(subject,'last_integration_failure',None),traceback=traceback.format_exc()))
                    raise
            return SimpleNamespace(fit_marginal=fit)
        namespace=dict(np=np,Model=factory,conditional_draw=original.conditional_draw,model=template,settings=settings,
            null_p=null_p,observed_lr=0.,p3_failed=[],p3_lr=[],p3_exceed=[],save=save)
        exec(compile(compiled_tree,'<unchanged_original_P3_draw>','exec'),namespace)
        namespace['original_draw'](a.draw,uniforms[a.draw])
        raw=read(out/f'p3_{a.draw:03d}.json'); result=raw['fit']
        require(not any(not x['known_numerical_refusal'] for x in caught),'Unqualified caught exception')
        save('null_polish_evidence',dict(trials=getattr(subject,'null_polish_evidence',[])))
        diagnostic=None
        if result is not None:
            ll,g=subject.marginal(np.asarray(result['parameters']),result['v']); point=np.r_[result['parameters'],result['v']]
            bounds=[b if keep else (result['null']['parameters'][i],result['null']['parameters'][i])
                    for i,(b,keep) in enumerate(zip(subject.bounds(),subject.informative()))]+[(0,None)]
            pg=g.copy()
            for i,(lo,hi) in enumerate(bounds):
                if (lo is not None and lo == hi) or (lo is not None and point[i]<=lo+1e-10 and g[i]<0) or (hi is not None and point[i]>=hi-1e-10 and g[i]>0): pg[i]=0.
            diagnostic=dict(loglik=ll,gradient=g,parameters=point,bounds=bounds,max_projected_gradient=float(np.max(np.abs(pg))),
                selected_null_converged=result['null']['converged'],selected_marginal_converged=result['converged'])
            save('returned_diagnostic',diagnostic)
        for path,digest in pins.items():require(sha(Path(path)) == digest,'Input changed during pilot')
        save('summary',dict(completed=True,draw=a.draw,failed_draw=a.draw in namespace['p3_failed'],LR=namespace['p3_lr'][0],
            exceed=namespace['p3_exceed'][0],caught_exceptions=caught,diagnostic=diagnostic,
            positive_marginal_calls=getattr(subject,'marginal_call_count',0),seconds=time.perf_counter()-started,
            no_existing_draw_reclassification=True,no_fifty_draw_agreement=True,no_calibration=True,no_human_validation=True))
    except Exception as exc:
        save('failure',dict(error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc(),
            integration_failure=getattr(subject,'last_integration_failure',None))); raise


if __name__ == '__main__':main()
