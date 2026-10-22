#!/usr/bin/env python3
"""Original synthetic observed-fit recipe in separate native/GL processes.

All fitter methods, starts, bounds, ties, null snaps and restarts are inherited.
Numerical refusals remain recorded; no automatic backend fallback is used.
No resampling, statistical calibration, human-data analysis or adoption.
"""
import argparse
import ast
import hashlib
import importlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import sys
import time
import traceback

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE=ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
WARM=REC/'bfix_fits_21984981/null.json'
BASELINE=REC/'p3_parallel_20261003T130323Z/baseline'
NATIVE=REC/'repair_p3_gamma_native_model_check_22079139'
GL=REC/'repair_p3_gamma_batched_full_check_22078827'
P3=(REC/'repair_p3_gamma_batched_P3_failed0gene2_22079995',
    REC/'repair_p3_gamma_batched_P3_positive18_22079996')
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
 'repair_p3_integral.py','repair_p3_local_gl.py','repair_p3_local_moments.py',
 'repair_p3_gamma_parallel_model.py','repair_p3_bb_gamma_difference_v3.py',
 'repair_p3_batched_conditional.py','repair_p3_prepared_conditional_v2.py',
 'repair_p3_bb_gamma_batched.py','repair_p3_gamma_batched_conditional.py',
 'repair_p3_gamma_batched_model.py','repair_p3_gamma_native_model.py')
CAP=32*1024**2
COUNT_SHA='7009be4adeec6157fa4bbf3f46e48ec24aa398a0bee4d7906c78c231b7d3fbad'


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--backend',choices=('native','batched_GL'),required=True)
    ap.add_argument('--pins',type=Path,required=True);ap.add_argument('--expected-pins-sha',required=True)
    ap.add_argument('--expected-runner-sha',required=True);ap.add_argument('--inspect',action='store_true');args=ap.parse_args()
    require(sha(Path(__file__))==args.expected_runner_sha and sha(args.pins)==args.expected_pins_sha,'Runner/manifest differs')
    pins={Path(p):v for p,v in json.loads(args.pins.read_text()).items()}
    for p,v in pins.items():
        require(p.is_absolute() and p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==v,'Pinned input differs: '+str(p))
        if p.suffix=='.py':ast.parse(p.read_text())
    required=[*(HERE/n for n in SOURCES),WARM,
        *(BASELINE/n for n in ('observed_null.json','observed_marginal.json')),
        *(FIXTURE/n for n in ('SHA256SUMS','globals.json','cohorts.tsv','individuals.tsv','rows.tsv')),
        *(NATIVE/n for n in ('protocol.json','summary.json','state_probes.json')),
        *(GL/n for n in ('protocol.json','summary.json','comparison.json'))]
    for d in P3:required += [d/'protocol.json',d/'summary.json']
    require(all(p in pins for p in required),'Required source/admission missing')
    read=lambda p:json.loads(p.read_text())
    ns=read(NATIVE/'summary.json');gs=read(GL/'summary.json')
    require(ns['completed'] and ns['byte_identical_integrals']==112 and ns['original_boundary_full_byte_calls']==16
        and ns['all_state_probes_passed'] and all(v is True for v in read(NATIVE/'state_probes.json').values()),
        'Native evaluator checks incomplete')
    require(gs['completed'] and gs['all_eight_gene_comparisons_passed'] and gs['full_comparison_passed']
        and gs['all_fixed_inputs_unchanged'] and read(GL/'comparison.json')['full']['passed'],'Full batched GL checks incomplete')
    for d in P3:
        s=read(d/'summary.json');require(all(s[k] is True for k in ('completed','returned','nominal_comparison_passed',
            'numerical_allowances_passed','all_fixed_inputs_unchanged')),'Batched P3 point incomplete')
    admitted_sources=set()
    for d in (NATIVE,GL,*P3):
        for path,digest in read(d/'protocol.json')['source_and_input_sha256'].items():
            if Path(path) in (HERE/n for n in SOURCES):
                require(pins[Path(path)]==digest,'Admitted source differs')
                admitted_sources.add(Path(path))
    require(admitted_sources=={HERE/n for n in SOURCES},'Incomplete admitted-source coverage')
    sums={name.removeprefix('./'):digest for digest,name in
          (line.split(maxsplit=1) for line in (FIXTURE/'SHA256SUMS').read_text().splitlines())}
    for n in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv'):
        require(pins[FIXTURE/n]==sums[n],'Original fixture differs')
    warm=read(WARM)['fit'];require(warm['converged'] is True,'Original warm null not converged')
    original_null=read(BASELINE/'observed_null.json');original_marginal=read(BASELINE/'observed_marginal.json')
    require(original_null['converged'] and original_marginal['converged'],'Original comparison fits incomplete')
    protocol=dict(backend=args.backend,source_and_input_sha256={str(p):v for p,v in pins.items()},
        fixture_seed=20260930,participants=60,rows=173,G=8,Z=7,P=28,count_vector_sha256=COUNT_SHA,
        recipe='model.fit(prior observed null parameters); model.fit_marginal(tau_max_fixture,new null parameters,null_fit=new null)',
        native_default_nodes=20,final_native_comparison_nodes=40,
        original_warm_cold_starts_ties_bounds_null_snap_restarts=True,original_projected_cutoff=1e-6,
        tau_max_is_cold_start_scale_not_upper_bound=True,no_automatic_backend_fallback=True,
        no_draws=True,no_calibration=True,no_biological_model_validation=True,no_actual_or_sealed_data=True,
        no_independent_full_Hessian=True,output_cap_bytes=CAP)
    if args.inspect:print(json.dumps(protocol));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
            'Threads not fixed before imports')
    out=REC/('repair_p3_gamma_observed_'+args.backend+'_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in SOURCES)):
        dest=archive/p.name;dest.write_bytes(p.read_bytes())
        require(sha(dest)==(args.expected_runner_sha if p==Path(__file__) else pins[p]),'Archive differs')
    import numpy as np
    import scipy
    def safe(x):
        if isinstance(x,np.ndarray):return safe(x.tolist())
        if isinstance(x,np.generic):return safe(x.item())
        if isinstance(x,float) and not math.isfinite(x):return {'nonfinite_float':x.hex()}
        if isinstance(x,dict):return {k:safe(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [safe(v) for v in x]
        return x
    def used():return sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
    def save(name,value):
        body=(json.dumps(safe(value),allow_nan=False,indent=2)+'\n').encode()
        require(used()+len(body)<CAP-65536,'Output cap exhausted')
        with (out/(name+'.json')).open('xb') as f:f.write(body);f.flush();os.fsync(f.fileno())
        print(json.dumps(dict(saved=name,seconds=time.perf_counter()-began)),flush=True)
    began=time.perf_counter();phase='initialization';model=None;events_seen=0;events_saved=0
    def event(fields):
        nonlocal events_seen,events_saved
        events_seen+=1;value=dict(phase=phase,elapsed_seconds=time.perf_counter()-began,**fields)
        body=(json.dumps(safe(value),allow_nan=False)+'\n').encode()
        if fields.get('event')=='marginal_evaluation' and used()+len(body)>4*1024**2:return
        require(used()+len(body)<CAP-8*1024**2,'Progress cap exhausted')
        with (out/'progress.jsonl').open('ab') as f:f.write(body);f.flush()
        events_saved+=1
        if fields.get('event')!='marginal_evaluation':print(json.dumps(safe(value)),flush=True)
    try:
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        require(environment==read(GL/'protocol.json')['environment'],'Environment differs')
        sys.path.insert(0,str(archive));loader=importlib.import_module('bfix_v1');importlib.import_module('tests_v1')
        for n in SOURCES[:4]:
            m=sys.modules[n[:-3]];require(Path(m.__file__).resolve().parent==archive.resolve()
                and sha(Path(m.__file__))==pins[HERE/n],'Original import escaped archive')
        design,settings=loader.load();private=load('observed_Gamma_original',archive/'likelihood_v1.py')
        gamma=load('observed_Gamma_scalar',archive/'repair_p3_bb_gamma_difference_v3.py')
        snapshots=load('observed_Gamma_snapshots',archive/'repair_p3_gamma_parallel_model.py')
        native=load('observed_Gamma_native',archive/'repair_p3_gamma_native_model.py')
        original_bb=private.bb_terms;base=private.Model(design.model().rows,gene_count=8)
        require(len(base.rows)==173 and base.Z==7 and base.size==28 and len(design.participants)==60
            and settings['seed']==20260930,'Synthetic scope differs')
        counts=lambda:hashlib.sha256(np.asarray([r.a for r in base.rows],dtype='<i8').tobytes()).hexdigest()
        require(counts()==COUNT_SHA,'Original counts differ');before=snapshots.row_signatures(base.rows,np)
        if args.backend=='native':
            Model=native.native_gamma_model_class(private,gamma,snapshots,event_sink=event)
        else:
            kernel=load('observed_Gamma_kernel',archive/'repair_p3_integral.py')
            gl=load('observed_Gamma_local_GL',archive/'repair_p3_local_gl.py')
            conditional=load('observed_Gamma_conditional',archive/'repair_p3_gamma_batched_conditional.py')
            batch=load('observed_Gamma_batched',archive/'repair_p3_gamma_batched_model.py')
            Model=batch.gamma_batched_model_class(private,kernel,gl,conditional,
                kernel_path=archive/'repair_p3_integral.py',event_sink=event)
        model=Model(base.rows,gene_count=8)
        require(all(getattr(Model,n) is getattr(private.Model,n) for n in native.INHERITED),'Original fit/design methods changed')
        protocol.update(environment=environment,tau_max_fixture=settings['tau_max_fixture'])
        save('protocol',protocol)
        def diagnostic(fit,marginal=False,null=None):
            p=np.asarray(fit['parameters'],dtype=float);tol=1e-10 if marginal else 1e-8
            if marginal:
                ll,g,h=model.marginal(p,float(fit['v']),with_hessian=True);point=np.r_[p,float(fit['v'])]
                bounds=[b if keep else (null['parameters'][j],null['parameters'][j])
                    for j,(b,keep) in enumerate(zip(model.bounds(),model.informative()))]+[(0,None)]
            else:
                ll,g,h=model.evaluate(p,need_hessian=True);point=p
                bounds=[b if keep else (p[j],p[j]) for j,(b,keep) in enumerate(zip(model.bounds(),model.informative()))]
            require(len(g)==len(point) and h.shape==(len(point),len(point)) and np.isfinite(ll)
                and np.isfinite(g).all() and np.isfinite(point).all(),'Final numerical arrays invalid')
            active=np.asarray([(lo is not None and lo==hi) or (lo is not None and point[j]<=lo+tol and g[j]<0)
                or (hi is not None and point[j]>=hi-tol and g[j]>0) for j,(lo,hi) in enumerate(bounds)])
            projected=g.copy();projected[active]=0.;maximum=float(np.max(np.abs(projected),initial=0))
            boundary=marginal and fit['v']==0.;core=h[:-1,:-1] if boundary else h
            require(np.isfinite(core).all(),'Nonfinite Hessian outside unavailable boundary variance entries')
            return dict(loglik=float(ll),gradient=g,hessian=h,hessian_finite_mask=np.isfinite(h),
                parameters=point,bounds=bounds,projection_bound_tolerance=tol,active=active,projected_gradient=projected,
                max_projected_gradient=maximum,projected_gradient_strictly_below_cutoff=maximum<1e-6,
                native_converged=bool(fit['converged']),variance_Hessian_unavailable_at_zero_boundary=boundary,
                fit_LL_recomputation_difference=float(ll-fit['loglik']),
                fit_gradient_recomputation_agrees=bool(np.allclose(g,fit['gradient'],atol=1e-8,rtol=1e-6)),
                positive_integration_error_allowances=getattr(model,'last_total_error_allowance',None),
                numerical_checks_are_not_global_optimum_or_biological_validation=True)
        phase='observed_null_fit';event(dict(event='fit_started',component='null'));start=time.perf_counter()
        null=model.fit(np.asarray(warm['parameters'],dtype=float));save('observed_null',dict(fit=null,seconds=time.perf_counter()-start))
        nd=diagnostic(null);save('observed_null_diagnostic',nd)
        require(null['converged'],'Original observed null stopping rule: nonconverged')
        phase='observed_marginal_fit';event(dict(event='fit_started',component='marginal'));start=time.perf_counter()
        fit=model.fit_marginal(settings['tau_max_fixture'],np.asarray(null['parameters'],dtype=float),null_fit=null)
        save('observed_marginal',dict(fit=fit,seconds=time.perf_counter()-start))
        phase='final_diagnostics';fd=diagnostic(fit,marginal=True,null=null);save('observed_marginal_diagnostic',fd)
        Comparison=native.native_gamma_model_class(private,gamma,snapshots)
        comparison_model=Comparison(base.rows,gene_count=8)
        try:
            cl,cg,ch=comparison_model.marginal(np.asarray(fit['parameters'],dtype=float),float(fit['v']),nodes=40,with_hessian=True)
        except (ValueError,ArithmeticError,np.linalg.LinAlgError) as exc:
            forty=dict(available=False,status='numerical_refusal',LL_score_agree=None,
                error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
        else:
            core=ch[:-1,:-1] if fit['v']==0 else ch;fh=fd['hessian'][:-1,:-1] if fit['v']==0 else fd['hessian']
            forty=dict(available=True,status='zero_boundary_same_original_evaluator' if fit['v']==0 else 'returned',
                order_refinement_applicable=fit['v']>0,loglik=cl,gradient=cg,hessian=ch,
                LL_difference=float(cl-fd['loglik']),score_max_difference=float(np.max(np.abs(cg-fd['gradient']))),
                nuisance_or_positive_Hessian_max_difference=float(np.max(np.abs(core-fh))),
                LL_score_agree=bool(np.isfinite(cl) and np.isfinite(cg).all()
                    and abs(cl-fd['loglik'])<=1e-8+1e-8*abs(cl)
                    and np.allclose(cg,fd['gradient'],atol=1e-8,rtol=1e-6)))
        save('final_native40_comparison',forty)
        lr=0. if fit['v']==0 else max(0.,2*(fit['loglik']-null['loglik']))
        oldlr=0. if original_marginal['v']==0 else max(0.,2*(original_marginal['loglik']-original_null['loglik']))
        comparison=dict(null_LL_difference=float(null['loglik']-original_null['loglik']),
            null_kappa_difference=float(null['parameters'][16]-original_null['parameters'][16]),
            marginal_LL_difference=float(fit['loglik']-original_marginal['loglik']),
            marginal_kappa_difference=float(fit['parameters'][16]-original_marginal['parameters'][16]),
            LR=lr,original_LR=oldlr,LR_difference=lr-oldlr,
            tau_difference=math.sqrt(fit['v'])-math.sqrt(original_marginal['v']))
        comparison['existing_LL_kappa_LR_or_tau_criteria_passed']=bool(abs(comparison['null_LL_difference'])<=1e-6
            and abs(comparison['marginal_LL_difference'])<=1e-6 and abs(comparison['null_kappa_difference'])<=1e-4
            and (abs(comparison['tau_difference'])<=1e-3 if max(lr,oldlr)>1 else abs(lr-oldlr)<=1e-6))
        save('previous_original_comparison',comparison)
        require(private.bb_terms is original_bb and counts()==COUNT_SHA and snapshots.row_signatures(base.rows,np)==before,
                'Original BB or synthetic inputs changed')
        for path,digest in pins.items():require(sha(path)==digest,'Pinned input changed during fit')
        save('summary',dict(completed=True,backend=args.backend,native_null_converged=bool(null['converged']),
            native_marginal_converged=bool(fit['converged']),null_max_projected_gradient=nd['max_projected_gradient'],
            marginal_max_projected_gradient=fd['max_projected_gradient'],
            both_projected_gradients_strictly_below_cutoff=bool(nd['max_projected_gradient']<1e-6 and fd['max_projected_gradient']<1e-6),
            final_native40_available=forty['available'],final_native40_status=forty['status'],
            final_native40_order_refinement_applicable=bool(fit['v']>0),
            final_native40_LL_score_agree=forty['LL_score_agree'],
            previous_original_fit_agreement=comparison['existing_LL_kappa_LR_or_tau_criteria_passed'],
            start_logliks=fit['start_logliks'],start_converged=fit['start_converged'],kept_start=fit['kept_start'],
            events_seen=events_seen,events_saved=events_saved,elapsed_seconds=time.perf_counter()-began,
            no_draws=True,no_calibration=True,no_independent_full_Hessian=True,no_biological_model_validation=True))
    except Exception as exc:
        save('failure',dict(completed=False,phase=phase,error_type=type(exc).__name__,message=str(exc),
            numerical_refusal=isinstance(exc,(ArithmeticError,np.linalg.LinAlgError))
                or (type(exc) is ValueError and str(exc)=='nonconcave mode in adaptive quadrature'),
            unclassified_ValueError=isinstance(exc,ValueError)
                and str(exc)!='nonconcave mode in adaptive quadrature',
            last_evaluated_point=getattr(model,'last_evaluated_point',None),
            integration_failure=getattr(model,'last_integration_failure',None),
            worker_records=getattr(model,'last_parallel_worker_records',None),traceback=traceback.format_exc(),
            no_draws=True,no_calibration=True,no_biological_model_validation=True))
        raise


if __name__=='__main__':main()
