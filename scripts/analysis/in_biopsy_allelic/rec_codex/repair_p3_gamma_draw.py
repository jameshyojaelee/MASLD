#!/usr/bin/env python3
"""One supplied P3 draw, with the original draw body and native Gamma model.

All fifty prespecified draws remain required. No held comparison outcomes,
calibration datasets, human outcomes, replacement draws or backend fallback.
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
OBSERVED=REC/'repair_p3_gamma_observed_native_22080548'
ORIGINAL_NULL=REC/'p3_parallel_20261003T130323Z/baseline/observed_null.json'
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
    'check_bfix_tests_v13.py','repair_p3_bb_gamma_difference_v3.py',
    'repair_p3_gamma_parallel_model.py','repair_p3_gamma_native_model.py')
CAP=32*1024**2


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m


def draw_function(source):
    loops=[n for n in ast.parse(source).body if isinstance(n,ast.For) and isinstance(n.target,ast.Tuple)
        and [getattr(x,'id',None) for x in n.target.elts]==['b','draw']]
    require(len(loops)==1,'Original P3 draw loop differs')
    function=ast.FunctionDef(name='original_draw',args=ast.arguments(posonlyargs=[],
        args=[ast.arg(arg='b'),ast.arg(arg='draw')],kwonlyargs=[],kw_defaults=[],defaults=[]),
        body=loops[0].body,decorator_list=[])
    module=ast.fix_missing_locations(ast.Module(body=[function],type_ignores=[]))
    require(ast.dump(ast.Module(body=module.body[0].body,type_ignores=[]),include_attributes=False)
        ==ast.dump(ast.Module(body=loops[0].body,type_ignores=[]),include_attributes=False),'Original draw body changed')
    return module


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--draw',type=int,choices=range(50),required=True)
    ap.add_argument('--pins',type=Path,required=True);ap.add_argument('--expected-pins-sha',required=True)
    ap.add_argument('--expected-runner-sha',required=True);ap.add_argument('--inspect',action='store_true');args=ap.parse_args()
    require(sha(Path(__file__))==args.expected_runner_sha and sha(args.pins)==args.expected_pins_sha,'Runner/manifest differs')
    pins={Path(p):h for p,h in json.loads(args.pins.read_text()).items()}
    for p,h in pins.items():
        require(p.is_absolute() and p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==h,'Pinned input differs')
        if p.suffix=='.py':ast.parse(p.read_text())
    required=[ORIGINAL_NULL,*(HERE/n for n in SOURCES),*(FIXTURE/n for n in ('SHA256SUMS','globals.json','cohorts.tsv',
        'individuals.tsv','rows.tsv','draws/p3_uniforms.npy')),
        *(OBSERVED/n for n in ('protocol.json','summary.json','observed_null.json','observed_marginal.json',
            'observed_null_diagnostic.json','observed_marginal_diagnostic.json'))]
    require(all(p in pins for p in required),'Missing source/input/observed-fit pin')
    require(pins[HERE/'check_bfix_tests_v13.py']=='a4ca55a52d28dfdbfeb0fc2611c2dcb3ec412987c8f7618f59002cbf72b08c2d',
            'Original supplied P3 recipe differs')
    read=lambda p:json.loads(p.read_text())
    summary=read(OBSERVED/'summary.json');protocol=read(OBSERVED/'protocol.json')
    require(all(summary[k] is True for k in ('completed','native_null_converged','native_marginal_converged',
        'both_projected_gradients_strictly_below_cutoff','previous_original_fit_agreement')),'Observed fit incomplete')
    for n in SOURCES:
        if n!='check_bfix_tests_v13.py':
            require(protocol['source_and_input_sha256'][str(HERE/n)]==pins[HERE/n],'Observed-fit numerical source differs')
    null=read(OBSERVED/'observed_null.json')['fit'];marginal=read(OBSERVED/'observed_marginal.json')['fit']
    require(null['converged'] and marginal['converged'] and marginal['v']==0.,'Observed scope differs')
    for n in ('null','marginal'):
        require(read(OBSERVED/('observed_'+n+'_diagnostic.json'))['max_projected_gradient']<1e-6,
                'Observed projected score failed')
    sums={n.removeprefix('./'):h for h,n in
        (line.split(maxsplit=1) for line in (FIXTURE/'SHA256SUMS').read_text().splitlines())}
    for p in required:
        if p.parent==FIXTURE or p==FIXTURE/'draws/p3_uniforms.npy':
            if p.name!='SHA256SUMS':require(sums[str(p.relative_to(FIXTURE))]==pins[p],'Supplied fixture differs')
    module=draw_function((HERE/'check_bfix_tests_v13.py').read_text())
    if args.inspect:print(json.dumps(dict(draw=args.draw,B=50,original_draw_body_unchanged=True)));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
            'Threads not fixed')
    out=REC/('repair_p3_gamma_draw_'+str(args.draw)+'_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    claim=args.pins.parent/f'draw_{args.draw:03d}.claim.json'
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in SOURCES)):
        q=archive/p.name;q.write_bytes(p.read_bytes())
        require(sha(q)==(args.expected_runner_sha if p==Path(__file__) else pins[p]),'Archive differs')
    import numpy as np
    import scipy
    def safe(x):
        if isinstance(x,np.ndarray):return safe(x.tolist())
        if isinstance(x,np.generic):return safe(x.item())
        if isinstance(x,float) and not math.isfinite(x):return {'nonfinite_float':x.hex()}
        if isinstance(x,dict):return {k:safe(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [safe(v) for v in x]
        return x
    def save(name,value):
        body=(json.dumps(safe(value),allow_nan=False,indent=2)+'\n').encode()
        require(sum(p.stat().st_size for p in out.rglob('*') if p.is_file())+len(body)<CAP-65536,'Output cap exceeded')
        with (out/(name+'.json')).open('xb') as f:f.write(body);f.flush();os.fsync(f.fileno())
    models=[];context=None;caught=[];start=time.perf_counter()
    def record_exception(where,exc):
        frames=[];tb=exc.__traceback__
        while tb is not None:
            frames.append(dict(file=tb.tb_frame.f_code.co_filename,line=tb.tb_lineno));tb=tb.tb_next
        known=isinstance(exc,(ArithmeticError,np.linalg.LinAlgError)) or (
            type(exc) is ValueError and str(exc)=='nonconcave mode in adaptive quadrature'
            and frames and Path(frames[-1]['file']).resolve()==(archive/'likelihood_v1.py').resolve())
        caught.append(dict(where=where,error_type=type(exc).__name__,message=str(exc),
            known_numerical_refusal=known,traceback=traceback.format_exc(),frames=frames))
    try:
        sys.path.insert(0,str(archive));loader=importlib.import_module('bfix_v1');importlib.import_module('tests_v1')
        for n in SOURCES[:4]:
            require(Path(sys.modules[n[:-3]].__file__).resolve()==(archive/n).resolve(),'Original import escaped archive')
        design,settings=loader.load();model=design.model();original=load('Gamma_draw_original',archive/'likelihood_v1.py')
        require(Path(loader.FIXTURE).resolve()==FIXTURE.resolve(),'Loader fixture path differs')
        gamma=load('Gamma_draw_arithmetic',archive/'repair_p3_bb_gamma_difference_v3.py')
        snapshots=load('Gamma_draw_snapshots',archive/'repair_p3_gamma_parallel_model.py')
        native=load('Gamma_draw_native',archive/'repair_p3_gamma_native_model.py')
        require(dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)==protocol['environment'],
                'Observed-fit environment differs')
        require(settings['B']==50 and len(model.rows)==173 and len(design.participants)==60
            and hashlib.sha256(np.asarray([r.a for r in model.rows],dtype='<i8').tobytes()).hexdigest()
            ==protocol['count_vector_sha256'],'Original observed fixture differs')
        ordered=[(r['individual_id'],r['gene_id']) for r in loader.table('rows.tsv')]
        require(ordered==[(r.individual,settings['gene_order'][r.gene]) for r in model.rows],'Supplied uniform row order differs')
        uniforms=np.load(FIXTURE/'draws/p3_uniforms.npy',mmap_mode='r',allow_pickle=False)
        require(uniforms.shape==(50,173,2) and np.isfinite(uniforms).all() and np.all((uniforms>=0)&(uniforms<1)),
                'Supplied uniforms invalid')
        original_bb=original.bb_terms;before=snapshots.row_signatures(model.rows,np)
        Model=native.native_gamma_model_class(original,gamma,snapshots)
        old_null=read(ORIGINAL_NULL);require(old_null['converged'] is True,'Original null reference incomplete')
        old_p=np.asarray(old_null['parameters'],dtype=float);expected=[]
        for row,(u1,u2) in zip(model.rows,uniforms[args.draw]):
            eta,r,_,_=model.predictors(row,old_p)
            expected.append(original.conditional_draw(row,eta,r,u1,u2).a)
        def recorded_draw(*a,**kw):
            try:return original.conditional_draw(*a,**kw)
            except Exception as exc:record_exception('conditional_draw',exc);raise
        def recorded_predictors(*a,**kw):
            try:return model.predictors(*a,**kw)
            except Exception as exc:record_exception('draw_predictors',exc);raise
        def factory(rows):
            from types import SimpleNamespace
            require(len(rows)==173,'Draw row count differs')
            require([r.a for r in rows]==expected,'Draw counts differ from original baseline null')
            try:m=Model(rows,gene_count=8)
            except Exception as exc:record_exception('draw_model_constructor',exc);raise
            models.append(m)
            def recorded_fit(*a,**kw):
                try:return m.fit_marginal(*a,**kw)
                except Exception as exc:record_exception('fit_marginal_including_null_fit',exc);raise
            return SimpleNamespace(fit_marginal=recorded_fit)
        null_p=np.asarray(null['parameters'],dtype=float)
        from types import SimpleNamespace
        draw_parent=SimpleNamespace(rows=model.rows,predictors=recorded_predictors)
        context=dict(np=np,Model=factory,conditional_draw=recorded_draw,model=draw_parent,
            settings=settings,null_p=null_p,observed_lr=0.,p3_failed=[],p3_lr=[],p3_exceed=[],save=save)
        exec(compile(module,'unchanged_original_Gamma_P3_draw','exec'),context)
        save('protocol',dict(draw=args.draw,B=50,source_and_input_sha256={str(p):h for p,h in pins.items()},
            observed_LR=0.,original_draw_body_unchanged=True,section4_from_null_unchanged=True,
            no_backend_fallback=True,no_held_results_read=True,no_calibration=True,no_human_data=True))
        with claim.open('x') as f:
            json.dump(dict(draw=args.draw,jobid=os.environ['SLURM_JOB_ID'],output=str(out)),f)
            f.flush();os.fsync(f.fileno())
        context['original_draw'](args.draw,uniforms[args.draw])
        raw=read(out/f'p3_{args.draw:03d}.json')
        require((raw['fit'] is None)==bool(caught),'Original draw refusal lacks a recorded exception')
        require(len(models)<=1,'More than one original draw model')
        subject=models[0] if models else None;counts=[r.a for r in subject.rows] if subject is not None else None
        save('draw_counts',dict(counts=counts,full_row_keys=ordered,
            count_vector_sha256=hashlib.sha256(np.asarray(counts,dtype='<i8').tobytes()).hexdigest() if counts is not None else None,
            original_baseline_count_vector_sha256=hashlib.sha256(np.asarray(expected,dtype='<i8').tobytes()).hexdigest(),
            exact_original_baseline_counts_equal=counts==expected if counts is not None else None))
        fitting_failure=getattr(subject,'last_integration_failure',None)
        require(original.bb_terms is original_bb and snapshots.row_signatures(model.rows,np)==before,
                'Original BB/observed rows changed')
        save('recipe_summary',dict(completed_original_draw_body=True,draw=args.draw,B=50,
            original_recipe_failed_draw=args.draw in context['p3_failed'],LR=context['p3_lr'][0],
            exceed=context['p3_exceed'][0],caught_exceptions=caught,fitting_integration_failure=fitting_failure,
            unknown_caught_exception_present=any(not x['known_numerical_refusal'] for x in caught),
            no_calibration=True,no_fifty_draw_agreement_claim=True))
        fit=raw['fit'];diagnostic=None
        diagnostic_error=None
        def returned_diagnostic():
            p=np.asarray(fit['parameters'],dtype=float);ll,g=subject.marginal(p,float(fit['v']))
            require(p.shape==(28,) and g.shape==(29,) and np.isfinite(p).all()
                and np.isfinite(g).all() and np.isfinite(ll),'Returned fit diagnostic invalid')
            bounds=[b if keep else (fit['null']['parameters'][j],fit['null']['parameters'][j])
                for j,(b,keep) in enumerate(zip(subject.bounds(),subject.informative()))]+[(0,None)]
            point=np.r_[p,float(fit['v'])];projected=g.copy()
            for j,(lo,hi) in enumerate(bounds):
                if (lo is not None and lo==hi) or (lo is not None and point[j]<=lo+1e-10 and g[j]<0) or (hi is not None and point[j]>=hi-1e-10 and g[j]>0):projected[j]=0.
            diagnostic=dict(loglik=ll,gradient=g,bounds=bounds,parameters=point,projected_gradient=projected,
                max_projected_gradient=float(np.max(np.abs(projected))),native_converged=fit['converged'])
            try:cl,cg=subject.marginal(p,float(fit['v']),nodes=40)
            except (ValueError,ArithmeticError,np.linalg.LinAlgError) as exc:
                diagnostic['native40']=dict(available=False,error_type=type(exc).__name__,message=str(exc))
            else:diagnostic['native40']=dict(available=True,order_refinement_applicable=fit['v']>0,
                loglik=cl,gradient=cg,LL_difference=cl-ll,score_max_difference=float(np.max(np.abs(cg-g))),
                LL_score_agree=bool(np.isfinite(cl) and np.isfinite(cg).all() and abs(cl-ll)<=1e-8+1e-8*abs(cl)
                    and np.allclose(cg,g,atol=1e-8,rtol=1e-6)))
            save('returned_fit_diagnostic',diagnostic)
            return diagnostic
        if fit is not None:
            try:diagnostic=returned_diagnostic()
            except Exception as exc:
                diagnostic_error=dict(error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
                save('diagnostic_failure',diagnostic_error)
        require(original.bb_terms is original_bb and snapshots.row_signatures(model.rows,np)==before,
                'Original BB/observed rows changed')
        for p,h in pins.items():require(sha(p)==h,'Pinned input changed during draw')
        save('summary',dict(completed=True,draw=args.draw,B=50,failed_draw=args.draw in context['p3_failed'],
            LR=context['p3_lr'][0],exceed=context['p3_exceed'][0],elapsed_seconds=time.perf_counter()-start,
            fitting_integration_failure=fitting_failure,caught_exceptions=caught,
            diagnostic_completed=fit is None or diagnostic is not None,diagnostic_failure=diagnostic_error,
            no_fifty_draw_agreement_claim=True,no_calibration=True,no_human_validation=True))
    except Exception as exc:
        save('failure',dict(completed=False,draw=args.draw,error_type=type(exc).__name__,message=str(exc),
            traceback=traceback.format_exc(),last_integration_failure=getattr(models[-1],'last_integration_failure',None) if models else None,
            no_calibration=True,no_human_validation=True));raise


if __name__=='__main__':main()
