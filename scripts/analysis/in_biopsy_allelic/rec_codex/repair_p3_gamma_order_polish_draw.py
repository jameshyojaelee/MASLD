#!/usr/bin/env python3
"""Bounded full-fit trials with the saved-point checked null polish.

These trials do not replace the existing fifty-draw family or its failures.
No held comparison outcomes, human outcomes, calibration or backend fallback.
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
OBSERVED=REC/'repair_p3_gamma_observed_batched_GL_22080549'
NATIVE_NULL=REC/'repair_p3_gamma_observed_native_22080548/observed_null.json'
ORIGINAL_NULL=REC/'p3_parallel_20261003T130323Z/baseline/observed_null.json'
FULL_ORDER=REC/'repair_p3_gamma_order_full_check_22083751'
ORDER_REVIEW=REC/'repair_p3_gamma_order_full_root_review_v2_20261007T065745Z/review.json'
REFERENCE=REC/'repair_p3_gamma_GL_draw_0_22082078'
REFERENCE_REVIEW=REC/'repair_p3_gamma_GL_draw0_root_review_20261007T064057Z/review.json'
FIT_PROOF=REC/'repair_p3_gamma_order_draw_0_22083807'
FIT_REVIEW=REC/'repair_p3_gamma_order_draw0_root_review_20261007T074750Z/review.json'
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
    'check_bfix_tests_v13.py','repair_p3_bb_gamma_difference_v3.py',
    'repair_p3_gamma_parallel_model.py','repair_p3_gamma_native_model.py',
    'repair_p3_integral.py','repair_p3_local_gl.py','repair_p3_local_moments.py',
    'repair_p3_batched_conditional.py','repair_p3_prepared_conditional_v2.py',
    'repair_p3_bb_gamma_batched.py','repair_p3_gamma_batched_conditional.py',
    'repair_p3_gamma_batched_model.py','repair_p3_local_gl_adaptive_order.py','repair_p3_gamma_order_model.py')
BASE_SOURCES=SOURCES
SOURCES=(*BASE_SOURCES,'repair_p3_score_integral_polish.py','repair_p3_gamma_order_polish_model.py')
POLISH_REVIEW=REC/'score_integral_polish_v3_root_review_20261007T084056Z/review.json'
POLISH_CHECKS=(REC/'score_integral_polish_check_22084403',REC/'score_integral_polish_check_22084404')
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
    ap.add_argument('--draw',type=int,choices=(0,11),required=True)
    ap.add_argument('--trial',choices=('intel0','amd0','failed11'),required=True)
    ap.add_argument('--pins',type=Path,required=True);ap.add_argument('--expected-pins-sha',required=True)
    ap.add_argument('--expected-runner-sha',required=True);ap.add_argument('--inspect',action='store_true');args=ap.parse_args()
    require(sha(Path(__file__))==args.expected_runner_sha and sha(args.pins)==args.expected_pins_sha,'Runner/manifest differs')
    require(args.pins.resolve().parent==Path(__file__).resolve().parent
            and args.pins.resolve().parent.name.startswith('repair_p3_gamma_order_polish_launch_'),
            'Dedicated bounded polish launch directory required')
    pins={Path(p):h for p,h in json.loads(args.pins.read_text()).items()}
    for p,h in pins.items():
        require(p.is_absolute() and p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==h,'Pinned input differs')
        if p.suffix=='.py':ast.parse(p.read_text())
    required=[POLISH_REVIEW,*(d/n for d in POLISH_CHECKS for n in ('summary.json','protocol.json',
        'proof0.json','repeat0.json','failed11.json')),ORIGINAL_NULL,NATIVE_NULL,FULL_ORDER/'summary.json',FULL_ORDER/'protocol.json',
        ORDER_REVIEW,REFERENCE_REVIEW,FIT_REVIEW,
        *(FIT_PROOF/n for n in ('summary.json','protocol.json','p3_000.json','comparison_to_GL32.json',
            'returned_fit_diagnostic.json')),
        *(REFERENCE/n for n in ('summary.json','protocol.json','p3_000.json')),
        *(HERE/n for n in SOURCES),*(FIXTURE/n for n in ('SHA256SUMS','globals.json','cohorts.tsv',
        'individuals.tsv','rows.tsv','draws/p3_uniforms.npy')),
        *(OBSERVED/n for n in ('protocol.json','summary.json','observed_null.json','observed_marginal.json',
            'observed_null_diagnostic.json','observed_marginal_diagnostic.json'))]
    require(all(p in pins for p in required),'Missing source/input/observed-fit pin')
    require(pins[HERE/'check_bfix_tests_v13.py']=='a4ca55a52d28dfdbfeb0fc2611c2dcb3ec412987c8f7618f59002cbf72b08c2d',
            'Original supplied P3 recipe differs')
    read=lambda p:json.loads(p.read_text())
    order_review=read(ORDER_REVIEW);reference_review=read(REFERENCE_REVIEW)
    fit_review=read(FIT_REVIEW);fit_summary=read(FIT_PROOF/'summary.json')
    fit_protocol=read(FIT_PROOF/'protocol.json');fit_proof=read(FIT_PROOF/'p3_000.json')['fit']
    old_fit=read(REFERENCE/'p3_000.json')['fit']
    require(fit_review['passed'] is True and fit_review['only_one_paired_full_fit'] is True
        and fit_summary['completed'] is True and fit_summary['failed_draw'] is False
        and read(FIT_PROOF/'comparison_to_GL32.json')['passed'] is True
        and read(FIT_PROOF/'returned_fit_diagnostic.json')['max_projected_gradient']<1e-6
        and fit_proof['converged'] and fit_proof['null']['converged']
        and fit_proof['null']==old_fit['null']
        and all(fit_proof[k]==old_fit[k] for k in ('parameters','loglik','gradient','v')),
        'Completed adaptive-order full-fit proof differs')
    require(all(fit_protocol['source_and_input_sha256'][str(HERE/n)]==pins[HERE/n] for n in BASE_SOURCES),
        'Completed adaptive-order fit helper sources differ')
    require(all(order_review['input_sha256'][str(FULL_ORDER/n)]==pins[FULL_ORDER/n]
        for n in ('summary.json','protocol.json')),'Full8gene reviewed outputs differ')
    require(all(reference_review['input_sha256'][str(REFERENCE/n)]==pins[REFERENCE/n]
        for n in ('summary.json','protocol.json','p3_000.json')),'Completed draw0 reference differs')
    require(order_review['full_eight_gene_fixed_point_verified'] is True
        and order_review['all14_preflight_pass'] is True
        and order_review['executed_rule_geometry_subset_verified'] is True
        and read(FULL_ORDER/'summary.json')['completed'] is True
        and reference_review['GL_original_recipe_converged'] is True
        and reference_review['draw']==0 and read(REFERENCE/'summary.json')['completed'] is True,
        'Completed full8gene and draw0 numerical prerequisites required')
    limits={'logmass':1e-9,'gradient':1e-8,'hessian':1e-5}
    require(set(order_review['full_differences_vs32'])==set(limits)
        and all(math.isfinite(order_review['full_differences_vs32'][k])
            and 0<=order_review['full_differences_vs32'][k]<=v for k,v in limits.items()),
        'Full8gene numerical comparisons differ')
    orderp=read(FULL_ORDER/'protocol.json')
    for n in BASE_SOURCES:
        if n not in ('check_bfix_tests_v13.py','repair_p3_gamma_native_model.py'):
            require(orderp['source_and_input_sha256'][str(HERE/n)]==pins[HERE/n],
                    'Checked adaptive-order source differs')
    summary=read(OBSERVED/'summary.json');protocol=read(OBSERVED/'protocol.json')
    require(all(summary[k] is True for k in ('completed','native_null_converged','native_marginal_converged',
        'both_projected_gradients_strictly_below_cutoff','previous_original_fit_agreement')),'Observed fit incomplete')
    for n in BASE_SOURCES:
        if n not in ('check_bfix_tests_v13.py','repair_p3_local_gl_adaptive_order.py','repair_p3_gamma_order_model.py'):
            require(protocol['source_and_input_sha256'][str(HERE/n)]==pins[HERE/n],'Observed-fit numerical source differs')
    null=read(OBSERVED/'observed_null.json')['fit'];marginal=read(OBSERVED/'observed_marginal.json')['fit']
    require(protocol['backend']=='batched_GL' and null['parameters']==read(NATIVE_NULL)['fit']['parameters'],
            'GL and native Gamma observed null parameters differ')
    require(null['converged'] and marginal['converged'] and marginal['v']==0.,'Observed scope differs')
    for n in ('null','marginal'):
        require(read(OBSERVED/('observed_'+n+'_diagnostic.json'))['max_projected_gradient']<1e-6,
                'Observed projected score failed')
    sums={n.removeprefix('./'):h for h,n in
        (line.split(maxsplit=1) for line in (FIXTURE/'SHA256SUMS').read_text().splitlines())}
    for p in required:
        if p.parent==FIXTURE or p==FIXTURE/'draws/p3_uniforms.npy':
            if p.name!='SHA256SUMS':require(sums[str(p.relative_to(FIXTURE))]==pins[p],'Supplied fixture differs')
    require(read(POLISH_REVIEW)['passed'] is True,'Saved-point polish comparison incomplete')
    for d in POLISH_CHECKS:
        require(read(d/'summary.json')['all_candidate_points_passed'] is True
            and read(d/'summary.json')['all_boundary_scores_compatible'] is True
            and read(d/'summary.json')['points']==3
            and read(d/'protocol.json')['source_and_input_sha256'][str(HERE/'repair_p3_score_integral_polish.py')]
                ==pins[HERE/'repair_p3_score_integral_polish.py'],'Polish check/source differs')
        require(all(read(POLISH_REVIEW)['input_sha256'][str(d/n)]==pins[d/n]
            for n in ('proof0.json','repeat0.json','failed11.json')),'Reviewed saved-point outputs differ')
    require((args.trial in ('intel0','amd0') and args.draw==0) or (args.trial=='failed11' and args.draw==11),
            'Trial/draw pairing differs')
    module=draw_function((HERE/'check_bfix_tests_v13.py').read_text())
    compiled=compile(module,'unchanged_original_Gamma_P3_draw','exec')
    if args.inspect:print(json.dumps(dict(draw=args.draw,B=50,original_draw_body_unchanged=True)));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    memory_mb=int(os.environ.get('SLURM_MEM_PER_NODE','0')) or (
        int(os.environ.get('SLURM_MEM_PER_CPU','0'))*int(os.environ.get('SLURM_CPUS_PER_TASK','0')))
    require(int(os.environ.get('SLURM_CPUS_PER_TASK','0'))>=8 and len(os.sched_getaffinity(0))>=8
            and memory_mb>=8192,'Eight allocated CPUs and eight GiB required before fitting')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
            'Threads not fixed')
    out=REC/('repair_p3_gamma_order_polish_draw_'+str(args.draw)+'_'+args.trial+'_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    claim=args.pins.parent/f'draw_{args.draw:03d}_{args.trial}.claim.json'
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
    models=[];context=None;caught=[];kernel=None;start=time.perf_counter()
    def event(fields):
        if fields.get('event')=='marginal_evaluation':
            reports=models[-1].last_integration_reports if models else []
            require(len(reports)==8,'Adaptive-order progress requires all8 gene reports')
            fields=dict(fields,selected_orders=[[r['selected_first_order'],r['selected_refined_order']] for r in reports],
                quadrature_nodes=sum(r['likelihood_evaluations'] for r in reports),
                order64_genes=sum(64 in (r['selected_first_order'],r['selected_refined_order']) for r in reports))
        body=(json.dumps(safe(dict(elapsed_seconds=time.perf_counter()-start,**fields)),allow_nan=False)+'\n').encode()
        progress=out/'progress.jsonl'
        if (progress.stat().st_size if progress.exists() else 0)+len(body)>4*1024**2:return
        require(sum(p.stat().st_size for p in out.rglob('*') if p.is_file())+len(body)<CAP-8*1024**2,
                'Progress reserve exhausted')
        with progress.open('ab') as f:f.write(body);f.flush()
    def record_exception(where,exc):
        frames=[];tb=exc.__traceback__
        while tb is not None:
            frames.append(dict(file=tb.tb_frame.f_code.co_filename,line=tb.tb_lineno));tb=tb.tb_next
        known=(kernel is not None and isinstance(exc,kernel.IntegrationFailure)) or isinstance(exc,(ArithmeticError,np.linalg.LinAlgError)) or (
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
        kernel=load('Gamma_GL_draw_kernel',archive/'repair_p3_integral.py')
        gl=load('Gamma_GL_draw_local',archive/'repair_p3_local_gl_adaptive_order.py')
        conditional=load('Gamma_GL_draw_conditional',archive/'repair_p3_gamma_batched_conditional.py')
        batch=load('Gamma_GL_draw_batch',archive/'repair_p3_gamma_order_model.py')
        Model=batch.gamma_batched_model_class(original,kernel,gl,conditional,
            kernel_path=archive/'repair_p3_integral.py',event_sink=event)
        require(all(getattr(Model,n) is getattr(original.Model,n) for n in native.INHERITED),
                'Original fit/design methods changed')
        polish=load('Gamma_null_polish',archive/'repair_p3_score_integral_polish.py')
        polished=load('Gamma_polish_model',archive/'repair_p3_gamma_order_polish_model.py')
        Model=polished.model_class(Model,original,gamma,polish,
            expected_polish_sha=pins[HERE/'repair_p3_score_integral_polish.py'],event_sink=event)
        require(all(getattr(Model,n) is getattr(original.Model,n) for n in native.INHERITED if n!='constrained_polish'),
                'Methods beyond constrained polish changed')
        old_null=read(ORIGINAL_NULL);require(old_null['converged'] is True,'Original null reference incomplete')
        old_p=np.asarray(old_null['parameters'],dtype=float);expected=[]
        for row,(u1,u2) in zip(model.rows,uniforms[args.draw]):
            eta,r,_,_=model.predictors(row,old_p)
            expected.append(original.conditional_draw(row,eta,r,u1,u2).a)
        null_p=np.asarray(null['parameters'],dtype=float);new_counts=[]
        for row,(u1,u2) in zip(model.rows,uniforms[args.draw]):
            eta,r,_,_=model.predictors(row,null_p)
            new_counts.append(original.conditional_draw(row,eta,r,u1,u2).a)
        require(new_counts==expected,'GL-null draw counts differ from original baseline before claim')
        if args.draw==0:
            require(hashlib.sha256(np.asarray(expected,dtype='<i8').tobytes()).hexdigest()
                ==reference_review['count_vector_sha256'],'Completed draw0 count vector differs')
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
                event(dict(event='fit_started',draw=args.draw))
                try:return m.fit_marginal(*a,**kw)
                except Exception as exc:record_exception('fit_marginal_including_null_fit',exc);raise
            return SimpleNamespace(fit_marginal=recorded_fit)
        from types import SimpleNamespace
        draw_parent=SimpleNamespace(rows=model.rows,predictors=recorded_predictors)
        context=dict(np=np,Model=factory,conditional_draw=recorded_draw,model=draw_parent,
            settings=settings,null_p=null_p,observed_lr=0.,p3_failed=[],p3_lr=[],p3_exceed=[],save=save)
        exec(compiled,context)
        save('protocol',dict(draw=args.draw,B=50,source_and_input_sha256={str(p):h for p,h in pins.items()},
            backend='adaptive_order_GL_with_null_score_integral_polish',trial=args.trial,
            separate_bounded_trial=True,existing_fifty_draws_unchanged=True,requested_node_batch_size=64,
            completed_fit_proof_directory=str(FIT_PROOF),observed_LR=0.,original_draw_body_unchanged=True,section4_from_null_unchanged=True,
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
        save('null_polish_evidence',dict(trials=getattr(subject,'null_polish_evidence',[]),
            only_constrained_polish_changed=True,original_cutoff=1e-6,existing_results_unchanged=True))
        fitting_failure=getattr(subject,'last_integration_failure',None)
        require(original.bb_terms is original_bb and snapshots.row_signatures(model.rows,np)==before,
                'Original BB/observed rows changed')
        save('recipe_summary',dict(completed_original_draw_body=True,draw=args.draw,B=50,
            original_recipe_failed_draw=args.draw in context['p3_failed'],LR=context['p3_lr'][0],
            exceed=context['p3_exceed'][0],caught_exceptions=caught,fitting_integration_failure=fitting_failure,
            unknown_caught_exception_present=any(not x['known_numerical_refusal'] for x in caught),
            no_calibration=True,no_fifty_draw_agreement_claim=True))
        require(not any(not x['known_numerical_refusal'] for x in caught),
                'Unknown caught exception is not an admitted numerical refusal')
        event(dict(event='original_recipe_completed',draw=args.draw))
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
                max_projected_gradient=float(np.max(np.abs(projected))),native_converged=fit['converged'],
                positive_integration_error_allowances=getattr(subject,'last_total_error_allowance',None),
                positive_integration_reports=getattr(subject,'last_integration_reports',None))
            save('GL_fit_diagnostic',diagnostic)
            Comparison=native.native_gamma_model_class(original,gamma,snapshots)
            comparison_model=Comparison(subject.rows,gene_count=8)
            try:cl,cg=comparison_model.marginal(p,float(fit['v']),nodes=40)
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
            positive_marginal_call_count=getattr(subject,'marginal_call_count',None),
            fitting_integration_failure=fitting_failure,caught_exceptions=caught,
            unknown_caught_exception_present=any(not x['known_numerical_refusal'] for x in caught),
            diagnostic_completed=fit is None or diagnostic is not None,diagnostic_failure=diagnostic_error,
            no_fifty_draw_agreement_claim=True,no_calibration=True,no_human_validation=True))
    except Exception as exc:
        save('failure',dict(completed=False,draw=args.draw,error_type=type(exc).__name__,message=str(exc),
            caught_exceptions=caught,
            traceback=traceback.format_exc(),last_integration_failure=getattr(models[-1],'last_integration_failure',None) if models else None,
            no_calibration=True,no_human_validation=True));raise


if __name__=='__main__':main()
