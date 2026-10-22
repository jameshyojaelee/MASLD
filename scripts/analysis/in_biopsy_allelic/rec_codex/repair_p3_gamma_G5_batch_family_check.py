#!/usr/bin/env python3
"""One of fifteen separate G5 batched fixed-point likelihood checks.

One fixed synthetic point; no fits, resampling or human outcomes.
The independent scalar results and their uncertainty are reused unchanged.
Also retains native GH20/GH40. Does not widen the existing G8 context or
independently validate the full gradient/Hessian.
"""
import argparse
import ast
from dataclasses import asdict
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time
import traceback

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
COLLECTION=REC/'repair_p3_five_gene_accuracy_launch_20261005T143525Z/five_gene_accuracy_collection.json'
NATIVE=REC/'repair_p3_gamma_five_gene_check_22081217'
DENSE=REC/'repair_p3_gamma_GL_five_gene_check_22082316'
PROBE=REC/'repair_p3_gamma_G5_batch_check_22083924'
PROBE_REVIEW=REC/'repair_p3_gamma_G5_probe_root_review_20261007T072503Z/review.json'
COUNTS=(5,25,75,200,500)
CAP=32*1024**2


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pins',type=Path,required=True);ap.add_argument('--expected-pins-sha',required=True)
    ap.add_argument('--expected-runner-sha',required=True);ap.add_argument('--inspect',action='store_true')
    ap.add_argument('--gene',type=int,choices=range(5),required=True)
    ap.add_argument('--tau-multiple',type=float,choices=(.25,1.,4.),required=True)
    args=ap.parse_args();require(sha(Path(__file__))==args.expected_runner_sha
        and sha(args.pins)==args.expected_pins_sha,'Runner/manifest differs')
    pins={Path(p):h for p,h in json.loads(args.pins.read_text()).items()}
    for p,h in pins.items():
        require(p.is_absolute() and p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==h,'Pinned input differs')
        if p.suffix=='.py':ast.parse(p.read_text())
    sources=('likelihood_v1.py','repair_p3_bb_gamma_difference_v3.py','check_quadrature_v1.py',
        'repair_p3_integral.py','repair_p3_composite_gl.py',
        'repair_p3_local_gl_adaptive_order.py','repair_p3_local_moments.py',
        'repair_p3_batched_conditional_g5.py','repair_p3_prepared_conditional_v2.py',
        'repair_p3_bb_gamma_batched.py','repair_p3_gamma_batched_conditional_g5.py',
        'repair_p3_gamma_parallel_model.py')
    require(all(HERE/n in pins for n in sources) and COLLECTION in pins,'Missing source/reference pin')
    require(PROBE_REVIEW in pins and json.loads(PROBE_REVIEW.read_text())['passed'] is True,
        'Completed independent probe review required')
    for name in ('protocol.json','summary.json','comparison_to_dense.json','conditional_probe.json',
                 'GL_returned_raw.json'):
        require(PROBE/name in pins,'Completed probe evidence not pinned')
    probe_protocol=json.loads((PROBE/'protocol.json').read_text())
    probe_summary=json.loads((PROBE/'summary.json').read_text())
    require(probe_summary['completed'] is True and probe_summary['selected_point_likelihood_criteria_passed'] is True
        and probe_summary['requested_gene']==1 and probe_summary['requested_tau_multiple']==1.,
        'Completed probe scope differs')
    for n in sources:require(probe_protocol['source_and_input_sha256'][str(HERE/n)]==pins[HERE/n],
        'Completed probe helper source differs')
    require(all(NATIVE/n in pins for n in ('protocol.json','summary.json')),'Missing checked Gamma family')
    require(all(DENSE/n in pins for n in ('protocol.json','summary.json','GL_returned_raw.json')),
        'Missing completed dense G5 comparison')
    densep=json.loads((DENSE/'protocol.json').read_text());denses=json.loads((DENSE/'summary.json').read_text())
    dense=json.loads((DENSE/'GL_returned_raw.json').read_text())
    require(denses['completed'] is True and denses['selected_point_likelihood_criteria_passed'] is True
        and densep['requested_gene']==1 and densep['requested_tau_multiple']==1.,'Dense G5 point differs')
    for n in sources[:5]:require(densep['source_and_input_sha256'][str(HERE/n)]==pins[HERE/n],
        'Checked dense G5 source differs')
    native_protocol=json.loads((NATIVE/'protocol.json').read_text())
    native_summary=json.loads((NATIVE/'summary.json').read_text())
    require(native_summary['completed'] and native_summary['all15_likelihood_criteria_passed']
        and native_protocol['source_and_input_sha256'][str(HERE/'repair_p3_bb_gamma_difference_v3.py')]
        ==pins[HERE/'repair_p3_bb_gamma_difference_v3.py'],'Checked Gamma arithmetic family differs')
    collection=json.loads(COLLECTION.read_text());entries=collection['results']
    require(collection['execution_complete'] is True and len(entries)==15,'Original family incomplete')
    require({(r['gene'],r['tau_multiple']) for r in entries}=={(g,m) for g in range(5) for m in (.25,1.,4.)},
            'Original point family differs')
    references=[];fixture_paths=[]
    for entry in entries:
        d=REC/('repair_p3_five_gene_accuracy_'+entry['job_id'])
        rp=d/f"gene{entry['gene']}_tau{entry['tau_multiple']}.json";fp=d/'fixture.json';pp=d/'protocol.json'
        require(all(p in pins for p in (rp,fp,pp)),'Raw reference/fixture/protocol not pinned')
        r=json.loads(rp.read_text());protocol=json.loads(pp.read_text())
        require(protocol['environment']==native_protocol['environment'],'Reference environments differ')
        require(r==entry['record'] and r['gene']==entry['gene'] and r['rows']==COUNTS[r['gene']]
            and r['tau']==.05*entry['tau_multiple'] and r['variance']==r['tau']**2,'Raw reference differs')
        original_sources={}
        for path,digest in protocol['source_sha256'].items():
            name=Path(path).name
            require(name not in original_sources,'Duplicate original source basename')
            original_sources[name]=digest
        require(protocol['seed']==20260923 and all(original_sources[n]==pins[HERE/n]
            for n in ('likelihood_v1.py','check_quadrature_v1.py','repair_p3_integral.py',
                'repair_p3_composite_gl.py')),'Original reference construction or GL rules differ')
        first=r['scalar']['first'];expanded=r['scalar']['expanded']
        require('scalar_refusal' not in r and all(math.isfinite(x['logmass']) for x in (first,expanded))
            and 0<=expanded['estimated_logmass_error']
            and abs(expanded['logmass']-first['logmass'])+expanded['estimated_logmass_error']<=1e-9
            and all(0<=x['relative_prior_mass_tail_allowance']<=1e-10 for x in (first,expanded))
            and expanded['z_bound']>first['z_bound'],'Scalar refinement unresolved')
        require(all(math.isfinite(r['original_GH'+str(n)]['logmass']) for n in (20,40)),
            'Original native reference unavailable')
        def finite_json(value):
            if isinstance(value,float):return math.isfinite(value)
            if isinstance(value,dict):return all(finite_json(v) for v in value.values())
            if isinstance(value,list):return all(finite_json(v) for v in value)
            return True
        require(finite_json(r['scalar']),'Nonfinite scalar diagnostic')
        references.append(r);fixture_paths.append(fp)
    require(len({sha(p) for p in fixture_paths})==1
        and sha(fixture_paths[0])==collection['shared_fixture_sha256'],'Original fixtures differ')
    selected=[r for r in references if r['gene']==args.gene and r['tau']==.05*args.tau_multiple]
    require(len(selected)==1,'Requested point is not unique')
    if args.inspect:print(json.dumps(dict(required_points=15,selected_points=1,gene=args.gene,
        tau_multiple=args.tau_multiple,orders=[20,40],references_checked=True)));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
            'Threads not fixed')
    out=REC/('repair_p3_gamma_G5_batch_family_check_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in sources)):
        q=archive/p.name;q.write_bytes(p.read_bytes())
        require(sha(q)==(args.expected_runner_sha if p==Path(__file__) else pins[p]),'Archive differs')
    import numpy as np
    import scipy
    require(dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__)==native_protocol['environment'],
            'Checked numerical environment differs')
    original=load('Gamma_five_original',archive/'likelihood_v1.py')
    gamma=load('Gamma_five_arithmetic',archive/'repair_p3_bb_gamma_difference_v3.py')
    gl=load('Gamma_five_GL',archive/'repair_p3_local_gl_adaptive_order.py')
    prepared=load('Gamma_five_prepared',archive/'repair_p3_gamma_batched_conditional_g5.py')
    snapshots=load('Gamma_five_snapshots',archive/'repair_p3_gamma_parallel_model.py')
    kernel=load('Gamma_five_kernel',archive/'repair_p3_integral.py');gl_settings=kernel.Settings()
    original_bb=original.bb_terms
    require(original.conditional_terms.__globals__['bb_terms'] is original_bb,'Conditional BB lookup differs')
    rows=[]
    for gene,count in enumerate(COUNTS):
        for j in range(count):
            rows.append(original.Row(gene,f'{j:04d}','fixture',20+j%17,10,.6,.3,.1,.01,.02,.05,
                (-1)**(gene+j),float(j%3-1),(float(j%5-2)/10,),float(j%7-3)/10,float(j%5-2)/10,10,2,10))
    model=original.Model(rows);p=model.initial();p[:5]=[-1.2,-.7,.4,.9,1.3];p[10]=.08;p[11]=.04
    rng=np.random.default_rng(20260923);drawn=[]
    for row in model.rows:
        eta,r,_,_=model.predictors(row,p)
        drawn.append(original.conditional_draw(row,eta,r,*rng.random(2)))
    fixture=json.loads(fixture_paths[0].read_text())
    require(np.array_equal(p,fixture['parameters']) and [r.a for r in drawn]==fixture['drawn_counts'],
            'Exact original parameters/counts differ')
    p.setflags(write=False)
    def save(name,value):
        def safe(x):
            if isinstance(x,np.ndarray):return safe(x.tolist())
            if isinstance(x,np.generic):return safe(x.item())
            if isinstance(x,float) and not math.isfinite(x):return {'nonfinite_float':x.hex()}
            if isinstance(x,dict):return {k:safe(v) for k,v in x.items()}
            if isinstance(x,(tuple,list)):return [safe(v) for v in x]
            return x
        body=(json.dumps(safe(value),allow_nan=False,indent=2)+'\n').encode()
        require(sum(f.stat().st_size for f in out.rglob('*') if f.is_file())+len(body)<CAP-65536,'Output cap exceeded')
        with (out/(name+'.json')).open('xb') as f:f.write(body);f.flush();os.fsync(f.fileno())
    save('protocol',dict(source_and_input_sha256={str(p):h for p,h in pins.items()},
        environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__),
        seed=20260923,rows_per_gene=COUNTS,required_points=15,requested_gene=args.gene,
        requested_tau_multiple=args.tau_multiple,scalar_reference_reused=True,
        GL_rules=gl.inspect_metadata(),executed_GL_settings=asdict(gl_settings),full_drawn_count_vector_sha256=hashlib.sha256(
            np.asarray([r.a for r in drawn],dtype='<i8').tobytes()).hexdigest(),
        adaptive_GL_order_with_unchanged_precision_and_caps=True,requested_node_batch_size=64,
        separate_G5_fixture_context=True,no_G8_scope_expansion=True,
        scalar_error_estimates_not_certified=True,no_fits=True,no_calibration=True,no_human_data=True))
    records=[]
    for reference in selected:
        gene=reference['gene'];rs=[r for r in drawn if r.gene==gene];model=original.Model(rs,gene_count=5)
        require(model.G==5 and model.Z==1 and model.size==16 and tuple(model.rows)==tuple(rs)
            and all(hasattr(r,'__dict__') for r in rs),'Exact G5 row/parameter context differs')
        result=dict(gene=gene,tau=reference['tau'],variance=reference['variance'],scalar=reference['scalar'],orders=[])
        for order in (20,40):
            start=time.perf_counter();before=p.tobytes()
            try:
                original.bb_terms=gamma.stable_bb_terms
                ll,g=model.marginal(p,reference['variance'],nodes=order)
            except (ValueError,ArithmeticError,np.linalg.LinAlgError) as exc:
                record=dict(order=order,returned=False,error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
            else:
                require(np.isfinite(ll) and g.shape==(17,) and np.isfinite(g).all(),'Returned arrays invalid')
                delta=float(ll-reference['scalar']['expanded']['logmass'])
                record=dict(order=order,returned=True,logmass=float(ll),gradient=g.tolist(),scalar_difference=delta,
                    scalar_relative_LL_pass=abs(delta)<=1e-8*abs(reference['scalar']['expanded']['logmass']),
                    original_GH_LL_difference=float(ll-reference['original_GH'+str(order)]['logmass']))
            finally:original.bb_terms=original_bb
            require(p.tobytes()==before and original.bb_terms is original_bb,'Parameter/BB changed')
            record['seconds']=time.perf_counter()-start;result['orders'].append(record)
        returned=all(r['returned'] for r in result['orders'])
        result['GH20_GH40_difference']=result['orders'][0]['logmass']-result['orders'][1]['logmass'] if returned else None
        result['GH20_GH40_absolute_pass']=abs(result['GH20_GH40_difference'])<=1e-8 if returned else None
        save('native_GH',result)
        Context=prepared.context_class(archive);probe=Context(model,p,rs)
        row_snapshot=snapshots.row_signatures(rs,np);component_maxima=dict(logmass=0.,gradient=0.,hessian=0.)
        try:
            original.bb_terms=gamma.stable_bb_terms
            for row in rs:
                pr=[model.predictors(row,p,u,need_hessian=False) for u in (0.,reference['tau'],-reference['tau'])]
                require(all(x[1]==pr[0][1] for x in pr),'Fixed dispersion varied across probe nodes')
                etas=np.asarray([x[0] for x in pr])
                for flag in (False,True):
                    actual=probe.callback(row,etas,pr[0][1],need_hessian=flag)
                    expected=[original.conditional_terms(row,float(eta),pr[0][1],need_hessian=flag) for eta in etas]
                    arrays=(np.asarray([x[0] for x in expected]),np.stack([x[1] for x in expected]),
                        np.stack([x[2] for x in expected]) if flag else None)
                    for key,a,b in zip(('logmass','gradient','hessian'),actual,arrays):
                        if key=='hessian' and not flag:continue
                        require(a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all(),
                            'Conditional probe shape/finiteness differs')
                        component_maxima[key]=max(component_maxima[key],float(np.max(np.abs(a-b))))
        finally:original.bb_terms=original_bb
        save('conditional_probe',dict(rows=len(rs),nodes_per_row=3,Hessian_modes=2,
            comparisons=len(rs)*3*2,maxima=component_maxima,Gamma_BB_batch_calls=probe.gamma_BB_calls))
        require(component_maxima['logmass']<=1e-11 and component_maxima['gradient']<=1e-11
            and component_maxima['hessian']<=1e-10,'G5 batched conditional differs from scalar Gamma')
        require(snapshots.row_signatures(rs,np)==row_snapshot,'G5 conditional probe changed rows')
        context=Context(model,p,rs)
        print(json.dumps(dict(event='GL_started',gene=gene,rows=len(rs),tau=reference['tau'])),flush=True)
        start=time.perf_counter();before=p.tobytes();row_snapshot=tuple(tuple(vars(r).items()) for r in rs)
        try:
            ll,g,h,report=gl.integrate_gene(model,p,reference['variance'],rs,context.callback,
                with_hessian=True,expected_kernel_sha=pins[HERE/'repair_p3_integral.py'],
                kernel_path=archive/'repair_p3_integral.py',settings=gl_settings,batch_size=64)
        except Exception as exc:
            checked_kernel=sys.modules.get('local_gl_reviewed_repair_p3_integral')
            known=(checked_kernel is not None and isinstance(exc,checked_kernel.IntegrationFailure)
                and str(exc) in ('Composite GL evaluation cap exhausted',
                    'GL order refinement exceeds unchanged precision allowance',
                    'Composite GL relative-tail domain cap exhausted','Composite GL expanded domain exceeds cap',
                    'Composite GL tail allowances remain unresolved','Composite GL final domain exceeds cap',
                    'Composite GL expanded domain/grid refinement disagrees','Nonfinite derivative moments','Invalid GL integral'))
            frames=[];tb=exc.__traceback__
            while tb is not None:
                frames.append(dict(file=tb.tb_frame.f_code.co_filename,line=tb.tb_lineno));tb=tb.tb_next
            known=known or (type(exc) is ValueError and str(exc) in ('Nonfinite rising product','Nonfinite candidate BB result')
                and any(Path(f['file']).resolve()==(archive/'repair_p3_bb_gamma_difference_v3.py').resolve() for f in frames))
            if not known:raise
            result['GL']=dict(returned=False,known_numerical_refusal=True,error_type=type(exc).__name__,
                message=str(exc),traceback=traceback.format_exc(),
                integration_reports=getattr(exc,'integration_reports',None),
                likelihood_evaluations=getattr(exc,'likelihood_evaluations',None))
        else:
            save('GL_returned_raw',dict(logmass=ll,gradient=g,hessian=h,report=report))
            require(np.isfinite(ll) and g.shape==(17,) and h.shape==(17,17)
                and np.isfinite(g).all() and np.isfinite(h).all(),'GL returned arrays invalid')
            delta=float(ll-reference['scalar']['expanded']['logmass'])
            result['GL']=dict(returned=True,logmass=float(ll),gradient=g,hessian=h,report=report,
                scalar_difference=delta,
                scalar_relative_LL_pass=abs(delta)<=1e-8*abs(reference['scalar']['expanded']['logmass']))
        finally:original.bb_terms=original_bb
        require(p.tobytes()==before and original.bb_terms is original_bb
            and tuple(tuple(vars(r).items()) for r in rs)==row_snapshot,'GL parameter/BB/rows changed')
        require(context.gamma_BB_calls>0,'GL did not execute batched Gamma BB arithmetic')
        result['GL']['Gamma_BB_batch_calls']=context.gamma_BB_calls
        result['GL']['conditional_chunk_ledger']=context.chunk_ledger()
        ledger=result['GL']['conditional_chunk_ledger']
        require(ledger and all(0<=x['row_ordinal']<len(rs) and x['calls']>0
            and x['need_hessian'] is True and x['all_requested_nodes_returned'] is True
            and 0<x['requested_nodes']<=64
            and x['effective_chunk_limit']==min(x['requested_nodes'],100000//(rs[x['row_ordinal']].n+1))
            and sum(x['chunk_lengths'])==x['requested_nodes']
            and all(0<n<=x['effective_chunk_limit'] for n in x['chunk_lengths']) for x in ledger),
            'G5 GL chunk dimensions or coverage differ')
        require(context.memory_estimate['prepared_count_axis_entries']<=100000
            and context.memory_estimate['estimated_maximum_batch_working_bytes']<=256*1024**2,
            'G5 conditional memory allowance exceeded')
        result['GL']['seconds']=time.perf_counter()-start
        save('gene'+str(gene)+'_tau'+str(reference['tau']),result);records.append(result)
        print(json.dumps(dict(gene=gene,tau=reference['tau'],completed_points=len(records))),flush=True)
        require(result['GL']['returned'],'G5 batch probe did not return an integral')
        if gene==1 and reference['tau']==.05:
            differences=dict(LL=abs(ll-dense['logmass']),gradient=float(np.max(np.abs(g-dense['gradient']))),
                Hessian=float(np.max(np.abs(h-dense['hessian']))))
            save('comparison_to_dense',dict(differences=differences,full_gradient_dimension=17,
                full_Hessian_dimension=[17,17],context_memory_estimate=context.memory_estimate))
            require(differences['LL']<=1e-9 and differences['gradient']<=1e-8 and differences['Hessian']<=1e-5,
                'G5 batched GL disagrees with completed dense GL')
        inactive=[i for i in range(10) if i not in (gene,5+gene)]
        require(np.all(g[inactive]==0) and np.all(h[inactive,:]==0) and np.all(h[:,inactive]==0),
            'G5 inactive derivatives differ from exact zero')
    for path,h in pins.items():require(sha(path)==h,'Input changed during evaluation')
    save('summary',dict(completed=True,points=len(records),returned_native_integrals=sum(r['returned'] for x in records for r in x['orders']),
        GL_returned=records[0]['GL']['returned'],
        requested_gene=args.gene,requested_tau_multiple=args.tau_multiple,
        GL_known_numerical_refusal=records[0]['GL'].get('known_numerical_refusal',False),
        GL_only_scalar_relative_LL_pass=records[0]['GL'].get('scalar_relative_LL_pass'),
        scalar_relative_LL_rule='abs(candidate-reference)<=1e-8*abs(reference)',
        required_family_points=15,selected_point_likelihood_criteria_passed=all(x['GL'].get('scalar_relative_LL_pass') is True
            and x['GH20_GH40_absolute_pass'] is True
            and all(r.get('scalar_relative_LL_pass') is True for r in x['orders']) for x in records),
        no_independent_gradient_or_Hessian_acceptance=True,no_fits=True,no_calibration=True,no_human_data=True))


if __name__=='__main__':main()
