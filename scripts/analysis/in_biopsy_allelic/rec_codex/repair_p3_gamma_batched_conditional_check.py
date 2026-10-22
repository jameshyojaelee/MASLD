#!/usr/bin/env python3
"""Synthetic fully batched Gamma conditional mixtures and local GL moments.

Compare every conditional call with its scalar counterpart before attempting
the saved failed point. No fit, calibration or independent Hessian claim.
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
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
POINT=REC/'repair_p3_observed_failure_anchor_22042612'
THETA=REC/'repair_p3_theta4_anchored_gl_gene2_22042347'
BASELINE=REC/'repair_p3_gamma_observed_point_GL_22076126'
PRIMITIVE=REC/'repair_p3_bb_gamma_batched_check_22078120'
COUNT_SHA='7009be4adeec6157fa4bbf3f46e48ec24aa398a0bee4d7906c78c231b7d3fbad'
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
         'repair_p3_bb_gamma_difference_v3.py','repair_p3_integral.py',
         'repair_p3_local_gl.py','repair_p3_local_moments.py','repair_p3_gamma_parallel_model.py',
         'repair_p3_batched_conditional.py','repair_p3_prepared_conditional_v2.py',
         'repair_p3_bb_gamma_batched.py','repair_p3_gamma_batched_conditional.py')
CAP=32*1024**2


def require(ok,message):
    if not ok:
        raise RuntimeError(message)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def load(name,p):
    spec=importlib.util.spec_from_file_location(name,p)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pins',type=Path,required=True)
    parser.add_argument('--expected-pins-sha',required=True)
    parser.add_argument('--expected-checker-sha',required=True)
    parser.add_argument('--inspect',action='store_true')
    args=parser.parse_args()
    require(sha(Path(__file__))==args.expected_checker_sha and sha(args.pins)==args.expected_pins_sha,
            'Executable/manifest differs')
    pins={Path(p):h for p,h in json.loads(args.pins.read_text()).items()}
    for p,h in pins.items():
        require(p.is_absolute() and p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==h,
                'Pinned source/input differs: '+str(p))
        if p.suffix=='.py':ast.parse(p.read_text())
    require(all(HERE/n in pins for n in SOURCES) and all(p in pins for p in
        (POINT/'protocol.json',THETA/'protocol.json',BASELINE/'protocol.json',BASELINE/'summary.json',
         BASELINE/'GL_candidate.json')),'Required input absent')
    require(PRIMITIVE/'summary.json' in pins and PRIMITIVE/'protocol.json' in pins,
            'Batched primitive admission absent')
    primitive=json.loads((PRIMITIVE/'summary.json').read_text())
    primitive_protocol=json.loads((PRIMITIVE/'protocol.json').read_text())
    require(all(primitive[k] is True for k in ('completed','independent_rising_passed',
            'independent_selected21_passed','all_scalar_family_passed','oversized_batch_refused'))
            and primitive_protocol['source_and_input_sha256'][str(HERE/'repair_p3_bb_gamma_batched.py')]
            ==pins[HERE/'repair_p3_bb_gamma_batched.py'],'Admitted batched primitive differs')
    point=json.loads((POINT/'protocol.json').read_text())
    theta=json.loads((THETA/'protocol.json').read_text())
    admitted=json.loads((BASELINE/'summary.json').read_text())
    source=json.loads((BASELINE/'protocol.json').read_text())
    reference=json.loads((BASELINE/'GL_candidate.json').read_text())
    require(point['count_vector_sha256']==COUNT_SHA and point['call']==12
            and point['fixed_v']==.0009636462009520182 and theta['fixed_v']==.0004,
            'Prespecified synthetic points differ')
    require(admitted['completed'] is True and admitted['whole_saved_score_integral_nominal_criteria_met'] is True
            and source['source_and_input_sha256'][str(HERE/'repair_p3_bb_gamma_difference_v3.py')]
            ==pins[HERE/'repair_p3_bb_gamma_difference_v3.py'],'Gamma whole-score prerequisite differs')
    require(source['gene']==2 and source['count_vector_sha256']==COUNT_SHA
            and source['fixed_p']==point['fixed_p'] and source['fixed_v']==point['fixed_v']
            and len(reference['gradient'])==29
            and all(math.isfinite(x) for x in reference['gradient'])
            and math.isfinite(reference['logmass']),'Whole score reference identity/dimensions differ')
    protocol=dict(source_and_input_sha256={str(p):h for p,h in pins.items()},
        synthetic_participants=60,rows=173,G=8,Z=7,P=28,count_vector_sha256=COUNT_SHA,
        component_vectors=['theta4','saved_call12'],component_u=['0','+sqrt(v)','-sqrt(v)'],
        component_comparisons=2076,component_byte_identity_required=False,
        component_LP_atol=1e-9,component_score_atol=1e-8,component_Hessian_atol=1e-5,
        wide_eta_points=64,wide_eta_range=[-30.,30.],wide_comparisons=22144,
        whole_gene=2,whole_parameters=point['fixed_p'],whole_variance=point['fixed_v'],
        whole_modes=['score_only','Hessian'],batch_size=32,Hessian_tolerance=1e-5/8,
        score_baseline_LL_atol=1e-9,score_baseline_gradient_atol=1e-8,
        unchanged_local_GL=True,batched_Gamma_native_conditional=True,
        no_fit=True,no_calibration=True,no_independent_full_Hessian=True,no_actual_or_sealed_data=True,
        output_cap_bytes=CAP)
    if args.inspect:print(json.dumps(protocol));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    out=REC/('repair_p3_gamma_batched_conditional_check_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in SOURCES)):
        dest=archive/p.name;dest.write_bytes(p.read_bytes())
        require(sha(dest)==(args.expected_checker_sha if p==Path(__file__) else pins[p]),'Executable archive differs')
    began=time.perf_counter()
    import numpy as np
    import scipy
    def serialize(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    def refusal_json(x):
        if isinstance(x,np.ndarray):return refusal_json(x.tolist())
        if isinstance(x,np.generic):return refusal_json(x.item())
        if isinstance(x,float) and not math.isfinite(x):return {'nonfinite_float':x.hex()}
        if isinstance(x,dict):return {k:refusal_json(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [refusal_json(v) for v in x]
        return x
    def save(name,value):
        body=(json.dumps(value,default=serialize,allow_nan=False,indent=2)+'\n').encode()
        require(sum(p.stat().st_size for p in out.rglob('*') if p.is_file())+len(body)<CAP-65536,
                'Output cap exceeded')
        with (out/(name+'.json')).open('xb') as f:f.write(body);f.flush();os.fsync(f.fileno())
        print(json.dumps(dict(saved=name,seconds=time.perf_counter()-began)),flush=True)
    try:
        require(dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
                ==point['environment'],'Numerical environment differs')
        save('protocol',protocol)
        sys.path.insert(0,str(archive))
        loader=importlib.import_module('bfix_v1');importlib.import_module('tests_v1')
        for n in SOURCES[:4]:
            module=sys.modules[n[:-3]]
            require(Path(module.__file__).resolve().parent==archive.resolve()
                    and sha(Path(module.__file__))==pins[HERE/n],'Original import escaped archive')
        design,settings=loader.load()
        private=load('gamma_local_original',archive/'likelihood_v1.py')
        bb=load('gamma_local_bb',archive/'repair_p3_bb_gamma_difference_v3.py')
        kernel=load('gamma_local_kernel',archive/'repair_p3_integral.py')
        local=load('gamma_local_GL',archive/'repair_p3_local_gl.py')
        state=load('gamma_local_state',archive/'repair_p3_gamma_parallel_model.py')
        model=private.Model(design.model().rows,gene_count=8)
        require(len(model.rows)==173 and model.size==28 and model.Z==7
                and len(design.participants)==60 and settings['seed']==20260930,'Synthetic design differs')
        counts=lambda:hashlib.sha256(np.asarray([r.a for r in model.rows],dtype='<i8').tobytes()).hexdigest()
        require(counts()==COUNT_SHA,'Synthetic counts differ')
        initial=state.row_signatures(model.rows,np)
        def row_array_flags(rows):
            return tuple((i,key,value.dtype.str,value.shape,bool(value.flags.writeable))
                         for i,row in enumerate(rows) for key,value in sorted(vars(row).items())
                         if isinstance(value,np.ndarray))
        initial_array_flags=row_array_flags(model.rows)
        require(all(fn.__globals__ is vars(private) for fn in
            (private.conditional_terms,private.het_inclusion,private.mixture_terms))
            and all('bb_terms' in fn.__code__.co_names for fn in (private.het_inclusion,private.mixture_terms)),
            'Native conditional BB lookup differs')
        original_bb=private.bb_terms;calls=[0]
        def counted(*a,**kw):
            calls[0]+=1
            return bb.stable_bb_terms(*a,**kw)
        private.bb_terms=counted
        assembly=load('gamma_batched_conditional',archive/'repair_p3_gamma_batched_conditional.py')
        Context=assembly.context_class(archive)
        contexts=[]
        def callback_for(p,rows):
            context=Context(model,p,rows);contexts.append(context)
            return context.callback
        records=[]
        try:
            for label,parent in (('theta4',theta),('saved_call12',point)):
                p=np.asarray(parent['fixed_p']);v=parent['fixed_v'];callback=callback_for(p,model.rows)
                for row in model.rows:
                    predictors=[model.predictors(row,p,u,False) for u in (0.,np.sqrt(v),-np.sqrt(v))]
                    etas=np.asarray([x[0] for x in predictors]);r=predictors[0][1]
                    require(all(x[1]==r for x in predictors),'Dispersion varied with u')
                    for with_hessian in (False,True):
                        actual=callback(row,etas,r,with_hessian)
                        expected=[private.conditional_terms(row,float(eta),r,need_hessian=with_hessian) for eta in etas]
                        delta=[float(np.max(np.abs(actual[h]-np.stack([x[h] for x in expected]))))
                               for h in range(3 if with_hessian else 2)]
                        exact=all(np.asarray(actual[h]).tobytes()==np.stack([x[h] for x in expected]).tobytes()
                                  for h in range(3 if with_hessian else 2))
                        passed=delta[0]<=1e-9 and delta[1]<=1e-8 and (not with_hessian or delta[2]<=1e-5)
                        records.append(dict(vector=label,individual=row.individual,gene=row.gene,
                            with_hessian=with_hessian,points=3,byte_identity=exact,max_difference=delta,passed=bool(passed)))
            require(sum(r['points'] for r in records)==2076 and calls[0]>0,'Component family incomplete')
            save('components',dict(records=records,comparisons=2076,native_Gamma_BB_calls=calls[0],
                batched_Gamma_BB_calls=sum(c.gamma_BB_calls for c in contexts),
                all_nominal_criteria_met=all(x['passed'] for x in records)))
            require(all(x['passed'] for x in records) and all(c.gamma_BB_calls>0 for c in contexts),
                    'Conditional component criteria or backend call proof failed')
            wide=[];p=np.asarray(point['fixed_p']);callback=callback_for(p,model.rows)
            etas=np.linspace(-30.,30.,64)
            for row in model.rows:
                r=model.predictors(row,p,0.,False)[1]
                for flag in (False,True):
                    actual=callback(row,etas,r,flag)
                    expected=[private.conditional_terms(row,float(eta),r,need_hessian=flag) for eta in etas]
                    delta=[float(np.max(np.abs(actual[h]-np.stack([x[h] for x in expected]))))
                           for h in range(3 if flag else 2)]
                    passed=delta[0]<=1e-9 and delta[1]<=1e-8 and (not flag or delta[2]<=1e-5)
                    wide.append(dict(individual=row.individual,gene=row.gene,Hessian=flag,
                                     points=64,max_difference=delta,passed=bool(passed)))
            save('wide_components',dict(records=wide,comparisons=sum(x['points'] for x in wide),
                support_calls=contexts[-1].support_calls,chunk_ledger=contexts[-1].chunk_ledger()))
            require(sum(x['points'] for x in wide)==22144 and all(x['passed'] for x in wide),
                    'Wide-node conditional family failed')
            require(all(contexts[-1].support_calls[k]>0 for k in ('observed','tails','selected')),
                    'Not all native integer-selection support branches were exercised')
            require(any(len(x['chunk_lengths'])>1 for x in contexts[-1].chunk_ledger())
                    and all(x['all_requested_nodes_returned'] for x in contexts[-1].chunk_ledger()),
                    'Bounded node splitting not exercised or nodes dropped')
            from dataclasses import replace
            probe_p=np.asarray(point['fixed_p']).copy()
            probe_row=replace(model.rows[0],z=np.asarray(model.rows[0].z,dtype=float).copy())
            clone_array_flags=row_array_flags((probe_row,))
            probe=Context(model,probe_p,(probe_row,))
            require(row_array_flags((probe_row,))==clone_array_flags,
                    'Preparing constants changed caller row-array flags')
            eta,r,*_=model.predictors(probe_row,probe_p,0.,False)
            def accepted_probe():
                actual=probe.callback(probe_row,np.asarray([eta]),r,False)
                expected=private.conditional_terms(probe_row,eta,r,need_hessian=False)
                require(np.isfinite(actual[0]).all() and np.isfinite(actual[1]).all()
                        and abs(actual[0][0]-expected[0])<=1e-9
                        and np.max(np.abs(actual[1][0]-expected[1]))<=1e-8,
                        'Unmodified/restored probe context did not return a valid conditional')
            accepted_probe()
            original_alpha=probe_p[0];probe_p[0]+=1.
            try:probe.callback(probe_row,np.asarray([eta]),r,False)
            except RuntimeError as exc:parameter_refused='Fixed parameters changed' in str(exc)
            else:parameter_refused=False
            probe_p[0]=original_alpha
            accepted_probe()
            d=probe.prepared.rows[id(probe_row)];original_hom=d['hom'];d['hom']=float(original_hom)+1.
            try:probe.callback(probe_row,np.asarray([eta]),r,False)
            except RuntimeError as exc:constant_refused='Prepared constants changed' in str(exc)
            else:constant_refused=False
            d['hom']=original_hom;accepted_probe()
            original_z=probe_row.z[0];probe_row.z[0]+=1.
            try:probe.callback(probe_row,np.asarray([eta]),r,False)
            except RuntimeError as exc:row_refused='Fixed row changed' in str(exc)
            else:row_refused=False
            probe_row.z[0]=original_z;accepted_probe()
            save('state_probes',dict(changed_parameters_refused=parameter_refused,
                changed_derived_constants_refused=constant_refused,changed_row_array_refused=row_refused,
                input_row_array_flags_preserved=row_array_flags((probe_row,))==clone_array_flags,
                accepted_unmodified_and_after_all_restores=True,only_cloned_synthetic_row_mutated=True))
            require(parameter_refused and constant_refused and row_refused,'Stale context accepted changed inputs')
            p=np.asarray(point['fixed_p']);v=point['fixed_v'];rows=tuple(r for r in model.rows if r.gene==2)
            callback=callback_for(p,rows)
            from dataclasses import replace
            per_gene=replace(kernel.Settings(),hessian_tolerance=1e-5/8)
            outcomes=[]
            for flag in (False,True):
                start=time.perf_counter();before=contexts[-1].gamma_BB_calls
                try:
                    ll,g,h,report=local.integrate_gene(model,p,v,rows,callback,with_hessian=flag,
                        expected_kernel_sha=state.KERNEL_SHA,kernel_path=archive/'repair_p3_integral.py',
                        settings=per_gene,batch_size=32)
                    result=dict(returned=True,logmass=ll,gradient=g,hessian=h,report=report)
                    require(np.isfinite(ll) and g.shape==(29,) and np.isfinite(g).all()
                            and (not flag or (h.shape==(29,29) and np.isfinite(h).all())),
                            'Returned numerical values are not finite with original dimensions')
                    # Preserve rejected-ladder diagnostics without sanitizing
                    # the accepted likelihood, scores or Hessian.
                    result['report']=refusal_json(report)
                    save('Hessian_integral' if flag else 'score_only_integral',result)
                    if not flag:
                        result['reference_keys']=list(reference)
                        result['LL_difference']=ll-reference['logmass']
                        result['score_max_difference']=float(np.max(np.abs(g-np.asarray(reference['gradient']))))
                        result['nominal_score_comparison_passed']=bool(abs(result['LL_difference'])<=1e-9
                            and result['score_max_difference']<=1e-8)
                except Exception as exc:
                    numerical=(isinstance(exc,ValueError) and hasattr(exc,'integration_reports')
                               and hasattr(exc,'likelihood_evaluations'))
                    result=dict(returned=False,error_type=type(exc).__name__,message=str(exc),
                        numerical_refusal=numerical,traceback=traceback.format_exc(),
                        refused_gene_ladders=getattr(exc,'integration_reports',None),
                        refused_gene_evaluations=getattr(exc,'likelihood_evaluations',None))
                    if not numerical:
                        save('integrity_failure',refusal_json(result));raise
                result.update(with_hessian=flag,seconds=time.perf_counter()-start,gamma_BB_calls=contexts[-1].gamma_BB_calls-before)
                save('Hessian' if flag else 'score_only',result if result['returned'] else refusal_json(result))
                outcomes.append(result)
        finally:
            private.bb_terms=original_bb
        require(counts()==COUNT_SHA and state.row_signatures(model.rows,np)==initial
                and row_array_flags(model.rows)==initial_array_flags
                and private.bb_terms is original_bb,'Fixed inputs/parent callable changed')
        for p,h in pins.items():require(sha(p)==h,'Source/input changed while running')
        save('context_ledgers',[dict(context=i,support_calls=c.support_calls,chunk_ledger=c.chunk_ledger())
                                for i,c in enumerate(contexts)])
        save('summary',dict(completed=True,component_comparisons=2076,all_component_nominal_criteria_met=True,
            wide_comparisons=22144,all_wide_nominal_criteria_met=True,
            score_returned=outcomes[0]['returned'],Hessian_returned=outcomes[1]['returned'],
            nominal_score_comparison_passed=outcomes[0].get('nominal_score_comparison_passed',False),
            elapsed_seconds=time.perf_counter()-began,no_fit=True,no_calibration=True,
            no_independent_full_Hessian=True))
    except Exception as exc:
        save('failure',refusal_json(dict(error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())))
        raise


if __name__=='__main__':main()
