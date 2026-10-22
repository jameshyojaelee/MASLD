#!/usr/bin/env python3
"""Synthetic native-Gamma conditional batches with existing local GL moments.

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
COUNT_SHA='7009be4adeec6157fa4bbf3f46e48ec24aa398a0bee4d7906c78c231b7d3fbad'
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
         'repair_p3_bb_gamma_difference_v3.py','repair_p3_integral.py',
         'repair_p3_local_gl.py','repair_p3_local_moments.py','repair_p3_gamma_parallel_model.py')
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
        component_comparisons=2076,component_byte_identity_required=True,
        whole_gene=2,whole_parameters=point['fixed_p'],whole_variance=point['fixed_v'],
        whole_modes=['score_only','Hessian'],batch_size=32,Hessian_tolerance=1e-5/8,
        score_baseline_LL_atol=1e-9,score_baseline_gradient_atol=1e-8,
        unchanged_local_GL=True,scalar_native_conditional_per_node=True,
        no_fit=True,no_calibration=True,no_independent_full_Hessian=True,no_actual_or_sealed_data=True,
        output_cap_bytes=CAP)
    if args.inspect:print(json.dumps(protocol));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    out=REC/('repair_p3_gamma_local_check_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
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
        require(all(fn.__globals__ is vars(private) for fn in
            (private.conditional_terms,private.het_inclusion,private.mixture_terms))
            and all('bb_terms' in fn.__code__.co_names for fn in (private.het_inclusion,private.mixture_terms)),
            'Native conditional BB lookup differs')
        original_bb=private.bb_terms;calls=[0]
        def counted(*a,**kw):
            calls[0]+=1
            return bb.stable_bb_terms(*a,**kw)
        private.bb_terms=counted
        def callback_for(p,rows):
            fixed=np.asarray(p).copy();fixed.setflags(write=False)
            contexts={id(r):(r,state.row_signatures((r,),np),float(model.predictors(r,fixed,0.,False)[1]))
                      for r in rows}
            def callback(row,etas,dispersion,need_hessian=True):
                context=contexts.get(id(row))
                require(context is not None and context[0] is row
                        and state.row_signatures((row,),np)==context[1]
                        and float(dispersion).hex()==context[2].hex(),'Fixed row/dispersion changed')
                eta=np.asarray(etas)
                require(eta.ndim==1 and 0<len(eta)<=64 and np.isfinite(eta).all()
                        and row.n<=5000,'Bounded finite node batch required')
                answers=[private.conditional_terms(row,float(x),dispersion,need_hessian=need_hessian) for x in eta]
                return (np.asarray([x[0] for x in answers]),np.stack([x[1] for x in answers]),
                        np.stack([x[2] for x in answers]) if need_hessian else None)
            return callback
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
                        for j,answer in enumerate(expected):
                            require(all(np.asarray(a[j]).dtype==np.asarray(b).dtype
                                and np.asarray(a[j]).shape==np.asarray(b).shape
                                and np.asarray(a[j]).tobytes()==np.asarray(b).tobytes()
                                for a,b in zip(actual[:3 if with_hessian else 2],answer)),
                                'Batched callback changed native conditional bytes')
                        records.append(dict(vector=label,individual=row.individual,gene=row.gene,
                                            with_hessian=with_hessian,points=3,byte_identity=True))
            require(sum(r['points'] for r in records)==2076 and calls[0]>0,'Component family incomplete')
            save('components',dict(records=records,comparisons=2076,gamma_BB_calls=calls[0],all_byte_identical=True))
            p=np.asarray(point['fixed_p']);v=point['fixed_v'];rows=tuple(r for r in model.rows if r.gene==2)
            callback=callback_for(p,rows)
            from dataclasses import replace
            per_gene=replace(kernel.Settings(),hessian_tolerance=1e-5/8)
            outcomes=[]
            for flag in (False,True):
                start=time.perf_counter();before=calls[0]
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
                result.update(with_hessian=flag,seconds=time.perf_counter()-start,gamma_BB_calls=calls[0]-before)
                save('Hessian' if flag else 'score_only',result if result['returned'] else refusal_json(result))
                outcomes.append(result)
        finally:
            private.bb_terms=original_bb
        require(counts()==COUNT_SHA and state.row_signatures(model.rows,np)==initial
                and private.bb_terms is original_bb,'Fixed inputs/parent callable changed')
        for p,h in pins.items():require(sha(p)==h,'Source/input changed while running')
        save('summary',dict(completed=True,component_comparisons=2076,all_component_bytes_identical=True,
            score_returned=outcomes[0]['returned'],Hessian_returned=outcomes[1]['returned'],
            nominal_score_comparison_passed=outcomes[0].get('nominal_score_comparison_passed',False),
            elapsed_seconds=time.perf_counter()-began,no_fit=True,no_calibration=True,
            no_independent_full_Hessian=True))
    except Exception as exc:
        save('failure',refusal_json(dict(error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())))
        raise


if __name__=='__main__':main()
