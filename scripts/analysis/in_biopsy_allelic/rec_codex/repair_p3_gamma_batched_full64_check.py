#!/usr/bin/env python3
"""Compare64-node batching with the checked32-node full eight-gene result.

No optimizer, calibration, independent full Hessian or human-data claim.
"""
import argparse
import ast
import copy
import hashlib
import importlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import resource
import sys
import time
import traceback

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
REFERENCE=REC/'repair_p3_gamma_batched_full_check_22078827'
PAIR=REC/'repair_p3_gamma_batch64_pair_review_20261007T060216Z/review.json'
ADMISSION=REC/'repair_p3_gamma_batched_conditional_check_22078218'
POINT=REC/'repair_p3_observed_failure_anchor_22042612'
COUNT_SHA='7009be4adeec6157fa4bbf3f46e48ec24aa398a0bee4d7906c78c231b7d3fbad'
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
 'repair_p3_integral.py','repair_p3_local_gl.py','repair_p3_local_moments.py',
 'repair_p3_gamma_parallel_model.py','repair_p3_bb_gamma_difference_v3.py',
 'repair_p3_batched_conditional.py','repair_p3_prepared_conditional_v2.py',
 'repair_p3_bb_gamma_batched.py','repair_p3_gamma_batched_conditional.py',
 'repair_p3_gamma_batched_model.py','repair_p3_gamma_batched_model64.py')
CAP=32*1024**2


def require(ok,message):
    if not ok:raise RuntimeError(message)


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
            'Checker/manifest differs')
    pins={Path(p):h for p,h in json.loads(args.pins.read_text()).items()}
    for p,h in pins.items():
        require(p.is_absolute() and p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==h,
                'Pinned source/input differs: '+str(p))
        if p.suffix=='.py':ast.parse(p.read_text())
    required=[*(HERE/n for n in SOURCES),POINT/'protocol.json',PAIR,
        *(REFERENCE/n for n in ('summary.json','protocol.json','answer.json','workers.json')),
        *(ADMISSION/n for n in ('summary.json','protocol.json','Hessian.json','state_probes.json'))]
    require(all(p in pins for p in required),'Required input absent')
    read=lambda p:json.loads(p.read_text())
    pair=read(PAIR)
    require(pair['both_fixed_point_numerical_returns_verified'] and len(pair['records'])==2
        and {r['case'] for r in pair['records']}=={'failed0gene2','positive18'}
        and all(r['nominal_fixed_point_checks_passed'] and r['all24_rows_obey_capacity_rule']
            and r['node_path_and_settings_identical'] for r in pair['records']),
        'Both supported64-node fixed points required')
    point=read(POINT/'protocol.json');refp=read(REFERENCE/'protocol.json')
    refs=read(REFERENCE/'summary.json');ref=read(REFERENCE/'answer.json')
    refworkers=read(REFERENCE/'workers.json');admit=read(ADMISSION/'summary.json')
    admitp=read(ADMISSION/'protocol.json')
    require(point['count_vector_sha256']==COUNT_SHA and point['call']==12
        and point['fixed_v']==.0009636462009520182
        and refp['parameters']==point['fixed_p'] and refp['variance']==point['fixed_v']
        and refp['count_vector_sha256']==COUNT_SHA and refs['completed'] and refs['returned']
        and refs['all_fixed_inputs_unchanged'] and ref['returned']
        and [x['gene'] for x in refworkers]==list(range(8))
        and all(x['success'] and x['present'] for x in refworkers),'Reference incomplete or different')
    require(admit['completed'] and admit['score_returned'] and admit['Hessian_returned']
        and admit['all_component_nominal_criteria_met'] and admit['all_wide_nominal_criteria_met']
        and all(read(ADMISSION/'state_probes.json').values())
        and admitp['whole_parameters']==point['fixed_p'] and admitp['whole_variance']==point['fixed_v']
        and admitp['count_vector_sha256']==COUNT_SHA,'Batched conditional admission incomplete')
    for name in SOURCES:
        path=HERE/name
        if str(path) in admitp['source_and_input_sha256']:
            require(admitp['source_and_input_sha256'][str(path)]==pins[path],
                    'Batched admitted source differs: '+name)
        if str(path) in refp['source_and_input_sha256']:
            require(refp['source_and_input_sha256'][str(path)]==pins[path],
                    'Scalar reference source differs: '+name)
    protocol=dict(source_and_input_sha256={str(p):h for p,h in pins.items()},
        parameters=point['fixed_p'],variance=point['fixed_v'],count_vector_sha256=COUNT_SHA,
        participants=60,rows=173,G=8,Z=7,P=28,fixture_seed=20260930,
        with_hessian=True,node_batch_size=64,reference_node_batch_size=32,comparison_LL_atol=1e-9,comparison_score_atol=1e-8,
        comparison_Hessian_atol=1e-5,all_eight_genes_required=True,
        global_LL_allowance=1e-8,global_score_allowance=1e-7,global_Hessian_allowance=1e-5,
        same_rule_geometry_required=True,local_errors_embedded_before_sum=True,
        no_fit=True,no_calibration=True,no_independent_full_Hessian=True,no_actual_or_sealed_outcomes=True)
    if args.inspect:print(json.dumps(protocol));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS',
            'MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),'Threads not fixed before imports')
    out=REC/('repair_p3_gamma_batched_full64_check_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in SOURCES)):
        dest=archive/p.name;dest.write_bytes(p.read_bytes())
        require(sha(dest)==(args.expected_checker_sha if p==Path(__file__) else pins[p]),'Archive differs')
    began=time.perf_counter()
    import numpy as np
    import scipy
    def serialize(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    def failed_json(x):
        if isinstance(x,np.ndarray):return failed_json(x.tolist())
        if isinstance(x,np.generic):return failed_json(x.item())
        if isinstance(x,float) and not math.isfinite(x):return {'nonfinite_float':x.hex()}
        if isinstance(x,dict):return {k:failed_json(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [failed_json(v) for v in x]
        return x
    def save(name,value):
        body=(json.dumps(value,default=serialize,allow_nan=False,indent=2)+'\n').encode()
        require(sum(p.stat().st_size for p in out.rglob('*') if p.is_file())+len(body)<CAP-65536,
                'Output cap exhausted')
        with (out/(name+'.json')).open('xb') as f:f.write(body);f.flush();os.fsync(f.fileno())
        print(json.dumps(dict(saved=name,seconds=time.perf_counter()-began)),flush=True)
    try:
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        require(environment==point['environment']==refp['environment'],'Numerical environment differs')
        protocol['environment']=environment;save('protocol',protocol)
        sys.path.insert(0,str(archive));loader=importlib.import_module('bfix_v1')
        importlib.import_module('tests_v1')
        for name in SOURCES[:4]:
            m=sys.modules[name[:-3]]
            require(Path(m.__file__).resolve().parent==archive.resolve() and sha(Path(m.__file__))==pins[HERE/name],
                    'Original import escaped archive')
        design,settings=loader.load()
        private=load('batched_full_original',archive/'likelihood_v1.py')
        kernel=load('batched_full_kernel',archive/'repair_p3_integral.py')
        gl=load('batched_full_gl',archive/'repair_p3_local_gl.py')
        conditional=load('batched_full_conditional',archive/'repair_p3_gamma_batched_conditional.py')
        wrapper=load('batched_full_wrapper',archive/'repair_p3_gamma_batched_model64.py')
        base=private.Model(design.model().rows,gene_count=8)
        require(len(base.rows)==173 and base.Z==7 and base.size==28 and len(design.participants)==60
                and settings['seed']==20260930,'Synthetic dimensions differ')
        counts=lambda:hashlib.sha256(np.asarray([r.a for r in base.rows],dtype='<i8').tobytes()).hexdigest()
        require(counts()==COUNT_SHA,'Synthetic count vector differs')
        p=np.asarray(point['fixed_p'],dtype=float);p.setflags(write=False)
        signatures=wrapper.row_signatures(base.rows,np);original_bb=private.bb_terms
        Model=wrapper.gamma_batched_model_class(private,kernel,gl,conditional,
            kernel_path=archive/'repair_p3_integral.py')
        model=Model(base.rows,gene_count=8)
        preflight={}
        for h in (False,True):
            old=base.marginal(p,0.,with_hessian=h);new=model.marginal(p,0.,with_hessian=h)
            preflight['boundary_'+str(h)]=all(np.asarray(a).tobytes()==np.asarray(b).tobytes()
                for a,b in zip(old,new))
        report=read(ADMISSION/'Hessian.json')['report']
        errors=wrapper.embed_report_errors(report,model,2,np,True)
        indices=report['augmented_global_indices'];inactive=[i for i in range(29) if i not in indices]
        preflight['embedding_inactive_exact_zero']=bool(np.all(errors['gradient'][inactive]==0)
            and np.all(errors['hessian'][inactive,:]==0) and np.all(errors['hessian'][:,inactive]==0))
        for label,edit in [('wrong_coordinates',lambda r:r['augmented_global_indices'].reverse()),
            ('negative_allowance',lambda r:r['refined']['gradient_error_allowance'].__setitem__(0,-1.)),
            ('wrong_shape',lambda r:r['refined']['hessian_error_allowance'].pop())]:
            bad=copy.deepcopy(report);edit(bad)
            try:wrapper.embed_report_errors(bad,model,2,np,True)
            except wrapper.ModelIntegrityError:preflight[label+'_refused']=True
            else:preflight[label+'_refused']=False
        preflight['inherited_methods_unchanged']=all(getattr(Model,n) is getattr(private.Model,n) for n in wrapper.INHERITED)
        try:kernel.require(False,'Nonfinite derivative moments')
        except Exception as exc:preflight['known_arithmetic_refusal_classified']=wrapper.numerical_worker_refusal(exc)
        try:kernel.require(False,'Actual-node predictor invariance failed')
        except Exception as exc:
            preflight['state_refusal_not_numeric']=not wrapper.numerical_worker_refusal(exc)
            exc.integration_reports=[];exc.likelihood_evaluations=0
            preflight['state_refusal_with_ladder_not_numeric']=not wrapper.numerical_worker_refusal(exc)
        preflight['unrelated_value_error_not_numeric']=not wrapper.numerical_worker_refusal(ValueError('unexpected'))
        numeric=dict(gene=1,success=False,numerical_integration_refusal=True,linalg_error=False)
        integrity=dict(gene=5,success=False,numerical_integration_refusal=False,linalg_error=False)
        preflight['later_integrity_failure_has_priority']=wrapper.first_failed_worker([numeric,integrity]) is integrity
        preflight['numeric_only_failure_retained']=wrapper.first_failed_worker([numeric]) is numeric
        preflight['all_success_no_failure']=wrapper.first_failed_worker([dict(success=True)]) is None
        require(all(preflight.values()),'Boundary or local error mapping failed');save('preflight',preflight)
        terms={};start=time.perf_counter()
        try:
            ll,g,h=model.marginal(p,float(point['fixed_v']),with_hessian=True,gene_terms=terms)
        except Exception as exc:
            save('workers',failed_json(getattr(model,'last_parallel_worker_records',[])))
            save('refusal',failed_json(dict(error_type=type(exc).__name__,message=str(exc),
                failure=getattr(model,'last_integration_failure',None),gene_terms=terms)))
            raise
        elapsed=time.perf_counter()-start;workers=model.last_parallel_worker_records
        answer=dict(returned=True,logmass=ll,gradient=g,hessian=h,gene_terms=terms,
            error_allowances=model.last_total_error_allowance,seconds=elapsed)
        save('workers',workers);save('answer',answer)
        require([x['gene'] for x in workers]==list(range(8)) and all(x['present'] and x['success'] for x in workers),
                'All eight genes required')
        def compare(actual,expected):
            a,b=np.asarray(actual),np.asarray(expected)
            require(a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all(),'Array shape/finiteness differs')
            return float(np.max(np.abs(a-b)))
        records=[]
        geometry=lambda r:[tuple(t[k] for k in ('z_bound','grid_spacing','panels','order','node_evaluations',
            'cumulative_evaluations')) for t in r['ladder']]
        for w,old in zip(workers,refworkers):
            gene=w['gene'];d=dict(gene=gene,LL=abs(w['logmass']-old['logmass']),
                score=compare(w['gradient'],old['gradient']),Hessian=compare(w['hessian'],old['hessian']),
                same_geometry=geometry(w['report'])==geometry(old['report']))
            require(w['report']['likelihood_evaluations']==old['report']['likelihood_evaluations'],
                    'Reference32 and64 logical node counts differ')
            rows=tuple(r for r in base.rows if r.gene==gene);ledger=w['conditional_chunk_ledger']
            require({x['row_ordinal'] for x in ledger}==set(range(len(rows))),'Worker row ledger differs')
            for entry in ledger:
                requested=entry['requested_nodes'];n=rows[entry['row_ordinal']].n
                effective=min(requested,100000//(n+1))
                require(0<requested<=64 and entry['effective_chunk_limit']==effective
                    and entry['chunk_lengths']==[min(effective,requested-start) for start in range(0,requested,effective)]
                    and entry['all_requested_nodes_returned'] and entry['need_hessian'],
                    'Worker violated unchanged count-axis chunk rule')
            require(all(any(x['row_ordinal']==i and x['requested_nodes']==64 for x in ledger)
                for i in range(len(rows))),'Requested64-node batch missing for a worker row')
            d['count_axis_chunk_rule_verified']=True
            d['Gamma_BB_calls']=w['conditional_BB_calls']
            require(w['conditional_BB_calls']>0,'Gamma arithmetic did not execute')
            d['reference32_arrays_exactly_equal']=bool(w['logmass']==old['logmass']
                and np.array_equal(w['gradient'],old['gradient']) and np.array_equal(w['hessian'],old['hessian']))
            indices=w['report']['augmented_global_indices'];inactive=[i for i in range(29) if i not in indices]
            d['inactive_exact_zero']=bool(np.all(np.asarray(w['gradient'])[inactive]==0)
                and np.all(np.asarray(w['hessian'])[inactive,:]==0) and np.all(np.asarray(w['hessian'])[:,inactive]==0))
            mapped=wrapper.embed_report_errors(w['report'],model,gene,np,True)
            d['errors_rederived']=mapped['logmass']==w['error_allowances']['logmass'] and all(
                compare(mapped[k],w['error_allowances'][k])==0 for k in ('gradient','hessian'))
            d['passed']=d['LL']<=1e-9 and d['score']<=1e-8 and d['Hessian']<=1e-5 and all(
                d[k] for k in ('same_geometry','inactive_exact_zero','errors_rederived'))
            records.append(d)
        full=dict(LL=abs(ll-ref['logmass']),score=compare(g,ref['gradient']),Hessian=compare(h,ref['hessian']))
        full['reference32_arrays_exactly_equal']=bool(ll==ref['logmass'] and np.array_equal(g,ref['gradient'])
            and np.array_equal(h,ref['hessian']))
        full['passed']=full['LL']<=1e-9 and full['score']<=1e-8 and full['Hessian']<=1e-5
        loop_ll=0.
        for w in workers:loop_ll+=w['logmass']
        full['ascending_gene_sum_exact']=ll==loop_ll and all(
            compare(a,sum((np.asarray(w[k]) for w in workers),np.zeros_like(a)))==0
            for a,k in ((g,'gradient'),(h,'hessian')))
        summed=model.last_total_error_allowance
        full['error_sum_exact']=summed['logmass']==sum(w['error_allowances']['logmass'] for w in workers) and all(
            compare(summed[k],sum((np.asarray(w['error_allowances'][k]) for w in workers),np.zeros_like(summed[k])))==0
            for k in ('gradient','hessian'))
        unchanged=counts()==COUNT_SHA and wrapper.row_signatures(base.rows,np)==signatures and private.bb_terms is original_bb
        save('comparison',dict(genes=records,full=full,all_fixed_inputs_unchanged=unchanged))
        require(all(x['passed'] for x in records) and all(full[k] for k in ('passed','ascending_gene_sum_exact','error_sum_exact'))
                and unchanged,'Full batched comparison failed')
        save('summary',dict(completed=True,returned=True,all_eight_gene_comparisons_passed=True,
            full_comparison_passed=True,all_fixed_inputs_unchanged=True,elapsed_seconds=time.perf_counter()-began,
            positive_evaluation_seconds=elapsed,measured_process_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            no_fit=True,no_calibration=True,no_independent_full_Hessian=True))
    except Exception as exc:
        save('failure',dict(completed=False,error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()))
        raise


if __name__=='__main__':main()
