#!/usr/bin/env python3
"""One guarded G8 synthetic point, sequential versus eight forked genes.

No fitting, simulation, numerical work on login, calibration or source adoption.
"""
import argparse
import ast
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys
import time
import traceback

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
CAMPAIGN=REC/'p3_parallel_20261003T130323Z'
FIXTURE=ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
HERE=Path(__file__).resolve().parent
PARALLEL=HERE/'repair_p3_parallel_gene_model.py'
SEQUENTIAL=HERE/'repair_p3_batched_model.py'
GL=HERE/'repair_p3_batched_gl.py'
CONDITIONAL=HERE/'repair_p3_batched_conditional.py'
V2=HERE/'repair_p3_prepared_conditional_v2.py'
KERNEL=REC/'repair_p3_launch_20261005T121740Z/repair_p3_integral.py'
NULL=REC/'repair_p3_batched_observed_fit_22040873/observed_null.json'
VARIANCE=REC/'repair_p3_batched_gl_positive18_22040259/protocol.json'
CAP=64*1024**2
GUARDS={
 CAMPAIGN/'campaign.json':'bb3e82833476f67b0b11d2388e2dc00358eb987367b7a6a666d11938c6baca94',
 SEQUENTIAL:'2f64fdaffb44e5d78fd327ca6caf8672242ff53e17e53ad89b7e249bd77c0c3c',
 GL:'8bffaf8057ace5a2593861f879068f7e4a3e3f1ef72d2fecc4199b94c8372c52',
 CONDITIONAL:'184b3ec43d205290e0c4e2d8df0237a9df770af38116d62762a8b1f9811f3d76',
 V2:'289e93b651411042093a2b987717d87cdb5ed156ab7c736e1ac5cf820670fd4b',
 KERNEL:'461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3',
 NULL:'edf9322d9bc5da1c843957496f7049824f5d89d665e2e15ca3d00c126b2b9a7b',
 VARIANCE:'4fff9a95d228ce2f12ab0ff08350df24d225cabda53044d00ab2a96582e43b56'}


def require(ok,message):
    if not ok:raise ValueError(message)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def guard(values):
    for p,d in values.items():require(sha(p)==d,'Reviewed source/synthetic input changed: '+str(p))
def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module);return module


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-checker-sha',required=True)
    parser.add_argument('--expected-parallel-sha',required=True)
    parser.add_argument('--inspect',action='store_true')
    args=parser.parse_args()
    guards={**GUARDS,Path(__file__).resolve():args.expected_checker_sha,PARALLEL:args.expected_parallel_sha}
    guard(guards)
    campaign=json.loads((CAMPAIGN/'campaign.json').read_text())
    for n in ('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py','check_bfix_tests_v13.py'):
        p=CAMPAIGN/'executed_sources'/n;guards[p]=campaign['files'][str(p)]
    guard(guards)
    for p in guards:
        if p.suffix=='.py':ast.parse(p.read_text())
    if args.inspect:
        print(json.dumps(dict(source_metadata_checks=True,numerical_execution=False,
            participants=60,rows=173,G=8,Z=7,P=28,point='saved converged observed null p; saved positive18 v',
            worker_count=8,start_method='Linux fork',node_batch=32,unchanged_batched_GL=True,
            fitting=False,calibration=False,biological_scale_validated=False,
            requested_allocation='8CPU/8GiB; generous runtime buffer needed; no submission here')))
        return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
        'Individual compute allocation required')
    require(int(os.environ.get('SLURM_CPUS_PER_TASK','0'))>=8 and len(os.sched_getaffinity(0))>=8,
        'Eight allocated/available CPUs required')
    require(all(os.environ.get(k)=='1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
        'Set numerical worker threads to1 before imports')
    for n in ('globals.json','individuals.tsv','rows.tsv','cohorts.tsv'):
        p=FIXTURE/n;guards[p]=campaign['files'][str(p)]
    guard(guards)
    output=REC/('repair_p3_parallel_gene_check_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(mode=0o750,exist_ok=False)
    archive=output/'executed_sources';archive.mkdir()
    archived={}
    for p,d in guards.items():
        if p.suffix=='.py':
            target=archive/p.name;require(not target.exists(),'Duplicate source basename')
            target.write_bytes(p.read_bytes());require(sha(target)==d,'Archived source differs');archived[target]=d
    np=None;phase='numerical_imports';started=time.perf_counter()
    def serial(value):
        if np is not None and isinstance(value,np.ndarray):return value.tolist()
        if np is not None and isinstance(value,np.generic):return value.item()
        raise TypeError(type(value).__name__)
    def save(n,value):
        body=json.dumps(value,indent=2,default=serial,allow_nan=False)+'\n'
        require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file())+len(body.encode())<=CAP-65536,
            'Comparison output cap exhausted')
        with (output/(n+'.json')).open('x') as h:h.write(body);h.flush();os.fsync(h.fileno())
    try:
        sys.path.insert(0,str(archive))
        require(not any(n in sys.modules for n in ('likelihood_v1','design_v1','bfix_v1','tests_v1')),
            'Original numerical modules already imported')
        import numpy as np
        import scipy
        loader=importlib.import_module('bfix_v1');original=importlib.import_module('likelihood_v1')
        importlib.import_module('tests_v1')
        for n in ('likelihood_v1','design_v1','bfix_v1','tests_v1'):
            require(Path(sys.modules[n].__file__).resolve().parent==archive.resolve(),'Original import escaped archive')
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        require(environment==dict(python='3.10.18',numpy='2.2.6',scipy='1.15.2'),'Frozen numerical environment differs')
        design,settings=loader.load();base=design.model()
        require(len(design.participants)==60 and len(base.rows)==173 and base.G==8 and base.Z==7 and base.size==28,
            'Exact synthetic fixture dimensions differ')
        null=json.loads(NULL.read_text())['fit'];require(null['converged'],'Saved observed null must converge')
        p=np.asarray(null['parameters'],float)
        variance=float(json.loads(VARIANCE.read_text())['fixed_v']);require(variance>0,'Saved positive variance required')
        counts=hashlib.sha256(np.asarray([r.a for r in base.rows],dtype='<i8').tobytes()).hexdigest()
        private=load_module('_repair_p3_parallel_private_original',archive/'likelihood_v1.py')
        kernel=load_module('_repair_p3_parallel_kernel',archive/KERNEL.name)
        gl=load_module('_repair_p3_parallel_GL',archive/GL.name)
        conditional=load_module('_repair_p3_parallel_conditional',archive/CONDITIONAL.name)
        sequential_source=load_module('_repair_p3_parallel_sequential_source',archive/SEQUENTIAL.name)
        parallel_source=load_module('_repair_p3_parallel_model_source',archive/PARALLEL.name)
        Sequential=sequential_source.batched_model_class(private,kernel,gl,conditional,kernel_path=archive/KERNEL.name,batch_size=32)
        Parallel=parallel_source.parallel_model_class(private,kernel,gl,conditional,sequential_source,
            kernel_path=archive/KERNEL.name,batch_size=32)
        seq=Sequential(base.rows,gene_count=8);par=Parallel(base.rows,gene_count=8)
        original_callables=(original.bb_terms,original.conditional_terms,original.Model.predictors,original.Model.evaluate)
        inherited={n:getattr(Parallel,n) is getattr(private.Model,n) for n in sequential_source.INHERITED}
        require(all(inherited.values()),'Original fitting/design methods changed')
        save('protocol',dict(source_guards={str(k):v for k,v in guards.items()},environment=environment,
            participants=60,rows=173,G=8,Z=7,P=28,fixed_p=p,fixed_v=variance,observed_count_vector_sha256=counts,
            null_source=str(NULL),variance_source=str(VARIANCE),point_predeclared=True,with_hessian=True,
            node_batch=32,workers=8,start_method='fork',inherited_methods=inherited,
            no_fits=True,no_new_draws=True,no_calibration=True,no_private_or_sealed_data=True,
            old20_node_adoption_not_asserted=True,local_GL_not_used=True,saved_output_cap_bytes=CAP))
        evaluations=[]
        for name,model in [('sequential',seq),('parallel',par)]:
            phase=name+'_full_point'
            began=time.perf_counter();terms={}
            try:
                result=model.marginal(p,variance,with_hessian=True,gene_terms=terms)
                record=dict(method=name,returned=True,seconds=time.perf_counter()-began,logmass=result[0],
                    gradient=result[1],hessian=result[2],gene_terms=terms,
                    reports=model.last_integration_reports,error_allowances=model.last_total_error_allowance)
            except Exception as exc:
                record=dict(method=name,returned=False,seconds=time.perf_counter()-began,
                    error_type=type(exc).__name__,message=str(exc),value_error=isinstance(exc,ValueError),
                    traceback=traceback.format_exc(),gene_terms=terms,
                    failure=getattr(model,'last_integration_failure',None),reports=model.last_integration_reports)
            save(name,record);evaluations.append(record)
        phase='fixed_point_comparison'
        a,b=evaluations;comparison=dict(both_returned=a['returned'] and b['returned'],
            original_failure_records_retained=True,no_forced_failure_injection=True)
        if comparison['both_returned']:
            difference=dict(logmass=abs(a['logmass']-b['logmass']),
                gradient_max=float(np.max(np.abs(a['gradient']-b['gradient']))),
                hessian_max=float(np.max(np.abs(a['hessian']-b['hessian']))))
            report_checks=[]
            def numerical_report(report):
                # Timing is observational, not part of method identity. All node,
                # arithmetic, allowance, callback memory and adaptive ledgers remain.
                result=dict(report);callback=dict(result['batched_callback']);callback.pop('setup_seconds')
                result['batched_callback']=callback
                return json.loads(json.dumps(result,default=serial,allow_nan=False))
            for x,y,w in zip(a['reports'],b['reports'],par.last_parallel_worker_records):
                gene=x['gene'];inactive=[j for j in range(28) if j not in (gene,8+gene) and j<16]
                g,h=w['gradient'],w['hessian']
                crosszeros=bool(np.all(g[inactive]==0) and np.all(h[inactive,:]==0) and np.all(h[:,inactive]==0))
                report_checks.append(dict(gene=gene,full_numerical_report_exact=numerical_report(x)==numerical_report(y),
                    crossgene_alpha_r_scores_and_hessians_exact_zero=crosszeros,
                    worker_pid=w['worker_pid'],worker_seconds=w['elapsed_seconds'],worker_fixed_inputs_unchanged=w['fixed_inputs_unchanged']))
            errors_exact=json.loads(json.dumps(a['error_allowances'],default=serial))==json.loads(json.dumps(b['error_allowances'],default=serial))
            passed=(difference['logmass']<=1e-9 and difference['gradient_max']<=1e-8 and difference['hessian_max']<=1e-5
                and len(report_checks)==8 and all(r['full_numerical_report_exact'] and r['crossgene_alpha_r_scores_and_hessians_exact_zero'] for r in report_checks)
                and errors_exact and a['gene_terms']==b['gene_terms'])
            comparison.update(differences=difference,pergene_checks=report_checks,error_allowances_exact=errors_exact,
                within_existing_fixed_point_budgets=bool(passed),gene_terms_exact=a['gene_terms']==b['gene_terms'],
                observed_sequential_to_parallel_time_ratio=a['seconds']/b['seconds'],
                parallel_capacity_estimate=par.last_parallel_capacity_estimate,
                runtime_single_point_not_fit_speedup=True,refusal_path_executed=False)
        else:
            af,bf=a.get('failure') or {},b.get('failure') or {}
            comparison.update(within_existing_fixed_point_budgets=False,
                same_return_disposition=a['returned']==b['returned'],same_first_failure_gene=af.get('gene')==bf.get('gene'),
                same_first_failure_message=a.get('message')==b.get('message'),
                no_success_claim_from_failure=True)
        guard(guards);guard(archived)
        require(original_callables==(original.bb_terms,original.conditional_terms,original.Model.predictors,original.Model.evaluate),
            'Unmodified original source reference changed')
        comparison['original_reference_callables_unchanged']=True
        comparison['source_input_postchecks_passed']=True
        save('comparison',comparison)
        save('summary',dict(comparison_records_completed=True,both_returned=comparison['both_returned'],
            within_existing_fixed_point_budgets=comparison['within_existing_fixed_point_budgets'],
            elapsed_seconds=time.perf_counter()-started,source_inputs_unchanged=True,
            no_fits=True,no_calibration=True,real558_person_scale_not_validated=True))
    except Exception as exc:
        save('failure',dict(phase=phase,error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()))
        raise


if __name__=='__main__':main()
