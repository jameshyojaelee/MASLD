#!/usr/bin/env python3
"""Allowed 64-node GL batches at either frozen P3 gene2 point.

Reuse only the exact reconstructed synthetic count vector. No native mode
search, optimizer, order selection or independent full-Hessian claim.
"""
import argparse
import ast
from dataclasses import asdict,replace
import hashlib
import importlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import struct
import sys
import time
import traceback

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE=ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
DRAW=REC/'repair_p3_gamma_native_P3_check_22079693'
ADMISSION=REC/'repair_p3_gamma_batched_full_check_22078827'
CASES={'failed0gene2':22039200,'positive18':22039201}
REFERENCE32={'failed0gene2':REC/'repair_p3_gamma_batched_P3_failed0gene2_22079995',
    'positive18':REC/'repair_p3_gamma_batched_P3_positive18_22079996'}
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
 'repair_p3_integral.py','repair_p3_local_gl.py','repair_p3_local_moments.py',
 'repair_p3_gamma_parallel_model.py','repair_p3_bb_gamma_difference_v3.py',
 'repair_p3_batched_conditional.py','repair_p3_prepared_conditional_v2.py',
 'repair_p3_bb_gamma_batched.py','repair_p3_gamma_batched_conditional.py','repair_p3_gamma_batched_model.py')


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case',choices=tuple(CASES),required=True)
    ap.add_argument('--pins',type=Path,required=True);ap.add_argument('--expected-pins-sha',required=True)
    ap.add_argument('--expected-checker-sha',required=True);ap.add_argument('--inspect',action='store_true');args=ap.parse_args()
    require(sha(Path(__file__))==args.expected_checker_sha and sha(args.pins)==args.expected_pins_sha,'Checker/manifest differs')
    pins={Path(p):v for p,v in json.loads(args.pins.read_text()).items()}
    for p,v in pins.items():
        require(p.is_absolute() and p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==v,'Pinned input differs')
        if p.suffix=='.py':ast.parse(p.read_text())
    parentdir=REC/('repair_p3_fixed_point_'+args.case+'_'+str(CASES[args.case]))
    required=[*(REFERENCE32[args.case]/n for n in ('protocol.json','summary.json','answer.json','settings.json')),
        *(HERE/n for n in SOURCES),*(DRAW/n for n in ('protocol.json','draw_reconstruction.json')),
        *(ADMISSION/n for n in ('protocol.json','summary.json','comparison.json')),
        *(parentdir/n for n in ('protocol.json','candidate_baseline.json')),
        *(FIXTURE/n for n in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv'))]
    require(all(p in pins for p in required),'Required input absent')
    read=lambda p:json.loads(p.read_text())
    admitted=read(ADMISSION/'summary.json');admitp=read(ADMISSION/'protocol.json')
    require(admitted['completed'] and admitted['all_eight_gene_comparisons_passed']
        and admitted['full_comparison_passed'] and admitted['all_fixed_inputs_unchanged']
        and read(ADMISSION/'comparison.json')['full']['passed'],'Batched full G8 admission incomplete')
    for n in SOURCES:require(admitp['source_and_input_sha256'][str(HERE/n)]==pins[HERE/n],'Admitted source differs: '+n)
    priorp=read(REFERENCE32[args.case]/'protocol.json');priors=read(REFERENCE32[args.case]/'summary.json')
    prior=read(REFERENCE32[args.case]/'answer.json')
    require(all(priors[k] is True for k in ('completed','returned','nominal_comparison_passed',
        'numerical_allowances_passed','all_fixed_inputs_unchanged')),'Checked batch32 reference incomplete')
    for n in SOURCES:require(priorp['source_and_input_sha256'][str(HERE/n)]==pins[HERE/n],'Batch32 source differs')
    require(prior['returned'] and len(prior['gradient'])==29 and len(prior['hessian'])==29
        and all(len(row)==29 for row in prior['hessian']) and all(math.isfinite(x) for x in
            [prior['logmass'],*prior['gradient'],*(x for row in prior['hessian'] for x in row)]),
        'Batch32 reference arrays invalid')
    drawp=read(DRAW/'protocol.json');receipt=read(DRAW/'draw_reconstruction.json')[args.case]
    parent=read(parentdir/'protocol.json');reference=read(parentdir/'candidate_baseline.json')
    for name in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv'):
        path=FIXTURE/name
        require(str(path) in drawp['source_and_input_sha256'] and path in pins
            and pins[path]==drawp['source_and_input_sha256'][str(path)],'Reconstructed fixture differs')
    require(len(receipt['full_counts'])==173 and len(receipt['full_row_keys'])==173 and receipt['target_rows']==24
        and receipt['original_BB_identity_preserved'] is True
        and hashlib.sha256(struct.pack('<173q',*receipt['full_counts'])).hexdigest()==parent['count_vector_sha256']
        ==receipt['count_vector_sha256']==drawp['cases'][args.case]['count_vector_sha256'],'Exact P3 count vector differs')
    require(parent['fixed_p']==drawp['cases'][args.case]['fixed_p'] and parent['fixed_v']==drawp['cases'][args.case]['fixed_v']
        and len(parent['fixed_p'])==28 and all(math.isfinite(v) for v in parent['fixed_p'])
        and math.isfinite(parent['fixed_v']) and parent['fixed_v']>0 and parent['fixed_gene']==2
        and parent['fixed_gene_rows']==24,'Frozen P3 point differs')
    require(len(reference['gradient'])==29 and len(reference['hessian'])==29
        and all(len(row)==29 for row in reference['hessian'])
        and all(math.isfinite(v) for v in [reference['logmass'],*reference['gradient'],
            *(v for row in reference['hessian'] for v in row)]),'Reference derivative arrays invalid')
    require(priorp['fixed_p']==parent['fixed_p'] and priorp['fixed_v']==parent['fixed_v']
        and priorp['count_vector_sha256']==parent['count_vector_sha256'] and priorp['node_batch_size']==32,
        'Batch32 point or batch scope differs')
    protocol=dict(case=args.case,source_and_input_sha256={str(p):v for p,v in pins.items()},
        fixed_p=parent['fixed_p'],fixed_v=parent['fixed_v'],count_vector_sha256=parent['count_vector_sha256'],
        fixture_seed=20260930,participants=60,full_rows=173,target_gene=2,target_rows=24,G=8,Z=7,P=28,
        with_hessian=True,node_batch_size=64,LL_tolerance=1e-9,score_tolerance=1e-8,Hessian_tolerance=1e-5/8,
        reference_LL_atol=1e-9,reference_score_atol=1e-8,reference_Hessian_atol=1e-5,
        unchanged_GL_rules_and_caps=True,no_native_mode_search=True,mixed_entries=[[10,24],[10,28]],
        no_fit=True,no_calibration=True,no_independent_full_Hessian=True,no_actual_or_sealed_data=True)
    if args.inspect:print(json.dumps(protocol));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
            'Threads not fixed before imports')
    out=REC/('repair_p3_gamma_batched_P3_batch64_'+args.case+'_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in SOURCES)):
        dest=archive/p.name;dest.write_bytes(p.read_bytes())
        require(sha(dest)==(args.expected_checker_sha if p==Path(__file__) else pins[p]),'Archive differs')
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
        if isinstance(x,(tuple,list)):return [failed_json(v) for v in x]
        return x
    def save(name,value):
        body=(json.dumps(value,default=serialize,allow_nan=False,indent=2)+'\n').encode()
        require(sum(p.stat().st_size for p in out.rglob('*') if p.is_file())+len(body)<32*1024**2-65536,'Output cap exhausted')
        with (out/(name+'.json')).open('xb') as f:f.write(body);f.flush();os.fsync(f.fileno())
        print(json.dumps(dict(saved=name)),flush=True)
    began=time.perf_counter()
    try:
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        require(environment==parent['environment']==priorp['environment'],'Environment differs')
        protocol.update(environment=environment,hostname=platform.node(),machine=platform.machine())
        with Path('/proc/cpuinfo').open() as f:
            protocol['cpu_model']=next((line.split(':',1)[1].strip() for line in f if line.startswith('model name')),None)
        save('protocol',protocol)
        sys.path.insert(0,str(archive));loader=importlib.import_module('bfix_v1');importlib.import_module('tests_v1')
        imported_original_bb=sys.modules['likelihood_v1'].bb_terms
        for n in SOURCES[:4]:
            m=sys.modules[n[:-3]];require(Path(m.__file__).resolve().parent==archive.resolve() and sha(Path(m.__file__))==pins[HERE/n],
                                        'Original import escaped archive')
        design,settings=loader.load();private=load('P3_batch_original',archive/'likelihood_v1.py');original_bb=private.bb_terms
        base=private.Model(design.model().rows,gene_count=8)
        require(len(base.rows)==173 and len(design.participants)==60 and settings['seed']==20260930
            and [[r.individual,r.gene] for r in base.rows]==receipt['full_row_keys'],'Original full row mapping differs')
        drawn=[replace(row,a=a) for row,a in zip(base.rows,receipt['full_counts'])]
        require(all(0<=r.a<=r.n for r in drawn),'Count outside native support')
        rows=tuple(r for r in drawn if r.gene==2);model=private.Model(rows,gene_count=8)
        require(len(rows)==24 and model.Z==7 and model.size==28,'Target dimensions differ')
        kernel=load('P3_batch_kernel',archive/'repair_p3_integral.py')
        gl=load('P3_batch_local_GL',archive/'repair_p3_local_gl.py')
        conditional=load('P3_batch_conditional',archive/'repair_p3_gamma_batched_conditional.py')
        mapping=load('P3_batch_mapping',archive/'repair_p3_gamma_batched_model.py')
        require(mapping.KERNEL_SHA==pins[HERE/'repair_p3_integral.py'],'Integration kernel pin differs')
        p=np.asarray(parent['fixed_p'],dtype=float);p.setflags(write=False);before=mapping.row_signatures(model.rows,np)
        Context=conditional.context_class(archive);context=Context(model,p,model.rows)
        settings=replace(kernel.Settings(),logmass_tolerance=1e-9,gradient_tolerance=1e-8,hessian_tolerance=1e-5/8)
        defaults=asdict(kernel.Settings());actual=asdict(settings)
        prior_settings=read(REFERENCE32[args.case]/'settings.json')
        require(defaults==prior_settings['defaults'] and actual==prior_settings['actual'],'Batch32 settings differ')
        require(all(actual[k]==v for k,v in defaults.items() if k!='hessian_tolerance')
            and actual['hessian_tolerance']==1e-5/8<=defaults['hessian_tolerance'],
            'GL rule/cap/LL/score settings differ from defaults')
        save('settings',dict(defaults=defaults,actual=actual,only_Hessian_budget_tightened=True))
        try:
            ll,g,h,report=gl.integrate_gene(model,p,parent['fixed_v'],model.rows,context.callback,with_hessian=True,
                expected_kernel_sha=mapping.KERNEL_SHA,kernel_path=archive/'repair_p3_integral.py',settings=settings,batch_size=64)
        except Exception as exc:
            save('refusal',failed_json(dict(error_type=type(exc).__name__,message=str(exc),
                integration_reports=getattr(exc,'integration_reports',None),likelihood_evaluations=getattr(exc,'likelihood_evaluations',None),
                gamma_BB_calls=context.gamma_BB_calls,traceback=traceback.format_exc())))
            raise
        require(g.shape==(29,) and h.shape==(29,29) and np.isfinite(ll) and np.isfinite(g).all() and np.isfinite(h).all(),
                'Returned arrays invalid')
        chunk_ledger=context.chunk_ledger()
        save('raw_integral_return',dict(logmass=ll,gradient=g,hessian=h,report=failed_json(report),
            conditional_chunk_ledger=chunk_ledger,gamma_BB_calls=context.gamma_BB_calls))
        def valid_chunk(x):
            requested=x['requested_nodes'];n=model.rows[x['row_ordinal']].n
            effective=min(requested,100000//(n+1))
            return (0<requested<=64 and x['effective_chunk_limit']==effective
                and x['chunk_lengths']==[min(effective,requested-start) for start in range(0,requested,effective)]
                and x['all_requested_nodes_returned'] and x['need_hessian'])
        require({x['row_ordinal'] for x in chunk_ledger}==set(range(24)) and all(
            any(x['row_ordinal']==i and x['requested_nodes']==64 for x in chunk_ledger) for i in range(24))
            and all(valid_chunk(x) for x in chunk_ledger),
            'Actual64-node requests violate the unchanged count-axis chunk rule')
        path_fields=('z_bound','grid_spacing','panels','order','node_evaluations','cumulative_evaluations')
        path=lambda r:[[x[k] for k in path_fields] for x in r['ladder']]
        require(path(report)==path(prior['report'])
            and report['likelihood_evaluations']==prior['report']['likelihood_evaluations'],
            'Batch32 and batch64 integration paths differ')
        errors=mapping.embed_report_errors(report,model,2,np,True)
        differences=dict(LL=abs(ll-reference['logmass']),score=float(np.max(np.abs(g-reference['gradient']))),
                         Hessian=float(np.max(np.abs(h-reference['hessian']))))
        batch32_differences=dict(LL=abs(ll-prior['logmass']),
            score=float(np.max(np.abs(g-prior['gradient']))),Hessian=float(np.max(np.abs(h-prior['hessian']))))
        batch32_passed=(batch32_differences['LL']<=1e-9 and batch32_differences['score']<=1e-8
            and batch32_differences['Hessian']<=1e-5)
        batch32_exact=bool(ll==prior['logmass'] and np.array_equal(g,prior['gradient']) and np.array_equal(h,prior['hessian']))
        require(context.gamma_BB_calls>0,'Gamma arithmetic did not execute')
        allowances=errors['logmass']<=1e-9 and np.max(errors['gradient'])<=1e-8 and np.max(errors['hessian'])<=1e-5/8
        passed=differences['LL']<=1e-9 and differences['score']<=1e-8 and differences['Hessian']<=1e-5
        unchanged=(private.bb_terms is original_bb and mapping.row_signatures(model.rows,np)==before
            and sys.modules['likelihood_v1'].bb_terms is imported_original_bb)
        save('answer',dict(returned=True,logmass=ll,gradient=g,hessian=h,report=failed_json(report),error_allowances=errors,
            differences=differences,batch32_differences=batch32_differences,
            batch32_arrays_exactly_equal=batch32_exact,identical_batch32_integration_path=True,
            all24_rows_received_requested64_nodes=True,
            rows_requiring_count_axis_splits=sum(r.n+1>100000//64 for r in model.rows),
            batch32_comparison_passed=bool(batch32_passed),nominal_comparison_passed=bool(passed),numerical_allowances_passed=bool(allowances),
            gamma_BB_calls=context.gamma_BB_calls,conditional_chunk_ledger=context.chunk_ledger(),
            mixed_entries=[dict(indices=[i,j],native=float(h[i,j]),reference=reference['hessian'][i][j],
                difference=float(h[i,j]-reference['hessian'][i][j])) for i,j in protocol['mixed_entries']]))
        require(unchanged,'Fixed inputs or original BB changed')
        for path,digest in pins.items():require(sha(path)==digest,'Pinned input changed during check')
        save('summary',dict(completed=True,returned=True,nominal_comparison_passed=bool(passed),
            batch32_comparison_passed=bool(batch32_passed),
            batch32_arrays_exactly_equal=batch32_exact,
            numerical_allowances_passed=bool(allowances),all_fixed_inputs_unchanged=True,
            elapsed_seconds=time.perf_counter()-began,no_fit=True,no_calibration=True,no_independent_full_Hessian=True))
    except Exception as exc:
        save('failure',dict(completed=False,error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()));raise


if __name__=='__main__':main()
