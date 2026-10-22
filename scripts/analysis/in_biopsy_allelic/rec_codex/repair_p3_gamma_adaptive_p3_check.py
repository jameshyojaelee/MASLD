#!/usr/bin/env python3
"""Unchanged adaptive GK with checked scalar Gamma at two frozen P3 points.

Reuse exact reconstructed synthetic counts and original adaptive settings.
No optimizer, fallback, calibration or independent full-Hessian claim.
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
GL_REFERENCES={'failed0gene2':22079995,'positive18':22079996}
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
    reference_dir=REC/('repair_p3_gamma_batched_P3_'+args.case+'_'+str(GL_REFERENCES[args.case]))
    required=[*(HERE/n for n in SOURCES),*(DRAW/n for n in ('protocol.json','draw_reconstruction.json')),
        *(ADMISSION/n for n in ('protocol.json','summary.json','comparison.json')),
        *(parentdir/n for n in ('protocol.json','candidate_baseline.json')),
        *(reference_dir/n for n in ('protocol.json','summary.json','answer.json')),
        *(FIXTURE/n for n in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv'))]
    require(all(p in pins for p in required),'Required input absent')
    read=lambda p:json.loads(p.read_text())
    admitted=read(ADMISSION/'summary.json');admitp=read(ADMISSION/'protocol.json')
    require(admitted['completed'] and admitted['all_eight_gene_comparisons_passed']
        and admitted['full_comparison_passed'] and admitted['all_fixed_inputs_unchanged']
        and read(ADMISSION/'comparison.json')['full']['passed'],'Batched full G8 admission incomplete')
    for n in SOURCES:require(admitp['source_and_input_sha256'][str(HERE/n)]==pins[HERE/n],'Admitted source differs: '+n)
    drawp=read(DRAW/'protocol.json');receipt=read(DRAW/'draw_reconstruction.json')[args.case]
    parent=read(parentdir/'protocol.json');reference=read(reference_dir/'answer.json')
    refp=read(reference_dir/'protocol.json');refs=read(reference_dir/'summary.json')
    require(refs['completed'] is True and refs['returned'] is True
        and refs['nominal_comparison_passed'] is True and refs['numerical_allowances_passed'] is True
        and refp['fixed_p']==parent['fixed_p'] and refp['fixed_v']==parent['fixed_v']
        and refp['count_vector_sha256']==parent['count_vector_sha256'],
        'Checked Gamma GL reference differs')
    for n in SOURCES:
        require(refp['source_and_input_sha256'][str(HERE/n)]==pins[HERE/n],
                'Checked Gamma GL source differs: '+n)
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
    protocol=dict(case=args.case,source_and_input_sha256={str(p):v for p,v in pins.items()},
        fixed_p=parent['fixed_p'],fixed_v=parent['fixed_v'],count_vector_sha256=parent['count_vector_sha256'],
        fixture_seed=20260930,participants=60,full_rows=173,target_gene=2,target_rows=24,G=8,Z=7,P=28,
        with_hessian=True,method='unchanged adaptive Gauss-Kronrod with scalar Gamma3',
        LL_tolerance=1e-9,score_tolerance=1e-8,Hessian_tolerance=1e-5/8,
        reference_LL_atol=1e-9,reference_score_atol=1e-8,reference_Hessian_atol=1e-5,
        unchanged_adaptive_rules_and_caps=True,no_native_mode_search=True,mixed_entries=[[10,24],[10,28]],
        no_fit=True,no_calibration=True,no_independent_full_Hessian=True,no_actual_or_sealed_data=True)
    if args.inspect:print(json.dumps(protocol));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
            'Threads not fixed before imports')
    out=REC/('repair_p3_gamma_adaptive_P3_'+args.case+'_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
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
        require(environment==parent['environment'],'Environment differs');protocol['environment']=environment;save('protocol',protocol)
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
        mapping=load('P3_batch_mapping',archive/'repair_p3_gamma_batched_model.py')
        gamma=load('P3_adaptive_Gamma',archive/'repair_p3_bb_gamma_difference_v3.py')
        require(mapping.KERNEL_SHA==pins[HERE/'repair_p3_integral.py'],'Integration kernel pin differs')
        p=np.asarray(parent['fixed_p'],dtype=float);p.setflags(write=False);before=mapping.row_signatures(model.rows,np)
        require(all(fn.__globals__ is vars(private) for fn in
            (private.conditional_terms,private.het_inclusion,private.mixture_terms))
            and all('bb_terms' in fn.__code__.co_names for fn in (private.het_inclusion,private.mixture_terms)),
            'Native conditional BB lookup differs')
        gamma_calls=[0];conditional_calls=[0];first_row_calls=[0]
        def checked_bb(*a,**kw):
            gamma_calls[0]+=1
            return gamma.stable_bb_terms(*a,**kw)
        def conditional(row,*a,**kw):
            conditional_calls[0]+=1
            if row is model.rows[0]:first_row_calls[0]+=1
            return private.conditional_terms(row,*a,**kw)
        settings=replace(kernel.Settings(),logmass_tolerance=1e-9,gradient_tolerance=1e-8,hessian_tolerance=1e-5/8)
        defaults=asdict(kernel.Settings());actual=asdict(settings)
        require(all(actual[k]==v for k,v in defaults.items() if k!='hessian_tolerance')
            and actual['hessian_tolerance']==1e-5/8<=defaults['hessian_tolerance'],
            'Adaptive rule/cap/LL/score settings differ from defaults')
        save('settings',dict(defaults=defaults,actual=actual,only_Hessian_budget_tightened=True))
        private.bb_terms=checked_bb
        try:
            ll,g,h,report=kernel.integrate_gene(model,p,parent['fixed_v'],model.rows,conditional,
                with_hessian=True,settings=settings)
        except Exception as exc:
            save('refusal',failed_json(dict(error_type=type(exc).__name__,message=str(exc),
                integration_reports=getattr(exc,'integration_reports',None),likelihood_evaluations=getattr(exc,'likelihood_evaluations',None),
                gamma_BB_calls=gamma_calls[0],conditional_calls=conditional_calls[0],
                first_row_calls=first_row_calls[0],known_adaptive_refusal=isinstance(exc,kernel.IntegrationFailure),
                traceback=traceback.format_exc())))
            raise
        finally:
            private.bb_terms=original_bb
        save('raw_integral_return',dict(logmass=ll,gradient=g,hessian=h,report=failed_json(report),
            gamma_BB_calls=gamma_calls[0],conditional_calls=conditional_calls[0],
            first_row_calls=first_row_calls[0]))
        require(g.shape==(29,) and h.shape==(29,29) and np.isfinite(ll) and np.isfinite(g).all() and np.isfinite(h).all(),
                'Returned arrays invalid')
        errors=dict(logmass=report['refined']['logmass_error_allowance']+report['logmass_refinement_difference'],
            gradient=np.asarray(report['refined']['gradient_error_allowance'])+report['gradient_refinement_difference'],
            hessian=np.asarray(report['refined']['hessian_error_allowance'])+report['hessian_refinement_difference'],
            empirical_estimates_not_certified=True)
        require(errors['gradient'].shape==(29,) and errors['hessian'].shape==(29,29)
            and all(np.isfinite(errors[k]).all() and np.all(np.asarray(errors[k])>=0)
                for k in ('logmass','gradient','hessian')),'Invalid dense empirical error allowances')
        inactive=[i for i in range(16) if i not in (2,10)]
        require(np.all(g[inactive]==0) and np.all(h[inactive,:]==0) and np.all(h[:,inactive]==0),
                'Inactive gene derivatives differ from exact zero')
        require(gamma_calls[0]>0 and conditional_calls[0]==24*first_row_calls[0]
            and first_row_calls[0]==report['unique_likelihood_evaluations'],
            'Adaptive native Gamma conditional did not execute for every row/node')
        differences=dict(LL=abs(ll-reference['logmass']),score=float(np.max(np.abs(g-reference['gradient']))),
                         Hessian=float(np.max(np.abs(h-reference['hessian']))))
        allowances=errors['logmass']<=1e-9 and np.max(errors['gradient'])<=1e-8 and np.max(errors['hessian'])<=1e-5/8
        passed=differences['LL']<=1e-9 and differences['score']<=1e-8 and differences['Hessian']<=1e-5
        unchanged=(private.bb_terms is original_bb and mapping.row_signatures(model.rows,np)==before
            and sys.modules['likelihood_v1'].bb_terms is imported_original_bb)
        save('answer',dict(returned=True,logmass=ll,gradient=g,hessian=h,report=failed_json(report),error_allowances=errors,
            differences=differences,nominal_comparison_passed=bool(passed),numerical_allowances_passed=bool(allowances),
            gamma_BB_calls=gamma_calls[0],conditional_calls=conditional_calls[0],
            first_row_calls=first_row_calls[0],
            mixed_entries=[dict(indices=[i,j],native=float(h[i,j]),reference=reference['hessian'][i][j],
                difference=float(h[i,j]-reference['hessian'][i][j])) for i,j in protocol['mixed_entries']]))
        require(unchanged,'Fixed inputs or original BB changed')
        for path,digest in pins.items():require(sha(path)==digest,'Pinned input changed during check')
        require(passed and allowances,'Adaptive return disagrees with checked GL or exceeds unchanged budgets')
        save('summary',dict(completed=True,returned=True,nominal_comparison_passed=bool(passed),
            numerical_allowances_passed=bool(allowances),all_fixed_inputs_unchanged=True,
            elapsed_seconds=time.perf_counter()-began,no_fit=True,no_calibration=True,no_independent_full_Hessian=True))
    except Exception as exc:
        save('failure',dict(completed=False,error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()));raise


if __name__=='__main__':main()
