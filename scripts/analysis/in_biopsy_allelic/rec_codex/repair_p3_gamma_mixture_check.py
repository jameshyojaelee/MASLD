#!/usr/bin/env python3
"""Recheck the unchanged six mixture-GH policies with checked Gamma arithmetic.

Both frozen synthetic P3 points and every policy remain required. No policy
selection, fit, replacement draw, calibration or human outcome access.
"""
import argparse
import ast
from dataclasses import replace,asdict
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
DECIMAL=REC/'repair_p3_decimal_GL_parallel_collection_20261007T033640Z/collection.json'
DECIMAL_LAUNCH=REC/'repair_p3_decimal_GL_parallel_launch_20261006T235003Z'
CASES={'failed0gene2':(22039200,22079995),'positive18':(22039201,22079996)}
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py','repair_p3_integral.py',
    'repair_p3_gamma_parallel_model.py','repair_p3_bb_gamma_difference_v3.py','repair_p3_mixture_gh.py')
MIXTURE_SHA='8a21be3a5d4f7c731b13cdfe5a01165f640be36556fc8227ec6a8c67441828a0'
GAMMA_SHA='546e7ea2f1db21e6fcdd6aed72204d3bbe48c2d818996ab21de5f5ab5215bcf3'


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--case',choices=tuple(CASES),required=True)
    ap.add_argument('--policy',choices=('equal','scale_sqrt','scale_precision'),required=True)
    ap.add_argument('--wide',type=int,choices=(0,1),required=True);ap.add_argument('--pins',type=Path,required=True)
    ap.add_argument('--expected-pins-sha',required=True);ap.add_argument('--expected-checker-sha',required=True)
    ap.add_argument('--inspect',action='store_true');args=ap.parse_args()
    require(sha(Path(__file__))==args.expected_checker_sha and sha(args.pins)==args.expected_pins_sha,'Checker/manifest differs')
    pins={Path(p):h for p,h in json.loads(args.pins.read_text()).items()}
    for p,h in pins.items():
        require(p.is_absolute() and p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==h,'Pinned input differs')
        if p.suffix=='.py':ast.parse(p.read_text())
    oldjob,refjob=CASES[args.case];point=REC/f'repair_p3_fixed_point_{args.case}_{oldjob}'
    ref=REC/f'repair_p3_gamma_batched_P3_{args.case}_{refjob}'
    required=[*(HERE/n for n in SOURCES),DECIMAL,DECIMAL_LAUNCH/'pins.json',DECIMAL_LAUNCH/'jobs.json',
        DRAW/'protocol.json',DRAW/'draw_reconstruction.json',
        point/'protocol.json',ref/'protocol.json',ref/'summary.json',ref/'answer.json',
        *(FIXTURE/n for n in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv'))]
    require(all(p in pins for p in required),'Required input missing')
    require(pins[HERE/'repair_p3_mixture_gh.py']==MIXTURE_SHA
        and pins[HERE/'repair_p3_bb_gamma_difference_v3.py']==GAMMA_SHA,'Frozen mixture/Gamma source differs')
    read=lambda p:json.loads(p.read_text())
    parent=read(point/'protocol.json');drawp=read(DRAW/'protocol.json');receipt=read(DRAW/'draw_reconstruction.json')[args.case]
    reference=read(ref/'answer.json');admission=read(ref/'summary.json');refp=read(ref/'protocol.json');decimal=read(DECIMAL)
    require(len(reference['gradient'])==29 and len(reference['hessian'])==29
        and all(len(r)==29 for r in reference['hessian'])
        and all(math.isfinite(x) for x in [reference['logmass'],*reference['gradient'],
            *(x for r in reference['hessian'] for x in r)]),'Reference arrays invalid')
    require(all(admission[k] is True for k in ('completed','returned','nominal_comparison_passed',
        'numerical_allowances_passed','all_fixed_inputs_unchanged')),'Gamma reference incomplete')
    for n in SOURCES:
        if n!='repair_p3_mixture_gh.py':require(refp['source_and_input_sha256'][str(HERE/n)]==pins[HERE/n],'Reference source differs')
    for n in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv'):
        require(drawp['source_and_input_sha256'][str(FIXTURE/n)]==pins[FIXTURE/n],'Reconstructed fixture differs')
    require(len(receipt['full_counts'])==173 and len(receipt['full_row_keys'])==173 and receipt['target_rows']==24
        and hashlib.sha256(struct.pack('<173q',*receipt['full_counts'])).hexdigest()==receipt['count_vector_sha256']
        ==parent['count_vector_sha256']==refp['count_vector_sha256'],'Count vector differs')
    require(parent['fixed_p']==refp['fixed_p']==drawp['cases'][args.case]['fixed_p']
        and parent['fixed_v']==refp['fixed_v']==drawp['cases'][args.case]['fixed_v']
        and parent['fixed_gene']==refp['target_gene']==2 and parent['fixed_gene_rows']==refp['target_rows']==24
        and len(parent['fixed_p'])==28 and all(math.isfinite(x) for x in parent['fixed_p'])
        and math.isfinite(parent['fixed_v']) and parent['fixed_v']>0,'Frozen point differs')
    require(decimal['complete_family'] is True and decimal['completed_chunks']==96
        and decimal['families'][args.case]['precision_within_1e_minus30'] is True
        and decimal['families'][args.case]['refinement_within_1e_minus8'] is True,'Independent mixed reference incomplete')
    decimal_pins=read(DECIMAL_LAUNCH/'pins.json')
    require(decimal_pins[str(point/'protocol.json')]==pins[point/'protocol.json'],'Independent point binding differs')
    decimal_jobs=[j for j in read(DECIMAL_LAUNCH/'jobs.json') if j['case']==args.case]
    require(len(decimal_jobs)==48 and {(j['rule'],j['chunk']) for j in decimal_jobs}
        =={(r,c) for r in range(3) for c in range(16)},'Independent chunk family differs')
    for j in decimal_jobs:
        d=Path(j['output_directory']);require(d.resolve().is_relative_to(REC.resolve()),'Independent chunk path escaped')
        cp=read(d/'protocol.json')
        require(cp['case']==args.case and cp['fixed_p']==parent['fixed_p'] and cp['fixed_v']==parent['fixed_v']
            and cp['count_vector_sha256']==parent['count_vector_sha256'] and cp['source_and_input_sha256']==decimal_pins
            and cp['target_mixed_entries']==[[10,24],[10,28]],'Independent raw point identity differs')
        for digits in (80,110):
            p=d/f'moments_{digits}.json'
            require(sha(p)==decimal['input_result_sha256'][str(p)],'Independent raw moments changed')
    finest=[r for r in decimal['families'][args.case]['all_rule_and_precision_results'] if r['rule']==2 and r['digits']==110]
    require(len(finest)==1,'Independent finest rule is not unique')
    protocol=dict(case=args.case,policy=args.policy,wide=bool(args.wide),source_and_input_sha256={str(p):h for p,h in pins.items()},
        fixed_p=parent['fixed_p'],fixed_v=parent['fixed_v'],count_vector_sha256=parent['count_vector_sha256'],
        policies=['equal','scale_sqrt','scale_precision'],wide_options=[False,True],required_case_policy_points=12,
        unchanged_mixture_rules=True,only_BB_arithmetic_changed=True,LL_atol=1e-9,score_atol=1e-8,Hessian_atol=1e-5,
        no_policy_selection=True,no_fit=True,no_calibration=True,no_independent_full_Hessian=True,no_human_data=True)
    if args.inspect:print(json.dumps(protocol));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
            'Threads not fixed')
    out=REC/f"repair_p3_gamma_mixture_{args.case}_{args.policy}_wide{args.wide}_{os.environ['SLURM_JOB_ID']}";out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in SOURCES)):
        q=archive/p.name;q.write_bytes(p.read_bytes())
        require(sha(q)==(args.expected_checker_sha if p==Path(__file__) else pins[p]),'Archive differs')
    import numpy as np
    import scipy
    def safe(x):
        if isinstance(x,np.ndarray):return safe(x.tolist())
        if isinstance(x,np.generic):return safe(x.item())
        if isinstance(x,float) and not math.isfinite(x):return {'nonfinite_float':x.hex()}
        if isinstance(x,dict):return {k:safe(v) for k,v in x.items()}
        if isinstance(x,(tuple,list)):return [safe(v) for v in x]
        return x
    def save(name,value):
        body=(json.dumps(safe(value),allow_nan=False,indent=2)+'\n').encode()
        require(sum(p.stat().st_size for p in out.rglob('*') if p.is_file())+len(body)<32*1024**2-65536,'Output cap exceeded')
        with (out/(name+'.json')).open('xb') as f:f.write(body);f.flush();os.fsync(f.fileno())
    began=time.perf_counter();private=None;original_bb=None
    try:
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        require(environment==parent['environment'],'Numerical environment differs');protocol['environment']=environment
        sys.path.insert(0,str(archive));loader=importlib.import_module('bfix_v1');importlib.import_module('tests_v1')
        require(Path(loader.FIXTURE).resolve()==FIXTURE.resolve(),'Loader fixture differs')
        for n in SOURCES[:4]:require(Path(sys.modules[n[:-3]].__file__).resolve()==(archive/n).resolve(),'Original import escaped archive')
        design,settings=loader.load();private=load('Gamma_mixture_original',archive/'likelihood_v1.py');original_bb=private.bb_terms
        base=private.Model(design.model().rows,gene_count=8)
        require(len(base.rows)==173 and len(design.participants)==60 and settings['seed']==20260930
            and [[r.individual,r.gene] for r in base.rows]==receipt['full_row_keys'],'Original row mapping differs')
        drawn=[replace(r,a=a) for r,a in zip(base.rows,receipt['full_counts'])]
        rows=tuple(r for r in drawn if r.gene==2);model=private.Model(rows,gene_count=8)
        require(len(rows)==24 and model.Z==7 and model.size==28 and all(r.a in r.allowed for r in rows),'Target support differs')
        require((model.G+2,2*model.G+1+model.Z,model.size)==(10,24,28),'Mixed coordinate layout differs')
        gamma=load('Gamma_mixture_arithmetic',archive/'repair_p3_bb_gamma_difference_v3.py')
        snapshots=load('Gamma_mixture_snapshots',archive/'repair_p3_gamma_parallel_model.py')
        kernel=load('Gamma_mixture_kernel',archive/'repair_p3_integral.py');mixture=load('Gamma_mixture_rule',archive/'repair_p3_mixture_gh.py')
        p=np.asarray(parent['fixed_p'],dtype=float);p.setflags(write=False);pbytes=p.tobytes();before=snapshots.row_signatures(rows,np)
        numerical_settings=kernel.Settings()
        protocol.update(settings=asdict(numerical_settings),mixture_metadata=mixture.metadata())
        require(private.conditional_terms.__globals__ is vars(private) and private.mixture_terms.__globals__ is vars(private)
            and 'mixture_terms' in private.conditional_terms.__code__.co_names
            and 'bb_terms' in private.mixture_terms.__code__.co_names,'Conditional BB lookup differs')
        gamma_calls=0;gamma_bb=gamma.stable_bb_terms
        def counted_gamma(*a,**kw):
            nonlocal gamma_calls
            gamma_calls+=1;return gamma_bb(*a,**kw)
        save('protocol',protocol);returned=False
        try:
            private.bb_terms=counted_gamma
            ll,g,h,report=mixture.integrate_gene(model,p,parent['fixed_v'],rows,private.conditional_terms,
                with_hessian=True,expected_kernel_sha=pins[HERE/'repair_p3_integral.py'],
                kernel_path=archive/'repair_p3_integral.py',settings=numerical_settings,
                prior_policy=args.policy,wide_components=bool(args.wide))
        except Exception as exc:
            live_kernel=sys.modules.get('mixture_gh_reviewed_kernel')
            known=(live_kernel is not None and Path(live_kernel.__file__).resolve()==(archive/'repair_p3_integral.py').resolve()
                and isinstance(exc,live_kernel.IntegrationFailure))
            tb=exc.__traceback__
            while tb is not None and tb.tb_next is not None:tb=tb.tb_next
            known=known or (type(exc) is ValueError and str(exc) in ('Nonfinite rising product','Nonfinite candidate BB result')
                and tb is not None and Path(tb.tb_frame.f_code.co_filename).resolve()==(archive/'repair_p3_bb_gamma_difference_v3.py').resolve())
            if not known:raise
            save('refusal',dict(error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc(),
                integration_reports=getattr(exc,'integration_reports',None),geometry_reports=getattr(exc,'geometry_reports',None),
                likelihood_evaluations=getattr(exc,'likelihood_evaluations',None),Gamma_BB_calls=gamma_calls))
        else:
            require(g.shape==(29,) and h.shape==(29,29) and np.isfinite(ll) and np.isfinite(g).all() and np.isfinite(h).all(),
                    'Invalid returned arrays')
            differences=dict(LL=abs(ll-reference['logmass']),score=float(np.max(np.abs(g-reference['gradient']))),
                Hessian=float(np.max(np.abs(h-reference['hessian']))))
            require(gamma_calls>0,'Gamma arithmetic was not exercised')
            dref=finest[0]['mixed_curvature_r_omega_and_r_variance']
            mixed=[float(h[10,j])-float(x) for j,x in zip((24,28),dref)]
            passed=(differences['LL']<=1e-9 and differences['score']<=1e-8 and differences['Hessian']<=1e-5
                and all(abs(x)<=1e-5 for x in mixed))
            save('answer',dict(logmass=ll,gradient=g,hessian=h,report=report,differences=differences,
                mixed_minus_independent_Decimal=mixed,nominal_comparison_passed=bool(passed),Gamma_BB_calls=gamma_calls,
                G8_Hessian_budget_postcheck_passed=bool(np.max(report['empirical_refinement_plus_tail_allowances'][2])<=1e-5/8),
                unchanged_default_Hessian_budget=True))
            returned=True
        finally:private.bb_terms=original_bb
        require(private.bb_terms is original_bb and p.tobytes()==pbytes and not p.flags.writeable
            and snapshots.row_signatures(rows,np)==before,'Fixed state changed')
        for path,digest in pins.items():require(sha(path)==digest,'Pinned input changed during check')
        save('summary',dict(completed=True,returned=returned,elapsed_seconds=time.perf_counter()-began,
            nominal_comparison_passed=bool(passed) if returned else None,all_fixed_inputs_unchanged=True,
            Gamma_BB_calls=gamma_calls,required_case_policy_points=12,no_policy_selection=True,
            no_fit=True,no_calibration=True,no_independent_full_Hessian=True))
    except Exception as exc:
        if private is not None and original_bb is not None:private.bb_terms=original_bb
        save('failure',dict(completed=False,error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()));raise


if __name__=='__main__':main()
