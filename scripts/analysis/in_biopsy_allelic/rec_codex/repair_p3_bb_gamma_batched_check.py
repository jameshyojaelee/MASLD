#!/usr/bin/env python3
"""Independent rising products and selected BB checks of node batching only."""
import argparse
import ast
from decimal import Decimal,localcontext
import hashlib
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
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
CASES=REC/'repair_p3_bb_mode_anchor_batch_fast_launch_20261005T172617Z/cases.json'
ADMISSION=REC/'repair_p3_bb_gamma_difference_v3_check_22076101'
SOURCES=('repair_p3_bb_gamma_difference_v3.py','repair_p3_bb_gamma_batched.py')
CAP=32*1024**2


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load(name,p):
    spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pins',type=Path,required=True)
    parser.add_argument('--expected-pins-sha',required=True)
    parser.add_argument('--expected-checker-sha',required=True)
    parser.add_argument('--inspect',action='store_true')
    args=parser.parse_args()
    require(sha(Path(__file__))==args.expected_checker_sha and sha(args.pins)==args.expected_pins_sha,
            'Executable/manifest changed')
    pins={Path(p):h for p,h in json.loads(args.pins.read_text()).items()}
    for p,h in pins.items():
        require(p.is_absolute() and p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==h,
                'Source/reference changed: '+str(p))
        if p.suffix=='.py':ast.parse(p.read_text())
    require(all(p in pins for p in (CASES,ADMISSION/'protocol.json',ADMISSION/'summary.json',
                                  *(HERE/n for n in SOURCES))),'Required reference absent')
    admitted=json.loads((ADMISSION/'summary.json').read_text())
    require(all(admitted[k] is True for k in ('completed','selected21_criterion_met',
        'independent_rising_criterion_met','conditional_nominal_LP_score_H_criteria_met')),
        'Scalar Gamma prerequisite incomplete')
    cases=json.loads(CASES.read_text())
    require(len(cases['cases'])==21,'Independent selected BB family differs')
    for name,h in cases['input_sha256'].items():
        require(REC/name in pins and pins[REC/name]==h,'Independent reference receipt absent')
    protocol=dict(source_and_input_sha256={str(p):h for p,h in pins.items()},
        independent_BB_points=21,independent_rising_digits=80,
        logshapes=['-1200','log(1e-8)','0','log(31.999)','log(32)','log(32.001)','log(1000)','log(1e6)'],
        rising_counts=[0,1,2,999,1000,1001,2500,5000],
        rising_LP_atol=1e-9,rising_derivative_atol=1e-10,rising_derivative_rtol=1e-10,
        BB_LP_atol=1e-9,BB_score_atol=1e-8,BB_Hessian_atol=1e-5,
        scalar_family_n=[0,5,100,3000,5000],scalar_family_r=[-15.,-3.,5.],
        scalar_family_eta_nodes=64,scalar_family_eta_range=[-30.,30.],
        timing_n=3000,timing_nodes=32,timing_tail_counts=300,timing_repeats=3,
        no_conditional_mixture=True,no_integral=True,no_fit=True,no_calibration=True,
        maximum_entries=100000,maximum_nodes=64,output_cap_bytes=CAP)
    if args.inspect:print(json.dumps(protocol));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    out=REC/('repair_p3_bb_gamma_batched_check_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in SOURCES)):
        dest=archive/p.name;dest.write_bytes(p.read_bytes())
        require(sha(dest)==(args.expected_checker_sha if p==Path(__file__) else pins[p]),'Archive differs')
    began=time.perf_counter()
    import numpy as np
    import scipy
    def serial(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    def save(name,value):
        body=(json.dumps(value,default=serial,allow_nan=False,indent=2)+'\n').encode()
        require(sum(p.stat().st_size for p in out.rglob('*') if p.is_file())+len(body)<CAP-65536,
                'Output capacity exceeded')
        with (out/(name+'.json')).open('xb') as f:f.write(body);f.flush();os.fsync(f.fileno())
        print(json.dumps(dict(saved=name,seconds=time.perf_counter()-began)),flush=True)
    try:
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        require(environment==json.loads((ADMISSION/'protocol.json').read_text())['environment'],
                'Numerical environment differs')
        protocol['environment']=environment;save('protocol',protocol)
        scalar=load('batched_gamma_scalar',archive/SOURCES[0]);batch=load('batched_gamma_candidate',archive/SOURCES[1])
        logshapes=np.asarray([-1200.,math.log(1e-8),0.,math.log(31.999),math.log(32.),
                              math.log(32.001),math.log(1000.),math.log(1e6)])
        counts=np.asarray(protocol['rising_counts']);got=batch.rising(counts,logshapes,True)
        rising=[]
        with localcontext() as ctx:
            ctx.prec=80
            for i,ls in enumerate(logshapes):
                logshape=Decimal.from_float(float(ls));a=logshape.exp();lp,d,q=Decimal(0),Decimal(0),Decimal(0)
                expected={0:(lp,d,q)}
                for j in range(5000):
                    if j==0:lp,d,q=logshape,Decimal(1),Decimal(0)
                    else:
                        denominator=a+Decimal(j)
                        lp+=denominator.ln();d+=a/denominator;q+=a*Decimal(j)/(denominator*denominator)
                    if j+1 in counts:expected[j+1]=(lp,d,q)
                differences=[];passed=True
                for k_index,k in enumerate(counts):
                    truth=expected[int(k)]
                    delta=[Decimal.from_float(float(got[h][i,k_index]))-truth[h] for h in range(3)]
                    valid=abs(delta[0])<=Decimal('1e-9') and all(abs(x)<=Decimal('1e-10') for x in delta[1:])
                    relative=[None if abs(truth[h])<=Decimal('1e-250') else str(abs(delta[h]/truth[h])) for h in (1,2)]
                    valid=valid and all(x is None or Decimal(x)<=Decimal('1e-10') for x in relative)
                    differences.append(dict(count=int(k),expected=[str(x) for x in truth],
                        batch_minus_Decimal=[str(x) for x in delta],relative_derivative_errors=relative,passed=bool(valid)))
                    passed=passed and valid
                rising.append(dict(logshape=float(ls),records=differences,passed=bool(passed)))
        save('independent_rising',rising)
        require(all(x['passed'] for x in rising),'Independent rising-product criterion failed')
        no_second=batch.rising(counts,logshapes,False)
        require(no_second[2] is None and np.array_equal(no_second[0],got[0])
                and np.array_equal(no_second[1],got[1]),'Score-only rising branch changed values')
        expanded=batch.rising(counts,np.tile(logshapes,4)[::-1],True)
        permutation_max=[float(np.max(np.abs(expanded[h]-np.tile(got[h],(4,1))[::-1]))) for h in range(3)]
        require(permutation_max[0]<=1e-9 and max(permutation_max[1:])<=1e-10,'Mixed-shape batch mapping changed')
        selected=[]
        for case in cases['cases']:
            actual=batch.bb_terms([case['k']],case['n'],[case['eta']],case['r'],scalar,True)
            raw=[actual[0][0,0],*actual[1][0,0],*actual[2][0,0].ravel()]
            ref=next(x for x in case['reference_ladder'] if x['digits']==110 and x['step']=='1e-20')
            expected=[ref['logpmf'],*ref['score'],*(v for row in ref['hessian'] for v in row)]
            with localcontext() as ctx:
                ctx.prec=110
                delta=[Decimal.from_float(float(x))-Decimal(y) for x,y in zip(raw,expected)]
                passed=abs(delta[0])<=Decimal('1e-9') and all(abs(x)<=Decimal('1e-8') for x in delta[1:3])
                passed=passed and all(abs(x)<=Decimal('1e-5') for x in delta[3:])
            selected.append(dict(origin=case['origin'],n=case['n'],k=case['k'],eta=case['eta'],r=case['r'],
                                 batch_minus_Decimal=[str(x) for x in delta],passed=bool(passed)))
        save('independent_selected',selected)
        require(all(x['passed'] for x in selected),'Independent selected BB criterion failed')
        family=[]
        for n in protocol['scalar_family_n']:
            k=np.unique(np.clip([0,1,2,n//4,n//2,3*n//4,n-1,n],0,n))
            eta=np.linspace(-30.,30.,64)
            for r in protocol['scalar_family_r']:
                for flag in (False,True):
                    actual=batch.bb_terms(k,n,eta,r,scalar,flag)
                    references=[scalar.stable_bb_terms(k,n,float(x),r,flag) for x in eta]
                    maximum=[float(np.max(np.abs(actual[h]-np.stack([v[h] for v in references]))))
                             for h in range(3 if flag else 2)]
                    passed=maximum[0]<=1e-9 and maximum[1]<=1e-8 and (not flag or maximum[2]<=1e-5)
                    family.append(dict(n=n,r=r,Hessian=flag,nodes=64,counts=k,max_difference=maximum,passed=bool(passed)))
        save('scalar_family',family);require(all(x['passed'] for x in family),'Scalar BB batch family failed')
        refused=False
        try:batch.bb_terms(np.arange(5001),5000,np.zeros(64),-3.,scalar,True)
        except RuntimeError:refused=True
        require(refused,'Oversized count-axis allocation not refused')
        timings=[];k=np.r_[np.arange(150),np.arange(2851,3001)];eta=np.linspace(-2.,2.,32)
        for flag in (False,True):
            for repeat in range(3):
                times={};answers={}
                for mode in (('scalar','batch') if repeat%2==0 else ('batch','scalar')):
                    start=time.perf_counter()
                    if mode=='batch':answer=batch.bb_terms(k,3000,eta,-15.,scalar,flag)
                    else:answer=[scalar.stable_bb_terms(k,3000,float(x),-15.,flag) for x in eta]
                    times[mode]=time.perf_counter()-start
                    answers[mode]=answer
                maximum=[float(np.max(np.abs(answers['batch'][h]-np.stack([v[h] for v in answers['scalar']]))))
                         for h in range(3 if flag else 2)]
                passed=maximum[0]<=1e-9 and maximum[1]<=1e-8 and (not flag or maximum[2]<=1e-5)
                timings.append(dict(Hessian=flag,repeat=repeat,seconds=times,max_difference=maximum,
                                    passed=bool(passed),scalar_to_batch=times['scalar']/times['batch']))
                require(passed,'Timed workload changed numerical result')
        save('timing',timings)
        for p,h in pins.items():require(sha(p)==h,'Source/input changed while running')
        save('summary',dict(completed=True,independent_rising_passed=True,independent_selected21_passed=True,
            all_scalar_family_passed=True,scalar_family_records=len(family),permutation_max=permutation_max,
            oversized_batch_refused=True,measured_process_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            elapsed_seconds=time.perf_counter()-began,no_conditional_mixture=True,no_integral=True,no_fit=True,no_calibration=True))
    except Exception as exc:
        save('failure',dict(error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()));raise


if __name__=='__main__':main()
