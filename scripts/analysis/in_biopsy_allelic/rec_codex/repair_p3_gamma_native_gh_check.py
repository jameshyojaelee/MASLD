#!/usr/bin/env python3
"""All four native adaptive GH orders at the frozen synthetic failed point.

Only BB arithmetic changes. No fitting, choice of a passing order or admission.
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
import resource
import sys
import time
import traceback

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
REF=REC/'repair_p3_gamma_full_parallel_22077927'
POINT=REC/'repair_p3_observed_failure_anchor_22042612'
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
         'repair_p3_bb_gamma_difference_v3.py')
COUNT_SHA='7009be4adeec6157fa4bbf3f46e48ec24aa398a0bee4d7906c78c231b7d3fbad'


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pins',type=Path,required=True)
    ap.add_argument('--expected-pins-sha',required=True)
    ap.add_argument('--expected-checker-sha',required=True)
    ap.add_argument('--inspect',action='store_true');args=ap.parse_args()
    require(sha(Path(__file__))==args.expected_checker_sha and sha(args.pins)==args.expected_pins_sha,
            'Checker/manifest differs')
    pins={Path(p):v for p,v in json.loads(args.pins.read_text()).items()}
    for p,v in pins.items():
        require(p.is_absolute() and p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==v,'Pinned input differs')
        if p.suffix=='.py':ast.parse(p.read_text())
    require(all(p in pins for p in [*(HERE/n for n in SOURCES),POINT/'protocol.json',
        *(REF/n for n in ('summary.json','protocol.json','answer.json','workers.json'))]),'Required pin absent')
    read=lambda p:json.loads(p.read_text())
    point=read(POINT/'protocol.json');refp=read(REF/'protocol.json');refs=read(REF/'summary.json')
    reference=read(REF/'answer.json');workers=read(REF/'workers.json')
    require(point['count_vector_sha256']==COUNT_SHA and refp['count_vector_sha256']==COUNT_SHA
        and point['fixed_p']==refp['parameters'] and point['fixed_v']==refp['variance']
        and refs['completed'] and refs['returned'] and refs['all_fixed_inputs_unchanged']
        and reference['returned'] and [w['gene'] for w in workers]==list(range(8))
        and all(w['success'] and w['present'] for w in workers),'Reference identity differs')
    for n in SOURCES:
        require(refp['source_and_input_sha256'][str(HERE/n)]==pins[HERE/n],'Reference used different original/Gamma source')
    protocol=dict(source_and_input_sha256={str(p):v for p,v in pins.items()},
        parameters=point['fixed_p'],variance=point['fixed_v'],count_vector_sha256=COUNT_SHA,
        native_adaptive_GH_orders=[20,40,80,160],all_orders_retained=True,
        no_posthoc_order_selection=True,original_mode_search_unchanged=True,
        strict_LL_atol=1e-9,strict_score_atol=1e-8,strict_Hessian_atol=1e-5,
        owner_score_atol=1e-8,owner_score_rtol=1e-6,owner_LL_atol=1e-6,
        no_fit=True,no_calibration=True,no_independent_full_Hessian=True,no_actual_or_sealed_data=True)
    if args.inspect:print(json.dumps(protocol));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
            'Threads not fixed before imports')
    out=REC/('repair_p3_gamma_native_GH_check_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in SOURCES)):
        q=archive/p.name;q.write_bytes(p.read_bytes())
        require(sha(q)==(args.expected_checker_sha if p==Path(__file__) else pins[p]),'Archive differs')
    import numpy as np
    import scipy
    def serialize(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    def save(name,x):
        body=(json.dumps(x,default=serialize,allow_nan=False,indent=2)+'\n').encode()
        require(sum(p.stat().st_size for p in out.rglob('*') if p.is_file())+len(body)<32*1024**2-65536,'Output cap exhausted')
        with (out/(name+'.json')).open('xb') as f:f.write(body);f.flush();os.fsync(f.fileno())
        print(json.dumps(dict(saved=name)),flush=True)
    began=time.perf_counter()
    try:
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        require(environment==point['environment']==refp['environment'],'Environment differs')
        protocol['environment']=environment;save('protocol',protocol)
        sys.path.insert(0,str(archive));loader=importlib.import_module('bfix_v1');importlib.import_module('tests_v1')
        for n in SOURCES[:4]:
            m=sys.modules[n[:-3]];require(Path(m.__file__).resolve().parent==archive.resolve() and sha(Path(m.__file__))==pins[HERE/n],
                                        'Original import escaped archive')
        design,settings=loader.load();private=load('native_GH_original',archive/'likelihood_v1.py')
        gamma=load('native_GH_gamma',archive/'repair_p3_bb_gamma_difference_v3.py')
        model=private.Model(design.model().rows,gene_count=8)
        require(len(model.rows)==173 and model.G==8 and model.Z==7 and model.size==28
                and settings['seed']==20260930 and len(design.participants)==60,'Synthetic dimensions differ')
        counts=lambda:hashlib.sha256(np.asarray([r.a for r in model.rows],dtype='<i8').tobytes()).hexdigest()
        require(counts()==COUNT_SHA,'Counts differ')
        p=np.asarray(point['fixed_p'],dtype=float);pbytes=p.tobytes();p.setflags(write=False)
        original=private.bb_terms;records=[]
        try:
            private.bb_terms=gamma.stable_bb_terms
            for order in protocol['native_adaptive_GH_orders']:
                start=time.perf_counter();terms={}
                try:
                    ll,g,h=model.marginal(p,point['fixed_v'],nodes=order,with_hessian=True,gene_terms=terms)
                    require(g.shape==(29,) and h.shape==(29,29),'Unexpected native GH output shape')
                    if not (np.isfinite(ll) and np.isfinite(g).all() and np.isfinite(h).all()):
                        raise ValueError('Nonfinite native GH answer')
                    dg=np.abs(g-np.asarray(reference['gradient']));dh=np.abs(h-np.asarray(reference['hessian']))
                    dll=abs(ll-reference['logmass'])
                    record=dict(order=order,returned=True,logmass=ll,gradient=g,hessian=h,gene_terms=terms,
                        LL_absolute_difference=dll,score_max_absolute_difference=float(dg.max()),
                        Hessian_max_absolute_difference=float(dh.max()),
                        all_eight_gene_LL_differences=[terms[w['gene']]-w['logmass'] for w in workers],
                        strict_criteria_passed=bool(dll<=1e-9 and np.all(dg<=1e-8) and np.all(dh<=1e-5)),
                        owner_LL_score_criteria_passed=bool(dll<=1e-6 and np.all(dg<=1e-8+1e-6*np.abs(reference['gradient']))))
                except (ValueError,np.linalg.LinAlgError) as exc:
                    record=dict(order=order,returned=False,error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
                record['seconds']=time.perf_counter()-start;records.append(record);save('order'+str(order),record)
        finally:private.bb_terms=original
        require(private.bb_terms is original and counts()==COUNT_SHA and p.tobytes()==pbytes,'Fixed state changed')
        save('summary',dict(completed=True,orders=[r['order'] for r in records],returned_orders=[r['order'] for r in records if r['returned']],
            all_orders_retained=True,elapsed_seconds=time.perf_counter()-began,
            measured_process_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            no_order_admitted=True,no_fit=True,no_calibration=True,no_independent_full_Hessian=True))
    except Exception as exc:
        save('failure',dict(completed=False,error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()));raise


if __name__=='__main__':main()
