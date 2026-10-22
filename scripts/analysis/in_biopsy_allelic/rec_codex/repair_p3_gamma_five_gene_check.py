#!/usr/bin/env python3
"""Recheck native GH20/GH40 with Gamma arithmetic against saved scalar masses.

Fifteen fixed synthetic points; no fits, resampling or human outcomes.
The independent scalar results and their uncertainty are reused unchanged.
This tests native numerical methods directly, not the G8-only model wrapper.
"""
import argparse
import ast
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
    args=ap.parse_args();require(sha(Path(__file__))==args.expected_runner_sha
        and sha(args.pins)==args.expected_pins_sha,'Runner/manifest differs')
    pins={Path(p):h for p,h in json.loads(args.pins.read_text()).items()}
    for p,h in pins.items():
        require(p.is_absolute() and p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==h,'Pinned input differs')
        if p.suffix=='.py':ast.parse(p.read_text())
    sources=('likelihood_v1.py','repair_p3_bb_gamma_difference_v3.py','check_quadrature_v1.py')
    require(all(HERE/n in pins for n in sources) and COLLECTION in pins,'Missing source/reference pin')
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
        require(r==entry['record'] and r['gene']==entry['gene'] and r['rows']==COUNTS[r['gene']]
            and r['tau']==.05*entry['tau_multiple'] and r['variance']==r['tau']**2,'Raw reference differs')
        original_sources={}
        for path,digest in protocol['source_sha256'].items():
            name=Path(path).name
            require(name not in original_sources,'Duplicate original source basename')
            original_sources[name]=digest
        require(protocol['seed']==20260923 and all(original_sources[n]==pins[HERE/n]
            for n in ('likelihood_v1.py','check_quadrature_v1.py')),'Original reference construction differs')
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
    if args.inspect:print(json.dumps(dict(points=15,orders=[20,40],references_checked=True)));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
            'Threads not fixed')
    out=REC/('repair_p3_gamma_five_gene_check_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in sources)):
        q=archive/p.name;q.write_bytes(p.read_bytes())
        require(sha(q)==(args.expected_runner_sha if p==Path(__file__) else pins[p]),'Archive differs')
    import numpy as np
    import scipy
    original=load('Gamma_five_original',archive/'likelihood_v1.py')
    gamma=load('Gamma_five_arithmetic',archive/'repair_p3_bb_gamma_difference_v3.py')
    original_bb=original.bb_terms
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
    def save(name,value):
        body=(json.dumps(value,allow_nan=False,indent=2)+'\n').encode()
        require(sum(f.stat().st_size for f in out.rglob('*') if f.is_file())+len(body)<CAP-65536,'Output cap exceeded')
        with (out/(name+'.json')).open('xb') as f:f.write(body);f.flush();os.fsync(f.fileno())
    save('protocol',dict(source_and_input_sha256={str(p):h for p,h in pins.items()},
        environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__),
        seed=20260923,rows_per_gene=COUNTS,scalar_reference_reused=True,
        scalar_error_estimates_not_certified=True,no_fits=True,no_calibration=True,no_human_data=True))
    records=[]
    for reference in references:
        gene=reference['gene'];rs=[r for r in drawn if r.gene==gene];model=original.Model(rs,gene_count=5)
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
        save('gene'+str(gene)+'_tau'+str(reference['tau']),result);records.append(result)
        print(json.dumps(dict(gene=gene,tau=reference['tau'],completed_points=len(records))),flush=True)
    for path,h in pins.items():require(sha(path)==h,'Input changed during evaluation')
    save('summary',dict(completed=True,points=len(records),returned_integrals=sum(r['returned'] for x in records for r in x['orders']),
        all15_likelihood_criteria_passed=all(x['GH20_GH40_absolute_pass'] is True
            and all(r.get('scalar_relative_LL_pass') is True for r in x['orders']) for x in records),
        no_independent_gradient_or_Hessian_acceptance=True,no_fits=True,no_calibration=True,no_human_data=True))


if __name__=='__main__':main()
