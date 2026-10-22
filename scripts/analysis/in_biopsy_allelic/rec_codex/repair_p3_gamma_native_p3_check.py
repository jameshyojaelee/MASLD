#!/usr/bin/env python3
"""All four native-Gamma GH orders at the two frozen synthetic P3 points.

Original BB generates the exact supplied count draws before any Gamma patch.
The existing whole GL references and mixed Hessian entries are retained.
No optimizer, calibration, order selection or independent-Hessian claim.
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
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE=ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
NULL=REC/'p3_parallel_20261003T130323Z/baseline/observed_null.json'
ADMISSION=REC/'repair_p3_gamma_native_model_check_22079139'
CASES={'failed0gene2':(0,22039200,'5f4a1aa455cd09b0d5f31889021ec88b145dcbb50c8e5b600343ea9a3a992691'),
       'positive18':(18,22039201,'6370f80717f29f477686c2df3eb0b7151ec907ebd8055d8baec3c314d79c47ca')}
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
 'repair_p3_bb_gamma_difference_v3.py','repair_p3_gamma_parallel_model.py','repair_p3_gamma_native_model.py')


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pins',type=Path,required=True);ap.add_argument('--expected-pins-sha',required=True)
    ap.add_argument('--expected-checker-sha',required=True);ap.add_argument('--inspect',action='store_true');args=ap.parse_args()
    require(sha(Path(__file__))==args.expected_checker_sha and sha(args.pins)==args.expected_pins_sha,'Checker/manifest differs')
    pins={Path(p):v for p,v in json.loads(args.pins.read_text()).items()}
    for p,v in pins.items():
        require(p.is_absolute() and p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==v,'Pinned input differs')
        if p.suffix=='.py':ast.parse(p.read_text())
    required=[*(HERE/n for n in SOURCES),NULL,*(ADMISSION/n for n in ('summary.json','protocol.json','state_probes.json')),
        *(FIXTURE/n for n in ('SHA256SUMS','globals.json','cohorts.tsv','individuals.tsv','rows.tsv','draws/p3_uniforms.npy'))]
    for case,(draw,job,digest) in CASES.items():
        required += [REC/('repair_p3_fixed_point_'+case+'_'+str(job))/n for n in ('protocol.json','candidate_baseline.json')]
    require(all(p in pins for p in required),'Required input absent')
    read=lambda p:json.loads(p.read_text())
    admission=read(ADMISSION/'summary.json');admitp=read(ADMISSION/'protocol.json')
    require(admission['completed'] and admission['byte_identical_integrals']==112
        and admission['original_boundary_full_byte_calls']==16 and admission['all_state_probes_passed']
        and all(v is True for v in read(ADMISSION/'state_probes.json').values()),'Scoped model checks incomplete')
    for n in SOURCES:require(admitp['source_and_input_sha256'][str(HERE/n)]==pins[HERE/n],'Admitted source differs')
    checksums={name.removeprefix('./'):h for h,name in
        (line.split(maxsplit=1) for line in (FIXTURE/'SHA256SUMS').read_text().splitlines())}
    for p in (FIXTURE/n for n in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv','draws/p3_uniforms.npy')):
        require(pins[p]==checksums[str(p.relative_to(FIXTURE))],'Original fixture differs')
    parents={};refs={}
    for case,(draw,job,digest) in CASES.items():
        d=REC/('repair_p3_fixed_point_'+case+'_'+str(job));parent=read(d/'protocol.json');ref=read(d/'candidate_baseline.json')
        require(parent['draw']==draw and parent['count_vector_sha256']==digest and parent['fixed_p'][10]==-15.
            and parent['fixed_gene']==2 and parent['fixed_gene_rows']==24,'Saved P3 identity differs')
        require(len(parent['fixed_p'])==28 and all(isinstance(x,(int,float)) and math.isfinite(x) for x in parent['fixed_p'])
            and isinstance(parent['fixed_v'],(int,float)) and math.isfinite(parent['fixed_v']) and parent['fixed_v']>0,
            'Finite original P3 parameters/variance required')
        require(len(ref['gradient'])==29 and len(ref['hessian'])==29
            and all(len(row)==29 for row in ref['hessian']) and all(math.isfinite(v) for v in
                [ref['logmass'],*ref['gradient'],*(v for row in ref['hessian'] for v in row)]),'P3 reference dimensions differ')
        parents[case]=parent;refs[case]=ref
    protocol=dict(source_and_input_sha256={str(p):v for p,v in pins.items()},
        cases={case:dict(draw=draw,count_vector_sha256=digest,fixed_p=parents[case]['fixed_p'],
             fixed_v=parents[case]['fixed_v']) for case,(draw,job,digest) in CASES.items()},
        orders=[20,40,80,160],gene=2,rows=24,full_draw_rows=173,participants=60,G=8,Z=7,P=28,fixture_seed=20260930,
        comparison_LL_atol=1e-9,comparison_score_atol=1e-8,comparison_Hessian_atol=1e-5,
        mixed_entries=[[10,24],[10,28]],all_orders_retained=True,no_order_admitted=True,
        original_BB_draw_generation=True,no_fit=True,no_calibration=True,no_independent_full_Hessian=True,
        no_actual_or_sealed_data=True)
    if args.inspect:print(json.dumps(protocol));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
            'Threads not fixed before imports')
    out=REC/('repair_p3_gamma_native_P3_check_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
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
    def nonfinite_json(x):
        if isinstance(x,np.ndarray):return nonfinite_json(x.tolist())
        if isinstance(x,np.generic):return nonfinite_json(x.item())
        if isinstance(x,float) and not math.isfinite(x):return {'nonfinite_float':x.hex()}
        if isinstance(x,(tuple,list)):return [nonfinite_json(v) for v in x]
        return x
    def save(name,value):
        body=(json.dumps(value,default=serialize,allow_nan=False,indent=2)+'\n').encode()
        require(sum(p.stat().st_size for p in out.rglob('*') if p.is_file())+len(body)<32*1024**2-65536,'Output cap exhausted')
        with (out/(name+'.json')).open('xb') as f:f.write(body);f.flush();os.fsync(f.fileno())
        print(json.dumps(dict(saved=name)),flush=True)
    began=time.perf_counter()
    try:
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        require(all(environment==p['environment'] for p in parents.values()),'Environment differs')
        protocol['environment']=environment;save('protocol',protocol)
        sys.path.insert(0,str(archive));loader=importlib.import_module('bfix_v1');importlib.import_module('tests_v1')
        for n in SOURCES[:4]:
            m=sys.modules[n[:-3]];require(Path(m.__file__).resolve().parent==archive.resolve() and sha(Path(m.__file__))==pins[HERE/n],
                                        'Original import escaped archive')
        design,settings=loader.load();private=load('P3_native_original',archive/'likelihood_v1.py')
        original_bb=private.bb_terms
        require(original_bb.__name__=='bb_terms' and Path(original_bb.__code__.co_filename).resolve()
            ==(archive/'likelihood_v1.py').resolve(),'Original BB source differs')
        base=private.Model(design.model().rows,gene_count=8)
        null=read(NULL);require(null['converged'] and settings['B']==50 and settings['seed']==20260930
            and len(base.rows)==173,'Generating null/design differs')
        nullp=np.asarray(null['parameters'],dtype=float)
        require(nullp.shape==(28,) and np.isfinite(nullp).all(),'Generating null parameters differ')
        uniforms=np.load(FIXTURE/'draws/p3_uniforms.npy',mmap_mode='r',allow_pickle=False)
        require(uniforms.shape==(50,173,2) and [(r['individual_id'],r['gene_id']) for r in loader.table('rows.tsv')]
            ==[(r.individual,settings['gene_order'][r.gene]) for r in base.rows],'Uniform row order differs')
        drawn_rows={};draw_receipts={}
        for case,(draw,job,digest) in CASES.items():
            require(private.bb_terms is original_bb,'Original BB required for exact draw generation')
            drawn=[]
            for row,(uc,uk) in zip(base.rows,uniforms[draw]):
                eta,r=base.predictors(row,nullp)[:2];drawn.append(private.conditional_draw(row,eta,r,uc,uk))
            actual_digest=hashlib.sha256(np.asarray([r.a for r in drawn],dtype='<i8').tobytes()).hexdigest()
            require(actual_digest==digest,'Exact count draw differs')
            drawn_rows[case]=tuple(r for r in private.Model(drawn,gene_count=8).rows if r.gene==2)
            draw_receipts[case]=dict(count_vector_sha256=actual_digest,full_counts=[r.a for r in drawn],
                full_row_keys=[[r.individual,r.gene] for r in drawn],target_rows=len(drawn_rows[case]),
                original_BB_identity_preserved=private.bb_terms is original_bb)
        save('draw_reconstruction',draw_receipts)
        # Both draws are complete before importing/building the Gamma evaluator.
        gamma=load('P3_native_gamma',archive/'repair_p3_bb_gamma_difference_v3.py')
        snapshots=load('P3_native_snapshots',archive/'repair_p3_gamma_parallel_model.py')
        wrapper=load('P3_native_wrapper',archive/'repair_p3_gamma_native_model.py')
        Model=wrapper.native_gamma_model_class(private,gamma,snapshots)
        require(private.bb_terms is original_bb,'Gamma setup changed original BB')
        allrecords=[];refinements=[]
        for case,(draw,job,digest) in CASES.items():
            rows=drawn_rows[case]
            model=Model(rows,gene_count=8);require(len(rows)==24 and model.Z==7 and model.size==28,'Target dimensions differ')
            before=snapshots.row_signatures(model.rows,np);p=np.asarray(parents[case]['fixed_p'],dtype=float);p.setflags(write=False)
            v=parents[case]['fixed_v'];reference=refs[case];records=[]
            for order in protocol['orders']:
                start=time.perf_counter()
                try:ll,g,h=model.marginal(p,v,nodes=order,with_hessian=True)
                except (ValueError,ArithmeticError,np.linalg.LinAlgError) as exc:
                    record=dict(case=case,order=order,returned=False,error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
                else:
                    require(g.shape==(29,) and h.shape==(29,29),'Unexpected native P3 output shape')
                    if np.isfinite(ll) and np.isfinite(g).all() and np.isfinite(h).all():
                        dl=abs(ll-reference['logmass']);dg=float(np.max(np.abs(g-reference['gradient'])));dh=float(np.max(np.abs(h-reference['hessian'])))
                        record=dict(case=case,order=order,returned=True,logmass=ll,gradient=g,hessian=h,
                            differences=dict(LL=dl,score=dg,Hessian=dh),nominal_comparison_passed=bool(dl<=1e-9 and dg<=1e-8 and dh<=1e-5),
                            mixed_entries=[dict(indices=[i,j],native=float(h[i,j]),reference=reference['hessian'][i][j],
                                difference=float(h[i,j]-reference['hessian'][i][j])) for i,j in protocol['mixed_entries']])
                    else:record=dict(case=case,order=order,returned=False,error_type='nonfinite_original_return',
                        logmass=nonfinite_json(ll),gradient=nonfinite_json(g),hessian=nonfinite_json(h),
                        nonfinite_entries=int(not np.isfinite(ll))+int((~np.isfinite(g)).sum())+int((~np.isfinite(h)).sum()))
                record['seconds']=time.perf_counter()-start;records.append(record);allrecords.append(record);save(case+'_order'+str(order),record)
                require(private.bb_terms is original_bb and snapshots.row_signatures(model.rows,np)==before,'Native evaluator state changed')
            for first,second in zip(records,records[1:]):
                q=dict(case=case,orders=[first['order'],second['order']],available=first['returned'] and second['returned'])
                if q['available']:
                    dl=abs(first['logmass']-second['logmass']);dg=np.abs(first['gradient']-second['gradient'])
                    q.update(LL_difference=dl,score_max_difference=float(dg.max()),
                        Hessian_max_difference=float(np.max(np.abs(first['hessian']-second['hessian']))),
                        owner_LL_score_passed=bool(dl<=1e-8+1e-8*abs(second['logmass']) and np.all(dg<=1e-8+1e-6*np.abs(second['gradient']))))
                else:q['owner_LL_score_passed']=False
                refinements.append(q)
        save('refinements',refinements)
        for path,digest in pins.items():require(sha(path)==digest,'Pinned input changed during check')
        save('summary',dict(completed=True,records=len(allrecords),returned_records=sum(r['returned'] for r in allrecords),
            all_nominal_comparisons_passed=all(r.get('nominal_comparison_passed',False) for r in allrecords),
            all_owner_refinements_passed=all(r['owner_LL_score_passed'] for r in refinements),
            elapsed_seconds=time.perf_counter()-began,no_order_admitted=True,no_fit=True,no_calibration=True,
            no_independent_full_Hessian=True,no_actual_or_sealed_data=True))
    except Exception as exc:
        save('failure',dict(completed=False,error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()));raise


if __name__=='__main__':main()
