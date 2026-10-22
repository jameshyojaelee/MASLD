#!/usr/bin/env python3
"""Native adaptive GH at every frozen B-FIX vector/gene and failed call12.

All four orders are retained. Boundary calls use original BB arithmetic.
Independent scalar likelihood, quadrature refinement and per-gene derivative
comparisons are distinct; no optimizer, calibration or order admission.
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
import resource
import sys
import time
import traceback

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE=ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
REF=REC/'repair_p3_gamma_full_parallel_22077927'
ADMISSION=REC/'repair_p3_gamma_native_GH_check_22078807'
THETAS=('theta1_start_like','theta2_truth_v0','theta3_truth_vpos','theta4_perturbed')
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
         'repair_p3_bb_gamma_difference_v3.py','repair_p3_gamma_parallel_model.py')
COUNT_SHA='7009be4adeec6157fa4bbf3f46e48ec24aa398a0bee4d7906c78c231b7d3fbad'
ORDERS=(20,40,80,160)


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module;spec.loader.exec_module(module);return module


def reference_paths():
    return [(theta,gene,REC/('repair_p3_bfix_vectors_'+str(22040192+i*8+gene)))
            for i,theta in enumerate(THETAS) for gene in range(8)]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pins',type=Path,required=True)
    ap.add_argument('--expected-pins-sha',required=True)
    ap.add_argument('--expected-checker-sha',required=True)
    ap.add_argument('--inspect',action='store_true');args=ap.parse_args()
    require(sha(Path(__file__))==args.expected_checker_sha and sha(args.pins)==args.expected_pins_sha,'Checker/manifest differs')
    pins={Path(p):v for p,v in json.loads(args.pins.read_text()).items()}
    for p,v in pins.items():
        require(p.is_absolute() and p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==v,'Pinned input differs: '+str(p))
        if p.suffix=='.py':ast.parse(p.read_text())
    required=[*(HERE/n for n in SOURCES),*(REF/n for n in ('protocol.json','summary.json','answer.json','workers.json')),
        *(ADMISSION/n for n in ('summary.json','protocol.json')),
        *(FIXTURE/n for n in ('SHA256SUMS','globals.json','cohorts.tsv','individuals.tsv','rows.tsv')),
        *(FIXTURE/'params'/(t+'.json') for t in THETAS)]
    for theta,gene,directory in reference_paths():
        required += [directory/'protocol.json',directory/(theta+'_gene'+str(gene)+'.json')]
    require(all(p in pins for p in required),'Required input absent')
    read=lambda p:json.loads(p.read_text())
    refp=read(REF/'protocol.json');refs=read(REF/'summary.json');workers=read(REF/'workers.json');refanswer=read(REF/'answer.json')
    def finite_vector(value,width):
        return isinstance(value,list) and len(value)==width and all(isinstance(x,(int,float)) and math.isfinite(x) for x in value)
    def check_reference_arrays(value,width):
        require(finite_vector(value['gradient'],29) and isinstance(value['hessian'],list)
            and len(value['hessian'])==width and all(finite_vector(row,width) for row in value['hessian']),
            'Reference derivative dimensions/finiteness differ')
    require(refs['completed'] and refs['returned'] and refs['all_fixed_inputs_unchanged']
        and refp['count_vector_sha256']==COUNT_SHA and [x['gene'] for x in workers]==list(range(8))
        and all(x['success'] and x['present'] for x in workers),'Call12 reference incomplete')
    require(refanswer['returned'],'Call12 full reference did not return')
    check_reference_arrays(refanswer,29)
    for worker in workers:check_reference_arrays(worker,29)
    admit=read(ADMISSION/'summary.json');admitp=read(ADMISSION/'protocol.json')
    require(admit['completed'] and admit['orders']==list(ORDERS) and admit['returned_orders']==list(ORDERS)
        and admitp['parameters']==refp['parameters'] and admitp['variance']==refp['variance']
        and admitp['count_vector_sha256']==COUNT_SHA,'Native GH prerequisite differs')
    for name in SOURCES:
        require(refp['source_and_input_sha256'][str(HERE/name)]==pins[HERE/name],'Reference source differs: '+name)
    checksums={name.removeprefix('./'):h for h,name in
               (line.split(maxsplit=1) for line in (FIXTURE/'SHA256SUMS').read_text().splitlines())}
    for p in [FIXTURE/n for n in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv')]+[
            FIXTURE/'params'/(t+'.json') for t in THETAS]:
        require(pins[p]==checksums[str(p.relative_to(FIXTURE))],'Original fixture checksum differs')
    reference={}
    for theta,gene,directory in reference_paths():
        old=read(directory/(theta+'_gene'+str(gene)+'.json'));oldp=read(directory/'protocol.json')
        require(old['theta']==theta and old['gene']==gene and oldp['requested_points']==[[theta,gene]]
            and isinstance(old.get('independent_scalar'),dict)
            and math.isfinite(old['independent_scalar']['logmass']),'Original scalar reference incomplete')
        for p in [FIXTURE/n for n in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv')]+[
                FIXTURE/'params'/(t+'.json') for t in THETAS]:
            require(oldp['source_and_input_sha256'][str(p)]==pins[p],'Reference fixture differs')
        reference[theta,gene]=old
        if old.get('candidate') is not None:
            check_reference_arrays(old['candidate'],28 if old['v']==0 else 29)
    protocol=dict(source_and_input_sha256={str(p):v for p,v in pins.items()},
        count_vector_sha256=COUNT_SHA,participants=60,rows=173,G=8,Z=7,P=28,
        cases=[*THETAS,'saved_call12'],all_eight_genes=True,positive_orders=list(ORDERS),
        required_gene_points=40,required_integrals=112,boundary_points=16,
        independent_scalar_LL_atol=1e-8,refinement_LL_atol=1e-8,refinement_LL_rtol=1e-8,
        refinement_score_atol=1e-8,refinement_score_rtol=1e-6,
        call12_LL_atol=1e-9,call12_score_atol=1e-8,call12_Hessian_atol=1e-5,
        all_orders_retained=True,no_order_admitted=True,no_fit=True,no_calibration=True,
        no_independent_full_Hessian=True,no_actual_or_sealed_data=True)
    if args.inspect:print(json.dumps(protocol));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
            'Threads not fixed before imports')
    out=REC/('repair_p3_gamma_native_GH_family_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
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
    def save(name,x):
        body=(json.dumps(x,default=serialize,allow_nan=False,indent=2)+'\n').encode()
        require(sum(p.stat().st_size for p in out.rglob('*') if p.is_file())+len(body)<32*1024**2-65536,'Output cap exhausted')
        with (out/(name+'.json')).open('xb') as f:f.write(body);f.flush();os.fsync(f.fileno())
        print(json.dumps(dict(saved=name)),flush=True)
    began=time.perf_counter()
    try:
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        require(environment==refp['environment'],'Environment differs');protocol['environment']=environment;save('protocol',protocol)
        sys.path.insert(0,str(archive));loader=importlib.import_module('bfix_v1');importlib.import_module('tests_v1')
        for n in SOURCES[:4]:
            module=sys.modules[n[:-3]];require(Path(module.__file__).resolve().parent==archive.resolve()
                and sha(Path(module.__file__))==pins[HERE/n],'Original import escaped archive')
        design,settings=loader.load();private=load('GH_family_original',archive/'likelihood_v1.py')
        gamma=load('GH_family_gamma',archive/'repair_p3_bb_gamma_difference_v3.py')
        snapshots=load('GH_family_snapshots',archive/'repair_p3_gamma_parallel_model.py')
        base=private.Model(design.model().rows,gene_count=8)
        require(len(base.rows)==173 and base.Z==7 and base.size==28 and len(design.participants)==60
                and settings['seed']==20260930 and all(0<=r.n<=5000 for r in base.rows),'Synthetic scope differs')
        counts=lambda:hashlib.sha256(np.asarray([r.a for r in base.rows],dtype='<i8').tobytes()).hexdigest()
        require(counts()==COUNT_SHA,'Counts differ')
        before=snapshots.row_signatures(base.rows,np);original_bb=private.bb_terms
        cases=[(theta,*loader.parameters(FIXTURE/'params'/(theta+'.json'))) for theta in THETAS]
        cases.append(('saved_call12',np.asarray(refp['parameters'],dtype=float),refp['variance']))
        allrecords=[];allpoints=[]
        for label,p,v in cases:
            require(p.shape==(28,) and np.isfinite(p).all(),'Parameter dimensions differ');p.setflags(write=False);pbytes=p.tobytes()
            for gene in range(8):
                model=private.Model([r for r in base.rows if r.gene==gene],gene_count=8)
                old=reference.get((label,gene));orders=(20,) if v==0 else ORDERS;records=[]
                if old is not None:require(old['p']==p.tolist() and old['v']==v and old['rows']==len(model.rows),'Reference point differs')
                for order in orders:
                    start=time.perf_counter()
                    try:
                        if v>0:private.bb_terms=gamma.stable_bb_terms
                        ll,g,h=model.marginal(p,v,nodes=order,with_hessian=True)
                        relevant_h=h[:-1,:-1] if v==0 else h
                        require(g.shape==(29,) and relevant_h.shape==((28,28) if v==0 else (29,29)),'Unexpected output shape')
                        if not (np.isfinite(ll) and np.isfinite(g).all() and np.isfinite(relevant_h).all()):
                            raise ValueError('Nonfinite native GH result')
                        inactive=[i for i in range(16) if i not in (gene,8+gene)]
                        require(np.all(g[inactive]==0) and np.all(relevant_h[inactive,:]==0)
                            and np.all(relevant_h[:,inactive]==0),'Inactive-gene derivatives nonzero')
                        record=dict(case=label,gene=gene,order=order,variance=float(v),returned=True,
                            logmass=ll,gradient=g,hessian=relevant_h,variance_Hessian_not_applicable=v==0)
                        if old is not None:
                            dl=abs(ll-old['independent_scalar']['logmass'])
                            record['independent_scalar_LL_difference']=dl
                            record['independent_scalar_LL_passed']=bool(dl<=1e-8)
                            if old.get('candidate') is not None:
                                c=old['candidate'];record['older_candidate_differences']=dict(LL=abs(ll-c['logmass']),
                                    score=float(np.max(np.abs(g-c['gradient']))),Hessian=float(np.max(np.abs(relevant_h-c['hessian']))))
                        else:
                            c=workers[gene];diff=dict(LL=abs(ll-c['logmass']),score=float(np.max(np.abs(g-c['gradient']))),
                                Hessian=float(np.max(np.abs(h-c['hessian']))))
                            record['call12_reference_differences']=diff
                            record['call12_reference_passed']=bool(diff['LL']<=1e-9 and diff['score']<=1e-8 and diff['Hessian']<=1e-5)
                    except (ValueError,np.linalg.LinAlgError) as exc:
                        record=dict(case=label,gene=gene,order=order,returned=False,error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
                    finally:private.bb_terms=original_bb
                    record['seconds']=time.perf_counter()-start;records.append(record);allrecords.append(record)
                refinements=[]
                for first,second in zip(records,records[1:]):
                    if not (first['returned'] and second['returned']):
                        refinements.append(dict(orders=[first['order'],second['order']],available=False));continue
                    dg=np.abs(np.asarray(first['gradient'])-second['gradient']);dl=abs(first['logmass']-second['logmass'])
                    refinements.append(dict(orders=[first['order'],second['order']],available=True,
                        LL_difference=dl,score_max_difference=float(dg.max()),
                        Hessian_max_difference=float(np.max(np.abs(np.asarray(first['hessian'])-second['hessian']))),
                        owner_LL_score_passed=bool(dl<=1e-8+1e-8*abs(second['logmass'])
                            and np.all(dg<=1e-8+1e-6*np.abs(second['gradient'])))))
                require(private.bb_terms is original_bb and p.tobytes()==pbytes and counts()==COUNT_SHA
                    and snapshots.row_signatures(base.rows,np)==before,'Fixed state changed')
                point=dict(case=label,gene=gene,records=records,refinements=refinements);allpoints.append(point)
                save(label+'_gene'+str(gene),point)
        require(len(allpoints)==40 and len(allrecords)==112,'Incomplete prescribed family')
        totals=[];full_refinements=[]
        for label,p,v in cases:
            previous=None
            for order in ((20,) if v==0 else ORDERS):
                pieces=[r for r in allrecords if r['case']==label and r['order']==order]
                require([r['gene'] for r in pieces]==list(range(8)),'Full sum gene coverage differs')
                if not all(r['returned'] for r in pieces):
                    total=dict(case=label,order=order,available=False);totals.append(total)
                    if previous is not None:
                        full_refinements.append(dict(case=label,orders=[previous['order'],order],
                            available=False,owner_LL_score_passed=False))
                    previous=total;continue
                ll=0.;g=np.zeros(29);width=28 if v==0 else 29;h=np.zeros((width,width))
                for r in pieces:ll+=r['logmass'];g+=r['gradient'];h+=r['hessian']
                total=dict(case=label,order=order,available=True,logmass=ll,gradient=g,hessian=h,
                           variance_Hessian_not_applicable=v==0)
                totals.append(total)
                if label=='saved_call12':
                    differences=dict(LL=abs(ll-refanswer['logmass']),
                        score=float(np.max(np.abs(g-refanswer['gradient']))),
                        Hessian=float(np.max(np.abs(h-refanswer['hessian']))))
                    total['call12_full_reference_differences']=differences
                    total['call12_full_reference_passed']=bool(differences['LL']<=1e-9
                        and differences['score']<=1e-8 and differences['Hessian']<=1e-5)
                if previous is not None:
                    if not previous['available']:
                        full_refinements.append(dict(case=label,orders=[previous['order'],order],
                            available=False,owner_LL_score_passed=False));previous=total;continue
                    dl=abs(previous['logmass']-ll);dg=np.abs(previous['gradient']-g)
                    full_refinements.append(dict(case=label,orders=[previous['order'],order],available=True,
                        LL_difference=dl,score_max_difference=float(dg.max()),
                        Hessian_max_difference=float(np.max(np.abs(previous['hessian']-h))),
                        owner_LL_score_passed=bool(dl<=1e-8+1e-8*abs(ll)
                            and np.all(dg<=1e-8+1e-6*np.abs(g)))))
                previous=total
        save('full_sums',dict(totals=totals,refinements=full_refinements))
        save('summary',dict(completed=True,points=len(allpoints),integrals=len(allrecords),
            returned_integrals=sum(r['returned'] for r in allrecords),
            independent_LL_failures=[dict(case=r['case'],gene=r['gene'],order=r['order']) for r in allrecords
                if r.get('independent_scalar_LL_passed') is False],
            call12_derivative_failures=[dict(case=r['case'],gene=r['gene'],order=r['order']) for r in allrecords
                if r.get('call12_reference_passed') is False],
            refinement_failures=[dict(case=p['case'],gene=p['gene'],orders=q['orders']) for p in allpoints
                for q in p['refinements'] if not q.get('owner_LL_score_passed',False)],
            full_refinement_failures=[dict(case=q['case'],orders=q['orders']) for q in full_refinements
                if not q['owner_LL_score_passed']],
            call12_full_reference_failures=[q['order'] for q in totals if q['case']=='saved_call12'
                and not q.get('call12_full_reference_passed',False)],
            elapsed_seconds=time.perf_counter()-began,
            measured_process_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            all_fixed_inputs_unchanged=True,no_order_admitted=True,no_fit=True,no_calibration=True,no_independent_full_Hessian=True))
    except Exception as exc:
        save('failure',dict(completed=False,error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()));raise


if __name__=='__main__':main()
