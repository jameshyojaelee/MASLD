#!/usr/bin/env python3
"""Check scoped Gamma evaluator against all frozen native-GH family answers.

Also verify original zero-boundary bytes, inheritance and failure restoration.
No fitting or calibration is performed.
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
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FAMILY=REC/'repair_p3_gamma_native_GH_family_22079062'
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
    require(all(p in pins for p in [*(HERE/n for n in SOURCES),FAMILY/'protocol.json',FAMILY/'summary.json']),
            'Required source/reference missing')
    source=json.loads((FAMILY/'protocol.json').read_text());summary=json.loads((FAMILY/'summary.json').read_text())
    require(summary['completed'] and summary['points']==40 and summary['integrals']==summary['returned_integrals']==112
        and summary['all_fixed_inputs_unchanged'] and all(not summary[k] for k in
        ('independent_LL_failures','call12_derivative_failures','refinement_failures','full_refinement_failures',
         'call12_full_reference_failures')),'Native GH family incomplete or failed')
    for path,digest in source['source_and_input_sha256'].items():
        require(Path(path) in pins and pins[Path(path)]==digest,'Native family input/source differs')
    cases=source['cases'];records=[]
    for case in cases:
        for gene in range(8):
            path=FAMILY/(case+'_gene'+str(gene)+'.json');require(path in pins,'Point reference missing')
            point=json.loads(path.read_text());require(point['case']==case and point['gene']==gene,'Point identity differs')
            records.extend(point['records'])
    require(len(records)==112 and all(r['returned'] for r in records),'Full family records required')
    protocol=dict(source_and_input_sha256={str(p):v for p,v in pins.items()},
        reference_points=40,reference_integrals=112,byte_identity_required=True,
        zero_boundary_full_NaN_Hessian_bytes_checked=True,all_inherited_methods_unchanged=True,
        numeric_and_integrity_failure_probes=True,no_fit=True,no_calibration=True,no_actual_or_sealed_data=True)
    if args.inspect:print(json.dumps(protocol));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
            'Threads not fixed before imports')
    out=REC/('repair_p3_gamma_native_model_check_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in SOURCES)):
        dest=archive/p.name;dest.write_bytes(p.read_bytes())
        require(sha(dest)==(args.expected_checker_sha if p==Path(__file__) else pins[p]),'Archive differs')
    import numpy as np
    import scipy
    def save(name,value):
        body=(json.dumps(value,allow_nan=False,indent=2)+'\n').encode()
        require(sum(p.stat().st_size for p in out.rglob('*') if p.is_file())+len(body)<32*1024**2-65536,'Output cap exhausted')
        with (out/(name+'.json')).open('xb') as f:f.write(body);f.flush();os.fsync(f.fileno())
        print(json.dumps(dict(saved=name)),flush=True)
    began=time.perf_counter()
    try:
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        require(environment==source['environment'],'Environment differs');protocol['environment']=environment;save('protocol',protocol)
        sys.path.insert(0,str(archive));loader=importlib.import_module('bfix_v1');importlib.import_module('tests_v1')
        for n in SOURCES[:4]:
            m=sys.modules[n[:-3]];require(Path(m.__file__).resolve().parent==archive.resolve() and sha(Path(m.__file__))==pins[HERE/n],
                                        'Original import escaped archive')
        design,_=loader.load();private=load('native_model_original',archive/'likelihood_v1.py')
        gamma=load('native_model_gamma',archive/'repair_p3_bb_gamma_difference_v3.py')
        snapshots=load('native_model_snapshots',archive/'repair_p3_gamma_parallel_model.py')
        wrapper=load('native_model_wrapper',archive/'repair_p3_gamma_native_model.py')
        Model=wrapper.native_gamma_model_class(private,gamma,snapshots)
        rows=design.model().rows;original_bb=private.bb_terms;original_gamma=gamma.stable_bb_terms
        original_conditional=private.conditional_terms;before=snapshots.row_signatures(rows,np)
        points={case:loader.parameters(loader.FIXTURE/'params'/(case+'.json')) for case in cases if case!='saved_call12'}
        refp=json.loads((REC/'repair_p3_gamma_full_parallel_22077927/protocol.json').read_text())
        points['saved_call12']=(np.asarray(refp['parameters'],dtype=float),refp['variance'])
        checked=[];boundary_calls=0
        for r in records:
            p,v=points[r['case']];model=Model([row for row in rows if row.gene==r['gene']],gene_count=8)
            answer=model.marginal(p,v,nodes=r['order'],with_hessian=True)
            h=answer[2][:-1,:-1] if v==0 else answer[2]
            identical=(float(answer[0]).hex()==float(r['logmass']).hex()
                and answer[1].tobytes()==np.asarray(r['gradient'],dtype=float).tobytes()
                and h.tobytes()==np.asarray(r['hessian'],dtype=float).tobytes())
            if v==0:
                base=private.Model(model.rows,gene_count=8).marginal(p,v,nodes=r['order'],with_hessian=True)
                require(all(np.asarray(a).tobytes()==np.asarray(b).tobytes() for a,b in zip(answer,base)),
                        'Full original zero-boundary bytes differ');boundary_calls+=1
            require(identical and private.bb_terms is original_bb and snapshots.row_signatures(rows,np)==before,
                    'Scoped evaluator differs from checked native family or changes state')
            checked.append(dict(case=r['case'],gene=r['gene'],order=r['order'],byte_identical=True))
            if len(checked)%8==0:print(json.dumps(dict(checked_integrals=len(checked))),flush=True)
        save('byte_comparisons',checked)
        model=Model([row for row in rows if row.gene==4],gene_count=8)
        p,v=points['saved_call12'];probes={}
        def numeric_failure(*args,**kwargs):raise ValueError('injected numeric refusal')
        original_hermgauss=private.np.polynomial.hermite.hermgauss
        try:
            private.np.polynomial.hermite.hermgauss=numeric_failure
            try:model.marginal(p,v,with_hessian=True)
            except ValueError as exc:
                probes['numeric_refusal_preserved']=str(exc)=='injected numeric refusal'
                probes['numeric_state_intact']=model.last_integration_failure['fixed_state_unchanged'] is True
            else:probes['numeric_refusal_preserved']=False
        finally:private.np.polynomial.hermite.hermgauss=original_hermgauss
        probes['BB_restored_after_numeric_failure']=private.bb_terms is original_bb
        def corrupt_failure(*args,**kwargs):
            private.conditional_terms=numeric_failure
            raise ValueError('injected state corruption')
        try:
            private.np.polynomial.hermite.hermgauss=corrupt_failure
            try:model.marginal(p,v,with_hessian=True)
            except wrapper.ModelIntegrityError:
                probes['state_corruption_overrides_numeric_refusal']=model.last_integration_failure['fixed_state_unchanged'] is False
            else:probes['state_corruption_overrides_numeric_refusal']=False
        finally:
            private.np.polynomial.hermite.hermgauss=original_hermgauss;private.conditional_terms=original_conditional
        probes['BB_restored_after_integrity_failure']=private.bb_terms is original_bb
        recovered=model.marginal(p,v,with_hessian=True)
        probes['recovered_positive_call_finite']=all(np.isfinite(a).all() for a in recovered)
        probes['failure_metadata_cleared_on_recovery']=model.last_integration_failure is None
        try:
            gamma.stable_bb_terms=numeric_failure
            try:model.marginal(p,v,with_hessian=True)
            except wrapper.ModelIntegrityError:probes['changed_Gamma_callable_refused']=True
            else:probes['changed_Gamma_callable_refused']=False
        finally:gamma.stable_bb_terms=original_gamma
        try:
            gamma.stable_bb_terms=numeric_failure
            base=private.Model(model.rows,gene_count=8)
            for mode in (False,True):
                expected=base.marginal(p,0.,with_hessian=mode);actual=model.marginal(p,0.,with_hessian=mode)
                probes['boundary_bypasses_Gamma_'+str(mode)]=all(np.asarray(a).tobytes()==np.asarray(b).tobytes()
                    for a,b in zip(actual,expected))
        finally:gamma.stable_bb_terms=original_gamma
        def failing_sink(event):raise OSError('injected output failure')
        SinkModel=wrapper.native_gamma_model_class(private,gamma,snapshots,event_sink=failing_sink)
        try:SinkModel(model.rows,gene_count=8).marginal(p,v,with_hessian=False)
        except wrapper.ModelIntegrityError as exc:probes['event_failure_is_integrity_error']='Event callback failed' in str(exc)
        else:probes['event_failure_is_integrity_error']=False
        probes['BB_restored_after_event_failure']=private.bb_terms is original_bb
        score_only=model.marginal(p,v,with_hessian=False)
        probes['score_only_original_two_item_interface']=len(score_only)==2 and bool(np.isfinite(score_only[0]))
        probes['score_only_agrees_with_Hessian_call']=all(np.asarray(a).tobytes()==np.asarray(b).tobytes()
            for a,b in zip(score_only,recovered[:2]))
        probes['all_inherited_methods_unchanged']=all(getattr(Model,n) is getattr(private.Model,n) for n in wrapper.INHERITED)
        probes['final_original_BB_and_rows_unchanged']=private.bb_terms is original_bb and snapshots.row_signatures(rows,np)==before
        save('state_probes',probes);require(all(probes.values()),'State/failure probe failed')
        for path,digest in pins.items():require(sha(path)==digest,'Pinned input changed during check')
        save('summary',dict(completed=True,byte_identical_integrals=len(checked),original_boundary_full_byte_calls=boundary_calls,
            all_state_probes_passed=True,elapsed_seconds=time.perf_counter()-began,
            no_fit=True,no_calibration=True,no_actual_or_sealed_data=True))
    except Exception as exc:
        save('failure',dict(completed=False,error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()));raise


if __name__=='__main__':main()
