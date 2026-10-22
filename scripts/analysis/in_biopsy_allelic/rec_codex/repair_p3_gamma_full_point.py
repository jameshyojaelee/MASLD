#!/usr/bin/env python3
"""Fixed failed synthetic point: independent sequential/parallel GL jobs.

No optimization, calibration, actual outcomes or scientific adoption.
"""
import argparse
import ast
from dataclasses import replace
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

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE = ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
BOUNDARY = REC/'repair_p3_gamma_boundary_check_22077922'
POINT = REC/'repair_p3_observed_failure_anchor_22042612'
COUNT_SHA = '7009be4adeec6157fa4bbf3f46e48ec24aa398a0bee4d7906c78c231b7d3fbad'
SOURCES = ('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
           'repair_p3_integral.py','repair_p3_composite_gl.py',
           'repair_p3_bb_gamma_difference_v3.py','repair_p3_gamma_parallel_model.py')
CAP = 32*1024**2


def require(ok,message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name,path):
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def numerical_refusal(exc):
    # GL loads a separate kernel module, so class identity is not portable.
    return (isinstance(exc,ValueError)
            and isinstance(getattr(exc,'integration_reports',None),list)
            and isinstance(getattr(exc,'likelihood_evaluations',None),int))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('sequential','parallel'),required=True)
    parser.add_argument('--pins',type=Path,required=True)
    parser.add_argument('--expected-pins-sha',required=True)
    parser.add_argument('--expected-checker-sha',required=True)
    parser.add_argument('--inspect',action='store_true')
    args = parser.parse_args()
    require(sha(Path(__file__)) == args.expected_checker_sha
            and sha(args.pins) == args.expected_pins_sha,'Executable or manifest changed')
    pins = {Path(p):h for p,h in json.loads(args.pins.read_text()).items()}
    for path,digest in pins.items():
        require(path.is_absolute() and path.resolve().is_relative_to(ROOT.resolve())
                and sha(path) == digest,'Pinned input differs: '+str(path))
        if path.suffix == '.py':
            ast.parse(path.read_text())
    require(all(HERE/n in pins for n in SOURCES)
            and all(BOUNDARY/n in pins for n in ('protocol.json','summary.json','preflight_integrity.json'))
            and POINT/'protocol.json' in pins,'Required inputs missing')
    admitted = json.loads((BOUNDARY/'summary.json').read_text())
    boundary = json.loads((BOUNDARY/'protocol.json').read_text())
    preflight = json.loads((BOUNDARY/'preflight_integrity.json').read_text())
    require(all(admitted[k] is True for k in ('completed','all_boundary_nominal_criteria_met',
        'all_parent_zero_variance_calls_byte_identical','original_parent_BB_unchanged'))
        and admitted['records'] == 9 and all(x is True for x in preflight.values()),
        'Corrected source boundary admission incomplete')
    require(boundary['source_and_input_sha256'][str(HERE/'repair_p3_gamma_parallel_model.py')]
            == pins[HERE/'repair_p3_gamma_parallel_model.py'],'Boundary admitted a different model')
    point = json.loads((POINT/'protocol.json').read_text())
    require(point['count_vector_sha256'] == COUNT_SHA and point['call'] == 12
            and point['fixed_v'] == .0009636462009520182,'Saved failed point differs')
    protocol = dict(mode=args.mode,source_and_input_sha256={str(p):h for p,h in pins.items()},
        parameters=point['fixed_p'],variance=point['fixed_v'],count_vector_sha256=COUNT_SHA,
        participants=60,rows=173,G=8,Z=7,P=28,fixture_seed=20260930,with_hessian=True,
        settings='original scalar GL; per-gene LL min(1e-9,1e-8/8), score min(1e-8,1e-7/8), H min(1e-5,1e-5/8)',
        global_LL_allowance=1e-8,global_score_allowance=1e-7,global_Hessian_allowance=1e-5,
        comparison_LL_atol=1e-9,comparison_score_atol=1e-8,comparison_Hessian_atol=1e-5,
        no_fit=True,no_calibration=True,no_independent_full_Hessian=True,no_actual_or_sealed_outcomes=True,
        output_cap_bytes=CAP,retain_numerical_refusals=True)
    if args.inspect:
        print(json.dumps(protocol))
        return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    require(all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS',
                'MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),'Threads must be fixed before imports')
    out = REC/('repair_p3_gamma_full_'+args.mode+'_'+os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)
    archive = out/'executed_sources'
    archive.mkdir()
    for path in (Path(__file__),*(HERE/n for n in SOURCES)):
        target = archive/path.name
        target.write_bytes(path.read_bytes())
        require(sha(target) == (args.expected_checker_sha if path == Path(__file__) else pins[path]),
                'Executable archive differs')
    began = time.perf_counter()
    import numpy as np
    import scipy
    def serialize(x):
        if isinstance(x,np.ndarray):
            return x.tolist()
        if isinstance(x,np.generic):
            return x.item()
        raise TypeError(type(x).__name__)
    def refusal_json(x):
        # Preserve failed-domain nonfinite diagnostics explicitly. Successful
        # results still use strict finite JSON and are never sanitized.
        if isinstance(x,np.ndarray):
            return refusal_json(x.tolist())
        if isinstance(x,np.generic):
            return refusal_json(x.item())
        if isinstance(x,float) and not math.isfinite(x):
            return {'nonfinite_float':x.hex()}
        if isinstance(x,dict):
            return {k:refusal_json(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):
            return [refusal_json(v) for v in x]
        return x
    def save(name,value):
        body = (json.dumps(value,default=serialize,allow_nan=False,indent=2)+'\n').encode()
        require(sum(p.stat().st_size for p in out.rglob('*') if p.is_file())+len(body)<CAP-65536,
                'Output cap exhausted; refusal reserve retained')
        with (out/(name+'.json')).open('xb') as f:
            f.write(body);f.flush();os.fsync(f.fileno())
        print(json.dumps(dict(saved=name,seconds=time.perf_counter()-began)),flush=True)
    try:
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        require(environment == point['environment'],'Numerical environment differs')
        protocol['environment'] = environment
        save('protocol',protocol)
        sys.path.insert(0,str(archive))
        require(not any(n[:-3] in sys.modules for n in SOURCES[:4]),'Original modules already imported')
        loader = importlib.import_module('bfix_v1')
        importlib.import_module('tests_v1')
        for n in SOURCES[:4]:
            module=sys.modules[n[:-3]]
            require(Path(module.__file__).resolve().parent==archive.resolve()
                    and sha(Path(module.__file__))==pins[HERE/n],'Original import escaped archive')
        design,settings = loader.load()
        private = load('full_point_original',archive/'likelihood_v1.py')
        kernel = load('full_point_kernel',archive/'repair_p3_integral.py')
        gl = load('full_point_gl',archive/'repair_p3_composite_gl.py')
        bb = load('full_point_bb',archive/'repair_p3_bb_gamma_difference_v3.py')
        wrapper = load('full_point_wrapper',archive/'repair_p3_gamma_parallel_model.py')
        base = private.Model(design.model().rows,gene_count=8)
        require(len(base.rows)==173 and base.Z==7 and base.size==28 and len(design.participants)==60
                and settings['seed']==20260930,'Original synthetic dimensions differ')
        counts=lambda:hashlib.sha256(np.asarray([r.a for r in base.rows],dtype='<i8').tobytes()).hexdigest()
        require(counts()==COUNT_SHA,'Original synthetic count vector differs')
        signatures=wrapper.row_signatures(base.rows,np)
        original_bb=private.bb_terms
        p=np.asarray(point['fixed_p'],dtype=float);p.setflags(write=False)
        variance=float(point['fixed_v'])
        if args.mode=='parallel':
            Model=wrapper.gamma_parallel_model_class(private,kernel,gl,bb,kernel_path=archive/'repair_p3_integral.py')
            model=Model(base.rows,gene_count=8)
            terms={}
            try:
                ll,g,h=model.marginal(p,variance,with_hessian=True,gene_terms=terms)
                answer=dict(returned=True,logmass=ll,gradient=g,hessian=h,gene_terms=terms,
                    reports=model.last_integration_reports,error_allowances=model.last_total_error_allowance)
            except Exception as exc:
                answer=dict(returned=False,error_type=type(exc).__name__,message=str(exc),
                    numerical_refusal=numerical_refusal(exc),traceback=traceback.format_exc(),
                    integrity_failure=not numerical_refusal(exc),
                    failure=getattr(model,'last_integration_failure',None),gene_terms=terms)
            worker_records=getattr(model,'last_parallel_worker_records',[])
            save('workers',worker_records if answer['returned'] else refusal_json(worker_records))
        else:
            per_gene=replace(kernel.Settings(),logmass_tolerance=1e-9,
                             gradient_tolerance=1e-8,hessian_tolerance=1e-5/8)
            records=[]
            private.bb_terms=bb.stable_bb_terms
            try:
                for gene in range(8):
                    rows=tuple(row for row in base.rows if row.gene==gene)
                    start=time.perf_counter()
                    try:
                        ll,g,h,report=gl.integrate_gene(base,p,variance,rows,private.conditional_terms,
                            with_hessian=True,expected_kernel_sha=wrapper.KERNEL_SHA,
                            kernel_path=archive/'repair_p3_integral.py',settings=per_gene)
                        record=dict(gene=gene,success=True,logmass=ll,gradient=g,hessian=h,
                                    report=dict(gene=gene,**report))
                    except Exception as exc:
                        record=dict(gene=gene,success=False,error_type=type(exc).__name__,message=str(exc),
                            numerical_refusal=numerical_refusal(exc),
                            refused_gene_ladders=getattr(exc,'integration_reports',None),
                            refused_gene_evaluations=getattr(exc,'likelihood_evaluations',None),
                            traceback=traceback.format_exc())
                    record['seconds']=time.perf_counter()-start
                    records.append(record)
                    save('gene'+str(gene),record if record['success'] else refusal_json(record))
            finally:
                private.bb_terms=original_bb
            failed=[r for r in records if not r['success']]
            if failed:
                answer=dict(returned=False,first_failed_gene=failed[0]['gene'],all_genes_retained=True,
                    integrity_failure=any(not r['numerical_refusal'] for r in failed))
            else:
                reports=[r['report'] for r in records]
                errors=dict(logmass=sum(r['refined']['logmass_error_allowance']+r['logmass_refinement_difference'] for r in reports),
                    gradient=sum((r['refined']['gradient_error_allowance']+r['gradient_refinement_difference'] for r in reports),np.zeros(29)),
                    hessian=sum((r['refined']['hessian_error_allowance']+r['hessian_refinement_difference'] for r in reports),np.zeros((29,29))))
                valid=all(np.isfinite(x).all() and np.all(np.asarray(x)>=0) for x in errors.values())
                passed=valid and errors['logmass']<=1e-8 and np.max(errors['gradient'])<=1e-7 and np.max(errors['hessian'])<=1e-5
                answer=dict(returned=bool(passed),all_genes_retained=True,error_allowances=errors,
                    reports=reports,full_sum_formed_but_not_admitted=not bool(passed),
                    logmass=sum(r['logmass'] for r in records),
                    gradient=sum((r['gradient'] for r in records),np.zeros(29)),
                    hessian=sum((r['hessian'] for r in records),np.zeros((29,29))),
                    gene_terms={r['gene']:r['logmass'] for r in records})
        save('answer',answer if answer['returned'] else refusal_json(answer))
        require(not answer.get('integrity_failure',False),
                'Evaluator integrity failure retained; do not count as a numerical refusal')
        require(private.bb_terms is original_bb and counts()==COUNT_SHA
                and wrapper.row_signatures(base.rows,np)==signatures,'Original fixed inputs changed')
        for path,digest in pins.items():
            require(sha(path)==digest,'Input/source changed while running: '+str(path))
        save('summary',dict(completed=True,mode=args.mode,returned=answer['returned'],
            all_fixed_inputs_unchanged=True,elapsed_seconds=time.perf_counter()-began,
            measured_process_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            no_fit=True,no_calibration=True,no_independent_full_Hessian=True))
    except Exception as exc:
        save('failure',dict(error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()))
        raise


if __name__=='__main__':
    main()
