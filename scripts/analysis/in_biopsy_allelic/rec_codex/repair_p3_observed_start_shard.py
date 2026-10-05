#!/usr/bin/env python3
"""Future synthetic warm/cold starts; no submission and no completed-fit claim.

Reuse the pinned observed runner's exact setup and inherited single-start fitter.
Only its outer fit calls are replaced. Completed shared null packets are required.
"""
import argparse
import ast
import hashlib
import json
import math
import os
from pathlib import Path
import sys

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
BASE = REC/'repair_p3_batched_observed_fit_launch_20261005T152152Z/repair_p3_batched_observed_fit.py'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('self', 'base', 'model', 'null', 'protocol'):
        parser.add_argument('--expected-'+name+'-sha', required=True)
    parser.add_argument('--null-dir', type=Path, required=True)
    parser.add_argument('--start', choices=('warm', 'cold'), required=True)
    parser.add_argument('--batch-size', type=int, choices=(32, 64), required=True)
    parser.add_argument('--inspect', action='store_true')
    args = parser.parse_args()
    own = Path(__file__).resolve()
    null_dir = args.null_dir.resolve()
    require(null_dir.is_relative_to(REC.resolve()), 'Shared null must be inside synthetic REC')
    null_path, protocol_path = null_dir/'observed_null.json', null_dir/'protocol.json'
    pins = {own:args.expected_self_sha, BASE:args.expected_base_sha,
            null_path:args.expected_null_sha, protocol_path:args.expected_protocol_sha}
    for path, digest in pins.items():
        require(sha(path)==digest, 'Source/shared-null pin differs: '+str(path))
    null = json.loads(null_path.read_text())['fit']
    protocol = json.loads(protocol_path.read_text())
    require(null['converged'] is True and len(null['parameters'])==28
            and all(isinstance(x, (int, float)) and math.isfinite(x) for x in null['parameters']),
            'Completed converged finite28-coordinate null required')
    require(protocol['only_synthetic_observed_fit'] is True
            and protocol['participants']==60 and protocol['rows']==173 and protocol['genes']==8
            and protocol['requested_batch_size']==args.batch_size,
            'Shared-null synthetic protocol/batch differs')
    ast.parse(BASE.read_text())
    import importlib.util
    import platform
    import time
    spec = importlib.util.spec_from_file_location('_observed_start_setup', BASE)
    setup = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(setup)
    setup_args = argparse.Namespace(expected_runner_sha=args.expected_base_sha,
                                   expected_model_sha=args.expected_model_sha)
    campaign, guards = setup.metadata(setup_args)
    for name in ('globals.json','individuals.tsv','rows.tsv','cohorts.tsv'):
        path = setup.FIXTURE/name
        guards[path] = campaign['files'][str(path)]
    setup.guard(guards)
    for path,digest in guards.items():
        if path.suffix=='.py' and path!=BASE:
            require(protocol['source_guards'].get(str(path))==digest,
                    'Shared null source differs: '+str(path))
    if args.inspect:
        print(json.dumps(dict(pins_checked=True,completed_shared_null_available=True,
            start=args.start,numerical_execution=False,whole_fit_complete=False,
            deferred_boundary_phases=True)))
        return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required before numerical imports')
    output = REC/('repair_p3_observed_start_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(mode=0o750,exist_ok=False)
    archive = output/'executed_sources'; archive.mkdir()
    for path,digest in {**guards,**pins}.items():
        destination = archive/path.name if path.suffix=='.py' else output/('shared_'+digest[:12]+'_'+path.name)
        require(not destination.exists(),'Duplicate archive basename')
        destination.write_bytes(path.read_bytes())
        require(sha(destination)==digest,'Archive copy differs')
    np = None
    def safe(value):
        if np is not None and isinstance(value,np.ndarray): return safe(value.tolist())
        if np is not None and isinstance(value,np.generic): return safe(value.item())
        if isinstance(value,float) and not math.isfinite(value):
            return {'nonfinite_float':str(value)}
        if isinstance(value,dict): return {str(k):safe(v) for k,v in value.items()}
        if isinstance(value,(tuple,list)): return [safe(v) for v in value]
        return value
    def save(name,value):
        body = (json.dumps(safe(value),allow_nan=False,indent=2)+'\n').encode()
        require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file())+len(body)
                <=32*1024**2,'Saved-output32MiB cap exceeded')
        with (output/(name+'.json')).open('xb') as handle:
            handle.write(body); handle.flush(); os.fsync(handle.fileno())
    save('start_protocol',dict(start=args.start,
        pins={str(p):h for p,h in {**guards,**pins}.items()},environment=protocol['environment'],
        shared_null_sha=args.expected_null_sha,shared_protocol_sha=args.expected_protocol_sha,
        base_sha=args.expected_base_sha,model_sha=args.expected_model_sha,batch_size=args.batch_size,
        observed_count_vector_sha256=protocol['observed_count_vector_sha256'],
        kappa_attempt_order=[0.,-.1,.1],fixed_v=None,whole_fit_complete=False,
        deferred=['warm/cold selection','conditional boundary restart','null snap','final diagnostics']))
    model = None
    try:
        sys.path.insert(0,str(archive))
        require(not any(n in sys.modules for n in ('likelihood_v1','bfix_v1','design_v1','tests_v1')),
                'Numerical fixture modules already imported')
        import numpy as np
        import scipy
        import importlib
        loader = importlib.import_module('bfix_v1')
        original = importlib.import_module('likelihood_v1')
        importlib.import_module('tests_v1')
        for name in ('likelihood_v1','bfix_v1','design_v1','tests_v1'):
            require(Path(sys.modules[name].__file__).resolve().parent==archive.resolve(),
                    'Original import escaped archive')
        environment = dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        require(environment==dict(python='3.10.18',numpy='2.2.6',scipy='1.15.2'),
                'Original numerical environment changed')
        require(environment==protocol['environment'],'Shared-null numerical environment differs')
        design,settings = loader.load()
        base = design.model()
        require(len(design.participants)==60 and len(base.rows)==173 and base.G==8
                and base.Z==7 and base.size==28,'Synthetic fixture dimensions changed')
        old_bb = original.bb_terms
        private = setup.load_module('_start_private',archive/'likelihood_v1.py')
        kernel = setup.load_module('_start_kernel',archive/setup.KERNEL.name)
        gl = setup.load_module('_start_gl',archive/setup.GL.name)
        conditional = setup.load_module('_start_conditional',archive/setup.CONDITIONAL.name)
        model_source = setup.load_module('_start_model',archive/setup.MODEL.name)
        Candidate = model_source.batched_model_class(private,kernel,gl,conditional,
            kernel_path=archive/setup.KERNEL.name,batch_size=args.batch_size)
        model = Candidate(base.rows,gene_count=base.G)
        require(original.bb_terms is old_bb,'Private candidate leaked into reference')
        require(all(getattr(Candidate,n) is getattr(private.Model,n) for n in model_source.INHERITED),
                'Inherited fitter method changed')
        counts = hashlib.sha256(np.asarray([row.a for row in model.rows],dtype='<i8').tobytes()).hexdigest()
        require(counts==protocol['observed_count_vector_sha256'],'Shared null counts differ')
        p0 = np.asarray(null['parameters'],float)
        start_vector = np.r_[p0,0.] if args.start=='warm' else np.r_[model.cold_start(),settings['tau_max_fixture']**2]
        if args.start=='cold': start_vector[:-1][~model.informative()] = p0[~model.informative()]
        save('start_input',dict(start=args.start,start_vector=start_vector,
                               whole_fit_complete=False))
        began = time.perf_counter()
        fit = model.fit_marginal_single(start_vector,None)
        save('start_result',dict(start=args.start,fit=fit,seconds=time.perf_counter()-began,
                                whole_fit_complete=False))
    except BaseException as exc:
        save('failure',dict(start=args.start,error_type=type(exc).__name__,message=str(exc),
            whole_fit_complete=False,last_evaluated_point=getattr(model,'last_evaluated_point',None),
            last_integration_failure=getattr(model,'last_integration_failure',None)))
        raise


if __name__=='__main__':
    main()
