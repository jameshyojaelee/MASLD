#!/usr/bin/env python3
"""Compute-only matched-domain diagnostic for three preserved BFIX refusals."""
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

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
FROZEN = REC/'repair_p3_bfix_vectors_launch_20261005T144634Z'
HERE = Path(__file__).resolve().parent
THETA = 'theta4_perturbed'
CAP = 32*1024**2
OLD = {
    1:(22040217,15.,'63581cfe0ffd63d4e6f5445705fd50b99cd99ba9e3669dde54bbd039b403d2f2'),
    2:(22040218,16.,'f4a4e82a99d719ecf695e6fceb67ce3ea14aeb3c4f895c73a4d76da6be7842dc'),
    3:(22040219,16.,'d3d2fb0579ef0c6f4bb95c03daf4e1fc9fd15c2ef7d5ad2989475485200b4fac'),
}
KERNEL_SHA = '461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'
SOURCE_MANIFEST_SHA = 'b6366a979cbe250dedcb5d5c443f4068cb0528fcc5ffa398ffda2f9b3b1c4520'
PARENT_PROTOCOL_SHA = {
    1:'aa474755d353c90a4b897effcd915e0b300015e9ebae8dcc62ba8ebcb938a463',
    2:'bc35c4c185310b0895464071cc49e0af860bc83d842319d78b9129d1dac72512',
    3:'ff8b0e4bef2728088841cd27a138ff46156233a4cec3a578c5df281dc8804732',
}
MODULES = ('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py')


def require(ok,message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name,path):
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gene',type=int,choices=(1,2,3),required=True)
    parser.add_argument('--expected-checker-sha',required=True)
    parser.add_argument('--expected-diagnostic-sha',required=True)
    parser.add_argument('--inspect',action='store_true')
    args = parser.parse_args()
    job,bound,result_sha = OLD[args.gene]
    parent = REC/f'repair_p3_bfix_vectors_{job}'
    old_path = parent/f'{THETA}_gene{args.gene}.json'
    require(sha(old_path) == result_sha,'Preserved refusal result changed')
    require(sha(parent/'protocol.json') == PARENT_PROTOCOL_SHA[args.gene], 'Preserved parent protocol changed')
    require(sha(FROZEN/'source_hashes.json') == SOURCE_MANIFEST_SHA, 'Frozen source manifest changed')
    old = json.loads(old_path.read_text())
    parent_protocol = json.loads((parent/'protocol.json').read_text())
    require(old['theta'] == THETA and old['gene'] == args.gene and 'candidate_refusal' in old and
            old['candidate_refusal']['partial_reports'][-1]['z_bound'] == bound,
            'Wrong original refused point/domain')
    pins = {old_path:result_sha,Path(__file__):args.expected_checker_sha,
            HERE/'repair_p3_gl_order_diagnostic.py':args.expected_diagnostic_sha,
            parent/'protocol.json':PARENT_PROTOCOL_SHA[args.gene],
            FROZEN/'source_hashes.json':SOURCE_MANIFEST_SHA}
    frozen_sources = json.loads((FROZEN/'source_hashes.json').read_text())
    require(frozen_sources['repair_p3_integral.py'] == KERNEL_SHA,'Frozen kernel identity differs')
    pins[FROZEN/'repair_p3_integral.py'] = KERNEL_SHA
    pins[FROZEN/'repair_p3_composite_gl.py'] = frozen_sources['repair_p3_composite_gl.py']
    pins[FROZEN/'repair_p3_bfix_vectors.py'] = frozen_sources['repair_p3_bfix_vectors.py']
    for name in MODULES:
        pins[FROZEN/name] = frozen_sources[name]
    # All exact small synthetic fixture inputs from the completed parent.
    for path,digest in parent_protocol['source_and_input_sha256'].items():
        path = Path(path)
        if path.resolve().is_relative_to(FIXTURE.resolve()):
            pins[path] = digest
    for path,digest in pins.items():
        require(sha(path) == digest,'Reviewed source/fixture bytes differ: '+path.name)
        if path.suffix == '.py':
            ast.parse(path.read_text())
    protocol = dict(theta=THETA,gene=args.gene,matched_z_bound=bound,grid_spacing=.25,
        orders=[16,32,64,128],expected_likelihood_evaluations=int(2*bound/.25)*240,
        source_and_input_sha256={str(p):h for p,h in pins.items()},
        original_refusal_preserved=True,no_accuracy_acceptance=True,no_model_fit=True,no_calibration=True,
        scope='One fixed synthetic theta4/gene, matched retained expanded domain only',output_cap_bytes=CAP)
    if args.inspect:
        print(json.dumps(protocol,sort_keys=True))
        return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required before numerical imports')
    output = REC/f'repair_p3_theta4_order_gene{args.gene}_{os.environ["SLURM_JOB_ID"]}'
    output.mkdir(exist_ok=False)
    archive = output/'executed_sources'
    archive.mkdir()
    for path,digest in pins.items():
        if path.suffix == '.py':
            target = archive/path.name
            target.write_bytes(path.read_bytes())
            require(sha(target) == digest,'Archived source differs')
    sys.path.insert(0,str(archive))
    import numpy as np
    import scipy
    versions = dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    require(versions == parent_protocol['environment'],'Original environment differs')

    def diagnostic(value):
        if isinstance(value,np.ndarray):
            return diagnostic(value.tolist())
        if isinstance(value,np.generic):
            return diagnostic(value.item())
        if isinstance(value,dict):
            return {k:diagnostic(v) for k,v in value.items()}
        if isinstance(value,(list,tuple)):
            return [diagnostic(v) for v in value]
        if isinstance(value,float) and not np.isfinite(value):
            return dict(nonfinite_float='nan' if np.isnan(value) else 'positive_infinity' if value>0 else 'negative_infinity')
        return value

    def save(name,value):
        body = json.dumps(diagnostic(value),indent=2,allow_nan=False)+'\n'
        require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file())+len(body.encode()) <= CAP,
                'Diagnostic output cap exceeded')
        with (output/(name+'.json')).open('x') as handle:
            handle.write(body)

    began = time.monotonic()
    protocol['environment'] = versions
    protocol['fixture_seed'] = json.loads((FIXTURE/'globals.json').read_text())['seed']
    save('protocol',protocol)
    save('original_refusal',old)
    try:
        require(not any(n[:-3] in sys.modules for n in MODULES),'Original modules already imported')
        original = importlib.import_module('likelihood_v1')
        loader = importlib.import_module('bfix_v1')
        importlib.import_module('tests_v1')
        for name in MODULES:
            path = Path(sys.modules[name[:-3]].__file__)
            require(path.resolve().parent == archive.resolve() and sha(path) == frozen_sources[name],
                    'Original import escaped pinned archive')
        private = load('theta4_order_private_original',archive/'likelihood_v1.py')
        kernel = load('theta4_order_private_kernel',archive/'repair_p3_integral.py')
        helper = load('theta4_order_private_diagnostic',archive/'repair_p3_gl_order_diagnostic.py')
        design,settings = loader.load()
        model = design.model()
        p,v = loader.parameters(FIXTURE/'params'/(THETA+'.json'))
        require(model.G == 8 and model.Z == 7 and model.size == 28 and len(model.rows) == 173 and
                np.array_equal(p,np.asarray(old['p'])) and v == old['v'], 'Synthetic dimensions/point differ')
        rows = [r for r in model.rows if r.gene == args.gene]
        require(len(rows) == old['rows'],'Fixed gene rows differ')
        private.bb_terms = kernel.stable_bb_terms
        candidate = private.Model(rows,gene_count=8)
        report = helper.diagnostic_ladder(candidate,p,v,rows,private.conditional_terms,kernel,bound)
        require(report['likelihood_evaluations'] == protocol['expected_likelihood_evaluations'],
                'Matched-domain evaluation count differs')
        save('order_diagnostic',dict(p=p,v=v,rows=len(rows),report=report))
        for path,digest in pins.items():
            require(sha(path) == digest,'Source or fixture changed during diagnostic')
        save('summary',dict(diagnostic_completed=True,original_refusal_preserved=True,
            elapsed_seconds=time.monotonic()-began,likelihood_evaluations=report['likelihood_evaluations'],
            no_accuracy_acceptance=True,independent_derivative_accuracy_verified=False,
            no_domain_or_grid_refinement=True,no_model_fit=True,no_calibration=True))
    except BaseException:
        save('failure',dict(diagnostic_completed=False,traceback=traceback.format_exc(),
                            original_refusal_preserved=True,no_accuracy_acceptance=True))
        raise


if __name__ == '__main__':
    main()
