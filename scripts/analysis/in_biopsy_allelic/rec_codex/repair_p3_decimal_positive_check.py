#!/usr/bin/env python3
"""Pinned synthetic BFIX theta4 genes 1/2/3 Decimal scalar/Hvv checks only."""
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
HERE = Path(__file__).resolve().parent
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
CAMPAIGN = REC/'p3_parallel_20261003T130323Z'
HELPER = HERE/'repair_p3_decimal_positive.py'
BB = HERE/'repair_p3_boundary_decimal_fast.py'
BB_SHA = 'd192f1aa311b54d3fe94af3744cf14a8c02856a95ffc5e3e590c7efa1432c017'
MODULES = ('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py')
CAP = 32*1024**2


def require(ok,message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


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
    parser.add_argument('--expected-helper-sha',required=True)
    parser.add_argument('--inspect',action='store_true')
    args = parser.parse_args()
    pins = {Path(__file__):args.expected_checker_sha,HELPER:args.expected_helper_sha,BB:BB_SHA,
        FIXTURE/'SHA256SUMS':'0d3d06fc32a8e588d5e1c54c943a582e5054b6d2ab565290a4d57ced545d2925',
        CAMPAIGN/'campaign.json':'bb3e82833476f67b0b11d2388e2dc00358eb987367b7a6a666d11938c6baca94'}
    for path,digest in pins.items():
        require(sha(path)==digest,'Reviewed source or metadata differs: '+path.name)
    campaign = json.loads((CAMPAIGN/'campaign.json').read_text())
    for name in MODULES:
        pins[HERE/name] = campaign['files'][str(CAMPAIGN/'executed_sources'/name)]
    published = {name.removeprefix('./'):digest for digest,name in
        (line.split(maxsplit=1) for line in (FIXTURE/'SHA256SUMS').read_text().splitlines())}
    for name in ('globals.json','individuals.tsv','rows.tsv','cohorts.tsv','genes.tsv','tags.tsv',
                 'params/theta4_perturbed.json'):
        pins[FIXTURE/name] = published[name]
    for path,digest in pins.items():
        require(sha(path)==digest,'Pinned source or synthetic input differs: '+path.name)
        if path.suffix=='.py':
            ast.parse(path.read_text())
    protocol = dict(theta='theta4_perturbed',gene=args.gene,orders=[20,40,80],decimal_digits=[80,110],
        maximum_probability_evaluations=100000,output_cap_bytes=CAP,
        source_and_input_sha256={str(p):d for p,d in pins.items()},
        no_model_fits=True,no_calibration=True,no_private_data=True,
        full_Hessian_accuracy_verified=False,independent_Hvv_accuracy_accepted=False)
    if args.inspect:
        print(json.dumps(protocol,sort_keys=True))
        return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required before numerical imports')
    output = REC/('repair_p3_decimal_positive_gene%d_%s' % (args.gene,os.environ['SLURM_JOB_ID']))
    output.mkdir(exist_ok=False)
    archive = output/'executed_sources'
    archive.mkdir()
    for source,digest in pins.items():
        if source.suffix=='.py':
            target = archive/source.name
            target.write_bytes(source.read_bytes())
            require(sha(target)==digest,'Archive copy differs')
    began = time.monotonic()

    def save(name,value):
        body = json.dumps(value,indent=2,allow_nan=False)+'\n'
        require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file())+len(body.encode())<=CAP-65536,
                'Diagnostic output cap exhausted; failure reserve retained')
        partial,final = output/(name+'.json.partial'),output/(name+'.json')
        with partial.open('x') as handle:
            handle.write(body)
            handle.flush()
            os.fsync(handle.fileno())
        require(not final.exists(),'Refusing result overwrite')
        partial.rename(final)
        print(json.dumps(dict(saved=name,elapsed_seconds=time.monotonic()-began)),flush=True)

    try:
        sys.path.insert(0,str(archive))
        require(not any(n[:-3] in sys.modules for n in MODULES),'Original numerical module already imported')
        import numpy as np
        import scipy
        loader = importlib.import_module('bfix_v1')
        importlib.import_module('tests_v1')
        for name in MODULES:
            module = sys.modules[name[:-3]]
            require(Path(module.__file__).resolve().parent==archive.resolve() and
                    sha(module.__file__)==pins[HERE/name],'Import escaped pinned archive')
        helper = load('decimal_positive_archived',archive/HELPER.name)
        bb = load('decimal_positive_BB_archived',archive/BB.name)
        require(helper.ORDERS==(20,40,80) and helper.DIGITS==(80,110) and
                helper.MAXIMUM_EVALUATIONS==100000,'Declared independent diagnostic family differs')
        versions = dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        require(versions==dict(python='3.10.18',numpy='2.2.6',scipy='1.15.2'),'Reviewed environment differs')
        design,settings = loader.load()
        model = design.model()
        require(len(design.participants)==60 and len(model.rows)==173 and model.G==8 and
                model.Z==7 and model.size==28,'Shared synthetic fixture dimensions differ')
        p,v = loader.parameters(FIXTURE/'params/theta4_perturbed.json')
        require(p.shape==(28,) and np.isfinite(p).all() and v==.0004,'Fixed theta4 parameters differ')
        rows = tuple(r for r in model.rows if r.gene==args.gene)
        require(rows,'Requested synthetic gene has no rows')
        count_sha = hashlib.sha256(np.asarray([r.a for r in model.rows],dtype='<i8').tobytes()).hexdigest()
        row_body = json.dumps([vars(r) for r in rows],sort_keys=True,allow_nan=False)
        protocol.update(environment=versions,fixture_seed=settings['seed'],
            fixed_p=p.tolist(),fixed_v=v,fixed_gene_rows=len(rows),
            full_fixture_count_vector_sha256=count_sha,
            selected_row_metadata_sha256=hashlib.sha256(row_body.encode()).hexdigest(),
            exact_binary_float_inputs=True,counts_are_original_synthetic_fixture_not_CRT_draw=True)
        save('protocol',protocol)
        # Double roots are initial Newton guesses only. Double weights are discarded.
        def seeds(order):
            return np.sqrt(2)*np.polynomial.hermite.hermgauss(order)[0]
        result = helper.check(rows,p,v,model.G,model.Z,bb.bb_vector,seeds,save)
        for path,digest in pins.items():
            require(sha(path)==digest,'Pinned source/input changed during diagnostic')
        result.update(case='theta4_perturbed',gene=args.gene,environment=versions,
            overall_elapsed_seconds=time.monotonic()-began,
            all_original_accuracy_requirements_remain_required=True)
        save('summary',result)
        require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file())<=CAP,'Final output cap exceeded')
    except BaseException:
        with (output/'failure.json').open('x') as handle:
            json.dump(dict(completed=False,traceback=traceback.format_exc(),partials_preserved=True,
                no_model_fit=True,no_calibration=True),handle,indent=2)
            handle.write('\n')
        raise


if __name__=='__main__':
    main()
