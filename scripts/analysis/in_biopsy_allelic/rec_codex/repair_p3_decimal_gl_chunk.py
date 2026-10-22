#!/usr/bin/env python3
"""Independent Decimal GL raw-moment chunk on an exact synthetic P3 point."""
import argparse
import ast
from decimal import Decimal, localcontext
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

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
HERE = ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
FIXTURE = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
CAMPAIGN = REC/'p3_parallel_20261003T130323Z'
CASES = {'failed0gene2':(0,22039200), 'positive18':(18,22039201)}
MODULES = ('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
           'repair_p3_boundary_decimal_fast.py','repair_p3_decimal_mixed.py',
           'repair_p3_decimal_gl_moments.py')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=CASES, required=True)
    parser.add_argument('--rule',type=int,choices=range(3),required=True)
    parser.add_argument('--chunk',type=int,choices=range(16),required=True)
    parser.add_argument('--pins', type=Path, required=True)
    parser.add_argument('--pins-sha', required=True)
    parser.add_argument('--self-sha', required=True)
    parser.add_argument('--inspect', action='store_true')
    args = parser.parse_args()
    require(sha(Path(__file__)) == args.self_sha and sha(args.pins) == args.pins_sha,
            'Reviewed checker/manifest differs')
    pins = {Path(p):h for p,h in json.loads(args.pins.read_text()).items()}
    for p,h in pins.items():
        require(p.is_absolute() and p.resolve().is_relative_to(ROOT) and sha(p) == h,
                'Reviewed source/input changed: '+str(p))
        if p.suffix == '.py':
            ast.parse(p.read_text())
    require(all(HERE/n in pins for n in MODULES), 'Executable dependency pin missing')
    draw, job = CASES[args.case]
    parent_dir = REC/f'repair_p3_fixed_point_{args.case}_{job}'
    require(parent_dir/'protocol.json' in pins and FIXTURE/'draws/p3_uniforms.npy' in pins
            and CAMPAIGN/'baseline/observed_null.json' in pins, 'Original fixture identity pin missing')
    parent = json.loads((parent_dir/'protocol.json').read_text())
    require(parent['fixed_p'][10] == -15. and parent['fixed_gene_rows'] == 24,
            'Target dispersion-bound point differs')
    for admission_case,admission_job in [('failed0gene2',22075551),('positive18',22075552)]:
        admission_dir = REC/f'repair_p3_decimal_global_chain_v2_{admission_case}_{admission_job}'
        require(admission_dir/'summary.json' in pins and admission_dir/'protocol.json' in pins,
                'Both complete strengthened admission results must be pinned')
        receipt = json.loads((admission_dir/'summary.json').read_text())
        require(receipt['completed'] is True and receipt['native_and_global_chain_criteria_met'] is True
                and receipt['target_rows'] == 24 and receipt['records'] == 6,
                'Strengthened native/global/per-row derivative admission incomplete')
        if admission_case == args.case:
            admitted = json.loads((admission_dir/'protocol.json').read_text())
            require(admitted['count_vector_sha256'] == parent['count_vector_sha256']
                    and admitted['fixed_p'] == parent['fixed_p'] and admitted['fixed_v'] == parent['fixed_v'],
                    'Admitted synthetic point differs')
    protocol = dict(case=args.case, draw=draw, fixed_p=parent['fixed_p'], fixed_v=parent['fixed_v'],
        count_vector_sha256=parent['count_vector_sha256'], source_and_input_sha256={str(p):h for p,h in pins.items()},
        rule_index=args.rule,chunk=args.chunk,chunks_per_rule=16,
        rule_family=[dict(bound=24,spacing='0.25',order=16),
                     dict(bound=26,spacing='0.125',order=16),
                     dict(bound=26,spacing='0.125',order=32)],
        decimal_digits=[80,110],fixed_log_rescaling='0',
        target_mixed_entries=[[10,24],[10,28]],
        raw_moments_only=True,no_chunk_covariance=True,
        no_certified_quadrature_error=True,no_integral_accuracy=True,no_fit=True,no_calibration=True)
    if args.inspect:
        print(json.dumps(protocol))
        return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required before numerical imports')
    out = REC/('repair_p3_decimal_GL_'+args.case+f'_r{args.rule}_c{args.chunk}_'+os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)
    archive = out/'executed_sources'
    archive.mkdir()
    for p in (Path(__file__), *(HERE/n for n in MODULES)):
        target = archive/p.name
        target.write_bytes(p.read_bytes())
        require(sha(target) == (args.self_sha if p == Path(__file__) else pins[p]), 'Executable archive differs')
    began = time.perf_counter()
    try:
        import numpy as np
        import scipy
        def save(name, value):
            body = json.dumps(value, indent=2, allow_nan=False)+'\n'
            require(sum(p.stat().st_size for p in out.rglob('*') if p.is_file())+len(body.encode())
                    < 32*1024**2-65536, 'Output cap exceeded')
            with (out/(name+'.json')).open('x') as handle:
                handle.write(body)
                handle.flush()
                os.fsync(handle.fileno())
            print(json.dumps(dict(saved=name, elapsed_seconds=time.perf_counter()-began)), flush=True)
        sys.path.insert(0,str(archive))
        require(not any(n[:-3] in sys.modules for n in MODULES[:4]), 'Original modules already imported')
        loader = importlib.import_module('bfix_v1')
        original = importlib.import_module('likelihood_v1')
        importlib.import_module('tests_v1')
        for name in MODULES[:4]:
            require(Path(sys.modules[name[:-3]].__file__).resolve().parent == archive.resolve()
                    and sha(Path(sys.modules[name[:-3]].__file__)) == pins[HERE/name], 'Original import escaped archive')
        mixed = load('mixed_check_Decimal', archive/'repair_p3_decimal_mixed.py')
        probability = load('mixed_check_probability', archive/'repair_p3_boundary_decimal_fast.py')
        versions = dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        require(versions == parent['environment'], 'Original numerical environment differs')
        design, settings = loader.load()
        model = design.model()
        null = json.loads((CAMPAIGN/'baseline/observed_null.json').read_text())
        require(null['converged'] is True and settings['B'] == 50 and len(model.rows) == 173,
                'Original generating null or synthetic design differs')
        null_p = np.asarray(null['parameters'])
        uniforms = np.load(FIXTURE/'draws/p3_uniforms.npy', mmap_mode='r', allow_pickle=False)
        require(uniforms.shape == (50,173,2) and
            [(r['individual_id'],r['gene_id']) for r in loader.table('rows.tsv')] ==
            [(r.individual,settings['gene_order'][r.gene]) for r in model.rows], 'Uniform-to-row order differs')
        drawn = []
        for row,(uc,uk) in zip(model.rows,uniforms[draw]):
            eta,r = model.predictors(row,null_p)[:2]
            drawn.append(original.conditional_draw(row,eta,r,uc,uk))
        digest = hashlib.sha256(np.asarray([r.a for r in drawn],dtype='<i8').tobytes()).hexdigest()
        require(digest == parent['count_vector_sha256'], 'Exact original P3 draw counts differ')
        rows = tuple(r for r in original.Model(drawn,gene_count=8).rows if r.gene == 2)
        require(len(rows) == 24 and model.G == 8 and model.Z == 7 and model.size == 28,
                'Original target dimensions differ')
        quadrature = load('independent_Decimal_GL', archive/'repair_p3_decimal_gl_moments.py')
        protocol.update(environment=versions)
        save('protocol',protocol)
        results = []
        for digits in (80,110):
            with localcontext() as context:
                context.prec = digits
                prepared = mixed.prepare(rows,parent['fixed_p'],8,7,probability.bb_vector)
                variance = Decimal.from_float(parent['fixed_v'])
                def progress(panel,calls):
                    print(json.dumps(dict(digits=digits,panel_completed=panel,node_calls=calls,
                        elapsed_seconds=time.perf_counter()-began)),flush=True)
                result = quadrature.integrate(prepared,variance,Decimal(0),args.rule,args.chunk,digits,
                                               mixed.conditional_moments,progress)
                result['analytic_absolute_tail_bounds_not_interval_evaluated'] = quadrature.whole_line_tail_bounds(
                    rows,variance,Decimal(0),result['bound'])
                results.append(result)
                save(f'moments_{digits}',result)
        for path,digest in pins.items():
            require(sha(path) == digest, 'Reviewed source/input changed during integral')
        save('summary',dict(completed=True,case=args.case,rule_index=args.rule,chunk=args.chunk,
            records=len(results),raw_moments_only=True,no_chunk_covariance=True,
            elapsed_seconds=time.perf_counter()-began,
            measured_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            no_certified_quadrature_error=True,no_integral_accuracy=True,no_fit=True,no_calibration=True))
    except BaseException:
        with (out/'failure.json').open('x') as handle:
            json.dump(dict(traceback=traceback.format_exc(),partials_preserved=True),handle,indent=2)
        raise


if __name__ == '__main__':
    main()

