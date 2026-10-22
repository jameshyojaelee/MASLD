#!/usr/bin/env python3
"""Selected synthetic Decimal mixed derivatives versus probability-only FD."""
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
           'repair_p3_boundary_decimal_fast.py','repair_p3_decimal_mixed.py')


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
    protocol = dict(case=args.case, draw=draw, fixed_p=parent['fixed_p'], fixed_v=parent['fixed_v'],
        count_vector_sha256=parent['count_vector_sha256'], source_and_input_sha256={str(p):h for p,h in pins.items()},
        selected_rows='deepest target-gene row plus deepest opposite-orientation row in full original draw',
        u_points=['0','+sqrt(v)','-sqrt(v)'],
        decimal_digits=[80,110], FD_steps=['1e-16','1e-20'],
        selected_row_derivative_absolute_tolerance='1e-18', precision_stability_tolerance='1e-50',
        independent_closed_native_BB_derivatives=True, probability_only_FD=True,
        no_quadrature=True, no_integral_accuracy=True, no_fit=True, no_calibration=True)
    if args.inspect:
        print(json.dumps(protocol))
        return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required before numerical imports')
    out = REC/('repair_p3_decimal_mixed_check_'+args.case+'_'+os.environ['SLURM_JOB_ID'])
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
        target = max(rows,key=lambda r:r.n)
        opposite = tuple(r for r in drawn if r.orientation == -target.orientation)
        require(opposite, 'Full original synthetic draw lacks opposite orientation')
        selected = (target,max(opposite,key=lambda r:r.n))
        protocol.update(environment=versions, selected_rows=[dict(gene=r.gene,individual=r.individual,n=r.n,
            observed=r.a,orientation=r.orientation) for r in selected])
        save('protocol',protocol)
        results = []
        for digits in (80,110):
            with localcontext() as context:
                context.prec = digits
                prepared = mixed.prepare(selected,parent['fixed_p'],8,7,probability.bb_vector)
                rootv = Decimal.from_float(parent['fixed_v']).sqrt()
                for index,d in enumerate(prepared):
                    for label,u in [('zero',Decimal(0)),('positive',rootv),('negative',-rootv)]:
                        eta,r = d['eta0']+u*d['stage'],d['r']
                        analytic = mixed.row_terms(d,eta,r)
                        baseline = mixed.row_logmass(d,eta,r,probability.bb_vector)
                        ladders = []
                        for step in ('1e-16','1e-20'):
                            h = Decimal(step)
                            evaluate = lambda de,dr:mixed.row_logmass(d,eta+de,r+dr,probability.bb_vector)
                            se = (evaluate(h,0)-evaluate(-h,0))/(2*h)
                            sr = (evaluate(0,h)-evaluate(0,-h))/(2*h)
                            her = (evaluate(h,h)-evaluate(h,-h)-evaluate(-h,h)+evaluate(-h,-h))/(4*h*h)
                            raw = [analytic[0]-baseline,analytic[1]-se,analytic[2]-sr,analytic[3]-her]
                            ladders.append(dict(step=step, probability_FD=list(map(str,(baseline,se,sr,her))),
                                analytic_minus_probability_FD=list(map(str,raw)),
                                within_selected_criterion=all(abs(x) <= Decimal('1e-18') for x in raw)))
                        record = dict(digits=digits,selected_row=index,u_label=label,eta=str(eta),r=str(r),
                            analytic=list(map(str,analytic)),FD_ladder=ladders,
                            selected_criterion_met=all(a['within_selected_criterion'] for a in ladders))
                        results.append(record)
                        save(f'row{index}_{label}_{digits}',record)
        with localcontext() as context:
            context.prec = 110
            stability = []
            for low,high in zip(results[:6],results[6:]):
                require((low['selected_row'],low['u_label']) == (high['selected_row'],high['u_label']),
                        'Precision family ordering differs')
                difference = max(abs(Decimal(a)-Decimal(b)) for a,b in zip(low['analytic'],high['analytic']))
                stability.append(str(difference))
        for p,h in pins.items():
            require(sha(p) == h, 'Reviewed source/input changed during check')
        passed = all(r['selected_criterion_met'] for r in results) and all(Decimal(x) <= Decimal('1e-50') for x in stability)
        save('summary',dict(completed=True,records=len(results),selected_row_FD_criterion_met=passed,
            analytic_precision_stabilities=stability,elapsed_seconds=time.perf_counter()-began,
            measured_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            no_integral_accuracy=True,no_fit=True,no_calibration=True))
        require(passed,'Selected Decimal mixed-derivative admission failed')
    except BaseException:
        with (out/'failure.json').open('x') as handle:
            json.dump(dict(traceback=traceback.format_exc(),partials_preserved=True),handle,indent=2)
        raise


if __name__ == '__main__':
    main()
