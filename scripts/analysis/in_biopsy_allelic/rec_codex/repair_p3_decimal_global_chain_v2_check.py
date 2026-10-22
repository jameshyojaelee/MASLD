#!/usr/bin/env python3
"""Native likelihood and all-row global parameter checks for Decimal derivatives."""
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
           'repair_p3_bb_mode_anchor.py')


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
        selected_rows='all 24 target-gene rows, fixed original draw',
        parameter_indices=dict(r_gene=10,omega_S=24),
        native_absolute_tolerances=['1e-9','1e-8','1e-8','1e-5'],
        u_points=['0','+sqrt(v)','-sqrt(v)'],
        decimal_digits=[80,110], FD_steps=['1e-16','1e-20'],
        per_row_and_global_derivative_absolute_tolerance='1e-18', precision_stability_tolerance='1e-50',
        independent_closed_native_BB_derivatives=True, probability_only_FD=True,
        no_quadrature=True, no_integral_accuracy=True, no_fit=True, no_calibration=True)
    if args.inspect:
        print(json.dumps(protocol))
        return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required before numerical imports')
    out = REC/('repair_p3_decimal_global_chain_v2_'+args.case+'_'+os.environ['SLURM_JOB_ID'])
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
        protocol.update(environment=versions)
        save('protocol',protocol)
        anchor = load('mixed_chain_native_anchor', archive/'repair_p3_bb_mode_anchor.py')
        original.bb_terms = anchor.stable_bb_terms
        require(any(row.stage != 0 for row in rows), 'All target stages are zero')
        p = parent['fixed_p']
        require(all(p[j] != 0 for j in (16,24,25,26,27)),
                'Native stage/dispersion mapping has an unexercised zero coefficient')
        results = []
        for digits in (80,110):
            with localcontext() as context:
                context.prec = digits
                prepared = mixed.prepare(rows,p,8,7,probability.bb_vector)
                q = [Decimal.from_float(float(x)) for x in p]
                rootv = Decimal.from_float(parent['fixed_v']).sqrt()
                def direct_predictors(row, shifted, u):
                    # Recompute the native parameter map independently of prepare
                    # and omega_slope. Retain Decimal perturbations throughout.
                    stage = Decimal.from_float(float(row.stage))
                    factor = Decimal(1)+shifted[16]*stage+sum(
                        (shifted[17+j]*Decimal.from_float(float(row.z[j])) for j in range(7)),Decimal(0))
                    eta = shifted[row.gene]*factor-int(row.orientation)*(
                        shifted[24]*stage+Decimal.from_float(float(row.omega)))+u*stage
                    rr = shifted[8+row.gene]+sum((shifted[25+j]*Decimal.from_float(float(x))
                        for j,x in enumerate((row.stage,row.background,row.duplication))),Decimal(0))
                    return eta,rr
                for label,u in [('zero',Decimal(0)),('positive',rootv),('negative',-rootv)]:
                    native_checks = []
                    row_analytic = []
                    for d in prepared:
                        row = d['row']
                        eta,rr = direct_predictors(row,q,u)
                        en,rn,jac,_ = model.predictors(row,np.asarray(p),float(u),False)
                        require(jac[1,10] == 1 and jac[0,24] == -row.orientation*row.stage
                                and jac[0,10] == 0 and jac[1,24] == 0,
                                'Native r/omega Jacobian differs')
                        mapped = (en-row.orientation*row.omega,rn)
                        mapping_difference = [abs(eta-Decimal.from_float(float(mapped[0]))),
                                              abs(rr-Decimal.from_float(float(mapped[1])))]
                        require(eta == d['eta0']+u*d['stage'] and rr == d['r'],
                                'Independent Decimal/native parameter layout differs')
                        terms = mixed.row_terms(d,eta,rr)
                        row_analytic.append((terms[0],terms[2],d['omega_slope']*terms[1],
                                             d['omega_slope']*terms[3]))
                        native = original.conditional_terms(row,en,rn,True)
                        diffs = [terms[0]-Decimal.from_float(float(native[0])),
                                 terms[1]-Decimal.from_float(float(native[1][0])),
                                 terms[2]-Decimal.from_float(float(native[1][1])),
                                 terms[3]-Decimal.from_float(float(native[2][0,1]))]
                        limits = list(map(Decimal,('1e-9','1e-8','1e-8','1e-5')))
                        a,b = mixed.shapes(eta,rr)[:2]
                        bb = probability.bb_vector(d['support'],d['n'],a,b,d['logchoose'])
                        het = [d['hetero_prior']*lp.exp() for lp in bb]
                        oi = d['observed_index']
                        weights = (het[oi]/(het[oi]+d['hom'][oi]),
                                   sum(het,Decimal(0))/(sum(het,Decimal(0))+sum(d['hom'],Decimal(0))))
                        native_checks.append(dict(individual=row.individual,n=row.n,a=row.a,
                            stage=row.stage,orientation=int(row.orientation),
                            predictor_absolute_differences=list(map(str,mapping_difference)),
                            Decimal_minus_native=list(map(str,diffs)),
                            posterior_heterozygote_weights_observed_and_support=list(map(str,weights)),
                            passed=all(x <= Decimal('1e-12') for x in mapping_difference)
                                and all(abs(x) <= t for x,t in zip(diffs,limits))))
                    analytic = mixed.conditional_moments(prepared,u)
                    native_global = model.evaluate(np.asarray(p),rows=rows,u=float(u),need_hessian=True)
                    native_global_values = (native_global[0],native_global[1][10],
                                            native_global[1][24],native_global[2][10,24])
                    global_native_differences = [a-Decimal.from_float(float(b))
                                                for a,b in zip(analytic,native_global_values)]
                    global_native_passed = all(abs(x) <= Decimal(t) for x,t in
                        zip(global_native_differences,('1e-9','1e-8','1e-8','1e-5')))
                    probability_rows = {}
                    def evaluate(dr,dw):
                        shifted = q.copy()
                        shifted[10] += dr
                        shifted[24] += dw
                        values = tuple(mixed.row_logmass(d,*direct_predictors(d['row'],shifted,u),
                            probability.bb_vector) for d in prepared)
                        probability_rows[(Decimal(dr),Decimal(dw))] = values
                        return sum(values,Decimal(0))
                    baseline = evaluate(Decimal(0),Decimal(0))
                    ladder = []
                    for step in ('1e-16','1e-20'):
                        h = Decimal(step)
                        sr = (evaluate(h,0)-evaluate(-h,0))/(2*h)
                        sw = (evaluate(0,h)-evaluate(0,-h))/(2*h)
                        hrw = (evaluate(h,h)-evaluate(h,-h)-evaluate(-h,h)+evaluate(-h,-h))/(4*h*h)
                        differences = [x-y for x,y in zip(analytic,(baseline,sr,sw,hrw))]
                        row_differences = []
                        for index,reference in enumerate(row_analytic):
                            val = lambda dr,dw:probability_rows[(Decimal(dr),Decimal(dw))][index]
                            row_fd = (val(0,0), (val(h,0)-val(-h,0))/(2*h),
                                (val(0,h)-val(0,-h))/(2*h),
                                (val(h,h)-val(h,-h)-val(-h,h)+val(-h,-h))/(4*h*h))
                            row_differences.append(tuple(a-b for a,b in zip(reference,row_fd)))
                        rows_passed = all(abs(x) <= Decimal('1e-18') for row_diff in row_differences for x in row_diff)
                        ladder.append(dict(step=step,probability_FD=list(map(str,(baseline,sr,sw,hrw))),
                            per_row_analytic_minus_FD=[list(map(str,diff)) for diff in row_differences],
                            per_row_criterion_met=rows_passed,
                            analytic_minus_FD=list(map(str,differences)),
                            passed=rows_passed and all(abs(x) <= Decimal('1e-18') for x in differences)))
                    record = dict(digits=digits,u_label=label,analytic=list(map(str,analytic)),
                        native_checks=native_checks,FD_ladder=ladder,
                        Decimal_minus_native_global=list(map(str,global_native_differences)),
                        native_global_criterion_met=global_native_passed,
                        passed=global_native_passed and all(x['passed'] for x in native_checks+ladder))
                    results.append(record)
                    save(f'chain_{label}_{digits}',record)
        with localcontext() as context:
            context.prec = 110
            stability = []
            for low,high in zip(results[:3],results[3:]):
                require(low['u_label'] == high['u_label'],'Precision ordering differs')
                stability.append(str(max(abs(Decimal(a)-Decimal(b))
                    for a,b in zip(low['analytic'],high['analytic']))))
        for path,digest in pins.items():
            require(sha(path) == digest,'Reviewed source/input changed during check')
        passed = all(r['passed'] for r in results) and all(Decimal(x) <= Decimal('1e-50') for x in stability)
        save('summary',dict(completed=True,records=len(results),target_rows=24,
            native_and_global_chain_criteria_met=passed,analytic_precision_stabilities=stability,
            elapsed_seconds=time.perf_counter()-began,
            measured_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            no_integral_accuracy=True,no_fit=True,no_calibration=True))
        require(passed,'Native/global Decimal chain comparison failed')
    except BaseException:
        with (out/'failure.json').open('x') as handle:
            json.dump(dict(traceback=traceback.format_exc(),partials_preserved=True),handle,indent=2)
        raise


if __name__ == '__main__':
    main()

