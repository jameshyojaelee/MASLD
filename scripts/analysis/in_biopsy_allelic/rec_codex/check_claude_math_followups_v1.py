#!/usr/bin/env python3
"""Bounded synthetic checks requested by the independent source reviewer.

Boundary11 uses a nonintegrated heat-identity score with independent native
long-double conditional derivatives. Uniform_refusal uses composite GL over
every grid interval and fixed-u Gaussian-density derivatives, without peaks.
No fitting, calibration, human data or reclassification of earlier results.
"""
import argparse
import ast
from dataclasses import replace
import hashlib
import importlib
import json
import multiprocessing
import os
from pathlib import Path
import sys
import time
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE = ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
BOUNDARY = REC/'adaptive_interval_pilot_draw11_22085228'
REFUSAL = REC/'standard_hermite_refusal_check_22085142'
REF = REC/'refusal_adaptive_intervals_22085179'
SOURCES = ('likelihood_v1.py', 'design_v1.py', 'bfix_v1.py', 'tests_v1.py',
    'repair_p3_local_moments.py', 'repair_p3_longdouble_conditional_v1.py',
    'repair_p3_local_gl_adaptive_order.py', 'repair_p3_integral.py')
WORK = None


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def worker(gene):
    import numpy as np
    model, p, variance, modules, kernel_sha = WORK
    rows = tuple(r for r in model.rows if r.gene == gene)
    context = modules['repair_p3_longdouble_conditional_v1'].LongDoubleConditional(model, p, rows)
    def batch(row, etas, dispersion, need_hessian=True):
        terms = [context.terms(row, float(eta), dispersion, need_hessian=need_hessian) for eta in etas]
        return (np.asarray([t[0] for t in terms]), np.asarray([t[1] for t in terms]),
                np.asarray([t[2] for t in terms]) if need_hessian else None)
    gl = modules['repair_p3_local_gl_adaptive_order']
    kernel = modules['repair_p3_integral']
    try:
        ll, g, h, report = gl.integrate_gene(model, p, variance, rows, batch, with_hessian=True,
            expected_kernel_sha=kernel_sha, kernel_path=Path(kernel.__file__),
            settings=kernel.Settings(initial_z_bound=16.), batch_size=32)
    except ValueError as exc:
        return dict(gene=gene, returned=False, message=str(exc), traceback=traceback.format_exc(),
            reports=getattr(exc, 'integration_reports', None))
    return dict(gene=gene, returned=True, logmass=ll, gradient=g, hessian=h, report=report)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case', choices=('boundary11', 'uniform_refusal'), required=True)
    ap.add_argument('--pins', type=Path, required=True)
    ap.add_argument('--expected-pins-sha', required=True)
    ap.add_argument('--expected-runner-sha', required=True)
    ap.add_argument('--inspect', action='store_true')
    a = ap.parse_args()
    require(sha(Path(__file__)) == a.expected_runner_sha and sha(a.pins) == a.expected_pins_sha, 'Runner/pins differ')
    pins = json.loads(a.pins.read_text())
    required = [*(HERE/n for n in SOURCES), *(FIXTURE/n for n in ('globals.json', 'cohorts.tsv', 'individuals.tsv', 'rows.tsv')),
        BOUNDARY/'p3_011.json', BOUNDARY/'draw_counts.json', REFUSAL/'refusal.json',
        REF/'summary.json', REC/'standard_hermite_pilot_draw0_22084966/draw_counts.json']
    require(all(str(p) in pins for p in required), 'Missing pin')
    for path, digest in pins.items():
        p = Path(path)
        require(p.resolve().is_relative_to(ROOT.resolve()) and sha(p) == digest, 'Input differs')
        if p.suffix == '.py':
            ast.parse(p.read_text())
    if a.inspect:
        print(json.dumps(dict(case=a.case, no_fit=True, no_peak_based_integration=True)))
        return
    cpus = 1 if a.case == 'boundary11' else 8
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID')
        and os.environ.get('SLURM_CPUS_PER_TASK') == str(cpus), 'Individual compute allocation differs')
    require(all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
        'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')), 'Threads differ')
    out = REC/('claude_math_'+a.case+'_'+os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)
    archive = out/'executed_sources'
    archive.mkdir()
    for p in (Path(__file__), *(HERE/n for n in SOURCES)):
        q = archive/p.name
        q.write_bytes(p.read_bytes())
        require(sha(q) == (a.expected_runner_sha if p == Path(__file__) else pins[str(p)]), 'Archive differs')
    import numpy as np
    import scipy
    sys.path.insert(0, str(archive))
    modules = {n[:-3]: importlib.import_module(n[:-3]) for n in SOURCES}
    require(all(Path(m.__file__).resolve().parent == archive.resolve() for m in modules.values()), 'Import escaped archive')
    loader = modules['bfix_v1']
    require(loader.FIXTURE.resolve() == FIXTURE.resolve(), 'Fixture differs')
    design, settings = loader.load()
    template = design.model()
    read = lambda path: json.loads(path.read_text())
    if a.case == 'boundary11':
        fit = read(BOUNDARY/'p3_011.json')['fit']
        require(fit['converged'] and fit['null']['converged'] and fit['v'] == 0, 'Boundary reference invalid')
        parameters, variance = fit['parameters'], 0.
        counts = read(BOUNDARY/'draw_counts.json')
    else:
        failure = read(REFUSAL/'refusal.json')['failure']
        parameters, variance = failure['parameters'], failure['variance']
        counts = read(REC/'standard_hermite_pilot_draw0_22084966/draw_counts.json')
    require(counts['full_row_keys'] == [[r.individual, settings['gene_order'][r.gene]] for r in template.rows]
        and len(counts['counts']) == len(template.rows) == 173
        and hashlib.sha256(np.asarray(counts['counts'], dtype='<i8').tobytes()).hexdigest() == counts['count_vector_sha256'], 'Counts differ')
    model = modules['likelihood_v1'].Model([replace(r, a=c) for r, c in zip(template.rows, counts['counts'])], gene_count=8)
    require((model.G, model.Z, model.size) == (8, 7, 28), 'Synthetic scope differs')
    require(all(type(r.percent) is int and type(r.minor) is int and type(r.depth) is int for r in model.rows), 'Integer support types differ')
    p = np.asarray(parameters)
    p.setflags(write=False)
    before, row_before = p.tobytes(), repr(model.rows)
    def safe(x):
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, np.generic):
            return x.item()
        if isinstance(x, dict):
            return {k: safe(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [safe(v) for v in x]
        return x
    def save(name, value):
        with (out/(name+'.json')).open('x') as f:
            json.dump(safe(value), f, indent=2, allow_nan=False)
            f.write('\n')
            f.flush()
            os.fsync(f.fileno())
    save('protocol', dict(case=a.case, parameters=p, variance=variance, source_and_input_sha256=pins,
        environment=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
        no_fit=True, no_calibration=True, no_human_data=True, no_independent_full_conditional_Hessian=True))
    started = time.perf_counter()
    if a.case == 'boundary11':
        original_ll, original_g = model.marginal(p, 0.)
        records = []
        for gene in range(8):
            rows = tuple(r for r in model.rows if r.gene == gene)
            context = modules['repair_p3_longdouble_conditional_v1'].LongDoubleConditional(model, p, rows)
            gu, huu, ll = 0., 0., 0.
            for row in rows:
                eta, dispersion = model.predictors(row, p, 0., need_hessian=False)[:2]
                value, g, h = context.terms(row, eta, dispersion)
                ll += value
                gu += row.stage*g[0]
                huu += row.stage**2*h[0, 0]
            records.append(dict(gene=gene, logmass=ll, score_u=gu, curvature_u=huu,
                variance_score=.5*(huu+gu*gu)))
        heat = sum(r['variance_score'] for r in records)
        ll = sum(r['logmass'] for r in records)
        differences = [abs(ll-original_ll), abs(heat-original_g[-1])]
        passed = differences[0] <= 1e-6 and differences[1] <= 1e-8+1e-6*abs(original_g[-1]) and heat <= 1e-6
        save('boundary', dict(records=records, original_variance_score=original_g[-1],
            heat_variance_score=heat, differences=differences, raw_heat_score_nonpositive=heat <= 0.,
            boundary_score_passes_original_tolerance=heat <= 1e-6,
            observed_orientations=sorted({r.orientation for r in model.rows}),
            nonzero_fixed_offsets=sum(r.omega != 0 for r in model.rows), integer_support_types_verified=True))
    else:
        global WORK
        WORK = model, p, variance, modules, pins[str(HERE/'repair_p3_integral.py')]
        with multiprocessing.get_context('fork').Pool(8) as pool:
            records = pool.map(worker, range(8))
        save('gene_integrals', dict(records=records))
        returned = all(r['returned'] for r in records)
        differences = None
        if returned:
            ll = sum(r['logmass'] for r in records)
            g = np.sum([r['gradient'] for r in records], axis=0)
            h = np.sum([r['hessian'] for r in records], axis=0)
            reference = read(REF/'summary.json')
            differences = [abs(ll-reference['logmass']), float(np.max(np.abs(g-reference['gradient']))),
                           float(np.max(np.abs(h-reference['hessian'])))]
            save('combined', dict(logmass=ll, gradient=g, hessian=h, differences=differences))
        passed = returned and all(x <= y for x, y in zip(differences, (1e-9, 1e-8, 1e-5)))
    require(p.tobytes() == before and repr(model.rows) == row_before, 'Fixed inputs changed')
    for path, digest in pins.items():
        require(sha(Path(path)) == digest, 'Input changed during check')
    save('summary', dict(completed=True, case=a.case, passed=bool(passed), differences=differences,
        seconds=time.perf_counter()-started, no_fit=True, no_calibration=True,
        no_independent_full_conditional_Hessian=True))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
