#!/usr/bin/env python3
"""Separate adaptive-interval reference at the saved synthetic refusal point.

Detected peaks supply interval breakpoints only. Integrate the original normal
measure on two domains with explicit derivative-tail and quadrature allowances.
No fitting, model selection, human data or reclassification of old draws.
"""
import argparse
import ast
from dataclasses import replace
import hashlib
import importlib
import json
import math
import multiprocessing
import os
from pathlib import Path
import sys
import time

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE = ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
POINT = REC/'standard_hermite_refusal_check_22085142/refusal.json'
COUNTS = REC/'standard_hermite_pilot_draw0_22084966/draw_counts.json'
PEAKS = REC/'centered_hermite_pathwise_refusal8_22085156/gene_integrals.json'
SOURCES = ('likelihood_v1.py', 'design_v1.py', 'bfix_v1.py', 'tests_v1.py',
           'repair_p3_local_moments.py', 'repair_p3_pathwise_moments_v1.py',
           'repair_p3_longdouble_conditional_v1.py')
WORK = None


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def worker(gene):
    import numpy as np
    from scipy.integrate import quad_vec
    model, p, variance, peaks, modules = WORK
    local = modules['repair_p3_local_moments'].LocalGeneMoments(
        model, p, tuple(r for r in model.rows if r.gene == gene), None)
    conditional = modules['repair_p3_longdouble_conditional_v1'].LongDoubleConditional(model, p, local.rows)
    moments = modules['repair_p3_pathwise_moments_v1'].PathwiseGeneMoments(local)
    shift = max(x['logheight'] for x in peaks[gene])-.5*math.log(2*math.pi)
    # The failed importance rule contributes only locations, not any integral.
    points = sorted({x['center']+multiple*x['scale'] for x in peaks[gene]
                     for multiple in (-8., -4., -2., -1., 0., 1., 2., 4., 8.)})
    def raw(z):
        row_values = []
        def callback(row, eta, r, need_hessian=True):
            value = conditional.terms(row, eta, r, need_hessian=need_hessian)
            row_values.append(value[0])
            return value
        _, g, _, second = moments.point(p, float(z), variance, callback)
        logmass = math.fsum(row_values)-z*z/2-.5*math.log(2*math.pi)
        return logmass, np.r_[1., g, second.ravel()]
    width = len(local.active)+1
    scales = np.ones(1+width+width*width)
    for z in (0., -1., 1., *points):
        scales = np.maximum(scales, np.abs(raw(z)[1]))
    def integrand(z):
        logmass, vector = raw(z)
        return math.exp(logmass-shift)*vector/scales
    records = []
    for bound in (24., 26.):
        values, error, info = quad_vec(integrand, -bound, bound,
            points=[z for z in points if -bound < z < bound], epsabs=1e-12,
            epsrel=1e-12, norm='max', cache_size=4*1024**2,
            limit=2048, workers=1, full_output=True)
        record = dict(bound=bound, success=bool(info.success), status=int(info.status),
            message=info.message, evaluations=info.neval, error_estimate=float(error),
            scaled_values=values, scales=scales, log_shift=shift)
        records.append(record)
        if not info.success or not np.isfinite(values).all() or not math.isfinite(error) or error < 0:
            return dict(gene=gene, returned=False, records=records, refusal='Adaptive quadrature refused')
        values = values*scales
        mass = float(values[0])
        mt, gt, st = moments.shifted_tail_bounds(bound, variance, shift)
        lower = mass-scales[0]*error-mt
        if lower <= 0:
            return dict(gene=gene, returned=False, records=records, refusal='Mass allowance unresolved')
        gradient = values[1:1+width]/mass
        second = values[1+width:].reshape(width, width)/mass
        hessian = second-np.outer(gradient, gradient)
        mass_error = scales[0]*error+mt
        ga = (scales[1:1+width]*error+gt+np.abs(gradient)*mass_error)/lower
        sa = (scales[1+width:].reshape(width, width)*error+st+np.abs(second)*mass_error)/lower
        ha = sa+np.outer(np.abs(gradient), ga)+np.outer(ga, np.abs(gradient))+np.outer(ga, ga)
        record.update(logmass=shift+math.log(mass), gradient=local.embed_vector(gradient),
            hessian=local.embed_matrix(hessian), shifted_moments=values,
            estimated_logmass_allowance=mass_error/lower,
            estimated_gradient_allowance=local.embed_vector(ga),
            estimated_hessian_allowance=local.embed_matrix(ha),
            shifted_tail_bounds=[mt, gt, st])
    a, b = records
    deltas = [abs(a['logmass']-b['logmass']), np.abs(a['gradient']-b['gradient']),
              np.abs(a['hessian']-b['hessian'])]
    allowed = [deltas[0]+a['estimated_logmass_allowance']+b['estimated_logmass_allowance'],
        float(np.max(deltas[1]+a['estimated_gradient_allowance']+b['estimated_gradient_allowance'])),
        float(np.max(deltas[2]+a['estimated_hessian_allowance']+b['estimated_hessian_allowance']))]
    return dict(gene=gene, returned=all(x <= y for x, y in zip(allowed, (1e-9, 1e-8, 1e-5))),
        records=records, differences_plus_allowances=allowed, logmass=b['logmass'],
        gradient=b['gradient'], hessian=b['hessian'], empirical_not_certified=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pins', type=Path, required=True)
    ap.add_argument('--expected-pins-sha', required=True)
    ap.add_argument('--expected-runner-sha', required=True)
    ap.add_argument('--inspect', action='store_true')
    a = ap.parse_args()
    require(sha(Path(__file__)) == a.expected_runner_sha and sha(a.pins) == a.expected_pins_sha, 'Runner/pins differ')
    pins = json.loads(a.pins.read_text())
    for path, digest in pins.items():
        p = Path(path)
        require(p.resolve().is_relative_to(ROOT.resolve()) and sha(p) == digest, 'Input differs')
        if p.suffix == '.py':
            ast.parse(p.read_text())
    require(all(str(p) in pins for p in (POINT, COUNTS, PEAKS, *(HERE/n for n in SOURCES),
        *(FIXTURE/n for n in ('globals.json', 'cohorts.tsv', 'individuals.tsv', 'rows.tsv')))), 'Missing pin')
    if a.inspect:
        print(json.dumps(dict(fixed_point_only=True, domains=[24., 26.], eight_genes=True)))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID')
        and os.environ.get('SLURM_CPUS_PER_TASK') == '8', 'Individual eight-CPU compute required')
    require(all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
        'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')), 'Threads differ')
    out = REC/('refusal_adaptive_intervals_'+os.environ['SLURM_JOB_ID'])
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
    counts = json.loads(COUNTS.read_text())
    require(counts['full_row_keys'] == [[r.individual, settings['gene_order'][r.gene]] for r in template.rows]
        and len(counts['counts']) == len(template.rows) == 173
        and hashlib.sha256(np.asarray(counts['counts'], dtype='<i8').tobytes()).hexdigest()
        == counts['count_vector_sha256'], 'Counts differ')
    model = modules['likelihood_v1'].Model([replace(r, a=c) for r, c in zip(template.rows, counts['counts'])], gene_count=8)
    require((model.G, model.Z, model.size) == (8, 7, 28), 'Synthetic scope differs')
    failure = json.loads(POINT.read_text())['failure']
    p = np.asarray(failure['parameters'])
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
    save('protocol', dict(parameters=p, variance=failure['variance'], source_and_input_sha256=pins,
        environment=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
        no_fit=True, no_calibration=True, no_human_data=True,
        different_quadrature_not_independent_conditional_derivatives=True))
    peaks = {g['gene']: g['peaks'] for g in json.loads(PEAKS.read_text())['records']}
    require(set(peaks) == set(range(8)), 'Peak roster differs')
    global WORK
    WORK = model, p, failure['variance'], peaks, modules
    started = time.perf_counter()
    with multiprocessing.get_context('fork').Pool(8) as pool:
        answers = pool.map(worker, range(8))
    save('gene_integrals', dict(records=answers))
    passed = all(g['returned'] for g in answers)
    require(p.tobytes() == before and repr(model.rows) == row_before, 'Inputs changed')
    for path, digest in pins.items():
        require(sha(Path(path)) == digest, 'Input changed during check')
    save('summary', dict(completed=True, passed=passed, seconds=time.perf_counter()-started,
        returned_genes=sum(g['returned'] for g in answers), no_fit=True, no_calibration=True,
        no_independent_full_Hessian_validation=True,
        logmass=sum(g['logmass'] for g in answers) if passed else None,
        gradient=np.sum([g['gradient'] for g in answers], axis=0) if passed else None,
        hessian=np.sum([g['hessian'] for g in answers], axis=0) if passed else None))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
