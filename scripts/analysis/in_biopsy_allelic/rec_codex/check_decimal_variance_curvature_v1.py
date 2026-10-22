#!/usr/bin/env python3
"""Independent probability-only Decimal variance score/curvature diagnostic.

For a normal random effect, I'/I=E[(z^2-1)/(2v)] and
I''/I=E[(z^4-6z^2+3)/(4v^2)]. Decimal arithmetic limits cancellation.
No candidate conditional derivatives are used. Whole-line Gauss-Hermite
orders32/64/128 are fixed before execution;64/128 provide the comparison.
This checks only one gene's variance coordinate, not the full Hessian.
"""
import argparse
import ast
from dataclasses import replace
from decimal import Decimal, localcontext
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
POINT = REC/'repair_p3_gamma_order_polish_draw_0_amd0_22084420'
CANDIDATE = {2: REC/'pathwise_longdouble_gene2_22084806', 5: REC/'pathwise_longdouble_gene5_22084793'}
SOURCES = ('likelihood_v1.py', 'design_v1.py', 'bfix_v1.py', 'tests_v1.py',
           'repair_p3_decimal_mixed.py', 'repair_p3_boundary_decimal_fast.py')
DIGITS = 80
ORDERS = (32, 64, 128)
WORK = None


def require(ok, message):
    if not ok: raise RuntimeError(message)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def hermite(n, z):
    previous, current = Decimal(1), z
    for k in range(2, n+1):
        previous, current = current, z*current-(k-1)*previous
    return current, n*previous


def rule(n, np):
    seeds = np.polynomial.hermite_e.hermegauss(n)[0]
    nodes = []
    tolerance = Decimal('1e-65')
    for seed in seeds:
        z = Decimal.from_float(float(seed))
        for _ in range(40):
            value, derivative = hermite(n, z)
            change = value/derivative
            z -= change
            if abs(change) < tolerance: break
        else: raise RuntimeError('Decimal Hermite root did not converge')
        _, derivative = hermite(n, z)
        weight = Decimal(math.factorial(n))/(derivative*derivative)
        require(weight > 0 and z.is_finite(), 'Invalid Hermite root/weight')
        nodes.append((z, weight))
    require(len({z for z, _ in nodes}) == n, 'Duplicate Hermite root')
    errors = []
    expected = Decimal(1)
    for degree in range(2*n):
        if degree > 0 and degree % 2 == 0: expected *= degree-1
        found = sum((w*z**degree for z, w in nodes), Decimal(0))
        target = Decimal(0) if degree % 2 else expected
        error = abs(found-target)/max(Decimal(1), abs(expected))
        require(error < Decimal('1e-55'), 'Decimal Hermite moment check failed')
        errors.append(str(error))
    return nodes, errors


def node_moments(task):
    order, index, z, weight = task
    prepared, v, shift, mixed, probability = WORK
    with localcontext() as ctx:
        ctx.prec = DIGITS
        u = v.sqrt()*z
        ell = sum((mixed.row_logmass(d, d['eta0']+u*d['stage'], d['r'], probability.bb_vector)
                   for d in prepared), Decimal(0))
        mass = weight*(ell-shift).exp()
        g = (z*z-1)/(2*v)
        second = (z**4-6*z*z+3)/(4*v*v)
        return dict(order=order, index=index, z=str(z), weight=str(weight),
                    moments=[str(mass), str(mass*g), str(mass*second)])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--gene', type=int, choices=(2, 5), required=True)
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
        if p.suffix == '.py': ast.parse(p.read_text())
    required = [*(HERE/n for n in SOURCES), POINT/'recipe_summary.json', POINT/'draw_counts.json',
        *(FIXTURE/n for n in ('globals.json', 'cohorts.tsv', 'individuals.tsv', 'rows.tsv')),
        *(CANDIDATE[a.gene]/n for n in ('protocol.json', 'summary.json', 'bound26.json'))]
    require(all(str(p) in pins for p in required), 'Missing input pin')
    read = lambda p: json.loads(p.read_text())
    cp = read(CANDIDATE[a.gene]/'protocol.json')
    candidate = read(CANDIDATE[a.gene]/'bound26.json')
    failure = read(POINT/'recipe_summary.json')['fitting_integration_failure']
    require(cp['gene'] == a.gene and cp['parameters'] == failure['parameters']
        and cp['variance'] == failure['variance'], 'Candidate fixed point differs')
    if a.inspect:
        print(json.dumps(dict(gene=a.gene, digits=DIGITS, orders=ORDERS, workers=8, no_fit=True)))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'), 'Individual compute required')
    require(os.environ.get('SLURM_CPUS_PER_TASK') == '8', 'Eight CPUs required')
    require(all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')), 'Threads differ')
    out = REC/('decimal_variance_curvature_gene'+str(a.gene)+'_'+os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)
    archive = out/'executed_sources'; archive.mkdir()
    for p in (Path(__file__), *(HERE/n for n in SOURCES)):
        q = archive/p.name; q.write_bytes(p.read_bytes())
        require(sha(q) == (a.expected_runner_sha if p == Path(__file__) else pins[str(p)]), 'Archive differs')
    import numpy as np
    import scipy
    sys.path.insert(0, str(archive))
    loader = importlib.import_module('bfix_v1')
    original = importlib.import_module('likelihood_v1')
    mixed = importlib.import_module('repair_p3_decimal_mixed')
    probability = importlib.import_module('repair_p3_boundary_decimal_fast')
    for name in SOURCES:
        require(Path(sys.modules[name[:-3]].__file__).resolve().parent == archive.resolve(), 'Import escaped archive')
    require(loader.FIXTURE.resolve() == FIXTURE.resolve(), 'Fixture loader differs')
    design, settings = loader.load(); template = design.model()
    counts = read(POINT/'draw_counts.json')
    require(counts['full_row_keys'] == [[r.individual, settings['gene_order'][r.gene]] for r in template.rows]
        and len(counts['counts']) == len(template.rows) == 173, 'Fixed count order differs')
    model = original.Model([replace(r, a=c) for r, c in zip(template.rows, counts['counts'])], gene_count=8)
    rows = tuple(r for r in model.rows if r.gene == a.gene)
    require((model.G, model.Z, model.size) == (8, 7, 28) and len(rows) == cp['rows'], 'Synthetic scope differs')
    p = np.asarray(failure['parameters']); before = p.tobytes(); row_before = repr(rows)
    def save(name, value):
        with (out/(name+'.json')).open('x') as f:
            json.dump(value, f, indent=2, allow_nan=False); f.write('\n'); f.flush(); os.fsync(f.fileno())
    save('protocol', dict(gene=a.gene, parameters=p.tolist(), variance=failure['variance'], rows=len(rows),
        digits=DIGITS, orders=ORDERS, workers=8, source_and_input_sha256=pins,
        environment=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
        node=os.environ.get('SLURMD_NODENAME'), no_candidate_conditional_derivatives=True,
        acceptance=dict(LL=1e-9, variance_score=1e-8, variance_curvature=1e-5),
        no_fit=True, no_calibration=True, no_human_data=True, empirical_not_certified=True))
    start = time.perf_counter()
    global WORK
    with localcontext() as ctx:
        ctx.prec = DIGITS
        v = Decimal.from_float(failure['variance'])
        shift = Decimal.from_float(candidate['logmass'])
        prepared = mixed.prepare(rows, p, model.G, model.Z, probability.bb_vector)
        WORK = prepared, v, shift, mixed, probability
        tasks = []; checks = {}
        for n in ORDERS:
            nodes, errors = rule(n, np); checks[str(n)] = errors
            tasks.extend((n, i, z, w) for i, (z, w) in enumerate(nodes))
        save('rules', dict(moment_relative_errors=checks, task_count=len(tasks)))
        records = []
        with multiprocessing.get_context('fork').Pool(8) as pool:
            for result in pool.imap(node_moments, tasks, chunksize=1):
                records.append(result)
                if len(records) % 16 == 0:
                    print(json.dumps(dict(completed_nodes=len(records), required_nodes=len(tasks))), flush=True)
        save('node_moments', dict(records=records))
        require([(r['order'], r['index']) for r in records] == [(t[0], t[1]) for t in tasks], 'Node coverage/order differs')
        results = []
        for n in ORDERS:
            raw = [sum((Decimal(r['moments'][j]) for r in records if r['order'] == n), Decimal(0)) for j in range(3)]
            mass, score, second = raw
            require(mass > 0, 'Invalid integral mass')
            gradient = score/mass
            curvature = second/mass-gradient*gradient
            results.append(dict(order=n, raw_moments=list(map(str, raw)), logmass=str(shift+mass.ln()),
                                variance_score=str(gradient), variance_curvature=str(curvature)))
        names = ('logmass', 'variance_score', 'variance_curvature')
        difference = {k: float(abs(Decimal(results[-1][k])-Decimal(results[-2][k]))) for k in names}
        candidate_difference = dict(logmass=float(Decimal(results[-1]['logmass'])-Decimal.from_float(candidate['logmass'])),
            variance_score=float(Decimal(results[-1]['variance_score'])-Decimal.from_float(candidate['gradient'][-1])),
            variance_curvature=float(Decimal(results[-1]['variance_curvature'])-Decimal.from_float(candidate['hessian'][-1][-1])))
    passed = all(difference[k] <= t and abs(candidate_difference[k]) <= t for k, t in zip(names, (1e-9, 1e-8, 1e-5)))
    require(p.tobytes() == before and repr(rows) == row_before, 'Inputs changed')
    for path, digest in pins.items(): require(sha(Path(path)) == digest, 'Input changed during check')
    save('summary', dict(completed=True, gene=a.gene, results=results, order64_128_differences=difference,
        candidate_differences=candidate_difference, passed=passed, seconds=time.perf_counter()-start,
        no_independent_full_Hessian_validation=True, no_fit=True, no_calibration=True, empirical_not_certified=True))
    if not passed: raise SystemExit(1)


if __name__ == '__main__':
    main()
