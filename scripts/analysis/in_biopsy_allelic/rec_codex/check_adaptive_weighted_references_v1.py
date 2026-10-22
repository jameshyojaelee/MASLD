#!/usr/bin/env python3
"""Four fixed points with weighted adaptive conditioning; no fitting."""
import argparse
import ast
from dataclasses import replace
from decimal import Decimal, localcontext
import hashlib
import importlib
import json
import multiprocessing
import os
from pathlib import Path
import sys
import time

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE = ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
FULL = REC/'repair_p3_gamma_batched_full_check_22078827'
POINT = REC/'repair_p3_gamma_order_polish_draw_0_amd0_22084420'
DECIMAL = {2: REC/'decimal_variance_curvature_gene2_22084836',
           5: REC/'decimal_variance_curvature_gene5_22084837'}
PEAKS = {'full8': REC/'centered_hermite_pathwise_full8_22085153',
         'variance2': REC/'centered_hermite_pathwise_variance2_22085154',
         'variance5': REC/'centered_hermite_pathwise_variance5_22085155',
         'refusal8': REC/'centered_hermite_pathwise_refusal8_22085156'}
SOURCES = ('likelihood_v1.py', 'design_v1.py', 'bfix_v1.py', 'tests_v1.py',
           'repair_p3_local_moments.py', 'repair_p3_pathwise_moments_v1.py',
           'repair_p3_longdouble_conditional_v1.py', 'repair_p3_adaptive_intervals_weighted_v1.py')


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case', choices=tuple(PEAKS), required=True)
    ap.add_argument('--pins', type=Path, required=True)
    ap.add_argument('--expected-pins-sha', required=True)
    ap.add_argument('--expected-runner-sha', required=True)
    ap.add_argument('--inspect', action='store_true')
    a = ap.parse_args()
    require(sha(Path(__file__)) == a.expected_runner_sha and sha(a.pins) == a.expected_pins_sha, 'Runner/pins differ')
    pins = json.loads(a.pins.read_text())
    required = [*(HERE/n for n in SOURCES),
        *(FIXTURE/n for n in ('globals.json', 'cohorts.tsv', 'individuals.tsv', 'rows.tsv')),
        *(FULL/n for n in ('protocol.json', 'answer.json', 'summary.json')),
        POINT/'draw_counts.json', POINT/'recipe_summary.json',
        REC/'standard_hermite_refusal_check_22085142/refusal.json',
        REC/'standard_hermite_pilot_draw0_22084966/draw_counts.json',
        REC/'refusal_adaptive_intervals_22085179/summary.json',
        *(d/n for d in DECIMAL.values() for n in ('protocol.json', 'summary.json')),
        *(d/n for d in PEAKS.values() for n in ('protocol.json', 'gene_integrals.json'))]
    require(all(str(p) in pins for p in required), 'Missing input pin')
    for path, digest in pins.items():
        p = Path(path)
        require(p.resolve().is_relative_to(ROOT.resolve()) and sha(p) == digest, 'Input differs')
        if p.suffix == '.py':
            ast.parse(p.read_text())
    read = lambda path: json.loads(path.read_text())
    require(read(FULL/'summary.json')['completed'] and read(FULL/'answer.json')['returned'], 'Full reference incomplete')
    require(all(read(d/'summary.json')['passed'] for d in DECIMAL.values()), 'Decimal reference failed')
    if a.inspect:
        print(json.dumps(dict(case=a.case, probability_weighted_conditioning=True, domains=[24., 26.], no_fit=True)))
        return
    cpus = 8 if a.case in ('full8', 'refusal8') else 1
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID')
        and os.environ.get('SLURM_CPUS_PER_TASK') == str(cpus), 'Individual compute allocation differs')
    require(all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
        'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')), 'Threads differ')
    out = REC/('adaptive_weighted_reference_'+a.case+'_'+os.environ['SLURM_JOB_ID'])
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
    model = design.model()
    if a.case == 'full8':
        reference = read(FULL/'protocol.json')
        parameters, variance = reference['parameters'], reference['variance']
        require(hashlib.sha256(np.asarray([r.a for r in model.rows], dtype='<i8').tobytes()).hexdigest()
            == reference['count_vector_sha256'], 'Full counts differ')
        genes = list(range(8))
    elif a.case == 'refusal8':
        counts = read(REC/'standard_hermite_pilot_draw0_22084966/draw_counts.json')
        require(counts['full_row_keys'] == [[r.individual, settings['gene_order'][r.gene]] for r in model.rows]
            and len(counts['counts']) == len(model.rows) == 173, 'Refusal count order differs')
        require(hashlib.sha256(np.asarray(counts['counts'], dtype='<i8').tobytes()).hexdigest()
            == counts['count_vector_sha256'], 'Refusal counts/hash differ')
        model = modules['likelihood_v1'].Model([replace(r, a=c) for r, c in zip(model.rows, counts['counts'])], gene_count=8)
        failure = read(REC/'standard_hermite_refusal_check_22085142/refusal.json')['failure']
        parameters, variance = failure['parameters'], failure['variance']
        genes = list(range(8))
    else:
        counts = read(POINT/'draw_counts.json')
        require(counts['full_row_keys'] == [[r.individual, settings['gene_order'][r.gene]] for r in model.rows]
            and len(counts['counts']) == len(model.rows) == 173, 'Count order differs')
        model = modules['likelihood_v1'].Model([replace(r, a=c) for r, c in zip(model.rows, counts['counts'])], gene_count=8)
        failure = read(POINT/'recipe_summary.json')['fitting_integration_failure']
        parameters, variance = failure['parameters'], failure['variance']
        genes = [int(a.case[-1])]
    require((model.G, model.Z, model.size, len(model.rows)) == (8, 7, 28, 173), 'Synthetic scope differs')
    peaks_protocol = read(PEAKS[a.case]/'protocol.json')
    require(peaks_protocol['parameters'] == parameters and peaks_protocol['variance'] == variance
        and peaks_protocol['genes'] == genes, 'Peak source point differs')
    peaks = {g['gene']: g['peaks'] for g in read(PEAKS[a.case]/'gene_integrals.json')['records']}
    require(set(peaks) == set(genes), 'Peak roster differs')
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
    save('protocol', dict(case=a.case, parameters=p, variance=variance, genes=genes,
        source_and_input_sha256=pins, environment=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
        probability_weighted_conditioning=True, no_fit=True, no_calibration=True, no_human_data=True))
    worker_module = modules['repair_p3_adaptive_intervals_weighted_v1']
    worker_module.WORK = model, p, variance, peaks, modules
    started = time.perf_counter()
    if cpus == 8:
        with multiprocessing.get_context('fork').Pool(8) as pool:
            answers = pool.map(worker_module.worker, genes)
    else:
        answers = [worker_module.worker(genes[0])]
    save('gene_integrals', dict(records=answers))
    numerical_pass = all(g['returned'] for g in answers)
    comparable = all('logmass' in g for g in answers)
    differences = None
    ll, gradient, hessian = None, None, None
    if comparable:
        ll = sum(g['logmass'] for g in answers)
        gradient = np.sum([g['gradient'] for g in answers], axis=0)
        hessian = np.sum([g['hessian'] for g in answers], axis=0)
        if a.case == 'full8':
            reference = read(FULL/'answer.json')
            differences = [abs(ll-reference['logmass']), float(np.max(np.abs(gradient-reference['gradient']))),
                           float(np.max(np.abs(hessian-reference['hessian'])))]
        elif a.case == 'refusal8':
            reference = read(REC/'refusal_adaptive_intervals_22085179/summary.json')
            require(reference['passed'], 'Refusal reference failed')
            differences = [abs(ll-reference['logmass']), float(np.max(np.abs(gradient-reference['gradient']))),
                           float(np.max(np.abs(hessian-reference['hessian'])))]
        else:
            reference = read(DECIMAL[genes[0]]/'summary.json')['results'][-1]
            dp = read(DECIMAL[genes[0]]/'protocol.json')
            require(dp['parameters'] == parameters and dp['variance'] == variance, 'Decimal point differs')
            with localcontext() as context:
                context.prec = 80
                differences = [float(abs(Decimal(reference[k])-Decimal.from_float(float(v)))) for k, v in
                    zip(('logmass', 'variance_score', 'variance_curvature'), (ll, gradient[-1], hessian[-1, -1]))]
    passed = numerical_pass and differences is not None and all(x <= y for x, y in zip(differences, (1e-9, 1e-8, 1e-5)))
    require(p.tobytes() == before and repr(model.rows) == row_before, 'Inputs changed')
    for path, digest in pins.items():
        require(sha(Path(path)) == digest, 'Input changed during check')
    save('summary', dict(completed=True, case=a.case, passed=passed, numerical_allowances_passed=numerical_pass,
        reference_differences=differences, logmass=ll, gradient=gradient, hessian=hessian,
        seconds=time.perf_counter()-started, no_fit=True, no_calibration=True,
        no_independent_full_Hessian_validation=True))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
