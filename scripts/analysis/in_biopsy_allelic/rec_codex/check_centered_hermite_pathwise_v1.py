#!/usr/bin/env python3
"""Centered importance Hermite fixed points; refusal case has no external reference."""
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
POINT = REC/'repair_p3_gamma_order_polish_draw_0_amd0_22084420'
FULL = REC/'repair_p3_gamma_batched_full_check_22078827'
REFERENCES = {2: (REC/'pathwise_longdouble_gene2_22084806', REC/'decimal_variance_curvature_gene2_22084836'),
              5: (REC/'pathwise_longdouble_gene5_22084793', REC/'decimal_variance_curvature_gene5_22084837')}
SOURCES = ('likelihood_v1.py', 'design_v1.py', 'bfix_v1.py', 'tests_v1.py',
           'repair_p3_local_moments.py', 'repair_p3_pathwise_moments_v1.py',
           'repair_p3_longdouble_conditional_v1.py', 'repair_p3_centered_hermite_pathwise_v1.py')
REFUSAL = REC/'standard_hermite_refusal_check_22085142'
WORK = None


def require(ok, message):
    if not ok: raise RuntimeError(message)


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def worker(gene):
    model, p, variance, local_module, path_module, conditional_module, integrator = WORK
    rows = tuple(r for r in model.rows if r.gene == gene)
    local = local_module.LocalGeneMoments(model, p, rows, None)
    context = conditional_module.LongDoubleConditional(model, p, rows)
    moments = path_module.PathwiseGeneMoments(local)
    began = time.perf_counter()
    answer = integrator.integrate(local, context, moments, variance)
    answer['seconds'] = time.perf_counter()-began
    return answer


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case', choices=('full8', 'variance2', 'variance5', 'refusal8'), required=True)
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
    required = [*(HERE/n for n in SOURCES), *(FIXTURE/n for n in ('globals.json', 'cohorts.tsv', 'individuals.tsv', 'rows.tsv')),
        *(FULL/n for n in ('protocol.json', 'answer.json', 'summary.json')),
        POINT/'draw_counts.json', POINT/'recipe_summary.json',
        REFUSAL/'refusal.json', REC/'standard_hermite_pilot_draw0_22084966/draw_counts.json',
        *(d/n for pair in REFERENCES.values() for d in pair for n in (('protocol.json', 'summary.json'))),
        *(pair[0]/'bound26.json' for pair in REFERENCES.values())]
    require(all(str(p) in pins for p in required), 'Missing source/reference pin')
    read = lambda p: json.loads(p.read_text())
    fullp, full = read(FULL/'protocol.json'), read(FULL/'answer.json')
    require(full['returned'] and read(FULL/'summary.json')['completed'], 'Full reference incomplete')
    for _, d in REFERENCES.values(): require(read(d/'summary.json')['passed'], 'Independent Decimal reference failed')
    if a.inspect:
        print(json.dumps(dict(case=a.case, fixed_orders=[32,64,128,256], minimum_final_order=128, no_fit=True)))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'), 'Individual compute required')
    cpus = 8 if a.case in ('full8', 'refusal8') else 1
    require(os.environ.get('SLURM_CPUS_PER_TASK') == str(cpus), 'CPU allocation differs')
    require(all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')), 'Threads differ')
    out = REC/('centered_hermite_pathwise_'+a.case+'_'+os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)
    archive = out/'executed_sources'; archive.mkdir()
    for p in (Path(__file__), *(HERE/n for n in SOURCES)):
        q = archive/p.name; q.write_bytes(p.read_bytes())
        require(sha(q) == (a.expected_runner_sha if p == Path(__file__) else pins[str(p)]), 'Archive differs')
    import numpy as np
    import scipy
    sys.path.insert(0, str(archive))
    modules = {n[:-3]: importlib.import_module(n[:-3]) for n in SOURCES}
    require(all(Path(m.__file__).resolve().parent == archive.resolve() for m in modules.values()), 'Import escaped archive')
    loader, original = modules['bfix_v1'], modules['likelihood_v1']
    require(loader.FIXTURE.resolve() == FIXTURE.resolve(), 'Fixture differs')
    design, settings = loader.load(); template = design.model()
    if a.case == 'full8':
        model = template
        parameters, variance = fullp['parameters'], fullp['variance']
        count_sha = hashlib.sha256(np.asarray([r.a for r in model.rows], dtype='<i8').tobytes()).hexdigest()
        require(count_sha == fullp['count_vector_sha256'] == '7009be4adeec6157fa4bbf3f46e48ec24aa398a0bee4d7906c78c231b7d3fbad', 'Full counts differ')
        genes = list(range(8))
    elif a.case == 'refusal8':
        counts = read(REC/'standard_hermite_pilot_draw0_22084966/draw_counts.json')
        require(counts['full_row_keys'] == [[r.individual, settings['gene_order'][r.gene]] for r in template.rows], 'Refusal count order differs')
        require(hashlib.sha256(np.asarray(counts['counts'], dtype='<i8').tobytes()).hexdigest() == counts['count_vector_sha256'], 'Refusal counts/hash differ')
        model = original.Model([replace(r, a=c) for r,c in zip(template.rows,counts['counts'])],gene_count=8)
        failure = read(REFUSAL/'refusal.json')['failure']
        parameters, variance = failure['parameters'], failure['variance']
        genes = list(range(8))
    else:
        gene = int(a.case[-1]); counts = read(POINT/'draw_counts.json')
        require(counts['full_row_keys'] == [[r.individual, settings['gene_order'][r.gene]] for r in template.rows]
            and len(counts['counts']) == len(template.rows) == 173, 'Fixed count order differs')
        model = original.Model([replace(r, a=c) for r, c in zip(template.rows, counts['counts'])], gene_count=8)
        failure = read(POINT/'recipe_summary.json')['fitting_integration_failure']
        parameters, variance = failure['parameters'], failure['variance']
        require(read(REFERENCES[gene][0]/'protocol.json')['parameters'] == parameters
            and read(REFERENCES[gene][0]/'protocol.json')['variance'] == variance, 'Variance reference differs')
        genes = [gene]
    require((model.G, model.Z, model.size, len(model.rows)) == (8,7,28,173), 'Synthetic scope differs')
    p = np.asarray(parameters); p.setflags(write=False)
    before, row_before = p.tobytes(), repr(model.rows)
    def safe(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, dict): return {k:safe(v) for k,v in x.items()}
        if isinstance(x, (tuple,list)): return [safe(v) for v in x]
        return x
    def save(name, value):
        with (out/(name+'.json')).open('x') as f:
            json.dump(safe(value), f, indent=2, allow_nan=False); f.write('\n'); f.flush(); os.fsync(f.fileno())
    save('protocol', dict(case=a.case, parameters=p, variance=variance, genes=genes, source_and_input_sha256=pins,
        environment=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
        node=os.environ.get('SLURMD_NODENAME'), no_fit=True, no_calibration=True,
        departure_from_original_adaptive20_approximation=True, no_human_data=True, empirical_not_certified=True))
    global WORK
    WORK = (model, p, variance, modules['repair_p3_local_moments'], modules['repair_p3_pathwise_moments_v1'],
            modules['repair_p3_longdouble_conditional_v1'], modules['repair_p3_centered_hermite_pathwise_v1'])
    start = time.perf_counter()
    if cpus == 8:
        with multiprocessing.get_context('fork').Pool(8) as pool: answers = pool.map(worker, genes)
    else: answers = [worker(genes[0])]
    save('gene_integrals', dict(records=answers))
    # Preserve complete numerical output even when a gene refuses order agreement.
    all_returned = all(x['returned'] for x in answers)
    ll = sum(x['logmass'] for x in answers)
    g = np.sum([x['gradient'] for x in answers], axis=0)
    h = np.sum([x['hessian'] for x in answers], axis=0)
    if a.case == 'full8':
        differences = [abs(ll-full['logmass']), float(np.max(np.abs(g-full['gradient']))), float(np.max(np.abs(h-full['hessian'])))]
        allowance = [sum(abs(x['orders'][-1]['logmass']-x['orders'][-2]['logmass']) for x in answers),
            float(np.max(np.sum([np.abs(x['orders'][-1]['gradient']-x['orders'][-2]['gradient']) for x in answers], axis=0))),
            float(np.max(np.sum([np.abs(x['orders'][-1]['hessian']-x['orders'][-2]['hessian']) for x in answers], axis=0)))]
        passed = all(a <= b for a,b in zip(differences,(1e-9,1e-8,1e-5))) and all(a <= b for a,b in zip(allowance,(1e-8,1e-7,1e-5)))
        detail = dict(full_GL_differences=differences, aggregate_order_allowance_maxima=allowance)
    elif a.case == 'refusal8':
        passed = all(x['returned'] for x in answers)
        detail = dict(order_check_only=True, independent_large_variance_reference_not_run=True)
    else:
        ref = read(REFERENCES[genes[0]][1]/'summary.json')['results'][-1]
        with localcontext() as ctx:
            ctx.prec = 80
            differences = [float(abs(Decimal(ref[k])-Decimal.from_float(float(v)))) for k,v in
                           zip(('logmass','variance_score','variance_curvature'),(ll,g[-1],h[-1,-1]))]
        passed = all(a <= b for a,b in zip(differences,(1e-9,1e-8,1e-5)))
        detail = dict(independent_Decimal_differences=differences)
    require(p.tobytes() == before and repr(model.rows) == row_before, 'Inputs changed')
    for path,digest in pins.items(): require(sha(Path(path)) == digest, 'Pinned input changed')
    save('summary', dict(completed=True, case=a.case, passed=bool(passed and all_returned), all_genes_returned=all_returned, logmass=ll, gradient=g, hessian=h,
        seconds=time.perf_counter()-start, node_evaluations=sum(x['node_evaluations'] for x in answers),
        **detail, no_fit=True, no_calibration=True, no_independent_full_Hessian_validation=True))
    if not passed or not all_returned: raise SystemExit(1)


if __name__ == '__main__': main()
