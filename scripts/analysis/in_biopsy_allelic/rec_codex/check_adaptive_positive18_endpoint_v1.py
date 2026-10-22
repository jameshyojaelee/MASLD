#!/usr/bin/env python3
"""Check a saved positive-LR Gamma endpoint without fitting or reclassification."""
import argparse
import ast
from dataclasses import replace
import hashlib
import importlib
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE = ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
POINT = REC/'repair_p3_gamma_draw_18_22081587'
SOURCES = ('likelihood_v1.py', 'design_v1.py', 'bfix_v1.py', 'tests_v1.py',
    'repair_p3_local_moments.py', 'repair_p3_pathwise_moments_v1.py',
    'repair_p3_longdouble_conditional_v1.py', 'repair_p3_adaptive_intervals_weighted_v1.py',
    'repair_p3_adaptive_peak_scan_v1.py', 'repair_p3_adaptive_interval_model_v1.py')


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pins', type=Path, required=True)
    ap.add_argument('--expected-pins-sha', required=True)
    ap.add_argument('--expected-runner-sha', required=True)
    ap.add_argument('--inspect', action='store_true')
    a = ap.parse_args()
    require(sha(Path(__file__)) == a.expected_runner_sha and sha(a.pins) == a.expected_pins_sha, 'Runner/pins differ')
    pins = json.loads(a.pins.read_text())
    require(all(str(p) in pins for p in [*(HERE/n for n in SOURCES),
        *(FIXTURE/n for n in ('globals.json', 'cohorts.tsv', 'individuals.tsv', 'rows.tsv')),
        *(POINT/n for n in ('p3_018.json', 'draw_counts.json', 'returned_fit_diagnostic.json', 'summary.json'))]), 'Missing pin')
    for path, digest in pins.items():
        p = Path(path)
        require(p.resolve().is_relative_to(ROOT.resolve()) and sha(p) == digest, 'Input differs')
        if p.suffix == '.py':
            ast.parse(p.read_text())
    read = lambda path: json.loads(path.read_text())
    raw = read(POINT/'p3_018.json')
    fit = raw['fit']
    require(fit['converged'] and fit['null']['converged'] and fit['v'] > 0 and raw['LR'] > 0, 'Positive fit unavailable')
    require(read(POINT/'summary.json')['diagnostic_completed'], 'Reference diagnostic incomplete')
    if a.inspect:
        print(json.dumps(dict(saved_positive_variance=fit['v'], saved_LR=raw['LR'], no_fit=True)))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID')
        and os.environ.get('SLURM_CPUS_PER_TASK') == '8', 'Individual eight-CPU compute required')
    require(all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
        'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')), 'Threads differ')
    out = REC/('adaptive_positive18_endpoint_'+os.environ['SLURM_JOB_ID'])
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
    counts = read(POINT/'draw_counts.json')
    require(counts['full_row_keys'] == [[r.individual, settings['gene_order'][r.gene]] for r in template.rows]
        and len(counts['counts']) == len(template.rows) == 173
        and hashlib.sha256(np.asarray(counts['counts'], dtype='<i8').tobytes()).hexdigest()
        == counts['count_vector_sha256'], 'Counts differ')
    rows = [replace(r, a=c) for r, c in zip(template.rows, counts['counts'])]
    wrapper = modules['repair_p3_adaptive_interval_model_v1']
    Model = wrapper.model_class(modules['likelihood_v1'], modules,
        modules['repair_p3_adaptive_intervals_weighted_v1'], modules['repair_p3_adaptive_peak_scan_v1'])
    model = Model(rows, gene_count=8)
    require((model.G, model.Z, model.size) == (8, 7, 28), 'Synthetic scope differs')
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
    save('protocol', dict(parameters=fit['parameters'], variance=fit['v'], source_and_input_sha256=pins,
        environment=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
        same_point_check_only=True, no_fit=True, no_calibration=True, no_human_data=True))
    started = time.perf_counter()
    try:
        ll, g, h = model.marginal(np.asarray(fit['parameters']), fit['v'], with_hessian=True)
    except wrapper.IntervalRefusal as exc:
        save('refusal', dict(message=str(exc), failure=model.last_integration_failure))
        raise
    save('gene_integrals', dict(records=model.last_interval_records))
    reference = read(POINT/'returned_fit_diagnostic.json')
    differences = [abs(ll-fit['loglik']), float(np.max(np.abs(g-fit['gradient'])))]
    agree = differences[0] <= 1e-6 and bool(np.all(np.abs(g-fit['gradient']) <= 1e-8+1e-6*np.abs(fit['gradient'])))
    point = np.r_[fit['parameters'], fit['v']]
    bounds = [b if keep else (point[i], point[i]) for i, (b, keep) in enumerate(zip(model.bounds(), model.informative()))]+[(0., None)]
    require(safe(bounds) == reference['bounds'], 'Projected-score bounds differ')
    pg = g.copy()
    for i, (lo, hi) in enumerate(bounds):
        if (lo == hi and lo is not None) or (lo is not None and point[i] <= lo+1e-10 and g[i] < 0) or (hi is not None and point[i] >= hi-1e-10 and g[i] > 0):
            pg[i] = 0.
    maximum = float(np.max(np.abs(pg)))
    null_ll, null_g = model.marginal(np.asarray(fit['null']['parameters']), 0.)
    lr = max(0., 2*(ll-null_ll))
    passed = agree and maximum < 1e-6 and abs(lr-raw['LR']) <= 1e-6
    for path, digest in pins.items():
        require(sha(Path(path)) == digest, 'Input changed during check')
    save('summary', dict(completed=True, passed=passed, loglik=ll, gradient=g, hessian=h,
        projected_gradient=pg, bounds=bounds, max_projected_gradient=maximum,
        reference_LL_gradient_differences=differences, same_point_null_loglik=null_ll,
        same_point_null_gradient=null_g, same_point_LR=lr, saved_LR=raw['LR'],
        seconds=time.perf_counter()-started, no_fit=True, no_full_refit_agreement=True,
        no_calibration=True, no_independent_full_Hessian_validation=True))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
