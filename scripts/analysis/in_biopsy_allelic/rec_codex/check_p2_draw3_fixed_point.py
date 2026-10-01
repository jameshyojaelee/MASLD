"""Compute-only conditional derivative check at the saved toy P2 draw3 point.

No fit, optimizer, marginal model, generated draw, or inference replacement.
The executed 21998518 source and its supplied bootstrap are the only inputs.
"""
import argparse
import csv
from dataclasses import asdict
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import sys

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
ORIGINAL = REC / 'bfix_p2_21998518'
SNAPSHOT = ORIGINAL / 'executed_sources'
FIXTURE = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
STABLE_PROBE = REC / 'codex_p3_gene2_integrand_21998880/executed_sources/codex_stable_bb_probe.py'
GUARDS = {
    ORIGINAL / 'receipt.json': '63b377f65b4f37c8b6d60ffca2745f639ffa44519a7bdd2a0a55b682402f8cd8',
    ORIGINAL / 'p2_003.json': 'bea38141306a6ddc45ed20cc537501767fb0461f5803fb6d80a47c6d72aa660c',
    ORIGINAL / 'observed_null.json': 'a617d331fd3fb93576e5bd8588581aaab8f77946936536401905490438660181',
    STABLE_PROBE: '723781113c17061b9a28341a91149bbd4ab71f195121a67f88d8cbc710a21488',
}
INPUT_NAMES = ('globals.json', 'individuals.tsv', 'rows.tsv', 'cohorts.tsv',
               'draws/p2_resample.tsv')
MODULE_NAMES = ('likelihood_v1', 'tests_v1', 'design_v1', 'bfix_v1')
COORDINATES = {7: 'alpha_G8', 9: 'r_G2_lower_bound', 12: 'r_G5_lower_bound',
               15: 'r_G8', 16: 'kappa'}
STEPS = (0.001, 0.0003, 0.0001, 0.00003)
OUTPUT_CAP = 32 << 20


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    require(path.is_file() and path.stat().st_size <= 1_000_000, 'Missing/bounded toy source: '+str(path))
    return hashlib.sha256(path.read_bytes()).hexdigest()


def guard(paths):
    for path, expected in paths.items():
        require(sha(path) == expected, 'Toy input/source changed: '+str(path))


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def projected_score(p, score, bounds):
    result = score.copy()
    for j, (lo, hi) in enumerate(bounds):
        if ((lo is not None and lo == hi) or
                (lo is not None and p[j] <= lo+1e-8 and score[j] < 0) or
                (hi is not None and p[j] >= hi-1e-8 and score[j] > 0)):
            result[j] = 0
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Numerical diagnostic requires compute SLURM job')
    require(args.output.parent.resolve() == REC.resolve()
            and args.output.name == 'p2_draw3_fixed_point_'+os.environ['SLURM_JOB_ID']
            and not args.output.exists(), 'Require fresh own REC/actual-job-ID output')
    guard(GUARDS)
    receipt = json.loads((ORIGINAL/'receipt.json').read_text())
    require(receipt['role'] == 'p2' and receipt['spec'] == 'v1.7' and receipt['B'] == 50
            and Path(receipt['fixture']) == FIXTURE, 'Different original toy role/fixture')
    inputs = {FIXTURE/name: receipt['fixture_sha256'][name] for name in INPUT_NAMES}
    sources = {SNAPSHOT/(name+'.py'): receipt['source_sha256'][name+'.py'] for name in MODULE_NAMES}
    guard(inputs); guard(sources)
    require(sources[SNAPSHOT/'likelihood_v1.py'] ==
            '41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6',
            'Original executed numerical source changed')
    # Import only the archived module definitions, not either executed runner.
    require(not any(name in sys.modules for name in MODULE_NAMES), 'Toy module already imported')
    sys.path.insert(0, str(SNAPSHOT))
    modules = {name: importlib.import_module(name) for name in MODULE_NAMES}
    for name, module in modules.items():
        require(Path(module.__file__).resolve() == (SNAPSHOT/(name+'.py')).resolve(),
                'Import did not use original archived module: '+name)
    import numpy as np
    import scipy
    design, settings = modules['bfix_v1'].load()
    require(settings['B'] == 50 and len(design.participants) == 60 and len(design.rows) == 173,
            'Supplied bootstrap/source population changed')
    groups = design.cohort_clusters()
    group_index = {c: {design.clusters[members[0]]: j for j, members in enumerate(v)}
                   for c, v in groups.items()}
    picks = {c: [] for c in groups}; choices = []
    with (FIXTURE/'draws/p2_resample.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if int(row['draw']) != 3:
                continue
            c = row['stratum']
            require(c in picks and int(row['slot']) == len(picks[c]), 'Original draw3 slot sequence differs')
            require(row['cluster_id'] in group_index[c], 'Original draw3 cluster missing')
            picks[c].append(group_index[c][row['cluster_id']]); choices.append(row)
    require(all(len(picks[c]) == len(groups[c]) for c in groups), 'Incomplete supplied draw3')
    require(len(choices) == 58, 'Original supplied cluster-slot count differs')
    sample = design.resample(picks); model = sample.model()
    saved = json.loads((ORIGINAL/'p2_003.json').read_text()); fitted = saved['fit']
    require(saved['b'] == 3 and fitted['kept_start'] == 'warm' and not fitted['converged']
            and fitted['start_converged'] == {'warm': False, 'cold': True}
            and np.isnan(saved['kappa']), 'Original failed warm slot differs')
    require(fitted['start_logliks']['cold']-fitted['start_logliks']['warm'] < 1e-6,
            'Original prescribed likelihood start selection differs')
    p = np.asarray(fitted['parameters'], dtype=float)
    warm = np.asarray(json.loads((ORIGINAL/'observed_null.json').read_text())['parameters'], dtype=float)
    require(model.G == 8 and model.Z == 7 and p.shape == warm.shape == (28,)
            and np.isfinite(p).all(), 'Original parameter schema differs')
    informative = model.informative()
    require(np.array_equal(p[~informative], warm[~informative]), 'No-information fixed coordinates differ')
    bounds = [b if informative[j] else (warm[j], warm[j]) for j, b in enumerate(model.bounds())]
    require(all((lo is None or p[j] >= lo) and (hi is None or p[j] <= hi)
                for j, (lo, hi) in enumerate(bounds)), 'Retained warm point violates original bounds')
    ll, score, hessian = model.evaluate(p)
    require(np.isfinite(ll) and np.isfinite(score).all() and np.isfinite(hessian).all(),
            'Nonfinite reconstructed conditional point')
    saved_score = np.asarray(fitted['gradient'], dtype=float)
    saved_hessian = np.asarray(fitted['hessian'], dtype=float)
    projected = projected_score(p, score, bounds)
    # A separate in-memory module replaces only BB arithmetic with the existing
    # guarded log-shape probe. Original executed modules and files stay intact.
    def load_named(name, path):
        require(name not in sys.modules, 'Diagnostic module already loaded')
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec); sys.modules[name] = module
        spec.loader.exec_module(module)
        return module
    stable_probe = load_named('p2_draw3_stable_probe', STABLE_PROBE)
    stable_likelihood = load_named('p2_draw3_stable_likelihood', SNAPSHOT/'likelihood_v1.py')
    def replacement(k, n, eta, linear_r, need_hessian=True):
        value, gradient, hh = stable_probe.stable_bb_terms(k, n, eta, linear_r)
        return value, gradient, hh if need_hessian else None
    stable_likelihood.bb_terms = replacement
    stable_model = stable_likelihood.Model(model.rows, gene_count=model.G)
    stable_ll, stable_score, stable_hessian = stable_model.evaluate(p)
    require(np.isfinite(stable_ll) and np.isfinite(stable_score).all()
            and np.isfinite(stable_hessian).all(), 'Nonfinite isolated stable conditional point')
    result = dict(scope='Synthetic conditional fixed point only; not refit, correction, or calibration',
        original_job=21998518, draw=3, kept_start='warm', original_failed_slot_retained=True,
        original_kappa_slot='NaN', original_lambda_pair=saved['lambda'],
        original_start_logliks=fitted['start_logliks'], original_start_converged=fitted['start_converged'],
        input_sha256={str(k):v for k,v in {**GUARDS, **inputs, **sources}.items()},
        script_sha256=sha(Path(__file__)),
        launcher_sha256=sha(Path(__file__).with_name('run_check_p2_draw3_fixed_point.sbatch')),
        python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
        bootstrap_choices=choices, bootstrap_cluster_indices=picks,
        original_toy_individual_n=60, original_toy_row_n=173, bootstrap_cluster_slots=58,
        sampled_participant_n=len(sample.participants), sampled_row_n=len(model.rows),
        model_row_sha256=canonical_hash([asdict(r) for r in model.rows]),
        count_vector_sha256=canonical_hash([[r.individual, r.gene, r.n, r.a] for r in model.rows]),
        original_row_hash_available=False,
        no_information_coordinates=np.flatnonzero(~informative).tolist(),
        original_effective_bounds=bounds, parameters=p.tolist(),
        recomputed_loglik=float(ll), saved_loglik=float(fitted['loglik']),
        saved_minus_recomputed_loglik=float(fitted['loglik']-ll),
        max_saved_score_difference=float(np.max(np.abs(saved_score-score))),
        max_saved_hessian_difference=float(np.max(np.abs(saved_hessian-hessian))),
        analytic_score=score.tolist(), projected_score=projected.tolist(),
        max_projected_score=float(np.max(np.abs(projected))),
        original_strict_tolerance=1e-6, recomputed_converged=bool(np.max(np.abs(projected)) < 1e-6),
        analytic_hessian_symmetry_error=float(np.max(np.abs(hessian-hessian.T))),
        isolated_stable=dict(loglik=float(stable_ll), score=stable_score.tolist(),
            projected_score=projected_score(p, stable_score, bounds).tolist(),
            original_minus_stable_loglik=float(ll-stable_ll),
            max_original_score_difference=float(np.max(np.abs(score-stable_score))),
            max_original_hessian_difference=float(np.max(np.abs(hessian-stable_hessian))),
            saved_minus_stable_loglik=float(fitted['loglik']-stable_ll),
            max_saved_score_difference=float(np.max(np.abs(saved_score-stable_score))),
            max_saved_hessian_difference=float(np.max(np.abs(saved_hessian-stable_hessian))),
            hessian_symmetry_error=float(np.max(np.abs(stable_hessian-stable_hessian.T)))),
        finite_difference=[])
    # Five-point score stencil and its Hessian-column analogue. At an
    # active bound use the admissible one-sided stencil; never step outside it.
    for representation, fixed_model, base_ll, base_score, base_hessian in (
            ('original', model, ll, score, hessian),
            ('isolated_stable_bb', stable_model, stable_ll, stable_score, stable_hessian)):
        for j, name in COORDINATES.items():
            if not informative[j]:
                result['finite_difference'].append(dict(representation=representation, coordinate=j,
                                                       name=name, fixed_no_information=True))
                continue
            lo, hi = bounds[j]
            for h in STEPS:
                if hi is not None and p[j]+2*h > hi:
                    offsets=(0, -1, -2, -3, -4); weights=(25, -48, 36, -16, 3)
                    second=(35, -104, 114, -56, 11); mode='backward_admissible'
                elif lo is not None and p[j]-2*h < lo:
                    offsets=(0, 1, 2, 3, 4); weights=(-25, 48, -36, 16, -3)
                    second=(35, -104, 114, -56, 11); mode='forward_admissible'
                else:
                    offsets=(-2, -1, 0, 1, 2); weights=(1, -8, 0, 8, -1)
                    second=(-1, 16, -30, 16, -1); mode='central'
                values=[]; gradients=[]
                for offset in offsets:
                    q=p.copy(); q[j]+=offset*h
                    require((lo is None or q[j] >= lo) and (hi is None or q[j] <= hi), 'FD left original bounds')
                    value, gradient, _ = fixed_model.evaluate(q, need_hessian=False)
                    values.append(value); gradients.append(gradient)
                # Center before weighted subtraction to avoid unnecessary loss
                # from repeatedly subtracting the common likelihood/score offset.
                centered_values=np.asarray(values)-base_ll
                fd_score=float(np.dot(weights, centered_values)/(12*h))
                fd_column=np.asarray(weights) @ (np.stack(gradients)-base_score)/(12*h)
                fd_diagonal=float(np.dot(second, centered_values)/(12*h*h))
                require(np.isfinite(fd_score) and np.isfinite(fd_column).all() and np.isfinite(fd_diagonal), 'Nonfinite FD')
                result['finite_difference'].append(dict(representation=representation, coordinate=j,
                    name=name, step=h, stencil=mode,
                    analytic_score=float(base_score[j]), finite_difference_loglik_score=fd_score,
                    score_difference=fd_score-float(base_score[j]),
                    analytic_hessian_column=base_hessian[:, j].tolist(), finite_difference_score_column=fd_column.tolist(),
                    max_hessian_column_difference=float(np.max(np.abs(fd_column-base_hessian[:, j]))),
                    analytic_hessian_diagonal=float(base_hessian[j, j]), finite_difference_loglik_diagonal=fd_diagonal,
                    diagonal_difference=fd_diagonal-float(base_hessian[j, j])))
    result['interpretation_limit']='Focused derivatives only, not whole-fit stationarity. Inspect stability across steps; cancellation/truncation may prevent a verdict. No start chosen, parameter updated, or failed draw replaced.'
    guard(GUARDS); guard(inputs); guard(sources)
    require(not args.output.exists(), 'Own diagnostic output already exists')
    args.output.mkdir(parents=False)
    archive=args.output/'executed_sources';archive.mkdir()
    for path in (*sources, STABLE_PROBE, Path(__file__), Path(__file__).with_name('run_check_p2_draw3_fixed_point.sbatch')):
        (archive/path.name).write_bytes(path.read_bytes())
    with (args.output/'summary.json').open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False);handle.write('\n')
    require(sum(f.stat().st_size for f in args.output.rglob('*') if f.is_file()) <= OUTPUT_CAP,
            'Own diagnostic output cap exceeded; preserve outputs')
    print(json.dumps(dict(draw=3, fitted=False, max_projected_score=result['max_projected_score'],
                          finite_difference_points=len(result['finite_difference'])), allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
