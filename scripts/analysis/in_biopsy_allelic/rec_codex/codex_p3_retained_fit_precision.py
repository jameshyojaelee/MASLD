"""Read-only synthetic retained-point precision check of the prescribed 20-node objective.

No optimizer is called. Logged AST differs only by a diagnostic callback;
immutable likelihood bytes and original bounds/start/statistic rules remain.
Finite-difference uncertainty is measured, not a certified error bound.
"""
import ast
import copy
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import signal
import sys
import time
import traceback
import warnings

import numpy as np
import scipy


ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
OUTER = REC / 'codex_p3_private_damping64_21998950'
REPLAY = OUTER / 'codex_replay_first_p3_stable_bb_21998950'
ROLE = REC / 'bfix_p3_21998519'
FIXTURE = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
HERE = Path(__file__).resolve().parent
EXPECTED_COUNT = '5f4a1aa455cd09b0d5f31889021ec88b145dcbb50c8e5b600343ea9a3a992691'
PROBE_SHA = '723781113c17061b9a28341a91149bbd4ab71f195121a67f88d8cbc710a21488'
SOURCES = {
    'likelihood_v1.py': '8880a8059ca18259a43b638ce9b0c60991ab7b555ec3db736d31f85d548f8ff8',
    'bfix_v1.py': '324594ce783ba2464921eaee5924065f579060eded6043a48bdf4b7e392bfa56',
    'design_v1.py': 'cb3163896f982f39ba94cdbf227a8e7d4c7f85a0a646b26542868e6fa14c9c5e',
    'tests_v1.py': '52cd14d43b649c20347330b64b467dbd81320f10252fe85ce3c538ceef2f092e',
    'check_bfix_tests_v13.py': 'a4ca55a52d28dfdbfeb0fc2611c2dcb3ec412987c8f7618f59002cbf72b08c2d',
}
GUARDS = {
    OUTER / 'protocol.json': 'd3a5ff66365d9ff828c9e0ce01533e65bf0fe8ab6ae436ee6fa4cdd43c7ae2b1',
    OUTER / 'fixed_gene_precheck.json': '85bb770af3bf814ff64b7035cff9bc106b47f165fcbbbe440310ea161159a1e8',
    REPLAY / 'receipt.json': '145c4e523616d4c1fe1eb54b144415d5eec403d768096de952925112c25d19dc',
    OUTER / 'executed_sources/codex_p3_private_damping64_replay.py': '4c0d22e8e05efae7d618123f419947043172e5a8e0e58f5b322327c4338c4315',
    REPLAY / 'executed_sources/codex_replay_first_p3_stable_bb.py': '53cbd331ce8d65d387c2034d1fd0a955197bc47c44447d68fa326b906774de2e',
    REPLAY / 'executed_sources/codex_stable_bb_probe.py': PROBE_SHA,
}
STEPS = (1e-3, 3e-4, 1e-4, 3e-5)
NODE_COUNT = 20
KKT_TOLERANCE = 1e-6
ACTIVE_DISTANCE = 1e-10
OUTPUT_CAP = 32 << 20
PHASE_SECONDS = {
    'completion_identity_and_reconstruction': 2*3600,
    'baseline_and_boundary': 2*3600,
    'focused_directions': 4*3600,
    'remaining_informative_coordinates': 76*3600,
    'summary': 3600,
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def serial(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def logged_marginal(source):
    """Compile a separate method with exactly one logging-only AST addition."""
    original = next(n for n in ast.walk(ast.parse(source))
                    if isinstance(n, ast.FunctionDef) and n.name == 'marginal')
    changed = copy.deepcopy(original)
    loop = next(n for n in changed.body if isinstance(n, ast.For))
    position = next(i for i, n in enumerate(loop.body) if isinstance(n, ast.Assign)
                    and any(isinstance(t, ast.Name) and t.id == 'curvature' for t in n.targets))
    callback = ast.parse('__precision_log_mode(gene,mode,curvature,hstats(mode),len(rs))').body[0]
    loop.body.insert(position+1, callback)
    stripped = copy.deepcopy(changed)
    stripped_loop = next(n for n in stripped.body if isinstance(n, ast.For))
    removed = stripped_loop.body.pop(position+1)
    require(ast.dump(removed, include_attributes=False) == ast.dump(callback, include_attributes=False)
            and ast.dump(stripped, include_attributes=False) == ast.dump(original, include_attributes=False),
            'Logging changes numerical AST')
    require(not any(isinstance(n, ast.Name) and n.id == '__precision_log_mode' for n in ast.walk(original)),
            'Logging name collides with source')
    program = ast.fix_missing_locations(ast.Module(body=[changed], type_ignores=[]))
    return compile(program, 'diagnostic_logged_marginal_ast', 'exec'), {
        'original_marginal_ast_sha256': hashlib.sha256(ast.dump(original, include_attributes=False).encode()).hexdigest(),
        'stripped_logged_ast_identical': True,
        'logging_statement': ast.unparse(callback),
        'logged_ast_sha256': hashlib.sha256(ast.dump(changed, include_attributes=False).encode()).hexdigest(),
    }


def main():
    require(os.environ.get('SLURM_JOB_ID', '').isdigit()
            and int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) >= 1, 'Compute allocation required')
    output = REC / ('codex_p3_retained_fit_precision_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(exist_ok=False)
    private = output / 'executed_sources'
    private.mkdir()
    started = time.monotonic()
    phase = 'completion_identity_and_reconstruction'
    phase_started = started
    objective_calls = 0
    current_modes = []

    def disk_bytes():
        return sum(path.stat().st_size for path in output.rglob('*') if path.is_file())

    def save(name, value):
        require(disk_bytes() < OUTPUT_CAP, 'Diagnostic saved-output allowance exceeded')
        with (output / (name+'.json')).open('x') as f:
            json.dump(value, f, default=serial, indent=2, allow_nan=False)
            f.write('\n')
        require(disk_bytes() <= OUTPUT_CAP, 'Diagnostic saved-output allowance exceeded')

    def event(kind, **fields):
        with (output / 'progress.jsonl').open('a') as f:
            f.write(json.dumps(dict(event=kind, phase=phase, elapsed_seconds=time.monotonic()-started,
                phase_elapsed_seconds=time.monotonic()-phase_started, objective_calls=objective_calls,
                **fields), default=serial, allow_nan=False)+'\n')
        require(disk_bytes() <= OUTPUT_CAP, 'Diagnostic saved-output allowance exceeded')
        print(kind, phase, flush=True)

    def phase_timeout(signum, frame):
        raise TimeoutError('Diagnostic phase budget exhausted: '+phase)

    signal.signal(signal.SIGALRM, phase_timeout)

    def begin(name):
        nonlocal phase, phase_started
        phase, phase_started = name, time.monotonic()
        signal.setitimer(signal.ITIMER_REAL, PHASE_SECONDS[name])
        event('phase_started', budget_seconds=PHASE_SECONDS[name])

    def bounded_json(path):
        require(path.is_file() and path.stat().st_size <= 8 << 20, 'Missing/incomplete bounded receipt: '+str(path))
        return json.loads(path.read_text())

    try:
        begin(phase)
        # Read completion first: this helper never fits an incomplete/failed replay.
        require(not (OUTER / 'failure.json').exists(), 'Replay wrapper recorded failure')
        outer_summary = bounded_json(OUTER / 'summary.json')
        result_summary = bounded_json(REPLAY / 'summary.json')
        result = bounded_json(REPLAY / 'p3_000.json')
        generation = bounded_json(REPLAY / 'generation.json')
        require(outer_summary.get('fixed_gene_guard_passed')
                and outer_summary.get('realized_original_draw_identity_confirmed')
                and Path(outer_summary['replay_output']) == REPLAY
                and outer_summary.get('owner_spec_amended') is False
                and outer_summary.get('statistical_law_changed') is False,
                'Outer replay completion/identity is not admitted')
        require(result_summary.get('replayed_draws') == 1 and result_summary.get('first_draw_okay') is True
                and result.get('b') == 0 and result.get('fit') is not None,
                'Replay draw incomplete or failed; no retained fit to check')
        fit = result['fit']
        require(fit.get('converged') is True and fit['null'].get('converged') is True
                and np.isfinite(result['LR']), 'Replay retained fit/null or statistic failed')
        require(json.dumps(result_summary['replay_result'], sort_keys=True) == json.dumps(result, sort_keys=True),
                'Nested final-result copies differ')
        require(generation.get('completed') is True and generation.get('synthetic_rows') == 173
                and generation['count_vector_sha256'] == EXPECTED_COUNT, 'Final replay draw identity differs')
        for path, digest in GUARDS.items():
            require(sha(path) == digest, 'Frozen replay guard differs: '+str(path))
        frozen = dict(GUARDS)
        for path in (OUTER / 'summary.json', REPLAY / 'summary.json', REPLAY / 'p3_000.json',
                     REPLAY / 'generation.json', REPLAY / 'observed_null.json'):
            frozen[path] = sha(path)
        receipt = bounded_json(REPLAY / 'receipt.json')
        require(receipt['source_sha256'] == SOURCES and receipt['stable_probe_sha256'] == PROBE_SHA
                and receipt['wrapper_sha256'] == GUARDS[REPLAY / 'executed_sources/codex_replay_first_p3_stable_bb.py']
                and receipt['draw_index'] == 0 and receipt['replayed_draws'] == 1
                and Path(receipt['fixture']) == FIXTURE and receipt.get('actual_genotypes_read') is False
                and receipt.get('receiving_outcomes_read') is False, 'Replay source/modality receipt differs')
        require(receipt['python'] == platform.python_version() and receipt['numpy'] == np.__version__
                and receipt['scipy'] == scipy.__version__, 'Replay numerical environment differs')
        for name, digest in SOURCES.items():
            path = REPLAY / 'executed_sources' / name
            require(sha(path) == digest and sha(OUTER / 'executed_sources' / name) == digest,
                    'Private executed source differs: '+name)
            frozen[path] = digest
            (private / name).write_bytes(path.read_bytes())
        source = (private / 'likelihood_v1.py').read_text()
        original_source = (ROLE / 'executed_sources/likelihood_v1.py').read_text()
        require(sha(ROLE / 'executed_sources/likelihood_v1.py')
                == '41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6'
                and source == original_source.replace('for damping in 2.**-np.arange(30):',
                                                       'for damping in 2.**-np.arange(64):'),
                'Private likelihood is not the frozen isolated damping copy')
        for name, digest in receipt['fixture_sha256'].items():
            path = FIXTURE / name
            require(sha(path) == digest, 'Synthetic fixture differs: '+name)
            frozen[path] = digest
        for path in (Path(__file__), HERE / 'run_codex_p3_retained_fit_precision.sbatch',
                     REPLAY / 'executed_sources/codex_stable_bb_probe.py'):
            (private / path.name).write_bytes(path.read_bytes())
        sys.path.insert(0, str(private))
        require(not any(name in sys.modules for name in ('likelihood_v1', 'bfix_v1', 'design_v1', 'tests_v1')),
                'Numerical module imported before private isolation')
        loader = importlib.import_module('bfix_v1')
        likelihood = importlib.import_module('likelihood_v1')
        for name in ('likelihood_v1', 'bfix_v1', 'design_v1', 'tests_v1'):
            require(Path(sys.modules[name].__file__).parent == private, 'Model import escaped immutable private copy')
        spec = importlib.util.spec_from_file_location('retained_point_stable_bb', private / 'codex_stable_bb_probe.py')
        stable = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = stable
        spec.loader.exec_module(stable)
        def replacement(k, n, eta, linear_r, need_hessian=True):
            value, score, hessian = stable.stable_bb_terms(k, n, eta, linear_r)
            return value, score, hessian if need_hessian else None
        likelihood.bb_terms = replacement
        design, settings = loader.load()
        source_model = design.model()
        require(len(source_model.rows) == 173 and settings['B'] == 50, 'Original synthetic dimensions differ')
        require([(r['individual_id'], r['gene_id']) for r in loader.table('rows.tsv')]
                == [(r.individual, settings['gene_order'][r.gene]) for r in source_model.rows], 'Synthetic row order differs')
        null_receipt = bounded_json(REPLAY / 'observed_null.json')
        require(null_receipt['converged'], 'Replay generation-null fit failed')
        generation_p = np.asarray(null_receipt['parameters'], dtype=float)
        uniforms = np.load(FIXTURE / 'draws/p3_uniforms.npy', allow_pickle=False, mmap_mode='r')
        require(uniforms.shape == (50, 173, 2), 'Frozen synthetic draw shape differs')
        rows = []
        for row, (u1, u2) in zip(source_model.rows, uniforms[0]):
            eta, linear_r, _, _ = source_model.predictors(row, generation_p)
            rows.append(likelihood.conditional_draw(row, eta, linear_r, u1, u2))
        actual_count = hashlib.sha256(np.asarray([r.a for r in rows], dtype='<i8').tobytes()).hexdigest()
        require(actual_count == EXPECTED_COUNT, 'Exact replay generation does not reproduce original count vector')
        model = likelihood.Model(rows)
        point = np.r_[np.asarray(fit['parameters'], dtype=float), float(fit['v'])]
        saved_gradient = np.asarray(fit['gradient'], dtype=float)
        require(point.shape == saved_gradient.shape == (model.size+1,)
                and np.isfinite(point).all() and np.isfinite(saved_gradient).all() and point[-1] >= 0,
                'Retained parameter/gradient shape or domain differs')
        require(fit['kept_start'] in ('warm', 'cold', 'boundary_restart')
                and set(fit['start_logliks']) == {'warm', 'cold'}
                and set(fit['start_converged']) == {'warm', 'cold'}, 'Retained original-start record differs')
        if fit['kept_start'] != 'boundary_restart':
            expected_start = 'cold' if fit['start_logliks']['cold']-fit['start_logliks']['warm'] >= 1e-6 else 'warm'
            require(fit['kept_start'] == expected_start, 'Original retained-start tie convention differs')
        informative = model.informative()
        null_point = np.asarray(fit['null']['parameters'])
        require(np.array_equal(point[:-1][~informative], null_point[~informative]), 'Fixed uninformative coordinates changed')
        bounds = [b if informative[j] else (point[j], point[j]) for j, b in enumerate(model.bounds())]+[(0., None)]
        for value, (lo, hi) in zip(point, bounds):
            require((lo is None or value >= lo) and (hi is None or value <= hi), 'Saved point lies outside original bounds')
        names = (['alpha_'+str(j) for j in range(model.G)]+['r_'+str(j) for j in range(model.G)]
                 +['kappa']+['delta_'+str(j) for j in range(model.Z)]+['omega_S', 'rho_S', 'rho_b', 'rho_d', 'v'])
        code, ast_metadata = logged_marginal(source)
        def log_mode(gene, mode, curvature, stats, gene_rows):
            require(np.isfinite(stats).all(), 'Nonfinite logged inner mode')
            current_modes.append(dict(gene=int(gene), rows=int(gene_rows), mode=float(mode),
                log_integrand=float(stats[0]), score=float(stats[1]), curvature=float(curvature)))
        namespace = {**likelihood.__dict__, '__precision_log_mode': log_mode}
        exec(code, namespace)
        logged = namespace['marginal']
        save('protocol', dict(replay=str(REPLAY), frozen_inputs={str(p): d for p, d in frozen.items()},
            count_vector_sha256=actual_count, source_sha256=sha(Path(__file__)), parameter_names=names,
            retained_point=point, retained_start=fit['kept_start'], start_logliks=fit['start_logliks'],
            start_converged=fit['start_converged'], bounds=bounds, informative=informative,
            nodes=NODE_COUNT, absolute_steps=STEPS, original_projected_kkt_threshold=KKT_TOLERANCE,
            original_active_bound_distance=ACTIVE_DISTANCE, logging=ast_metadata,
            phase_budgets_seconds=PHASE_SECONDS, output_cap_bytes=OUTPUT_CAP,
            model_fitted=False, true_integral_comparison=False, log_variance_gradient=False,
            production_source_changed=False, nodes_changed=False,
            environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
            executed_source_bytes=sum(p.stat().st_size for p in private.iterdir()),
            frozen_input_bytes={str(p):p.stat().st_size for p in frozen}))

        def evaluate(q, label):
            nonlocal objective_calls
            objective_calls += 1
            current_modes.clear()
            tick = time.monotonic()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                value, gradient = logged(model, q[:-1], q[-1], nodes=NODE_COUNT)
            require(np.isfinite(value) and np.isfinite(gradient).all(), 'Nonfinite diagnostic objective')
            record = dict(label=label, loglik=float(value), gradient=gradient,
                modes=copy.deepcopy(current_modes), seconds=time.monotonic()-tick,
                warnings=[str(w.message) for w in caught])
            with (output / 'objective_calls.jsonl').open('a') as f:
                f.write(json.dumps(record, default=serial, allow_nan=False)+'\n')
            require(disk_bytes() <= OUTPUT_CAP, 'Objective log exceeds diagnostic output allowance')
            return record

        def projected(gradient):
            result = np.asarray(gradient, dtype=float).copy()
            for j, (lo, hi) in enumerate(bounds):
                if ((lo is not None and lo == hi)
                    or (lo is not None and point[j] <= lo+ACTIVE_DISTANCE and result[j] < 0)
                    or (hi is not None and point[j] >= hi-ACTIVE_DISTANCE and result[j] > 0)):
                    result[j] = 0.
            return result

        begin('baseline_and_boundary')
        tick = time.monotonic()
        unlogged_value, unlogged_gradient = model.marginal(point[:-1], point[-1], nodes=NODE_COUNT)
        unlogged_seconds = time.monotonic()-tick
        baseline = evaluate(point, 'retained_point')
        require(baseline['loglik'] == unlogged_value
                and np.array_equal(baseline['gradient'], unlogged_gradient), 'Logging changed computed objective/score')
        require(abs(baseline['loglik']-fit['loglik']) <= 1e-6
                and np.max(np.abs(baseline['gradient']-saved_gradient)) <= 1e-6,
                'Saved retained point is not reproduced by exact replay source')
        supplied_projected = projected(baseline['gradient'])
        baseline.update(unlogged_seconds=unlogged_seconds, logged_unlogged_bitwise_equal=True,
            saved_retained_convergence=fit['converged'],
            saved_loglik_difference=baseline['loglik']-fit['loglik'],
            saved_gradient_difference=baseline['gradient']-saved_gradient,
            original_projected_gradient=supplied_projected,
            original_projected_maximum=float(np.max(np.abs(supplied_projected))),
            original_kkt_condition=bool(np.max(np.abs(supplied_projected)) < KKT_TOLERANCE))
        if point[-1] == 0.:
            ll, g, _ = model.evaluate(point[:-1], need_hessian=False)
            terms = []
            for gene in range(model.G):
                su, hu = 0., 0.
                for row in model.rows:
                    if row.gene != gene:
                        continue
                    eta, linear_r, _, _ = model.predictors(row, point[:-1], need_hessian=False)
                    _, score, hessian = likelihood.conditional_terms(row, eta, linear_r)
                    su += score[0]*row.stage
                    hu += hessian[0, 0]*row.stage**2
                terms.append(dict(gene=gene, u_score=su, u_curvature=hu, right_variance_score=.5*(su**2+hu)))
            right_score = float(sum(t['right_variance_score'] for t in terms))
            exact_gradient = np.r_[g, right_score]
            require(ll == baseline['loglik'] and np.array_equal(exact_gradient, baseline['gradient']),
                    'Exact zero-variance formula differs from source objective')
            save('zero_variance_boundary', dict(loglik=ll, exact_gradient=exact_gradient, gene_terms=terms,
                right_variance_score=right_score, variance_nonpositive_kkt=right_score <= 0.,
                original_variance_tolerance_condition=right_score <= KKT_TOLERANCE,
                note='Exact source right derivative; no log-v derivative or positive-v approximation replaces it'))
        def fixed_coordinate(j):
            lo, hi = bounds[j]
            return lo is not None and lo == hi
        directions = [2*model.G, 2*model.G+1+model.Z, model.size]
        directions = [j for j in directions if not fixed_coordinate(j)]
        remaining = [j for j in range(len(point)) if j not in directions and not fixed_coordinate(j)]
        upper_calls = 2+2*len(STEPS)*(len(directions)+len(remaining))
        save('baseline', baseline)
        event('baseline_complete', observed_unlogged_objective_seconds=unlogged_seconds,
            observed_logged_objective_seconds=baseline['seconds'], max_planned_objective_calls=upper_calls,
            conditional_runtime_estimate_seconds=upper_calls*baseline['seconds'],
            runtime_estimate_limit='Estimate assumes later mode searches cost the measured baseline; not an ETA or guarantee')

        def direction(j):
            lo, hi = bounds[j]
            records = []
            for nominal_step in STEPS:
                h = nominal_step
                if ((lo is None or point[j]-h >= lo) and (hi is None or point[j]+h <= hi)):
                    offsets, weights, scheme = (-h, h), (-1., 1.), 'central'
                    multiplier = 0.
                elif hi is None or hi-point[j] >= point[j]-(lo if lo is not None else -np.inf):
                    if hi is not None:
                        h = min(h, (hi-point[j])/2)
                    offsets, weights, scheme = (h, 2*h), (4., -1.), 'forward_three_point'
                    multiplier = -3.
                else:
                    if lo is not None:
                        h = min(h, (point[j]-lo)/2)
                    offsets, weights, scheme = (-h, -2*h), (-4., 1.), 'backward_three_point'
                    multiplier = 3.
                require(h > 0 and np.isfinite(h), 'No feasible nonzero finite-difference step')
                probes = []
                for k, offset in enumerate(offsets):
                    q = point.copy()
                    q[j] += offset
                    require((lo is None or q[j] >= lo) and (hi is None or q[j] <= hi), 'Probe crossed original bound')
                    try:
                        probes.append(evaluate(q, names[j]+'_h'+str(h)+'_probe'+str(k)))
                    except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                        if j != model.size or point[-1] != 0.:
                            raise
                        result = dict(coordinate=j, name=names[j], saved_parameter=0., original_bounds=[lo, hi],
                            supplied_score=baseline['gradient'][j], primary_boundary_score=baseline['gradient'][j],
                            refinements=records, positive_v_probes_are_secondary=True,
                            positive_v_secondary_probe_failure=dict(step=h, offset=offset, traceback=traceback.format_exc()),
                            finest_step_derivative=baseline['gradient'][j],
                            note='Positive-v secondary probe failed; exact source right derivative retained, no fallback or fit verdict')
                        save('direction_'+names[j], result)
                        event('secondary_variance_probe_unavailable', coordinate=j, actual_step=h)
                        return result
                derivative = (multiplier*baseline['loglik']+sum(w*q['loglik'] for w, q in zip(weights, probes)))/(2*h)
                rounding_scale = np.finfo(float).eps*(abs(multiplier*baseline['loglik'])
                    +sum(abs(w*q['loglik']) for w, q in zip(weights, probes)))/(2*h)
                g = baseline['gradient'].copy()
                g[j] = derivative
                records.append(dict(nominal_step=nominal_step, actual_step=h, scheme=scheme,
                    derivative=derivative, minus_supplied_score=derivative-baseline['gradient'][j],
                    projected_derivative=float(projected(g)[j]), floating_subtraction_scale_estimate=rounding_scale,
                    plus_or_minus_probe_records=probes,
                    mode_displacements=[[dict(gene=a['gene'], displacement=b['mode']-a['mode'])
                        for a, b in zip(baseline['modes'], probe['modes'])] for probe in probes],
                    mode_branch_continuity_certified=False))
                event('step_completed', coordinate=j, name=names[j], actual_step=h, scheme=scheme)
            changes = [records[k]['derivative']-records[k-1]['derivative'] for k in range(1, len(records))]
            uncertainty = abs(changes[-1])+records[-1]['floating_subtraction_scale_estimate']+records[-2]['floating_subtraction_scale_estimate']
            # This is a transparent refinement estimate, not a rigorous derivative bound.
            last = records[-1]['derivative']
            g = baseline['gradient'].copy()
            low, high = g.copy(), g.copy()
            low[j], high[j] = last-uncertainty, last+uncertainty
            projected_interval = sorted([float(projected(low)[j]), float(projected(high)[j])])
            result = dict(coordinate=j, name=names[j], saved_parameter=point[j], original_bounds=[lo, hi],
                supplied_score=baseline['gradient'][j], refinements=records,
                successive_derivative_changes=changes, finest_step_derivative=last,
                last_adjacent_refinement_plus_rounding_uncertainty_estimate=uncertainty,
                estimated_projected_interval=projected_interval,
                note='Smallest step is not automatically most accurate; arithmetic errors/branch continuity uncertified')
            if j == model.size and point[-1] == 0.:
                result['primary_boundary_score'] = baseline['gradient'][j]
                result['positive_v_probes_are_secondary'] = True
            save('direction_'+names[j], result)
            return result

        begin('focused_directions')
        results = [direction(j) for j in directions]
        # A whole projected-gradient assessment needs every informative coordinate,
        # including active bounds whose derivative sign might change under refinement.
        begin('remaining_informative_coordinates')
        results.extend(direction(j) for j in remaining)
        begin('summary')
        require({r['coordinate'] for r in results}
                == {j for j in range(len(point)) if not fixed_coordinate(j)}, 'Incomplete informative-coordinate check')
        fd_gradient = baseline['gradient'].copy()
        for r in results:
            fd_gradient[r['coordinate']] = r['finest_step_derivative']
        if point[-1] == 0.:
            fd_gradient[-1] = baseline['gradient'][-1] # Exact primary right-variance score.
        fd_projected = projected(fd_gradient)
        for path, digest in frozen.items():
            require(sha(path) == digest, 'Frozen replay/fixture changed during diagnostic: '+str(path))
        save('summary', dict(retained_start=fit['kept_start'], retained_point=point,
            full_informative_coordinate_check_completed=True, supplied_projected_gradient=supplied_projected,
            original_supplied_kkt_condition=baseline['original_kkt_condition'],
            finest_step_diagnostic_gradient=fd_gradient, finest_step_diagnostic_projected_gradient=fd_projected,
            finest_step_diagnostic_projected_maximum=np.max(np.abs(fd_projected)),
            original_kkt_threshold=KKT_TOLERANCE, coordinate_summaries=[{k:v for k,v in r.items() if k != 'refinements'} for r in results],
            zero_variance_exact_primary=point[-1] == 0., objective_calls=objective_calls,
            elapsed_seconds=time.monotonic()-started, saved_output_bytes_before_summary=disk_bytes(),
            model_fitted=False, nodes_changed=False, true_integral_comparison=False,
            fit_declared_bad=False, automatic_repair=False, statistical_calibration_established=False,
            limitations='Refinement/rounding estimates are not certified bounds. Inspect all steps and logged modes; a mode-switch secant is not local derivative evidence. This examines the finite 20-node objective, not the true integral. No inference about wrong parameters from abnormal termination alone.'))
        event('phase_completed', actual_saved_output_bytes=disk_bytes())
    except BaseException:
        signal.setitimer(signal.ITIMER_REAL, 0)
        save('failure', dict(phase=phase, objective_calls=objective_calls,
            elapsed_seconds=time.monotonic()-started, traceback=traceback.format_exc(),
            completed_fit_precision_not_established=True, fit_declared_bad=False,
            model_fitted=False, production_source_changed=False))
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == '__main__':
    main()
