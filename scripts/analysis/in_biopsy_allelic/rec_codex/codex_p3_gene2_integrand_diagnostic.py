"""Frozen synthetic first-draw gene2 integrand; no fit or model-law change."""
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys
import time
import warnings

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import logsumexp
from scipy.stats import norm


ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
DIAGNOSTIC = REC / 'p3_first_draw_diagnostic_21998696'
ROLE = REC / 'bfix_p3_21998519'
FIXTURE = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
HERE = Path(__file__).resolve().parent
PROBE_SHA = '723781113c17061b9a28341a91149bbd4ab71f195121a67f88d8cbc710a21488'
GUARDS = {
    DIAGNOSTIC / 'summary.json': 'b68b8fbd37264d6e93428e1368fce5518c07ea743e23b5fd5580baadff448bd0',
    DIAGNOSTIC / 'receipt.json': '3b313579039c403094bcea2464f4cc23c2ee0b11dec22189047abc5f6e4f3630',
    DIAGNOSTIC / 'executed_diagnostic.py': 'bd9544d4e0a74c859c059b4699d678d126310b156f6a3104fc9e533fff3f6f8f',
    HERE / 'codex_stable_bb_probe.py': PROBE_SHA,
}
INPUTS = ('globals.json', 'individuals.tsv', 'rows.tsv', 'cohorts.tsv', 'draws/p3_uniforms.npy')
G = 2
MIN_BOUND = 12.
BOUND_BUFFER = 2.
OPERATIONAL_BOUND_CAP = 32.
GRID_STEPS = (.01, .005)
FD_STEPS = (1e-3, 3e-4, 1e-4)
ROOT_XTOL = 1e-10
QUAD_ATOL = 1e-10
QUAD_RTOL = 1e-8


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


def main():
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Compute allocation required')
    start = time.monotonic()
    output = REC / ('codex_p3_gene2_integrand_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(exist_ok=False)
    private = output / 'executed_sources'
    private.mkdir()
    def save(name, value):
        (output / (name+'.json')).write_text(json.dumps(value, default=serial, indent=2)+'\n')
    for path, digest in GUARDS.items():
        require(sha(path) == digest, 'Frozen diagnostic/probe changed: '+str(path))
    receipt = json.loads((DIAGNOSTIC / 'receipt.json').read_text())
    report = json.loads((DIAGNOSTIC / 'summary.json').read_text())
    frame = next(f for f in report['archived_exception_frames'] if f['method'] == 'marginal')['locals']
    require(report['draw'] == 0 and report['phase'] == 'fit_marginal'
            and frame['gene'] == G and frame['mode'] == 0., 'Not the fixed failed call')
    require(sha(ROLE / 'receipt.json') == receipt['role_receipt_sha256'], 'Original role receipt differs')
    require(sha(ROLE / 'observed_null.json') == receipt['observed_null_sha256'], 'Original generation parameters differ')
    source = ROLE / 'executed_sources'
    for name, digest in receipt['source_sha256'].items():
        require(sha(source / name) == digest, 'Original executed source differs: '+name)
        (private / name).write_bytes((source / name).read_bytes())
    for name in INPUTS:
        require(sha(FIXTURE / name) == receipt['fixture_sha256'][name], 'Synthetic fixture differs: '+name)
    for path in (Path(__file__), HERE / 'codex_stable_bb_probe.py'):
        (private / path.name).write_bytes(path.read_bytes())
    # Original numerical environment is part of deterministic draw reconstruction.
    with (DIAGNOSTIC / 'progress.jsonl').open() as handle:
        environment = next(json.loads(line) for line in handle if json.loads(line)['event'] == 'imports_verified')
    require(environment['python'] == platform.python_version() and environment['numpy'] == np.__version__
            and environment['scipy'] == scipy.__version__, 'Cannot reproduce original draw under a different numerical environment')
    protocol = dict(gene=G, draw=0, fixed_call=report['latest_marginal'],
                    archived_curvature=frame['curvature'], guards={str(p): d for p, d in GUARDS.items()},
                    original_receipt=receipt, source_sha256=sha(Path(__file__)),
                    minimum_bound=MIN_BOUND, gaussian_envelope_buffer=BOUND_BUFFER,
                    operational_bound_cap=OPERATIONAL_BOUND_CAP, grid_steps=GRID_STEPS,
                    finite_difference_steps=FD_STEPS, root_xtol=ROOT_XTOL,
                    quad_epsabs_scaled_integrand=QUAD_ATOL, quad_epsrel=QUAD_RTOL,
                    quad_limit=200, nodes=20, model_fitted=False, law_changed=False,
                    actual_genotypes_read=False, receiving_outcomes_read=False,
                    note='Fixed synthetic scalar integrand diagnostic; no optimizer/model/test-law amendment')
    save('protocol', protocol)
    sys.path.insert(0, str(private))
    require(not any(name in sys.modules for name in ('likelihood_v1', 'bfix_v1', 'design_v1', 'tests_v1')),
            'Numerical module imported outside private source')
    loader = importlib.import_module('bfix_v1')
    likelihood = importlib.import_module('likelihood_v1')
    for name in ('likelihood_v1', 'bfix_v1', 'design_v1', 'tests_v1'):
        require(Path(sys.modules[name].__file__).parent == private, 'Import escaped private source')
    design, settings = loader.load()
    model = design.model()
    require(settings['B'] == 50 and len(model.rows) == 173, 'Original row/draw dimensions differ')
    observed = json.loads((ROLE / 'observed_null.json').read_text())
    require(observed['converged'], 'Original generation null not converged')
    null_p = np.asarray(observed['parameters'])
    require([(r['individual_id'], r['gene_id']) for r in loader.table('rows.tsv')]
            == [(r.individual, settings['gene_order'][r.gene]) for r in model.rows], 'Original row order differs')
    uniforms = np.load(FIXTURE / 'draws/p3_uniforms.npy', allow_pickle=False, mmap_mode='r')
    require(uniforms.shape == (50, 173, 2), 'Original uniforms differ')
    drawn = []
    # Reconstruct the originally executed generation, before replacing BB arithmetic.
    for row, (u1, u2) in zip(model.rows, uniforms[0]):
        eta, r, _, _ = model.predictors(row, null_p)
        drawn.append(likelihood.conditional_draw(row, eta, r, u1, u2))
    require(len(drawn) == 173 and all(row.a in row.allowed for row in drawn), 'Original conditional rows cannot be reproduced')
    drawn_model = likelihood.Model(drawn)
    p = np.asarray(report['latest_marginal']['parameters'])
    v = float(report['latest_marginal']['v'])
    require(len(p) == drawn_model.size and v == frame['v'] and v > 0, 'Frozen failed-call parameters differ')
    rows = [row for row in drawn if row.gene == G]
    require(rows, 'Frozen gene has no reproduced rows')
    save('generation', dict(reconstructed_from_original_frozen_recipe=True, synthetic_rows=173,
                           gene_rows=len(rows), archived_count_vector_available=False,
                           reconstructed_count_vector_sha256=hashlib.sha256(np.asarray([r.a for r in drawn], dtype='<i8').tobytes()).hexdigest()))

    def hstats(u):
        value, score, curvature = -u*u/(2*v), -u/v, -1/v
        for row in rows:
            eta, r, _, _ = drawn_model.predictors(row, p, u, need_hessian=False)
            ll, g, h = likelihood.conditional_terms(row, eta, r)
            value += ll
            score += g[0]*row.stage
            curvature += h[0, 0]*row.stage**2
        return np.array([value, score, curvature])

    old_zero = hstats(0.)
    require(np.isfinite(old_zero).all() and abs(old_zero[2]-frame['curvature']) < 1e-6,
            'Reconstructed rows/parameters do not reproduce archived curvature; stop diagnostic')
    stable_spec = importlib.util.spec_from_file_location('private_stable_bb', private / 'codex_stable_bb_probe.py')
    stable = importlib.util.module_from_spec(stable_spec)
    sys.modules[stable_spec.name] = stable
    stable_spec.loader.exec_module(stable)
    original_bb = likelihood.bb_terms
    def replacement(k, n, eta, linear_r, need_hessian=True):
        ll, g, h = stable.stable_bb_terms(k, n, eta, linear_r)
        return ll, g, h if need_hessian else None
    likelihood.bb_terms = replacement
    zero = hstats(0.)
    require(np.isfinite(zero).all() and zero[0] <= 1e-8, 'Invalid fixed proper-PMF log integrand')
    # Since conditional row probabilities <=1, h(u)<=-u²/(2v).
    # Outside this radius the integrand cannot exceed h(0).
    envelope_radius = np.sqrt(2*v*abs(zero[0]))
    requested_bound = max(MIN_BOUND, float(np.ceil(envelope_radius+BOUND_BUFFER)))
    bound = min(requested_bound, OPERATIONAL_BOUND_CAP)
    extent_verified = requested_bound <= OPERATIONAL_BOUND_CAP
    save('extent', dict(envelope_radius=envelope_radius, fixed_buffer=BOUND_BUFFER,
                        requested_bound=requested_bound, operational_bound=bound,
                        global_maximum_search_extent_covered=extent_verified,
                        limitation=None if extent_verified else 'UNVERIFIED: operational cap prevents envelope-plus-buffer coverage'))

    def evaluate_representation(u, representation, check_evaluate=False):
        previous = likelihood.bb_terms
        likelihood.bb_terms = representation
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
                    value = hstats(u)
                    direct = drawn_model.evaluate(p, rows, u=u, need_hessian=False)[0]-u*u/(2*v) if check_evaluate else None
            result = dict(value=value, finite=bool(np.isfinite(value).all()),
                          warnings=sorted({str(w.message) for w in caught}))
            if check_evaluate:
                result.update(evaluate_log_integrand=direct, absolute_evaluate_difference=float(abs(value[0]-direct)))
                require(np.isfinite(direct) and abs(value[0]-direct) <= 1e-8,
                        'Scalar hstats differs from original model gene likelihood')
            return result
        except Exception as error:
            return dict(value=None, finite=False, exception_type=type(error).__name__, exception=str(error))
        finally:
            likelihood.bb_terms = previous

    # Exactly the thirty proposals that the original rule would make at u0.
    # The terminal mode0 record does not reconstruct prior iteration history.
    original_step = old_zero[1]/max(-old_zero[2], 1e-8)
    damping_trials = []
    for damping in 2.**-np.arange(30):
        u = float(damping*original_step)
        old_trial = evaluate_representation(u, original_bb)
        stable_trial = evaluate_representation(u, replacement)
        old_trial['would_accept_at_zero'] = bool(old_trial['value'] is not None and old_trial['value'][0] >= old_zero[0]-1e-12)
        stable_trial['would_accept_at_zero'] = bool(stable_trial['value'] is not None and stable_trial['value'][0] >= zero[0]-1e-12)
        damping_trials.append(dict(damping=damping, u=u, original=old_trial, stable=stable_trial))
    smallest_displacement = abs(damping_trials[-1]['u'])
    probe_base = min(1e-3, smallest_displacement/100.) if smallest_displacement > 0 else 1e-3
    probes = []
    for u in (0., *(sign*probe_base*factor for factor in (1., .1, .01, .001) for sign in (-1., 1.))):
        old_trial = evaluate_representation(u, original_bb, check_evaluate=True)
        stable_trial = evaluate_representation(u, replacement, check_evaluate=True)
        require(stable_trial['finite'] and 'absolute_evaluate_difference' in stable_trial,
                'Stable scalar likelihood consistency check failed')
        probes.append(dict(u=u, displacement_over_smallest_original_trial=abs(u)/smallest_displacement if smallest_displacement else None,
                           original=old_trial, stable=stable_trial,
                           stable_log_gain=float(stable_trial['value'][0]-zero[0])))
    save('damping_and_ascent_probes', dict(original_zero_step=original_step,
                                          smallest_original_trial_displacement=smallest_displacement,
                                          thirty_original_zero_proposals=damping_trials,
                                          small_probes=probes,
                                          history_limit='These are proposals at the archived terminal u0; prior mode-iteration history is unavailable'))
    fd = []
    for eps in FD_STEPS:
        plus, minus = hstats(eps), hstats(-eps)
        fd.append(dict(step=eps, score_from_logmass=(plus[0]-minus[0])/(2*eps),
                       curvature_from_logmass=(plus[0]-2*zero[0]+minus[0])/eps**2,
                       curvature_from_score=(plus[1]-minus[1])/(2*eps)))
    grid_passes = []
    # Cache shared scalar grid points so half-spacing adds only new points.
    cached = {}
    for spacing in GRID_STEPS:
        grid = np.linspace(-bound, bound, int(round(2*bound/spacing))+1)
        for u in grid:
            key = round(float(u), 8)
            if key not in cached:
                cached[key] = hstats(u)
        values = np.stack([cached[round(float(u), 8)] for u in grid])
        require(np.isfinite(values).all(), 'Stable frozen integrand not finite')
        sign_brackets = [(float(grid[i]), float(grid[i+1])) for i in range(len(grid)-1)
                         if values[i, 1]*values[i+1, 1] < 0]
        grid_passes.append(dict(spacing=spacing, points=len(grid), sign_bracket_count=len(sign_brackets),
                                sign_brackets=sign_brackets, exact_zero_score_points=grid[values[:, 1] == 0.].tolist()))
    original_records = [evaluate_representation(u, original_bb) for u in grid]
    original_values = np.asarray([r['value'] if r['value'] is not None else [np.nan]*3 for r in original_records])
    finite_original = np.isfinite(original_values).all(axis=1)
    grid_comparison = dict(finite_original_points=int(finite_original.sum()), total_points=len(grid),
                           original_numerical_exceptions=sum('exception_type' in r for r in original_records),
                           maximum_absolute_difference_on_finite_points=
                           np.max(np.abs(original_values[finite_original]-values[finite_original]), axis=0).tolist())
    with (output / 'fixed_gene_grid.tsv').open('x') as handle:
        handle.write('u\tlog_integrand\tscore\tcurvature\toriginal_log_integrand\toriginal_score\toriginal_curvature\n')
        for u, row, original_row in zip(grid, values, original_values):
            handle.write('\t'.join(f'{x:.17g}' for x in (u, *row, *original_row))+'\n')
    extrema = []
    for i in range(len(grid)-1):
        a, b = grid[i:i+2]
        sa, sb = values[i:i+2, 1]
        if sa*sb < 0:
            point = brentq(lambda u: hstats(u)[1], a, b, xtol=ROOT_XTOL, rtol=1e-14)
            value, score, curvature = hstats(point)
            extrema.append(dict(u=point, log_integrand=value, score=score, curvature=curvature,
                                bracket=[a, b], classification='maximum' if curvature < 0 else 'minimum' if curvature > 0 else 'flat'))
    for i in np.flatnonzero(values[:, 1] == 0.):
        if not any(abs(e['u']-grid[i]) <= ROOT_XTOL for e in extrema):
            extrema.append(dict(u=grid[i], log_integrand=values[i, 0], score=0., curvature=values[i, 2],
                                bracket=None, classification='maximum' if values[i, 2] < 0 else 'minimum' if values[i, 2] > 0 else 'flat'))
    extrema.sort(key=lambda e: e['u'])
    shift = max([float(values[:, 0].max())]+[float(e['log_integrand']) for e in extrema])
    integral = quad(lambda u: np.exp(hstats(u)[0]-shift), -bound, bound,
                    points=[e['u'] for e in extrema], epsabs=QUAD_ATOL, epsrel=QUAD_RTOL,
                    limit=200, full_output=True)
    scaled_integral, scaled_error, integration_info = integral[:3]
    require(scaled_integral > 0 and np.isfinite(scaled_integral), 'Invalid bounded integral')
    log_integral = shift+np.log(scaled_integral)
    # Each conditional row is a proper integer pmf <=1, so the product <=1.
    # This gives an absolute Gaussian-prior tail bound, not a fitted tail model.
    log_tail_bound = .5*np.log(2*np.pi*v)+np.log(2.)+norm.logsf(bound/np.sqrt(v))
    roots, weights = np.polynomial.hermite.hermgauss(20)
    agh = []
    for e in extrema:
        if e['curvature'] >= 0:
            continue
        scale = (-e['curvature'])**-.5
        u = e['u']+np.sqrt(2)*scale*roots
        terms = np.log(weights)+roots*roots+np.asarray([hstats(x)[0] for x in u])
        agh.append(dict(center=e['u'], log_integral=float(np.log(np.sqrt(2)*scale)+logsumexp(terms))))
    comparison = dict(old_zero=old_zero, stable_zero=zero, difference=zero-old_zero,
                      archived_curvature_difference=float(old_zero[2]-frame['curvature']))
    summary = dict(fixed_gene=G, fixed_v=v, zero_comparison=comparison, finite_differences=fd,
                   bounded_grid=[-bound, bound], grid_passes=grid_passes,
                   global_maximum_search_extent_covered=extent_verified, gaussian_envelope_radius=envelope_radius,
                   requested_bound=requested_bound, operational_bound_cap=OPERATIONAL_BOUND_CAP,
                   original_grid_comparison=grid_comparison,
                   bracketed_extrema=extrema,
                   log_bounded_integral=float(log_integral), scaled_quad_error=float(scaled_error),
                   relative_quad_error=float(scaled_error/scaled_integral), quad_evaluations=integration_info['neval'],
                   quad_warning=integral[3] if len(integral)>3 else None,
                   log_absolute_tail_upper_bound=float(log_tail_bound),
                   log_tail_upper_bound_relative_to_bounded_integral=float(log_tail_bound-log_integral),
                   alternative_20node_centers=agh, model_fitted=False, law_changed=False,
                   elapsed_seconds=time.monotonic()-start,
                   limits='Sign-changing grid brackets can miss tangential roots/subgrid extrema; bounded integral and Gaussian bound do not constitute a new fit or calibration')
    likelihood.bb_terms = original_bb
    for path, digest in GUARDS.items():
        require(sha(path) == digest, 'Diagnostic source changed during calculation')
    for name in INPUTS:
        require(sha(FIXTURE / name) == receipt['fixture_sha256'][name], 'Synthetic fixture changed')
    save('summary', summary)
    print(json.dumps(dict(completed=True, output=str(output), bracketed_extrema=len(extrema))), flush=True)


if __name__ == '__main__':
    main()
