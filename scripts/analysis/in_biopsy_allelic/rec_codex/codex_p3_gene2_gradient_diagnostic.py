"""Compute-only fixed synthetic marginal-score check; no fit or fallback."""
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys
import time
import traceback
import warnings

import numpy as np
import scipy
from scipy.integrate import quad_vec
from scipy.optimize import brentq
from scipy.special import logsumexp
from scipy.stats import norm


ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
DIAGNOSTIC = REC / 'p3_first_draw_diagnostic_21998696'
MODE = REC / 'codex_p3_gene2_integrand_21998880'
ROLE = REC / 'bfix_p3_21998519'
FIXTURE = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
HERE = Path(__file__).resolve().parent
PROBE = MODE / 'executed_sources/codex_stable_bb_probe.py'
GUARDS = {
    DIAGNOSTIC / 'summary.json': 'b68b8fbd37264d6e93428e1368fce5518c07ea743e23b5fd5580baadff448bd0',
    DIAGNOSTIC / 'receipt.json': '3b313579039c403094bcea2464f4cc23c2ee0b11dec22189047abc5f6e4f3630',
    MODE / 'summary.json': '2d51e8d2eca5db5491ba379a51616c2dec25cd4711092606c892d724b8154792',
    MODE / 'protocol.json': '72481917c5adf771340624584d43f8d1a79eb6b431c6f3fe7ba3c8650871e256',
    MODE / 'generation.json': '23af2ebf77c6c0d52b6a1e59f8681bb0eb202f0db5172e2637f10c73682ed92c',
    PROBE: '723781113c17061b9a28341a91149bbd4ab71f195121a67f88d8cbc710a21488',
}
INPUTS = ('globals.json', 'individuals.tsv', 'rows.tsv', 'cohorts.tsv', 'draws/p3_uniforms.npy')
G = 2
BOUNDS = (16., 18.)
FD_STEPS = (1e-3, 3e-4, 1e-4, 3e-5)
NODES = (20, 40)
OUTPUT_CAP = 16 << 20


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
    require(os.environ.get('SLURM_JOB_ID', '').isdigit()
            and int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) >= 1, 'Compute allocation required')
    output = REC / ('codex_p3_gene2_gradient_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(exist_ok=False)
    private = output / 'executed_sources'
    private.mkdir()
    started = time.monotonic()
    phase = 'source_and_fixture_guards'
    evaluations = 0

    def saved_bytes():
        return sum(p.stat().st_size for p in output.rglob('*') if p.is_file())

    def save(name, value):
        require(saved_bytes() < OUTPUT_CAP, 'Diagnostic output cap exceeded')
        with (output / (name+'.json')).open('x') as f:
            json.dump(value, f, default=serial, indent=2, allow_nan=False)
            f.write('\n')
        require(saved_bytes() <= OUTPUT_CAP, 'Diagnostic output cap exceeded')

    def event(kind, **fields):
        with (output / 'progress.jsonl').open('a') as f:
            f.write(json.dumps(dict(event=kind, phase=phase, elapsed_seconds=time.monotonic()-started,
                                    likelihood_evaluations=evaluations, **fields), default=serial, allow_nan=False)+'\n')
        print(kind, phase, flush=True)

    try:
        event('phase_started')
        for path, digest in GUARDS.items():
            require(sha(path) == digest, 'Frozen diagnostic/mode source changed: '+str(path))
        receipt = json.loads((DIAGNOSTIC / 'receipt.json').read_text())
        failed = json.loads((DIAGNOSTIC / 'summary.json').read_text())
        modes = json.loads((MODE / 'summary.json').read_text())
        generation = json.loads((MODE / 'generation.json').read_text())
        mode_protocol = json.loads((MODE / 'protocol.json').read_text())
        frame = next(f for f in failed['archived_exception_frames'] if f['method'] == 'marginal')['locals']
        require(failed['draw'] == 0 and failed['phase'] == 'fit_marginal' and frame['gene'] == G
                and frame['mode'] == 0. and modes['fixed_gene'] == G, 'Wrong fixed synthetic call')
        require(sha(ROLE / 'receipt.json') == receipt['role_receipt_sha256']
                and sha(ROLE / 'observed_null.json') == receipt['observed_null_sha256'], 'Generation source differs')
        for name, digest in receipt['source_sha256'].items():
            path = ROLE / 'executed_sources' / name
            require(sha(path) == digest, 'Original executed source differs: '+name)
            require(sha(MODE / 'executed_sources' / name) == digest, 'Mode-check executed source differs: '+name)
            (private / name).write_bytes(path.read_bytes())
        require(sha(MODE / 'executed_sources/codex_p3_gene2_integrand_diagnostic.py')
                == mode_protocol['source_sha256'], 'Mode-check helper differs')
        for name in INPUTS:
            require(sha(FIXTURE / name) == receipt['fixture_sha256'][name], 'Synthetic fixture differs: '+name)
        for path in (Path(__file__), HERE / 'run_codex_p3_gene2_gradient_diagnostic.sbatch', PROBE):
            (private / path.name).write_bytes(path.read_bytes())
        with (DIAGNOSTIC / 'progress.jsonl').open() as f:
            environment = next(json.loads(line) for line in f if json.loads(line)['event'] == 'imports_verified')
        require(environment['python'] == platform.python_version() and environment['numpy'] == np.__version__
                and environment['scipy'] == scipy.__version__, 'Original numerical environment differs')
        sys.path.insert(0, str(private))
        require(not any(n in sys.modules for n in ('likelihood_v1', 'bfix_v1', 'design_v1', 'tests_v1')),
                'Numerical source already imported')
        loader = importlib.import_module('bfix_v1')
        likelihood = importlib.import_module('likelihood_v1')
        for name in ('likelihood_v1', 'bfix_v1', 'design_v1', 'tests_v1'):
            require(Path(sys.modules[name].__file__).parent == private, 'Import escaped private source')
        design, settings = loader.load()
        model = design.model()
        require(settings['B'] == 50 and len(model.rows) == 173, 'Original fixture dimensions differ')
        observed = json.loads((ROLE / 'observed_null.json').read_text())
        require(observed['converged'], 'Original synthetic generation null not converged')
        null_p = np.asarray(observed['parameters'])
        require([(r['individual_id'], r['gene_id']) for r in loader.table('rows.tsv')]
                == [(r.individual, settings['gene_order'][r.gene]) for r in model.rows], 'Original row order differs')
        uniforms = np.load(FIXTURE / 'draws/p3_uniforms.npy', allow_pickle=False, mmap_mode='r')
        require(uniforms.shape == (50, 173, 2), 'Original synthetic uniforms differ')
        drawn = []
        for row, (u1, u2) in zip(model.rows, uniforms[0]):
            eta, linear_r, _, _ = model.predictors(row, null_p)
            drawn.append(likelihood.conditional_draw(row, eta, linear_r, u1, u2))
        count_sha = hashlib.sha256(np.asarray([r.a for r in drawn], dtype='<i8').tobytes()).hexdigest()
        require(count_sha == generation['reconstructed_count_vector_sha256']
                and all(row.a in row.allowed for row in drawn), 'Synthetic count reconstruction differs')
        drawn_model = likelihood.Model(drawn)
        p = np.asarray(failed['latest_marginal']['parameters'], dtype=float)
        v = float(failed['latest_marginal']['v'])
        require(len(p) == drawn_model.size and v == frame['v'] == modes['fixed_v'] and v > 0,
                'Fixed failed-call parameters differ')
        rows = [r for r in drawn if r.gene == G]
        require(len(rows) == generation['gene_rows'] == 24, 'Fixed gene rows differ')
        spec = importlib.util.spec_from_file_location('private_gradient_bb', private / PROBE.name)
        stable = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = stable
        spec.loader.exec_module(stable)
        def replacement(k, n, eta, linear_r, need_hessian=True):
            ll, score, hessian = stable.stable_bb_terms(k, n, eta, linear_r)
            return ll, score, hessian if need_hessian else None
        likelihood.bb_terms = replacement
        extrema = modes['bracketed_extrema']
        dominant = max((e for e in extrema if e['curvature'] < 0), key=lambda e: e['log_integrand'])
        center, scale = dominant['u'], (-dominant['curvature'])**-.5
        fixed_shift = dominant['log_integrand']-.5*np.log(2*np.pi*v)
        names = (['alpha_'+str(j) for j in range(model.G)]+['r_'+str(j) for j in range(model.G)]
                 +['kappa']+['delta_'+str(j) for j in range(model.Z)]+['omega_S', 'rho_S', 'rho_b', 'rho_d'])
        require(len(names) == len(p), 'Parameter labels differ')
        save('protocol', dict(gene=G, draw=0, fixed_p=p, fixed_v=v, parameter_names=names,
            guards={str(path): digest for path, digest in GUARDS.items()}, original_receipt=receipt,
            reconstruction_count_sha256=count_sha, source_sha256=sha(Path(__file__)),
            python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
            bounds=BOUNDS, nodes=NODES, finite_difference_absolute_steps=FD_STEPS,
            finite_difference_direction_rule='alpha_gene2/r_gene2/kappa/omega_S/rho_S plus largest remaining absolute direct p score; v separately',
            normalized_prior=True, prior_log_normalizer='-0.5*log(2*pi*v)',
            variance_score='(u*u-v)/(2*v*v)', output_cap_bytes=OUTPUT_CAP,
            fixed_central_root_bracket=[center-.25, center+.25],
            source_snapshot_bytes={str(path): path.stat().st_size for path in GUARDS},
            fixture_bytes={name:(FIXTURE/name).stat().st_size for name in INPUTS},
            executed_source_bytes=sum(path.stat().st_size for path in private.iterdir()),
            model_fitted=False, live_sources_changed=False, actual_genotypes_read=False,
            receiving_outcomes_read=False, new_integration_fallback=False))

        def stats(q, variance, u, need_hessian=False):
            nonlocal evaluations
            evaluations += 1
            ll, gradient, _ = drawn_model.evaluate(q, rows, u=u, need_hessian=False)
            h = ll-u*u/(2*variance)-.5*np.log(2*np.pi*variance)
            score = np.r_[gradient, (u*u-variance)/(2*variance*variance)]
            require(np.isfinite(h) and np.isfinite(score).all() and ll <= 1e-8, 'Invalid proper-PMF normalized integrand')
            if not need_hessian:
                return h, score
            uscore, curvature = -u/variance, -1/variance
            for row in rows:
                eta, linear_r, _, _ = drawn_model.predictors(row, q, u, need_hessian=False)
                _, g, hh = likelihood.conditional_terms(row, eta, linear_r)
                uscore += g[0]*row.stage
                curvature += hh[0, 0]*row.stage**2
            return h, score, uscore, curvature

        # Uniform p-score bound: BB log-shape scores lie within +/-n,
        # mixing and selection normalization each contribute at most n.
        # Predictor Jacobians do not depend on u for this frozen model.
        def score_bound(q):
            result = np.zeros(len(q))
            for row in rows:
                jac = drawn_model.predictors(row, q, 0., need_hessian=False)[2]
                result += 2*row.n*np.abs(jac).sum(axis=0)
            return result

        points = sorted(set([e['u'] for e in extrema]+[center+t*scale for t in (-8., -4., -1., 0., 1., 4., 8.)]))

        def integrate(q, variance, bound, rtol, atol, scores=True):
            def function(u):
                h, score = stats(q, variance, u)
                weight = np.exp(h-fixed_shift)
                return weight*np.r_[1., score] if scores else np.array([weight])
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                integral, error, info = quad_vec(function, -bound, bound, epsabs=atol, epsrel=rtol,
                    norm='max', points=[x for x in points if -bound < x < bound], limit=400,
                    quadrature='gk21', full_output=True)
            require(info.success and np.isfinite(integral).all() and integral[0] > error >= 0,
                    'Direct quadrature did not resolve positive integral: '+info.message)
            mass = integral[0]
            logmass = fixed_shift+np.log(mass)
            logmass_lower_estimate = fixed_shift+np.log(mass-error)
            logtail = np.log(2.)+norm.logsf(bound/np.sqrt(variance))
            ratio = np.exp(logtail-logmass_lower_estimate)
            result = dict(bound=bound, log_normalized_integral=logmass, scaled_mass=mass,
                norm_max_absolute_error_estimate=error, relative_mass_error_estimate=error/mass,
                logmass_quad_error_estimate=-np.log1p(-error/mass),
                log_absolute_prior_tail_bound=logtail, tail_to_bounded_mass_upper_bound=ratio,
                tail_ratio_uses_quad_lower_mass_estimate=True,
                logmass_tail_error_upper_bound=np.log1p(ratio), evaluations=info.neval,
                scipy_status=int(info.status), scipy_message=info.message,
                warnings=[str(w.message) for w in caught], epsabs=atol, epsrel=rtol)
            if scores:
                gradient = integral[1:]/mass
                component_error = error*(1+np.abs(gradient))/(mass-error)
                p_tail_error = ratio*(score_bound(q)+np.abs(gradient[:-1])+component_error[:-1])
                a = bound/np.sqrt(variance)
                # Integral of prior*(u²+v)/(2v²) over |u|>bound.
                log_v_tail = logsumexp([np.log(a)+norm.logpdf(a), np.log(2.)+norm.logsf(a)])-np.log(variance)
                v_tail_error = np.exp(log_v_tail-logmass_lower_estimate)+(abs(gradient[-1])+component_error[-1])*ratio
                result.update(score=gradient, component_quad_error_estimates=component_error,
                    score_tail_error_upper_bounds=np.r_[p_tail_error, v_tail_error],
                    uniform_p_score_bounds=score_bound(q), log_absolute_v_score_tail_bound=log_v_tail)
            return result

        def agh(q, variance, nodes, recenter=False):
            c, s = center, scale
            if recenter:
                left, right = center-.25, center+.25
                require(stats(q, variance, left, True)[2] > 0 and stats(q, variance, right, True)[2] < 0,
                        'Perturbed central-root bracket no longer identifies central maximum')
                c = brentq(lambda u: stats(q, variance, u, True)[2], left, right, xtol=1e-12, rtol=1e-14)
                curvature = stats(q, variance, c, True)[3]
                require(curvature < 0, 'Perturbed central extremum is not a maximum')
                s = (-curvature)**-.5
            x, weights = np.polynomial.hermite.hermgauss(nodes)
            evaluations_at_nodes = [stats(q, variance, float(u)) for u in c+np.sqrt(2)*s*x]
            terms = np.log(weights)+x*x+np.array([value[0] for value in evaluations_at_nodes])
            logden = logsumexp(terms)
            score = np.exp(terms-logden)@np.stack([value[1] for value in evaluations_at_nodes])
            return dict(nodes=nodes, center=c, scale=s, log_normalized_integral=np.log(np.sqrt(2)*s)+logden,
                        score=score, recentered=recenter)

        phase = 'direct_scores_and_quadrature_refinement'
        event('phase_started')
        coarse = integrate(p, v, BOUNDS[0], 1e-8, 1e-10)
        direct = integrate(p, v, BOUNDS[0], 1e-10, 1e-12)
        expanded = integrate(p, v, BOUNDS[1], 1e-10, 1e-12)
        require(abs(direct['log_normalized_integral']+.5*np.log(2*np.pi*v)-modes['log_bounded_integral']) < 1e-6,
                'Normalized direct mass does not reproduce completed mode check')
        zero = stats(p, v, 0., True)
        require(abs(zero[3]-modes['zero_comparison']['stable_zero'][2]) < 1e-6,
                'Fixed synthetic curvature differs')
        node_checks = [agh(p, v, n) for n in NODES]
        save('integral_and_scores', dict(coarse=coarse, direct=direct, expanded=expanded,
            score_change_tight_minus_coarse=direct['score']-coarse['score'],
            score_change_expanded_minus_direct=expanded['score']-direct['score'],
            quadrature=node_checks, modecheck_logmass_normalizer=-.5*np.log(2*np.pi*v)))
        event('phase_completed')
        phase = 'selected_direction_finite_differences'
        event('phase_started')
        selected = [G, model.G+G, 2*model.G, 2*model.G+1+model.Z, len(p)-3]
        remaining = [j for j in range(len(p)) if j not in selected]
        selected.append(max(remaining, key=lambda j: abs(direct['score'][j])))
        direction_results = []
        for coordinate in selected+[len(p)]:
            name = names[coordinate] if coordinate < len(p) else 'v'
            refinements = []
            for step in FD_STEPS:
                plus, minus = p.copy(), p.copy()
                vp, vm = v, v
                if coordinate < len(p):
                    plus[coordinate] += step
                    minus[coordinate] -= step
                else:
                    vp += step
                    vm -= step
                require(vm > 0, 'Variance finite difference crossed zero')
                ip = integrate(plus, vp, BOUNDS[0], 1e-10, 1e-12, scores=False)
                im = integrate(minus, vm, BOUNDS[0], 1e-10, 1e-12, scores=False)
                fd = (ip['log_normalized_integral']-im['log_normalized_integral'])/(2*step)
                fd_error = sum(t['logmass_quad_error_estimate']+t['logmass_tail_error_upper_bound']
                               for t in (ip, im))/(2*step)
                approximations = []
                for nodes, base in zip(NODES, node_checks):
                    qp, qm = agh(plus, vp, nodes, True), agh(minus, vm, nodes, True)
                    qfd = (qp['log_normalized_integral']-qm['log_normalized_integral'])/(2*step)
                    approximations.append(dict(nodes=nodes, finite_difference=qfd,
                        minus_direct_score=qfd-direct['score'][coordinate],
                        minus_weighted_agh_score=qfd-base['score'][coordinate], plus=qp, minus=qm))
                refinements.append(dict(step=step, direct_finite_difference=fd,
                    minus_direct_score=fd-direct['score'][coordinate],
                    finite_difference_quad_plus_tail_error_estimate=fd_error,
                    plus_logmass=ip, minus_logmass=im, quadrature_finite_differences=approximations))
            direction_results.append(dict(coordinate=coordinate, name=name,
                direct_score=direct['score'][coordinate], refinements=refinements))
            save('finite_difference_'+name, direction_results[-1])
            event('direction_completed', coordinate=coordinate, name=name)
        phase = 'summary'
        event('phase_started')
        summaries = [dict(nodes=q['nodes'], score_minus_direct=q['score']-direct['score'],
            max_absolute_score_difference=np.max(np.abs(q['score']-direct['score'])),
            logmass_minus_direct=q['log_normalized_integral']-direct['log_normalized_integral']) for q in node_checks]
        save('summary', dict(fixed_gene=G, fixed_v=v, normalized_prior=True, parameter_names=names+['v'],
            direct_score=direct['score'], direct_logmass=direct['log_normalized_integral'],
            direct_component_quad_error_estimates=direct['component_quad_error_estimates'],
            direct_score_tail_error_upper_bounds=direct['score_tail_error_upper_bounds'],
            score_change_tight_minus_coarse=direct['score']-coarse['score'],
            score_change_expanded_minus_direct=expanded['score']-direct['score'], quadrature_comparisons=summaries,
            selected_finite_difference_directions=[r['name'] for r in direction_results],
            finite_difference_results=direction_results, likelihood_evaluations=evaluations,
            elapsed_seconds=time.monotonic()-started, saved_output_bytes_before_summary=saved_bytes(),
            executed_source_bytes=sum(path.stat().st_size for path in private.iterdir()),
            model_fitted=False, live_sources_changed=False, new_integration_fallback=False,
            limits='QUADPACK errors are estimates, not proofs; step refinements measure finite-difference stability, not a statistical margin. Fixed gene/draw/point only; no calibration, all-gene accuracy or outer convergence claim.'))
        event('phase_completed', actual_saved_output_bytes=saved_bytes())
    except BaseException:
        save('failure', dict(phase=phase, likelihood_evaluations=evaluations,
            elapsed_seconds=time.monotonic()-started, traceback=traceback.format_exc()))
        raise


if __name__ == '__main__':
    main()
