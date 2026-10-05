#!/usr/bin/env python3
"""Compute-only fixed synthetic likelihood checks; no fits or calibration."""
import argparse
import ast
from dataclasses import asdict
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

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
CAMPAIGN = REC / 'p3_parallel_20261003T130323Z'
PROTOCOL = REC / 'repair_p3_case_protocol_20261005.json'
HERE = Path(__file__).resolve().parent
KERNEL = HERE / 'repair_p3_integral.py'
FAST_REFERENCE = HERE / 'repair_p3_reference_fast.py'
FAST_REFERENCE_SHA = '2eaaf549739f64322d394f0674aadba98f0f2ae6bedda1df56db8e9a8a3949fd'
BOUNDARY_MP = HERE / 'repair_p3_boundary_mp.py'
BOUNDARY_DECIMAL = HERE / 'repair_p3_boundary_decimal.py'
BOUNDARY_DECIMAL_FAST = HERE / 'repair_p3_boundary_decimal_fast.py'
PREDICTOR_CACHE = HERE / 'repair_p3_cached_predictors.py'
COMPOSITE_GL = HERE / 'repair_p3_composite_gl.py'
BASELINE_REUSE = HERE / 'repair_p3_baseline_reuse.py'
CAP = 32 * 1024**2
CASES = {'failed0gene2': 0, 'positive18': 18, 'boundary19': 19}
G = 2
STEPS = (1e-3, 3e-4, 1e-4)
FIXED_AXES = (2,10,16,17,24,25,28)  # immutable gene2 / r2 / kappa / delta0 / omegaS / rhoS / v
GUARDS = {
    PROTOCOL: '8e7ddf8b0a7bb18afd9692b53c3560f4258fdc373c4b822aee36aa94395eabe9',
    CAMPAIGN / 'campaign.json': 'bb3e82833476f67b0b11d2388e2dc00358eb987367b7a6a666d11938c6baca94',
    HERE / 'codex_stable_bb_probe.py': '723781113c17061b9a28341a91149bbd4ab71f195121a67f88d8cbc710a21488',
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def metadata():
    for path, expected in GUARDS.items():
        require(sha(path) == expected, 'Pinned small protocol/source differs')
    protocol = json.loads(PROTOCOL.read_text())
    guards = dict(GUARDS)
    guards.update({ROOT / p: digest for p, digest in protocol['case_input_sha256'].items()})
    for path, expected in guards.items():
        require(sha(path) == expected, 'Pinned case metadata differs')
    campaign = json.loads((CAMPAIGN / 'campaign.json').read_text())
    # Inspect hashes only for source code here; numerical arrays are compute-only.
    for path, expected in campaign['files'].items():
        p = Path(path)
        if p.suffix == '.py':
            require(p.parent == CAMPAIGN / 'executed_sources' and sha(p) == expected, 'Immutable synthetic source differs')
            guards[p] = expected
    ast.parse(KERNEL.read_text())
    return protocol, campaign, guards


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=tuple(CASES), required=True)
    parser.add_argument('--inspect', action='store_true')
    parser.add_argument('--expected-helper-sha')
    parser.add_argument('--expected-kernel-sha')
    parser.add_argument('--baseline-only', action='store_true')
    parser.add_argument('--fast-reference', action='store_true')
    parser.add_argument('--boundary-mp-sha')
    parser.add_argument('--boundary-decimal-sha')
    parser.add_argument('--boundary-decimal-fast-sha')
    parser.add_argument('--predictor-points-sha')
    parser.add_argument('--composite-gl-sha')
    parser.add_argument('--reuse-baseline-dir', type=Path)
    parser.add_argument('--reuse-helper-sha')
    parser.add_argument('--reuse-hashes-json', type=Path)
    parser.add_argument('--reuse-hashes-sha')
    parser.add_argument('--axes', nargs='+', type=int, choices=FIXED_AXES)
    args = parser.parse_args()
    require(args.axes is None or len(set(args.axes)) == len(args.axes), 'Derivative axes must be a unique subset of the fixed seven')
    protocol, campaign, guards = metadata()
    if args.fast_reference:
        require(sha(FAST_REFERENCE) == FAST_REFERENCE_SHA, 'Reviewed fast reference source differs')
        ast.parse(FAST_REFERENCE.read_text())
        guards[FAST_REFERENCE] = FAST_REFERENCE_SHA
    require(sum(bool(s) for s in (args.boundary_mp_sha, args.boundary_decimal_sha, args.boundary_decimal_fast_sha)) <= 1,
            'Choose one high-precision arithmetic implementation')
    boundary_path = BOUNDARY_DECIMAL_FAST if args.boundary_decimal_fast_sha else BOUNDARY_DECIMAL if args.boundary_decimal_sha else BOUNDARY_MP
    boundary_sha = args.boundary_decimal_fast_sha or args.boundary_decimal_sha or args.boundary_mp_sha
    if boundary_sha:
        require(args.case == 'boundary19', 'High-precision diagnostic is restricted to the saved boundary case')
        require(sha(boundary_path) == boundary_sha, 'Reviewed high-precision boundary source differs')
        ast.parse(boundary_path.read_text())
        guards[boundary_path] = boundary_sha
    if args.predictor_points_sha:
        require(not boundary_sha, 'Predictor comparison and boundary diagnostic are separate checks')
        require(sha(PREDICTOR_CACHE) == args.predictor_points_sha, 'Reviewed predictor-cache source differs')
        ast.parse(PREDICTOR_CACHE.read_text())
        guards[PREDICTOR_CACHE] = args.predictor_points_sha
    if args.composite_gl_sha:
        require(args.case != 'boundary19' and not args.predictor_points_sha,
                'Composite GL comparison requires a saved positive-variance point')
        require(sha(COMPOSITE_GL) == args.composite_gl_sha, 'Reviewed composite GL source differs')
        ast.parse(COMPOSITE_GL.read_text())
        guards[COMPOSITE_GL] = args.composite_gl_sha
    reuse_options = (args.reuse_baseline_dir, args.reuse_helper_sha, args.reuse_hashes_json, args.reuse_hashes_sha)
    require(not any(reuse_options) or all(reuse_options), 'All baseline reuse pins are required together')
    reuse_hashes = None
    if args.reuse_baseline_dir:
        require(args.case != 'boundary19' and args.composite_gl_sha and args.fast_reference and
                not args.baseline_only and not args.predictor_points_sha and not boundary_sha,
                'Baseline reuse is for positive composite GL derivative checks only')
        require(sha(BASELINE_REUSE) == args.reuse_helper_sha, 'Reviewed baseline reuse helper differs')
        ast.parse(BASELINE_REUSE.read_text())
        guards[BASELINE_REUSE] = args.reuse_helper_sha
        receipt_dir = args.reuse_baseline_dir.resolve(strict=True)
        require(receipt_dir.parent == REC.resolve() and receipt_dir.name.startswith('repair_p3_fixed_point_'+args.case+'_') and
                receipt_dir.name.rsplit('_',1)[-1].isdigit(), 'Baseline receipt escaped the synthetic case directory')
        require(args.reuse_hashes_json.resolve(strict=True).is_relative_to(REC.resolve()) and
                sha(args.reuse_hashes_json) == args.reuse_hashes_sha, 'Reviewed baseline hash manifest differs')
        reuse_hashes = json.loads(args.reuse_hashes_json.read_text())
        require(set(reuse_hashes) == {'protocol.json','candidate_baseline.json','independent_baseline.json','summary.json'},
                'Exactly four parent baseline file hashes required')
        guards[args.reuse_hashes_json] = args.reuse_hashes_sha
        for name,digest in reuse_hashes.items():
            path = receipt_dir/name
            require(not path.is_symlink() and path.resolve(strict=True).parent == receipt_dir and sha(path) == digest,
                    'Reviewed parent baseline file differs')
            guards[path] = digest
    if args.inspect:
        print(json.dumps(dict(case=args.case, draw=CASES[args.case], fixed_gene=G,
              numerical_execution=False, model_fitting=False, calibration=False,
              inherited_outer_starts_bounds_ties_failures_unchanged=True,
              approximation_revision_not_frozen20_node_agreement=True,
              baseline_only=args.baseline_only, requested_derivative_axes=args.axes or FIXED_AXES,
              fast_reference=args.fast_reference,
              high_precision_boundary=bool(boundary_sha),
              predictor_points_only=bool(args.predictor_points_sha),
              composite_gl=bool(args.composite_gl_sha),
              baseline_reuse=bool(args.reuse_baseline_dir),
              kernel_sha256=sha(KERNEL), metadata_and_source_checks=True)))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'), 'Individual compute job required')
    require(args.expected_helper_sha == sha(Path(__file__)) and args.expected_kernel_sha == sha(KERNEL), 'Reviewed candidate/checker source drift')
    output = REC / ('repair_p3_fixed_point_' + args.case + '_' + os.environ['SLURM_JOB_ID'])
    output.mkdir(mode=0o750)
    archive = output / 'executed_sources'
    archive.mkdir()
    began, phase = time.monotonic(), 'archive_before_import'

    def serial(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(type(value).__name__)

    def save(name, value):
        body = json.dumps(value, indent=2, default=serial, allow_nan=False) + '\n'
        require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file()) + len(body.encode()) <= CAP-65536,
                'Check output cap reached; failure reserve retained')
        final = output / (name + '.json')
        partial = output / (name + '.json.partial')
        with partial.open('x') as handle:
            handle.write(body)
            handle.flush()
            os.fsync(handle.fileno())
        require(not final.exists(), 'Refusing result overwrite')
        partial.rename(final)

    def event(name, **fields):
        record = dict(phase=phase, event=name, elapsed_seconds=time.monotonic()-began, **fields)
        with (output / 'progress.jsonl').open('a') as handle:
            handle.write(json.dumps(record, default=serial, allow_nan=False)+'\n')
        print(json.dumps(record, default=serial), flush=True)

    try:
        own = {Path(__file__):args.expected_helper_sha, KERNEL:args.expected_kernel_sha}
        for source in (*own, *(p for p in guards if p.suffix == '.py')):
            (archive / source.name).write_bytes(source.read_bytes())
            require(sha(archive / source.name) == (own.get(source) or guards[source]), 'Archive copy differs')
        if reuse_hashes is not None:
            parent_archive = output/'reused_parent_baseline'
            parent_archive.mkdir()
            for name,digest in reuse_hashes.items():
                source = receipt_dir/name
                target = parent_archive/name
                target.write_bytes(source.read_bytes())
                require(sha(target) == digest, 'Archived parent baseline differs')
            manifest_copy = parent_archive/'hash_manifest.json'
            manifest_copy.write_bytes(args.reuse_hashes_json.read_bytes())
            require(sha(manifest_copy) == args.reuse_hashes_sha, 'Archived parent hash manifest differs')
        fixture = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
        fixture_names = ('globals.json','individuals.tsv','rows.tsv','cohorts.tsv','draws/p3_uniforms.npy')
        for name in fixture_names:
            path = fixture / name
            require(path.resolve().is_relative_to(fixture.resolve()), 'Fixture input escaped synthetic root')
            expected = campaign['files'][str(path)]
            require(sha(path) == expected, 'Synthetic fixture bytes differ')
            guards[path] = expected
        phase = 'archived_imports_and_synthetic_draw'
        sys.path.insert(0, str(archive))
        require(not any(n in sys.modules for n in ('likelihood_v1','bfix_v1','design_v1','tests_v1')), 'Numerical original already imported')
        global np
        import numpy as np
        import scipy
        from scipy.integrate import quad
        from scipy.optimize import minimize_scalar
        from scipy.special import expit, gammaln, logsumexp
        from scipy.stats import binom, norm
        loader = importlib.import_module('bfix_v1')
        original = importlib.import_module('likelihood_v1')
        importlib.import_module('tests_v1')  # explicit before the archived-module assertion
        for name in ('likelihood_v1','bfix_v1','design_v1','tests_v1'):
            require(Path(sys.modules[name].__file__).resolve().parent == archive.resolve() and
                    sha(Path(sys.modules[name].__file__)) == campaign['files'][str(CAMPAIGN/'executed_sources'/ (name+'.py'))],
                    'Original import escaped guarded archive')
        # Candidate BB substitution never touches this second module's globals.
        reference_original = load_module('repair_p3_unmodified_reference',archive/'likelihood_v1.py')
        require(reference_original is not original and reference_original.bb_terms is not original.bb_terms,
                'Independent original numerical module not isolated')
        kernel = load_module('repair_p3_private_kernel', archive / KERNEL.name)
        stable_probe = load_module('repair_p3_reviewed_stable_probe', archive / 'codex_stable_bb_probe.py')
        fast_reference = load_module('repair_p3_private_fast_reference', archive / FAST_REFERENCE.name) if args.fast_reference else None
        boundary_mp = load_module('repair_p3_private_boundary_high_precision', archive / boundary_path.name) if boundary_sha else None
        if args.composite_gl_sha:
            composite_gl = load_module('repair_p3_private_composite_gl', archive / COMPOSITE_GL.name)
            def integrate_gene(*call_args, **call_kwargs):
                return composite_gl.integrate_gene(*call_args, **call_kwargs,
                                                  expected_kernel_sha=args.expected_kernel_sha,
                                                  kernel_path=archive / KERNEL.name)
        else:
            integrate_gene = kernel.integrate_gene
        versions = dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)
        require(versions == protocol['environment'], 'Original numerical environment versions differ')
        design, settings = loader.load()
        model = design.model()
        require(settings['B'] == 50 and len(model.rows) == 173, 'Synthetic design shape differs')
        observed = json.loads((CAMPAIGN/'baseline/observed_null.json').read_text())
        require(observed['converged'], 'Saved generating null not converged')
        null_p = np.asarray(observed['parameters'])
        uniforms = np.load(fixture/'draws/p3_uniforms.npy', mmap_mode='r', allow_pickle=False)
        require(uniforms.shape == (50,173,2), 'Original supplied uniforms differ')
        require([(r['individual_id'],r['gene_id']) for r in loader.table('rows.tsv')] ==
                [(r.individual,settings['gene_order'][r.gene]) for r in model.rows], 'Original uniform-to-row order differs')
        draw = CASES[args.case]
        rows = []
        for row, (uc, uk) in zip(model.rows, uniforms[draw]):
            eta, r, _, _ = model.predictors(row, null_p)
            rows.append(original.conditional_draw(row, eta, r, uc, uk))
        drawn = original.Model(rows, gene_count=model.G)
        count_sha = hashlib.sha256(np.asarray([r.a for r in rows], dtype='<i8').tobytes()).hexdigest()
        if draw == 0:
            require(count_sha == '5f4a1aa455cd09b0d5f31889021ec88b145dcbb50c8e5b600343ea9a3a992691', 'Known draw0 reconstruction differs')
            failed = json.loads((REC/'p3_first_draw_diagnostic_21998696/summary.json').read_text())
            p, v = np.asarray(failed['latest_marginal']['parameters']), float(failed['latest_marginal']['v'])
            require(v == 0.9234818604998617, 'Wrong original failed point')
        else:
            record = json.loads((CAMPAIGN/f'baseline/p3_{draw:03d}.json').read_text())
            require(record['b'] == draw and record['fit']['converged'], 'Wrong original successful point')
            p, v = np.asarray(record['fit']['parameters']), float(record['fit']['v'])
            require(v == (0.00711885901713953 if draw == 18 else 0.), 'Wrong original successful variance')
        gene_rows = [r for r in drawn.rows if r.gene == G]
        require(len(gene_rows) == 24 and len(p) == drawn.size, 'Fixed gene/coordinate shape differs')
        reference_model = reference_original.Model(rows,gene_count=model.G)
        boundary_original = reference_model.marginal(p, 0.) if v == 0 else None
        # Exact representation equivalence to previously reviewed source, no broad rerun.
        bb_errors = []
        for eta, r in ((-6.,-15.),(0.,-3.),(6.,5.),(-1000.,-15.),(1000.,5.)):
            k = np.arange(21)
            a, b = kernel.stable_bb_terms(k,20,eta,r), stable_probe.stable_bb_terms(k,20,eta,r)
            errors = [float(np.max(np.abs(x-y))) for x,y in zip(a,b)]
            require(max(errors) <= 1e-10, 'Candidate deviates from reviewed exact BB representation')
            bb_errors.append(dict(eta=eta, linear_r=r, mass_score_hessian_max_errors=errors))
        Candidate = kernel.candidate_class(original)
        require(reference_original.bb_terms is not original.bb_terms, 'Candidate patch leaked into reference globals')
        candidate = Candidate(rows, gene_count=model.G)
        unchanged = {name:getattr(Candidate,name) is getattr(original.Model,name) for name in
                     ('fit','fit_single','bounds','initial','cold_start','informative','fit_marginal','fit_marginal_single')}
        require(all(unchanged.values()), 'Outer recipe method was changed')
        save('protocol', dict(case=args.case, draw=draw, fixed_gene=G, fixed_gene_rows=24, fixed_p=p, fixed_v=v,
             count_vector_sha256=count_sha, source_guards={str(k):d for k,d in guards.items()}, own_source_hashes={str(k):d for k,d in own.items()},
             environment=versions, candidate_settings=asdict(kernel.Settings()), unchanged_outer_methods=unchanged,
             finite_difference_steps=STEPS, no_model_fit=True, no_calibration=True,
             baseline_only=args.baseline_only, required_combined_derivative_family_axes=FIXED_AXES,
             reference_arithmetic='independent cumulative integer rising products' if args.fast_reference else 'independent literal integer rising products',
             high_precision_boundary=bool(boundary_sha),
             predictor_points_only=bool(args.predictor_points_sha),
             composite_gl=bool(args.composite_gl_sha),
             baseline_reuse=bool(args.reuse_baseline_dir),
             reused_parent_directory=str(receipt_dir) if reuse_hashes is not None else None,
             reused_parent_file_hashes=reuse_hashes,
             requested_derivative_axes=args.axes or FIXED_AXES, all_declared_axes_remain_required=True,
             initial_positive_variance_scope='Only fixed gene2; no all8-gene integral or optimizer yet.',
             original_records_preserved=True, approximation_revision_not_frozen20_node_agreement=True,
             stochastic_draw_source='Original supplied uniform row only; no new Monte Carlo seed.',
             seed_for_future_synthetic_extensions=20261005, output_cap_bytes=CAP,
             original_completed_agreement_fixture_B50_not_calibration=True,
             independent_reference='Separately written full-support mixture probabilities + scalar QUADPACK, not reused conditional_terms or analytic scores.'))
        save('stable_BB_representation_equivalence', bb_errors)
        if args.predictor_points_sha:
            phase='predictor_cache_fixed_point_comparison'
            cache_helper = load_module('repair_p3_private_predictor_cache', archive / PREDICTOR_CACHE.name)
            report = cache_helper.compare_cached_predictors(candidate, p, gene_rows)
            save('predictor_cache_comparison', report)
            for path, digest in {**guards, **own}.items():
                require(sha(path) == digest, 'Source or input changed during predictor comparison')
            save('summary', dict(predictor_comparison_completed=True,
                                fixed_point_check_completed=False, derivative_checks_completed=False,
                                no_integration=True, no_model_fits=True, no_statistical_calibration=True,
                                case=args.case, draw=draw, fixed_gene=G,
                                all_declared_accuracy_checks_remain_required=True))
            event('predictor_comparison_completed')
            return

        # Independent scalar probability code: literal support, explicit orientation,
        # three genotype classes and selected-support normalizer. No original evaluator.
        def rise(logshape, count):
            return 0. if count == 0 else logshape + np.logaddexp(logshape, np.log(np.arange(1,count,dtype=float))).sum()
        def bb(k, n, eta, rr):
            if fast_reference is not None:
                return fast_reference.fast_bb(k, n, eta, rr)
            rho = 1e-6 + (1-2e-6)*expit(rr)
            logt = np.log((1-rho)/rho)
            la, lb = logt-np.logaddexp(0.,-eta), logt-np.logaddexp(0.,eta)
            return np.array([gammaln(n+1)-gammaln(a+1)-gammaln(n-a+1)+rise(la,int(a))+rise(lb,int(n-a))-rise(logt,n) for a in k])
        def fixed(k, row, mean):
            if row.phi == 0:
                return binom.logpmf(k,row.n,mean)
            t = (1-row.phi)/row.phi
            return np.array([gammaln(row.n+1)-gammaln(a+1)-gammaln(row.n-a+1)+rise(np.log(mean*t),int(a))+
                     rise(np.log((1-mean)*t),int(row.n-a))-rise(np.log(t),row.n) for a in k])
        fixed_parts = {}
        for row in gene_rows:
            lo = max(row.minor,(row.percent*row.n+99)//100)
            support = np.arange(lo,row.n-lo+1)
            require(row.n >= row.depth and row.a in support, 'Selected synthetic support differs')
            prior = np.array([2*row.p_ref*row.p_alt,row.p_ref**2+2*row.p_ref*row.p_x,row.p_alt**2+2*row.p_alt*row.p_x])
            if row.orientation == -1:
                prior = prior[[0,2,1]]
            with np.errstate(divide='ignore'):
                lp = np.log(prior)
            fixed_parts[id(row)] = (support, lp, fixed(support,row,row.e), fixed(support,row,1-row.e),int(row.a-lo))
        def independent_loglik(q,u):
            value = 0.
            for row in gene_rows:
                eta = q[row.gene]*(1+q[2*model.G]*row.stage+np.dot(q[2*model.G+1:2*model.G+1+model.Z],row.z))
                eta += u*row.stage-row.orientation*q[2*model.G+1+model.Z]*row.stage-row.orientation*row.omega
                rr = q[model.G+row.gene]+np.dot(q[-3:],(row.stage,row.background,row.duplication))
                support, prior, homref, homalt, at = fixed_parts[id(row)]
                parts = np.stack((bb(support,row.n,eta,rr),homref,homalt))+prior[:,None]
                value += logsumexp(parts[:,at])-logsumexp(parts)
            require(np.isfinite(value) and value <= 1e-8, 'Independent conditional PMF invalid')
            return float(value)

        def reference(q, variance, bound=24., precision=1.):
            require(variance > 0, 'Positive variance scalar reference required')
            sqrtv = np.sqrt(variance)
            cache = {}
            def h(z):
                if float(z) not in cache:
                    require(len(cache) < 100000, 'Independent scalar evaluation cap reached')
                    cache[float(z)] = independent_loglik(q,sqrtv*z)-z*z/2
                return cache[float(z)]
            grid = np.linspace(-bound,bound,int(2*bound)+1)
            values = np.array([h(z) for z in grid])
            points = set(map(float,grid[1:-1]))
            centers = [float(grid[values.argmax()])]
            for i in range(1,len(grid)-1):
                if values[i] >= max(values[i-1],values[i+1]):
                    peak = minimize_scalar(lambda z:-h(z),bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-10})
                    require(peak.success,'Independent scalar center search failed')
                    centers.append(float(peak.x))
                    points.update(peak.x+a for a in (-.5,-.1,0.,.1,.5) if -bound<peak.x+a<bound)
            shift = max(h(c) for c in centers)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                result = quad(lambda z:np.exp(h(z)-shift),-bound,bound,points=sorted(points),limit=4000,
                              epsrel=2e-12*precision,epsabs=2e-12*precision,full_output=True)
            require(len(result)==3 and result[0]>result[1]>=0,'Independent scalar integration unresolved')
            mass,error,info=result
            logmass=shift+np.log(mass)-.5*np.log(2*np.pi)
            lower=shift+np.log(mass-error)-.5*np.log(2*np.pi)
            tail=np.exp(np.log(2.)+norm.logsf(bound)-lower)
            return dict(logmass=float(logmass), estimated_quad_plus_mass_tail_log_error=float(-np.log1p(-error/mass)+np.log1p(tail)),
                        relative_prior_mass_tail_allowance=float(tail),z_bound=bound,evaluations=int(info['neval']),
                        warnings=[str(w.message) for w in caught],reported_quad_error_not_certified=True,
                        probability_arithmetic_roundoff_not_included=True)

        phase='candidate_and_independent_baseline'
        baseline=None
        if v>0:
            if reuse_hashes is not None:
                reuse_helper = load_module('repair_p3_private_baseline_reuse', archive / BASELINE_REUSE.name)
                current_protocol = json.loads((output/'protocol.json').read_text())
                baseline, independent = reuse_helper.load_reused_baseline(receipt_dir, reuse_hashes, current_protocol)
                save('candidate_baseline',dict(logmass=baseline[0],gradient=baseline[1],hessian=baseline[2],integration=baseline[3],
                                              reused_parent_baseline=True))
                save('independent_baseline',independent)
                first,tight = independent['first'],independent['tighter_expanded']
                bound=max(24.,baseline[3]['refined']['z_bound']+2)
                event('parent_baselines_reused',receipt_directory=str(receipt_dir),
                      independent_derivative_accuracy_not_verified=True)
            else:
                try:
                    baseline=integrate_gene(candidate,p,v,gene_rows,original.conditional_terms,with_hessian=True)
                    save('candidate_baseline',dict(logmass=baseline[0],gradient=baseline[1],hessian=baseline[2],integration=baseline[3]))
                    event('candidate_integral_completed')
                except (ValueError,FloatingPointError,np.linalg.LinAlgError) as error:
                    save('candidate_refusal',dict(completed=False,exception_type=type(error).__name__,traceback=traceback.format_exc(),
                                                 integration_reports=getattr(error,'integration_reports',None),
                                                 likelihood_evaluations=getattr(error,'likelihood_evaluations',None),
                                                 not_a_replacement_failed_draw=True))
                    event('candidate_integral_refused',reason=str(error))
                bound=max(24.,baseline[3]['refined']['z_bound']+2) if baseline else 24.
                first=reference(p,v,bound)
                tight=reference(p,v,bound+2,.1)
                save('independent_baseline',dict(first=first,tighter_expanded=tight,refinement_difference=tight['logmass']-first['logmass'],
                        candidate_difference=baseline[0]-tight['logmass'] if baseline else None,
                        prior_known_draw0_normalized_logmass=-87.25979655580508 if draw==0 else None))
                event('independent_scalar_completed')
        else:
            ll,g,h=candidate.marginal(p,0.,with_hessian=True)
            su_hu=[]
            for gene in range(model.G):
                su=hu=0.
                for row in drawn.rows:
                    if row.gene!=gene: continue
                    eta,rr,_,_=drawn.predictors(row,p,0.,need_hessian=False)
                    _,s,hh=original.conditional_terms(row,eta,rr)
                    su+=s[0]*row.stage
                    hu+=hh[0,0]*row.stage**2
                su_hu.append(dict(gene=gene,su=su,hu=hu,right_v_score=.5*(su*su+hu)))
            gene_ll, gene_g, gene_h=candidate.evaluate(p,gene_rows,u=0.,need_hessian=True)
            gene_boundary_v=su_hu[G]['right_v_score']
            save('boundary_baseline',dict(full_logmass=ll,gradient=g,nuisance_hessian=h[:-1,:-1],
                     original_snapshot_logmass_difference=ll-boundary_original[0],
                     original_snapshot_gradient_difference=g-boundary_original[1],
                     exact_right_variance_score=g[-1],rederived_gene_right_scores=su_hu,
                     rederived_right_score_difference=sum(r['right_v_score'] for r in su_hu)-g[-1],
                     fixed_gene_logmass=gene_ll,independent_fixed_gene_logmass=independent_loglik(p,0.),
                     positive_v_discrepancy_not_claimed_at_retained_v0=True))
            event('exact_boundary_checked')
            if boundary_mp is not None:
                if args.boundary_decimal_fast_sha:
                    save('independent_high_precision_arithmetic_comparison', boundary_mp.compare_small())
                report = boundary_mp.boundary_check(gene_rows, p, model.G, model.Z)
                save('independent_high_precision_boundary', dict(
                    report=report, candidate_gene_right_variance_score=gene_boundary_v,
                    no_replacement_of_declared_finite_difference_records=True))
                event('independent_high_precision_boundary_completed')
            baseline=(gene_ll,np.r_[gene_g,gene_boundary_v],None,None)

        phase='fixed_finite_differences'
        coordinates=(G,model.G+G,2*model.G,2*model.G+1,2*model.G+1+model.Z,model.size-3,model.size)
        require(coordinates == FIXED_AXES, 'Immutable seven-coordinate derivative family differs')
        requested_axes=tuple(args.axes) if args.axes is not None else coordinates
        require(set(requested_axes)<=set(coordinates), 'Requested axes escape the fixed derivative family')
        executed_axes=() if args.baseline_only else requested_axes
        bounds=candidate.bounds()+[(0.,None)]
        fd_records=[]
        for axis in executed_axes:
            for step in STEPS:
                vector=np.r_[p,v]
                lo,hi=bounds[axis]
                central=(lo is None or vector[axis]-step>=lo) and (hi is None or vector[axis]+step<=hi)
                offsets=(-1.,1.) if central else (1.,2.) if lo is not None and vector[axis]-step<lo else (-1.,-2.)
                refs=[]
                for offset in offsets:
                    q=vector.copy();q[axis]+=offset*step
                    if q[-1]==0:
                        refs.append(dict(logmass=independent_loglik(q[:-1],0.),estimated_quad_plus_mass_tail_log_error=0.,
                                         arithmetic_roundoff_unbounded=True))
                    else:
                        refs.append(reference(q[:-1],q[-1],bound if v>0 else 24.,.1))
                center=tight['logmass'] if v>0 else independent_loglik(p,0.)
                center_error=tight['estimated_quad_plus_mass_tail_log_error'] if v>0 else 0.
                if central:
                    fd=(refs[1]['logmass']-refs[0]['logmass'])/(2*step)
                    estimated_error=sum(r['estimated_quad_plus_mass_tail_log_error'] for r in refs)/(2*step)
                else:
                    sign=1 if offsets[0]>0 else -1
                    fd=sign*(-3*center+4*refs[0]['logmass']-refs[1]['logmass'])/(2*step)
                    estimated_error=(3*center_error+4*refs[0]['estimated_quad_plus_mass_tail_log_error']+refs[1]['estimated_quad_plus_mass_tail_log_error'])/(2*step)
                row=dict(axis=axis,step=step,method='central' if central else 'three_point_one_sided',scalar_reference_gradient=fd,
                         estimated_integration_plus_tail_FD_error=estimated_error,truncation_and_arithmetic_error_not_bounded=True,
                         candidate_gradient_difference=fd-baseline[1][axis] if baseline else None)
                # Independent check of Hessian moments through score differentiation;
                # all moment coordinates retained, including p-v and v-v columns.
                if v>0 and baseline and step in STEPS[1:]:
                    scores=[];score_errors=[]
                    for offset in offsets:
                        q=vector.copy();q[axis]+=offset*step
                        try:
                            ll,gg,_,report=integrate_gene(candidate,q[:-1],q[-1],gene_rows,original.conditional_terms)
                            scores.append(gg)
                            score_errors.append(report['refined']['gradient_error_allowance']+report['gradient_refinement_difference'])
                        except (ValueError,FloatingPointError,np.linalg.LinAlgError):
                            row['hessian_score_difference_refused']=traceback.format_exc()
                            break
                    if len(scores)==2:
                        if central:
                            column=(scores[1]-scores[0])/(2*step)
                            estimate=(score_errors[0]+score_errors[1])/(2*step)
                        else:
                            sign=1 if offsets[0]>0 else -1
                            column=sign*(-3*baseline[1]+4*scores[0]-scores[1])/(2*step)
                            base_error=baseline[3]['refined']['gradient_error_allowance']+baseline[3]['gradient_refinement_difference']
                            estimate=(3*base_error+4*score_errors[0]+score_errors[1])/(2*step)
                        row.update(hessian_moment_column=baseline[2][:,axis],score_FD_hessian_column=column,
                                   Hessian_column_difference=column-baseline[2][:,axis],estimated_score_integration_FD_error=estimate)
                fd_records.append(row)
                save(f'finite_difference_{axis:02d}_{str(step).replace(".","_")}',row)
                event('finite_difference_completed',axis=axis,step=step)
        save('finite_differences',fd_records)
        phase='final_source_hashes'
        for path,expected in {**guards,**own}.items():
            require(sha(path)==expected,'Source/fixture changed during fixed-point check')
        save('summary',dict(fixed_point_check_completed=True,case=args.case,draw=draw,positive_variance_gene_scope=G,
             candidate_integral_returned=baseline is not None,no_model_fits=True,no_statistical_calibration=True,
             no_active20_node_source_or_statistical_disposition_changes=True,original_records_preserved=True,
             finite_difference_coordinates=coordinates,finite_difference_steps=STEPS,
             all_fixed_comparisons_retained=not args.baseline_only and set(executed_axes)==set(coordinates),
             all_executed_comparisons_retained=True,
             baseline_only=args.baseline_only,
             derivative_checks_completed=not args.baseline_only and set(executed_axes)==set(coordinates),
             requested_derivative_axes=requested_axes,executed_derivative_axes=executed_axes,
             required_combined_derivative_family_axes=coordinates,all_declared_axes_remain_required=True,
             candidate_derivative_comparison_available=baseline is not None,
             baseline_candidate_error_estimates_not_certified=True,independent_reference_error_estimates_not_certified=True,
             probability_arithmetic_and_FD_truncation_unbounded=True,
             numerical_file_hashes={str(Path(np.__file__).resolve()):sha(Path(np.__file__).resolve()),
                                   str(Path(scipy.__file__).resolve()):sha(Path(scipy.__file__).resolve()),
                                   str(Path(sys.executable).resolve()):sha(Path(sys.executable).resolve())},
             elapsed_seconds=time.monotonic()-began,output_cap_bytes=CAP,
             saved_bytes_before_summary=sum(q.stat().st_size for q in output.rglob('*') if q.is_file())))
        event('completed')
    except BaseException:
        with (output/'failure.json').open('x') as handle:
            json.dump(dict(completed=False,phase=phase,traceback=traceback.format_exc(),partials_preserved=True,
                           no_fit_or_calibration=True),handle,indent=2)
            handle.write('\n')
        raise


if __name__=='__main__':
    main()
