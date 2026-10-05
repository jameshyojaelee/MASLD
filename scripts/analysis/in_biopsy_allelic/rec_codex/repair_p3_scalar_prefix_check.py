#!/usr/bin/env python3
"""Compute-only scalar-prefix fixed points and full FD family; no candidate checks."""
import argparse
import ast
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
PREFIX = HERE / 'repair_p3_reference_prefix.py'
PREFIX_SHA = '3c4d0296612010a40e15f60d0d13c10489b9d82d86595f0da344c7299d0b2f66'
ANCESTOR_SHA = '89434fd7e32b06cb93d9d31f18b646eaa646736ec329ab48cc00fa4409561f72'
ANCESTOR = REC / 'repair_p3_launch_20261005T121740Z/repair_p3_fixed_point_check.py'
CAP = 32 * 1024**2
CASES = {'failed0gene2': 0, 'positive18': 18, 'boundary19': 19}
G = 2
STEPS = (1e-3, 3e-4, 1e-4)
FIXED_AXES = (2,10,16,17,24,25,28)
GUARDS = {
    PROTOCOL: '8e7ddf8b0a7bb18afd9692b53c3560f4258fdc373c4b822aee36aa94395eabe9',
    CAMPAIGN / 'campaign.json': 'bb3e82833476f67b0b11d2388e2dc00358eb987367b7a6a666d11938c6baca94',
    PREFIX: PREFIX_SHA,
    ANCESTOR: ANCESTOR_SHA,
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
    ast.parse(PREFIX.read_text())
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
    parser.add_argument('--expected-prefix-sha')
    args = parser.parse_args()
    protocol, campaign, guards = metadata()
    if args.inspect:
        print(json.dumps(dict(case=args.case, draw=CASES[args.case], fixed_gene=G,
              numerical_execution=False, model_fitting=False, calibration=False,
              candidate_integral_checks_completed=False, candidate_hessian_checks_completed=False,
              scalar_derivative_family_axes=FIXED_AXES, finite_difference_steps=STEPS,
              prefix_sha256=sha(PREFIX), ancestor_sha256=ANCESTOR_SHA, metadata_and_source_checks=True)))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'), 'Individual compute job required')
    require(args.expected_helper_sha == sha(Path(__file__)) and args.expected_prefix_sha == sha(PREFIX), 'Reviewed scalar/checker source drift')
    output = REC / ('repair_p3_scalar_prefix_' + args.case + '_' + os.environ['SLURM_JOB_ID'])
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
        own = {Path(__file__):args.expected_helper_sha, PREFIX:args.expected_prefix_sha}
        for source in (*own, *(p for p in guards if p.suffix == '.py')):
            (archive / source.name).write_bytes(source.read_bytes())
            require(sha(archive / source.name) == (own.get(source) or guards[source]), 'Archive copy differs')
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
        prefix = load_module('repair_p3_archived_scalar_prefix', archive / PREFIX.name)
        require(Path(prefix.__file__).resolve().parent == archive.resolve(), 'Prefix import escaped archive')
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
        save('protocol', dict(case=args.case, draw=draw, fixed_gene=G, fixed_gene_rows=24, fixed_p=p, fixed_v=v,
             count_vector_sha256=count_sha, source_guards={str(k):d for k,d in guards.items()}, own_source_hashes={str(k):d for k,d in own.items()},
             environment=versions, finite_difference_steps=STEPS, scalar_derivative_family_axes=FIXED_AXES,
             no_model_fit=True, no_calibration=True, ancestor_sha256=ANCESTOR_SHA,
             candidate_integral_checks_completed=False, candidate_hessian_checks_completed=False,
             candidate_analytical_boundary_checks_completed=False,
             initial_positive_variance_scope='Only fixed gene2; no all8-gene integral or optimizer.',
             original_records_preserved=True,
             stochastic_draw_source='Original supplied uniform row only; no new Monte Carlo seed.',
             seed_for_future_synthetic_extensions=20261005, output_cap_bytes=CAP,
             original_completed_agreement_fixture_B50_not_calibration=True,
             independent_reference='Literal mixture/selection assembly and scalar QUADPACK; only het BB replaced by reviewed prefix value module.',
             conditional_probability_roundoff_validation='Prefix equivalence check22039113; no candidate conditional_terms or derivatives reused.'))

        # Independent scalar probability assembly: literal selected support, original
        # raw three-class weights and orientation. Only BB arithmetic uses prefixes.
        def rise(logshape, count):
            return 0. if count == 0 else logshape + np.logaddexp(logshape, np.log(np.arange(1,count,dtype=float))).sum()
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
        def prepare_for_q(q):
            return {id(row):prefix.PreparedBB(row.n, q[model.G+row.gene]+np.dot(q[-3:],(row.stage,row.background,row.duplication)))
                    for row in gene_rows}

        def independent_loglik(q,u,prepared=None):
            prepared = prepare_for_q(q) if prepared is None else prepared
            value = 0.
            for row in gene_rows:
                eta = q[row.gene]*(1+q[2*model.G]*row.stage+np.dot(q[2*model.G+1:2*model.G+1+model.Z],row.z))
                eta += u*row.stage-row.orientation*q[2*model.G+1+model.Z]*row.stage-row.orientation*row.omega
                support, prior, homref, homalt, at = fixed_parts[id(row)]
                parts = np.stack((prepared[id(row)].logpmf(eta)[support],homref,homalt))+prior[:,None]
                value += logsumexp(parts[:,at])-logsumexp(parts)
            require(np.isfinite(value) and value <= 1e-8, 'Independent conditional PMF invalid')
            return float(value)

        def reference(q, variance, bound=24., precision=1.):
            require(variance > 0, 'Positive variance scalar reference required')
            sqrtv = np.sqrt(variance)
            prepared = prepare_for_q(q)  # fresh contexts for every FD parameter vector
            cache = {}
            def h(z):
                if float(z) not in cache:
                    require(len(cache) < 100000, 'Independent scalar evaluation cap reached')
                    cache[float(z)] = independent_loglik(q,sqrtv*z,prepared)-z*z/2
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
                        probability_arithmetic_roundoff_not_included=True,
                        mode_search_grid_step=1.,mode_search_not_exhaustive=True)

        phase='independent_scalar_baseline'
        bound=24.
        if v>0:
            first=reference(p,v,bound)
            save('independent_baseline_first',first)
            event('independent_scalar_first_completed',logmass=first['logmass'])
            tight=reference(p,v,bound+2,.1)
            save('independent_baseline',dict(first=first,tighter_expanded=tight,
                refinement_difference=tight['logmass']-first['logmass'],
                prior_known_draw0_normalized_logmass=-87.25979655580508 if draw==0 else None,
                candidate_comparison_not_run=True))
            event('independent_scalar_baseline_completed',logmass=tight['logmass'])
        else:
            first=tight=dict(logmass=independent_loglik(p,0.),estimated_quad_plus_mass_tail_log_error=0.,
                arithmetic_roundoff_unbounded=True,no_quadrature_at_exact_boundary=True)
            save('independent_boundary_baseline',tight)
            event('independent_exact_u0_completed',logmass=tight['logmass'])

        phase='fixed_finite_differences'
        coordinates=(G,model.G+G,2*model.G,2*model.G+1,2*model.G+1+model.Z,model.size-3,model.size)
        require(coordinates == FIXED_AXES, 'Immutable seven-coordinate derivative family differs')
        bounds=drawn.bounds()+[(0.,None)]
        fd_records=[]
        for axis in coordinates:
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
                         scalar_reference_nodes=refs,scalar_reference_center_logmass=center,
                         scalar_reference_center_estimated_error=center_error,
                         signed_offsets=offsets,candidate_comparison_not_run=True)
                fd_records.append(row)
                save(f'finite_difference_{axis:02d}_{str(step).replace(".","_")}',row)
                event('finite_difference_completed',axis=axis,step=step)
        save('finite_differences',fd_records)
        require(len(fd_records)==len(FIXED_AXES)*len(STEPS), 'Incomplete fixed scalar FD family')
        phase='final_source_hashes'
        for path,expected in {**guards,**own}.items():
            require(sha(path)==expected,'Source/fixture changed during fixed-point check')
        save('summary',dict(scalar_fixed_point_check_completed=True,case=args.case,draw=draw,positive_variance_gene_scope=G,
             candidate_integral_checks_completed=False,candidate_hessian_checks_completed=False,
             candidate_analytical_boundary_checks_completed=False,candidate_derivative_comparison_available=False,
             no_model_fits=True,no_statistical_calibration=True,original_records_preserved=True,
             finite_difference_coordinates=coordinates,finite_difference_steps=STEPS,scalar_fd_record_count=len(fd_records),
             all_fixed_scalar_comparisons_retained=True,derivative_accuracy_pass_not_asserted=True,
             independent_reference_error_estimates_not_certified=True,probability_arithmetic_and_FD_truncation_unbounded=True,
             full_model_stationarity_not_assessed=True,
             numerical_file_hashes={str(Path(np.__file__).resolve()):sha(Path(np.__file__).resolve()),
                                   str(Path(scipy.__file__).resolve()):sha(Path(scipy.__file__).resolve()),
                                   str(Path(sys.executable).resolve()):sha(Path(sys.executable).resolve())},
             elapsed_seconds=time.monotonic()-began,output_cap_bytes=CAP,
             saved_bytes_before_summary=sum(q.stat().st_size for q in output.rglob('*') if q.is_file())))
        event('completed')
        require(sum(q.stat().st_size for q in output.rglob('*') if q.is_file())<=CAP, 'Final inclusive output cap exceeded')
    except BaseException:
        with (output/'failure.json').open('x') as handle:
            json.dump(dict(completed=False,phase=phase,traceback=traceback.format_exc(),partials_preserved=True,
                           no_fit_or_calibration=True),handle,indent=2)
            handle.write('\n')
        raise


if __name__=='__main__':
    main()
