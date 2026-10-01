"""Fixed synthetic six-context integration diagnostic; no fits or calibration."""
import ast
import copy
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
from scipy.special import logsumexp
from scipy.stats import norm

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
ROLE = REC / 'bfix_p3_21998519'
DIAG = REC / 'p3_first_draw_diagnostic_21998696'
MODE = REC / 'codex_p3_gene2_integrand_21998880'
FIX = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
HERE = Path(__file__).resolve().parent
PROBE = MODE / 'executed_sources/codex_stable_bb_probe.py'
GUARDS = {
    DIAG / 'summary.json': 'b68b8fbd37264d6e93428e1368fce5518c07ea743e23b5fd5580baadff448bd0',
    DIAG / 'receipt.json': '3b313579039c403094bcea2464f4cc23c2ee0b11dec22189047abc5f6e4f3630',
    MODE / 'generation.json': '23af2ebf77c6c0d52b6a1e59f8681bb0eb202f0db5172e2637f10c73682ed92c',
    PROBE: '723781113c17061b9a28341a91149bbd4ab71f195121a67f88d8cbc710a21488',
}
INPUTS = ('globals.json', 'individuals.tsv', 'rows.tsv', 'cohorts.tsv', 'draws/p3_uniforms.npy')
VARIANCES = (.05, .5, 2.)
NODES = (20, 40, 80)
FD_STEPS = (3e-4, 1e-4)
CAP = 32 << 20


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def serial(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def mode_code(source):
    marginal = next(n for n in ast.walk(ast.parse(source))
                    if isinstance(n, ast.FunctionDef) and n.name == 'marginal')
    loop = next(n for n in marginal.body if isinstance(n, ast.For))
    start = next(i for i, n in enumerate(loop.body) if isinstance(n, ast.FunctionDef) and n.name == 'hstats')
    stop = next(i for i, n in enumerate(loop.body) if isinstance(n, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == 'scale' for t in n.targets))
    fn = ast.parse('def diagnostic_mode(self,p,v,rs):\n return None\n').body[0]
    fn.body = copy.deepcopy(loop.body[start:stop+1]) + ast.parse('return mode,scale,hstats(mode)').body
    tree = ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[]))
    return compile(tree, '<immutable-mode-extraction>', 'exec'), hashlib.sha256(ast.dump(tree).encode()).hexdigest()


def main():
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Compute allocation required')
    started = time.monotonic()
    out = REC / ('codex_positive_variance_grid_'+os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)
    private = out / 'executed_sources'
    private.mkdir()
    phase = 'guards'
    evaluations = 0
    cases = []

    def save(name, obj):
        (out / (name+'.json')).write_text(json.dumps(obj, default=serial, indent=2, allow_nan=False)+'\n')
        require(sum(p.stat().st_size for p in out.rglob('*') if p.is_file()) <= CAP, 'Output cap exceeded')

    def failure():
        return dict(phase=phase, exception=traceback.format_exc())

    try:
        guards = dict(GUARDS)
        for path, digest in guards.items():
            require(sha(path) == digest, 'Frozen metadata differs: '+str(path))
        receipt = json.loads((DIAG / 'receipt.json').read_text())
        failed = json.loads((DIAG / 'summary.json').read_text())
        guards[ROLE / 'receipt.json'] = receipt['role_receipt_sha256']
        guards[ROLE / 'observed_null.json'] = receipt['observed_null_sha256']
        guards.update({ROLE / 'executed_sources' / n: h for n, h in receipt['source_sha256'].items()})
        guards.update({FIX / n: receipt['fixture_sha256'][n] for n in INPUTS})
        guards[Path(__file__)] = sha(Path(__file__))
        launcher = HERE / 'run_codex_positive_variance_grid.sbatch'
        guards[launcher] = sha(launcher)
        for path, digest in guards.items():
            require(sha(path) == digest, 'Frozen source/input differs: '+str(path))
        for name in receipt['source_sha256']:
            (private / name).write_bytes((ROLE / 'executed_sources' / name).read_bytes())
        original = (private / 'likelihood_v1.py').read_text()
        require(original.count('2.**-np.arange(30)') == 1, 'Damping source occurrence differs')
        altered = original.replace('2.**-np.arange(30)', '2.**-np.arange(64)')
        # Reversal establishes there is exactly this one textual/AST change.
        require(ast.dump(ast.parse(altered.replace('2.**-np.arange(64)', '2.**-np.arange(30)')))
                == ast.dump(ast.parse(original)), 'Unexpected source AST change')
        (private / 'likelihood_v1.py').write_text(altered)
        node_source_hashes = {}
        for n in NODES:
            variant = altered.replace('nodes=20,', 'nodes='+str(n)+',')
            require(altered.count('nodes=20,') == 1, 'Node default occurrence differs')
            path = private / ('likelihood_nodes_'+str(n)+'.py')
            path.write_text(variant)
            node_source_hashes[str(n)] = sha(path)
        for path in (Path(__file__), launcher, PROBE):
            (private / path.name).write_bytes(path.read_bytes())
        private_guards = {path: sha(path) for path in private.iterdir() if path.is_file()}
        with (DIAG / 'progress.jsonl').open() as handle:
            env = next(json.loads(line) for line in handle if json.loads(line)['event'] == 'imports_verified')
        require((env['python'], env['numpy'], env['scipy']) ==
                (platform.python_version(), np.__version__, scipy.__version__), 'Draw environment differs')
        sys.path.insert(0, str(private))
        require(not any(n in sys.modules for n in ('likelihood_v1', 'bfix_v1', 'design_v1', 'tests_v1')), 'Already imported')
        loader = importlib.import_module('bfix_v1')
        likelihood = importlib.import_module('likelihood_v1')
        for name in ('likelihood_v1', 'bfix_v1', 'design_v1', 'tests_v1'):
            require(Path(sys.modules[name].__file__).parent == private, 'Import escaped private snapshot')
        design, settings = loader.load()
        original_model = design.model()
        observed = json.loads((ROLE / 'observed_null.json').read_text())
        null_p = np.asarray(observed['parameters'], float)
        require(observed['converged'] and settings['B'] == 50 and len(original_model.rows) == 173, 'Fixture differs')
        require([(r['individual_id'], r['gene_id']) for r in loader.table('rows.tsv')] ==
                [(r.individual, settings['gene_order'][r.gene]) for r in original_model.rows], 'Row identity/order differs')
        uniforms = np.load(FIX / 'draws/p3_uniforms.npy', allow_pickle=False, mmap_mode='r')
        require(uniforms.shape == (50, 173, 2), 'Uniforms shape differs')
        drawn = []
        for row, (u1, u2) in zip(original_model.rows, uniforms[0]):
            eta, r, _, _ = original_model.predictors(row, null_p)
            drawn.append(likelihood.conditional_draw(row, eta, r, u1, u2))
        count_sha = hashlib.sha256(np.asarray([r.a for r in drawn], dtype='<i8').tobytes()).hexdigest()
        generation = json.loads((MODE / 'generation.json').read_text())
        require(count_sha == generation['reconstructed_count_vector_sha256'] and len(drawn) == 173
                and all(r.a in r.allowed for r in drawn), 'Draw reproduction differs')
        model = likelihood.Model(drawn)
        proposal = np.asarray(failed['latest_marginal']['parameters'], float)
        require(model.G == 8 and len(null_p) == len(proposal) == model.size, 'Parameter/gene dimensions differ')
        spec = importlib.util.spec_from_file_location('isolated_stable_bb', private / PROBE.name)
        stable = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(stable)
        require(stable.EPS == likelihood.EPS and 0 < stable.EPS < .5, 'BB rho-floor contract differs')
        def replacement(k, n, eta, r, need_hessian=True):
            ll, s, h = stable.stable_bb_terms(k, n, eta, r)
            return ll, s, h if need_hessian else None
        likelihood.bb_terms = replacement
        code, mode_ast_sha = mode_code(altered)
        namespace = dict(likelihood.__dict__)
        exec(code, namespace)
        mode_function = namespace['diagnostic_mode']
        names = (['alpha_'+str(j) for j in range(model.G)]+['r_'+str(j) for j in range(model.G)]+['kappa']
                 +['delta_'+str(j) for j in range(model.Z)]+['omega_S', 'rho_S', 'rho_b', 'rho_d', 'v'])
        save('protocol', dict(variances=VARIANCES, nodes=NODES, all_genes=list(range(8)),
            parameter_vectors=dict(observed_null=null_p, failed_proposal=proposal), parameter_names=names,
            fd_genes=[0, 2], fd_directions='local alpha,r,kappa,v', fd_steps=FD_STEPS,
            source_guards={str(p): h for p, h in guards.items()}, count_sha256=count_sha,
            original_source_sha256=hashlib.sha256(original.encode()).hexdigest(),
            private64_source_sha256=sha(private / 'likelihood_v1.py'), mode_ast_sha256=mode_ast_sha,
            node_variant_source_sha256=node_source_hashes,
            private_source_guards={path.name: digest for path,digest in private_guards.items()},
            changes='stable BB arithmetic injection; sole private source change damping30 to64; explicit nodes20/40/80 diagnostic arguments',
            reference='normalized Gaussian-prior mass and Fisher-identity score vector; bounded adaptive vector quadrature',
            p_score_envelope='For each row, both unconditioned BB eta/r scores are bounded by n; fixed-prior mixtures retain that bound; conditioning subtracts another bounded mean, giving 2*n per predictor. Jacobian is independent of u. See uniform_p_score_bound source derivation.',
            reference_grid_spacing=.05, reference_bound_cap=64., gaussian_tail='proper conditional PMF <=1',
            environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform(),
                             threads={k:os.environ.get(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')}),
            no_fit=True, no_MC=True, deterministic=True, diagnostic_not_calibration=True, owner20_node_rule_unchanged=True))

        for label, p in (('observed_null', null_p), ('failed_proposal', proposal)):
            for v in VARIANCES:
                context = dict(parameter=label, v=v, genes=[])
                for gene in range(8):
                    phase = f'{label}_v{v}_gene{gene}'
                    rowset = [r for r in drawn if r.gene == gene]
                    result = dict(gene=gene, rows=len(rowset), production=[], reference=None, finite_differences=[])
                    def stats(q, variance, u):
                        nonlocal evaluations
                        evaluations += 1
                        ll, g, _ = model.evaluate(q, rowset, u=u, need_hessian=False)
                        require(np.isfinite(ll) and np.isfinite(g).all() and ll <= 1e-8, 'Improper/nonfinite mass')
                        return ll-u*u/(2*variance)-.5*np.log(2*np.pi*variance), np.r_[g, (u*u-variance)/(2*variance*variance)]
                    def uniform_p_score_bound(q):
                        """Exact-real-arithmetic bound for this frozen model.

                        For rising factorials, D_a=sum_{i<k} a/(a+i)
                        is in [0,k]; D_b in [0,n-k]; D_t in [0,n].
                        BB s_eta=(1-m)*D_a-m*D_b, so |s_eta|<=n.
                        BB s_r=(D_a+D_b-D_t)*d(log t)/dr. Both
                        D_a+D_b and D_t lie in [0,n]. Write c=1-2e,
                        x=logistic(r)*(1-logistic(r))<=1/4. Then
                        rho(1-rho)-rho'=e(1-e)-2*e*c*x>=e/2>0,
                        hence |d(log t)/dr|=rho'/[rho(1-rho)]<=1.
                        Thus |s_r|<=n. The frozen genotype prior and
                        homozygote PMFs have zero parameter derivatives;
                        a mixture score is posterior het weight times BB
                        score, still <=n. Inclusion score is its weighted
                        mean over allowed counts, still <=n. Conditional
                        score subtracts that mean: <=2n per predictor.
                        Predictor Jacobians do not depend on random u.
                        This does not bound floating-point arithmetic error.
                        """
                        return sum((2*r.n*np.abs(model.predictors(r,q,0.,need_hessian=False)[2]).sum(axis=0)
                                    for r in rowset), np.zeros(len(q)))
                    try:
                        mode, scale, mst = mode_function(model, p, v, rowset)
                        require(np.isfinite(mst).all(), 'Mode nonfinite')
                        for n in NODES:
                            x, w = np.polynomial.hermite.hermgauss(n)
                            values = [stats(p, v, float(u)) for u in mode+np.sqrt(2)*scale*x]
                            terms = np.log(w)+x*x+np.array([a[0] for a in values])
                            total = logsumexp(terms)
                            result['production'].append(dict(nodes=n, mode=mode, scale=scale, mode_stats=mst,
                                stationarity_within_1e_5=bool(abs(mst[1]) <= 1e-5),
                                ll=np.log(np.sqrt(2)*scale)+total, score=np.exp(terms-total)@np.stack([a[1] for a in values])))
                    except Exception:
                        result['production_failure'] = failure()
                    try:
                        zero = stats(p, v, 0.)[0]
                        requested = max(8*np.sqrt(v), np.sqrt(2*v*abs(zero))+4*np.sqrt(v))
                        bound = min(64., requested)
                        grid = np.linspace(-bound, bound, int(np.ceil(2*bound/.05))+1)
                        masses = np.array([stats(p, v, u)[0] for u in grid])
                        peaks = [i for i in range(1,len(grid)-1) if masses[i] >= masses[i-1] and masses[i] >= masses[i+1]]
                        points = [float(grid[i]) for i in peaks]
                        shift = float(masses.max())
                        def integrate(q, variance, scores=True, tight=True, extent=bound):
                            def function(u):
                                h, s = stats(q, variance, u)
                                weight = np.exp(h-shift)
                                return weight*np.r_[1., s] if scores else np.array([weight])
                            with warnings.catch_warnings(record=True) as caught:
                                warnings.simplefilter('always')
                                value, error, info = quad_vec(function, -extent, extent, points=points,
                                    epsabs=1e-12 if tight else 1e-10, epsrel=1e-10 if tight else 1e-8,
                                    norm='max', quadrature='gk21', limit=600, full_output=True)
                            require(info.success and np.isfinite(value).all() and value[0] > error >= 0, 'Reference unresolved: '+info.message)
                            mass = value[0]
                            lower = shift+np.log(mass-error)
                            logtail = np.log(2)+norm.logsf(extent/np.sqrt(variance))
                            ratio = np.exp(logtail-lower)
                            answer = dict(ll=shift+np.log(mass), bound=extent, error=error, relative_mass_error=error/mass,
                                logmass_error_estimate=-np.log1p(-error/mass), prior_tail_log_bound=logtail,
                                tail_to_mass_bound_using_estimated_lower_mass=ratio, logmass_tail_bound=np.log1p(ratio),
                                evaluations=info.neval, success=bool(info.success), warnings=[str(w.message) for w in caught])
                            if scores:
                                gradient = value[1:]/mass
                                ce = error*(1+np.abs(gradient))/(mass-error)
                                sb = uniform_p_score_bound(q)
                                a = extent/np.sqrt(variance)
                                vt = logsumexp([np.log(a)+norm.logpdf(a), np.log(2)+norm.logsf(a)])-np.log(variance)
                                tail = np.r_[ratio*(sb+np.abs(gradient[:-1])+ce[:-1]),
                                    np.exp(vt-lower)+ratio*(abs(gradient[-1])+ce[-1])]
                                answer.update(score=gradient, component_error_estimates=ce, score_tail_bounds=tail)
                            return answer
                        coarse = integrate(p,v,tight=False)
                        ref = integrate(p,v)
                        expanded = integrate(p,v,extent=min(64.,bound+2*np.sqrt(v)))
                        result['reference'] = dict(tight=ref, coarse=coarse, expanded=expanded,
                            requested_bound=requested, cap_truncated=requested>64, grid_spacing=float(grid[1]-grid[0]),
                            grid_peak_count=len(peaks), grid_peak_locations=points,
                            qualification='Grid need not locate all modes; QUADPACK errors/lower mass are estimates. Prior envelope bounds mass tails, not global optimization.')
                        for item in result['production']:
                            item.update(ll_minus_reference=item['ll']-ref['ll'], score_minus_reference=item['score']-ref['score'])
                        if gene in (0,2):
                            for coordinate in (gene, model.G+gene, 2*model.G, len(p)):
                                for step in FD_STEPS:
                                    qp,qm=p.copy(),p.copy();vp,vm=v,v
                                    if coordinate == len(p): vp+=step;vm-=step
                                    else: qp[coordinate]+=step;qm[coordinate]-=step
                                    plus,minus=integrate(qp,vp,scores=False),integrate(qm,vm,scores=False)
                                    fd=(plus['ll']-minus['ll'])/(2*step)
                                    err=sum(t['logmass_error_estimate']+t['logmass_tail_bound'] for t in (plus,minus))/(2*step)
                                    result['finite_differences'].append(dict(coordinate=coordinate,name=names[coordinate],step=step,
                                        score_fd=fd, minus_reference=fd-ref['score'][coordinate], quad_plus_tail_error_estimate=err))
                    except Exception:
                        result['reference_failure'] = failure()
                    context['genes'].append(result)
                    save(label+'_v'+str(v)+'_gene'+str(gene),result)
                    print(phase, 'saved', evaluations, flush=True)
                for n in NODES:
                    items=[next((x for x in g['production'] if x['nodes']==n),None) for g in context['genes']]
                    if all(x is not None for x in items):
                        context.setdefault('full_model_production',[]).append(dict(nodes=n,ll=sum(x['ll'] for x in items),score=sum((x['score'] for x in items),np.zeros(len(p)+1))))
                    try:
                        # Exercise the private production entry point, not just
                        # its extracted mode and independently assembled nodes.
                        ll, score = model.marginal(p, v, nodes=n)
                        check = dict(nodes=n, ll=ll, score=score)
                        if all(x is not None for x in items):
                            check.update(ll_minus_assembled=ll-sum(x['ll'] for x in items),
                                score_minus_assembled=score-sum((x['score'] for x in items),np.zeros(len(p)+1)))
                        context.setdefault('production_entry_point_checks',[]).append(check)
                    except Exception:
                        context.setdefault('production_entry_point_failures',[]).append(dict(nodes=n,**failure()))
                refs=[g['reference'] for g in context['genes']]
                if all(x is not None for x in refs):
                    context['full_model_reference']=dict(ll=sum(x['tight']['ll'] for x in refs),score=sum((np.asarray(x['tight']['score']) for x in refs),np.zeros(len(p)+1)))
                cases.append(context)
                save('contexts_so_far',cases)
        phase='post_guards'
        for path,digest in guards.items():
            require(sha(path)==digest,'Post-run frozen source/input differs: '+str(path))
        require(sha(private / 'likelihood_v1.py') == hashlib.sha256(altered.encode()).hexdigest(), 'Private source mutated')
        for n in NODES:
            require(sha(private / ('likelihood_nodes_'+str(n)+'.py')) == node_source_hashes[str(n)], 'Node source mutated')
        for path,digest in private_guards.items():
            require(sha(path)==digest,'Private copied source mutated: '+path.name)
        require(count_sha==hashlib.sha256(np.asarray([r.a for r in drawn],dtype='<i8').tobytes()).hexdigest(),'Rows mutated')
        save('summary',dict(contexts=cases,elapsed_seconds=time.monotonic()-started,evaluations=evaluations,
            all_cases_retained=True, numerical_diagnostic_only=True, no_fit=True, no_MC=True,
            limits='No change/adoption of owner20-node law. No global maximum, general accuracy, coverage or calibration claim. Reference never replaces production failure.'))
    except BaseException:
        save('failure',dict(**failure(),elapsed_seconds=time.monotonic()-started,evaluations=evaluations))
        raise


if __name__=='__main__':
    main()
