"""Isolated 30-to-64 mode-damping diagnostic; one supplied synthetic P3 draw.

An archived fixed-gene quadrature check must pass before any fitting. This
does not edit production sources, alter model/test law, or claim calibration.
"""
import ast
import copy
import difflib
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import traceback
import types

import numpy as np
from scipy.special import logsumexp


ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
HERE = Path(__file__).resolve().parent
ROLE = REC / 'bfix_p3_21998519'
DIAG = REC / 'p3_first_draw_diagnostic_21998696'
QUAD = REC / 'codex_p3_gene2_integrand_21998880/summary.json'
GENERATION = QUAD.with_name('generation.json')
GUARDS = {
    QUAD: '2d51e8d2eca5db5491ba379a51616c2dec25cd4711092606c892d724b8154792',
    GENERATION: '23af2ebf77c6c0d52b6a1e59f8681bb0eb202f0db5172e2637f10c73682ed92c',
    DIAG / 'summary.json': 'b68b8fbd37264d6e93428e1368fce5518c07ea743e23b5fd5580baadff448bd0',
    DIAG / 'receipt.json': '3b313579039c403094bcea2464f4cc23c2ee0b11dec22189047abc5f6e4f3630',
    HERE / 'codex_stable_bb_probe.py': '723781113c17061b9a28341a91149bbd4ab71f195121a67f88d8cbc710a21488',
    HERE / 'codex_replay_first_p3_stable_bb.py': '53cbd331ce8d65d387c2034d1fd0a955197bc47c44447d68fa326b906774de2e',
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


def exact_mode_function(source, filename):
    """Extract the executed hstats/mode statements, preserving their AST."""
    marginal = next(n for n in ast.walk(ast.parse(source))
                    if isinstance(n, ast.FunctionDef) and n.name == 'marginal')
    gene_loop = next(n for n in marginal.body if isinstance(n, ast.For))
    start = next(i for i, n in enumerate(gene_loop.body)
                 if isinstance(n, ast.FunctionDef) and n.name == 'hstats')
    stop = next(i for i, n in enumerate(gene_loop.body)
                if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'scale' for t in n.targets))
    template = ast.parse('def fixed_mode(self,p,v,rs):\n return None\n').body[0]
    template.body = copy.deepcopy(gene_loop.body[start:stop+1])
    template.body += ast.parse('return mode,scale,hstats(mode),hstats').body
    program = ast.fix_missing_locations(ast.Module(body=[template], type_ignores=[]))
    return compile(program, filename, 'exec'), hashlib.sha256(ast.dump(program).encode()).hexdigest()


def main():
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Compute allocation required')
    job = os.environ['SLURM_JOB_ID']
    output = REC / ('codex_p3_private_damping64_'+job)
    output.mkdir(exist_ok=False)
    private = output / 'executed_sources'
    private.mkdir()
    def save(name, value):
        (output / (name+'.json')).write_text(json.dumps(value, default=serial, indent=2)+'\n')
    for path, digest in GUARDS.items():
        require(sha(path) == digest, 'Reviewed source/reference changed: '+str(path))
    receipt = json.loads((DIAG / 'receipt.json').read_text())
    report = json.loads((DIAG / 'summary.json').read_text())
    reference = json.loads(QUAD.read_text())
    generation = json.loads(GENERATION.read_text())
    require(sha(ROLE / 'receipt.json') == receipt['role_receipt_sha256'], 'Original role receipt differs')
    require(sha(ROLE / 'observed_null.json') == receipt['observed_null_sha256'], 'Original generation fit differs')
    for name, digest in receipt['source_sha256'].items():
        path = ROLE / 'executed_sources' / name
        require(sha(path) == digest, 'Immutable executed source differs: '+name)
        (private / name).write_bytes(path.read_bytes())
    for path in (Path(__file__), HERE / 'codex_stable_bb_probe.py', HERE / 'codex_replay_first_p3_stable_bb.py'):
        (private / path.name).write_bytes(path.read_bytes())
    path = private / 'likelihood_v1.py'
    original = path.read_text()
    old = 'for damping in 2.**-np.arange(30):'
    new = 'for damping in 2.**-np.arange(64):'
    require(original.count(old) == 1, 'Require exactly one gene-mode damping occurrence')
    patched = original.replace(old, new)
    require(patched.replace(new, old) == original and patched.count(new) == 1, 'Source change is not isolated')
    path.write_text(patched)
    (output / 'private_source.diff').write_text(''.join(difflib.unified_diff(
        original.splitlines(keepends=True), patched.splitlines(keepends=True),
        fromfile='immutable/likelihood_v1.py', tofile='private/likelihood_v1.py')))
    protocol = dict(original_likelihood_sha256=receipt['source_sha256']['likelihood_v1.py'],
                    private_likelihood_sha256=sha(path), numerical_change='Only mode damping proposals30→64',
                    stable_probe_sha256=GUARDS[HERE / 'codex_stable_bb_probe.py'],
                    replay_driver_sha256=GUARDS[HERE / 'codex_replay_first_p3_stable_bb.py'],
                    reference_quad_sha256=GUARDS[QUAD], reference_log_integral=reference['log_bounded_integral'],
                    fixed_gene=2, nodes=20, score_tolerance=1e-6, integral_log_tolerance=1e-6,
                    root_absolute_tolerance=1e-6, original_failed_draw=str(ROLE / 'p3_000.json'),
                    statistical_law_changed=False, owner_spec_amended=False, calibration_claim=False,
                    nested_replay_labels_relative_to_already_patched_snapshot=True,
                    source_sha256=sha(Path(__file__)),
                    launcher_sha256=sha(HERE / 'run_codex_p3_private_damping64_replay.sbatch'))
    save('protocol', protocol)
    phase = 'fixed_parameter_precheck_no_fitting'
    try:
        sys.path.insert(0, str(private))
        require(not any(n in sys.modules for n in ('bfix_v1', 'design_v1', 'likelihood_v1', 'tests_v1')),
                'Source isolation violated before import')
        loader = importlib.import_module('bfix_v1')
        likelihood = importlib.import_module('likelihood_v1')
        design, settings = loader.load()
        model = design.model()
        require(len(model.rows) == 173 and settings['B'] == 50, 'Original fixture dimensions differ')
        for name in ('globals.json', 'individuals.tsv', 'rows.tsv', 'cohorts.tsv', 'draws/p3_uniforms.npy'):
            require(sha(loader.FIXTURE / name) == receipt['fixture_sha256'][name], 'Frozen synthetic input differs: '+name)
        observed = json.loads((ROLE / 'observed_null.json').read_text())
        null_p = np.asarray(observed['parameters'])
        uniforms = np.load(loader.FIXTURE / 'draws/p3_uniforms.npy', allow_pickle=False, mmap_mode='r')
        require(uniforms.shape == (50, 173, 2), 'Original supplied uniforms differ')
        drawn = []
        for row, (u1, u2) in zip(model.rows, uniforms[0]):
            eta, r, _, _ = model.predictors(row, null_p)
            drawn.append(likelihood.conditional_draw(row, eta, r, u1, u2))
        drawn_model = likelihood.Model(drawn)
        rows = [r for r in drawn if r.gene == 2]
        count_sha = hashlib.sha256(np.asarray([r.a for r in drawn], dtype='<i8').tobytes()).hexdigest()
        require(len(drawn) == generation['synthetic_rows'] == 173
                and len(rows) == generation['gene_rows'] == 24
                and count_sha == generation['reconstructed_count_vector_sha256'],
                'Reconstructed first draw differs from completed fixed-integrand diagnostic')
        p = np.asarray(report['latest_marginal']['parameters'])
        v = report['latest_marginal']['v']
        stable_spec = importlib.util.spec_from_file_location('private_damping64_bb', private / 'codex_stable_bb_probe.py')
        stable = importlib.util.module_from_spec(stable_spec)
        sys.modules[stable_spec.name] = stable
        stable_spec.loader.exec_module(stable)
        def replacement(k, n, eta, linear_r, need_hessian=True):
            ll, g, h = stable.stable_bb_terms(k, n, eta, linear_r)
            return ll, g, h if need_hessian else None
        likelihood.bb_terms = replacement
        code, ast_sha = exact_mode_function(patched, str(path))
        namespace = dict(np=np, conditional_terms=likelihood.conditional_terms)
        exec(code, namespace)
        mode, scale, stats, hstats = namespace['fixed_mode'](drawn_model, p, v, rows)
        roots, weights = np.polynomial.hermite.hermgauss(20)
        points = mode+np.sqrt(2)*scale*roots
        log_integral = float(np.log(np.sqrt(2)*scale)+logsumexp(
            np.log(weights)+roots*roots+np.asarray([hstats(u)[0] for u in points])))
        dominant = max((e for e in reference['bracketed_extrema'] if e['classification'] == 'maximum'),
                       key=lambda e: e['log_integrand'])
        precheck = dict(mode=mode, scale=scale, log_integrand=stats[0], score=stats[1], curvature=stats[2],
                        log_gene_integral=log_integral, reference_log_integral=reference['log_bounded_integral'],
                        absolute_log_integral_difference=abs(log_integral-reference['log_bounded_integral']),
                        reference_mode=dominant['u'], absolute_mode_difference=abs(mode-dominant['u']),
                        extracted_mode_ast_sha256=ast_sha, synthetic_rows=len(drawn),
                        count_vector_sha256=count_sha,
                        fitting_performed=False)
        save('fixed_gene_precheck', precheck)
        require(np.isfinite(stats).all() and np.isfinite(log_integral)
                and abs(stats[1]) <= 1e-6 and stats[2] < 0
                and abs(mode-dominant['u']) <= 1e-6
                and abs(log_integral-reference['log_bounded_integral']) <= 1e-6,
                '64-damping fixed-gene mode/integral guard failed; no fitting authorized by this helper')
        # Directly compare original30 arithmetic/source against private64 plus
        # stable BB at historical observed parameters, before either refit.
        old_spec = importlib.util.spec_from_file_location(
            'immutable_damping30_observed_check', ROLE / 'executed_sources/likelihood_v1.py')
        old_module = importlib.util.module_from_spec(old_spec)
        sys.modules[old_spec.name] = old_module
        old_spec.loader.exec_module(old_module)
        old_model = old_module.Model(model.rows)
        observed_equivalence = {}
        for name in ('observed_null', 'observed_marginal'):
            archived = json.loads((ROLE / (name+'.json')).read_text())
            observed_p = np.asarray(archived['parameters'])
            before, after = old_model.evaluate(observed_p), model.evaluate(observed_p)
            differences = [float(np.max(np.abs(a-b))) for a, b in zip(before, after)]
            require(all(np.isfinite(d) and d <= 1e-6 for d in differences),
                    'Original30/private64 observed LL/score/Hessian differs')
            observed_equivalence[name] = differences
            if name == 'observed_marginal':
                before = old_model.marginal(observed_p, archived['v'])
                after = model.marginal(observed_p, archived['v'])
                differences = [float(np.max(np.abs(a-b))) for a, b in zip(before, after)]
                require(all(np.isfinite(d) and d <= 1e-6 for d in differences),
                        'Original30/private64 observed marginal LL/score differs')
                observed_equivalence[name+'_marginal'] = differences
        save('original30_private64_observed_equivalence', observed_equivalence)
        # Invoke the immutable existing single-draw replay protocol, replacing
        # only its source location/hash and output parent with our private copy.
        for name in ('bfix_v1', 'design_v1', 'likelihood_v1', 'tests_v1'):
            sys.modules.pop(name, None)
        driver = types.ModuleType('private_damping64_replay_driver')
        driver.__file__ = str(private / 'codex_replay_first_p3_stable_bb.py')
        exec(compile(Path(driver.__file__).read_bytes(), driver.__file__, 'exec'), driver.__dict__)
        driver.SNAPSHOT = private
        driver.HASHES = {**driver.HASHES, 'likelihood_v1.py': sha(path)}
        driver.REC = output
        phase = 'unchanged_observed_prefix_and_only_first_draw_replay'
        driver.main()
        phase = 'replayed_count_vector_identity_postcheck'
        replay_output = output / ('codex_replay_first_p3_stable_bb_'+job)
        replay_generation = json.loads((replay_output / 'generation.json').read_text())
        require(replay_generation['count_vector_sha256'] == count_sha,
                'Refitted-null generation changed the realized supplied draw; no original-draw agreement claimed')
        for name, digest in receipt['source_sha256'].items():
            require(sha(ROLE / 'executed_sources' / name) == digest, 'Canonical snapshot changed')
        save('summary', dict(fixed_gene_guard_passed=True,
                             replay_output=str(replay_output), realized_original_draw_identity_confirmed=True,
                             owner_spec_amended=False, statistical_law_changed=False,
                             limitation='One frozen precheck and one synthetic draw; no global mode guarantee, full fixture or calibration pass'))
    except Exception:
        save('failure', dict(phase=phase, traceback=traceback.format_exc(),
                             additional_method_tuning=False, source_law_unchanged=True))
        raise


if __name__ == '__main__':
    main()
