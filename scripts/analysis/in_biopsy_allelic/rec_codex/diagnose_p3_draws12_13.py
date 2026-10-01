#!/usr/bin/env python3
"""Retain exact archived synthetic P3 exceptions for supplied draws12 and13."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import time
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
ROLE = REC / 'bfix_p3_21998519'
ARCHIVE = ROLE / 'executed_sources'
FIXTURE = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
MODULES = ('bfix_v1', 'design_v1', 'likelihood_v1', 'tests_v1')
CAP = 32 * 1024 * 1024
FIXED = {
    'receipt.json': '3a3a8c8f714ee69989409049a337485d61100b65b433701dacab4b67a60acd4a',
    'observed_null.json': 'a617d331fd3fb93576e5bd8588581aaab8f77946936536401905490438660181',
    'p3_012.json': '92e214b618f3e03e626cba50c48c8fa1e551de3f5f9f0dd07a6b497b81d20bf0',
    'p3_013.json': '85767c9d280f255fc653ae7baf262ace5ced06d9943557a1852001493c84f284',
}


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def metadata_guards():
    guards = {ROLE / name: expected for name, expected in FIXED.items()}
    require(all(digest(p) == h for p, h in guards.items()), 'fixed receipt/null/draw record differs')
    receipt = json.loads((ROLE / 'receipt.json').read_text())
    require(receipt['spec'] == 'v1.7' and receipt['role'] == 'p3' and receipt['B'] == 50
            and Path(receipt['fixture']).resolve() == FIXTURE.resolve(), 'wrong original role/fixture')
    require(receipt['source_sha256']['likelihood_v1.py'] == '41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6', 'wrong statistical source')
    guards.update({ARCHIVE / n: h for n, h in receipt['source_sha256'].items()})
    guards.update({FIXTURE / n: h for n, h in receipt['fixture_sha256'].items()})
    guards[ROLE / 'executed_role.py'] = receipt['executed_role_sha256']
    for path, expected in guards.items():
        require(path.resolve().is_relative_to(ROLE.resolve()) or path.resolve().is_relative_to(FIXTURE.resolve()), 'guard escaped synthetic roots')
        require(digest(path) == expected, 'source/fixture mismatch: ' + str(path))
    for b in (12, 13):
        record = json.loads((ROLE / f'p3_{b:03d}.json').read_text())
        require(record['b'] == b and record['fit'] is None and record['exceed'] is True, 'original exception record differs')
    return receipt, guards


def completion_guards(guards):
    result = json.loads((ROLE / 'result.json').read_text())
    p3 = json.loads((ROLE / 'p3_result.json').read_text())
    require(result['role'] == 'p3' and result['spec'] == 'v1.7' and result['B'] == 50
            and Path(result['fixture']).resolve() == FIXTURE.resolve(), 'original full role result missing/wrong')
    require(p3['computed'] is True and len(p3['draw_LR']) == 50
            and result['P3']['computed'] is True and len(result['P3']['draw_LR']) == 50
            and result['P3']['failed'] == p3['failed'] and {12, 13} <= set(p3['failed']), 'full50 P3 completion contract differs')
    for b in range(50):
        p = ROLE / f'p3_{b:03d}.json'
        require(json.loads(p.read_text())['b'] == b, 'original full50 draw sequence incomplete')
        guards[p] = digest(p)
    for name in ('result.json', 'p3_result.json'):
        guards[ROLE / name] = digest(ROLE / name)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true', help='metadata/hash inspection only; no numerical imports/arrays/fits')
    args = parser.parse_args()
    receipt, guards = metadata_guards()
    own = [Path(__file__), Path(__file__).with_name('run_diagnose_p3_draws12_13.sbatch')]
    initial = {p: digest(p) for p in own}
    record = dict(draw_indices=[12, 13], zero_based_original_indices=True, rows=173,
                  source_role=str(ROLE), input_sha256={str(p): h for p, h in guards.items()},
                  own_sha256={str(p): h for p, h in initial.items()},
                  fit_call='Model(rows).fit_marginal(tau_max_fixture, saved_null_parameters, section4_from_null=True)',
                  random_source='only supplied p3_uniforms[12] and [13]; no new seed/random draw',
                  output_cap_bytes=CAP, full50_completion_required_before_execution=True,
                  original_environment_versions_required=['python', 'numpy', 'scipy'],
                  interpretation='synthetic exception diagnostic only; no calibration/model change/adoption')
    if args.inspect:
        print(json.dumps(record, indent=2))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'compute-only numeric SLURM job required')
    out = REC / ('p3_draws12_13_diagnostic_' + os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)
    for p in own:
        shutil.copyfile(p, out / p.name)
    (out / 'receipt.json').write_text(json.dumps(record, indent=2) + '\n')
    phase = 'full50_completion'
    try:
        original_result = completion_guards(guards)
        record['input_sha256'] = {str(p): h for p, h in guards.items()}
        (out / 'completed_original_role_receipt.json').write_text(json.dumps(record, indent=2) + '\n')
        require(not any(name in sys.modules for name in MODULES), 'numerical module imported before source isolation')
        phase = 'archived_imports'
        sys.path.insert(0, str(ARCHIVE))
        import numpy as np
        import scipy
        from bfix_v1 import load, table, FIXTURE as imported_fixture
        from likelihood_v1 import Model, conditional_draw
        import tests_v1
        require(imported_fixture.resolve() == FIXTURE.resolve(), 'archived loader fixture differs')
        for name in MODULES:
            path = Path(sys.modules[name].__file__).resolve()
            require(path.parent == ARCHIVE.resolve() and digest(path) == receipt['source_sha256'][name + '.py'], 'archived import mismatch: ' + name)
        phase = 'original_environment_comparison'
        versions = dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)
        for name, value in versions.items():
            require(value == original_result[name], 'original environment version differs: ' + name)
        def serial(value):
            if isinstance(value, np.ndarray):
                return value.tolist()
            if isinstance(value, np.generic):
                return value.item()
            raise TypeError(type(value).__name__)
        def room(text):
            size = sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
            require(size + len(text.encode()) <= CAP - 65536, 'diagnostic output cap reached; failure reserve retained')
        def save(path, value):
            text = json.dumps(value, default=serial, indent=2) + '\n'
            room(text)
            with path.open('x') as handle:
                handle.write(text)
        save(out / 'environment.json', dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
             imported_paths={n: sys.modules[n].__file__ for n in MODULES}, threads={k: os.environ.get(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')}))
        phase = 'synthetic_setup'
        design, settings = load()
        model = design.model()
        require(settings['B'] == 50 and len(model.rows) == 173, 'synthetic dimensions differ')
        observed = json.loads((ROLE / 'observed_null.json').read_text())
        require(observed['converged'], 'saved observed null not converged')
        null_p = np.asarray(observed['parameters'])
        ordered = [(r['individual_id'], r['gene_id']) for r in table('rows.tsv')]
        require(ordered == [(r.individual, settings['gene_order'][r.gene]) for r in model.rows], 'original173 row order differs')
        uniforms = np.load(FIXTURE / 'draws/p3_uniforms.npy', allow_pickle=False)
        require(uniforms.shape == (50, 173, 2), 'supplied uniform shape differs')
        save(out / 'setup.json', dict(tau_max_fixture=settings['tau_max_fixture'], saved_null_parameters=null_p,
             row_order=ordered, supplied_uniform_shape=list(uniforms.shape), extra_fit_kwargs=dict(section4_from_null=True)))
        methods = {getattr(Model, n).__code__: n for n in ('fit','fit_single','fit_marginal','fit_marginal_single','marginal')}
        outcomes = []
        for b in (12, 13):
            directory = out / f'draw_{b:03d}'
            directory.mkdir()
            started = time.monotonic()
            latest = None
            current_method = None
            marginal_calls = 0
            row_index = None
            progress = (directory / 'progress.jsonl').open('x')
            def event(name, **payload):
                text = json.dumps(dict(event=name, elapsed_seconds=time.monotonic()-started, **payload), default=serial) + '\n'
                room(text)
                progress.write(text)
                progress.flush()
            def profile(frame, action, argument):
                nonlocal latest, current_method, marginal_calls
                if action != 'call' or frame.f_code not in methods:
                    return
                current_method = methods[frame.f_code]
                if current_method == 'marginal':
                    marginal_calls += 1
                    latest = dict(call=marginal_calls, v=float(frame.f_locals['v']), parameters=np.asarray(frame.f_locals['p']).copy())
                else:
                    values = dict(method=current_method)
                    if frame.f_locals.get('start') is not None:
                        values['start'] = np.asarray(frame.f_locals['start']).copy()
                    if 'fixed_v' in frame.f_locals:
                        values['fixed_v'] = frame.f_locals['fixed_v']
                    event('method_entry', **values)
            fit = None
            try:
                phase = 'conditional_draw_generation'
                event('phase', phase=phase, draw=b)
                rows = []
                for row_index, (row, (u1, u2)) in enumerate(zip(model.rows, uniforms[b])):
                    current_method = 'predictors'
                    eta, r, _, _ = model.predictors(row, null_p)
                    current_method = 'conditional_draw'
                    rows.append(conditional_draw(row, eta, r, u1, u2))
                phase = 'Model_rows_construction'
                current_method = 'Model.__init__'
                drawn_model = Model(rows)
                phase = 'fit_marginal'
                event('phase', phase=phase, draw=b)
                sys.setprofile(profile)
                try:
                    fit = drawn_model.fit_marginal(settings['tau_max_fixture'], null_p, section4_from_null=True)
                finally:
                    sys.setprofile(None)
                okay = fit['converged'] and fit['null']['converged']
                lr = 0. if fit['v'] == 0. else max(0., 2*(fit['loglik']-fit['null']['loglik']))
                save(directory / 'fit.json', fit)
                outcome = dict(draw=b, exception=False, converged=bool(okay), LR=lr if okay else None)
            except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
                sys.setprofile(None)
                tb = traceback.format_exc()
                frames = []
                trace = exc.__traceback__
                while trace is not None:
                    frame = trace.tb_frame
                    if Path(frame.f_code.co_filename).resolve().parent == ARCHIVE.resolve():
                        values = {}
                        for name in ('p','start','v','gene','mode','curvature','nodes','with_hessian'):
                            v = frame.f_locals.get(name)
                            if isinstance(v, (bool,int,float,str,np.generic,np.ndarray)):
                                values[name] = v.copy() if isinstance(v,np.ndarray) else v
                        frames.append(dict(method=frame.f_code.co_name, line=trace.tb_lineno, locals=values))
                    trace = trace.tb_next
                room(tb)
                (directory / 'traceback.txt').write_text(tb)
                outcome = dict(draw=b, exception=True, exception_type=type(exc).__name__, message=str(exc),
                               phase=phase, row_index=row_index, current_method=current_method,
                               latest_marginal=latest, archived_frames=frames, traceback=tb)
            finally:
                sys.setprofile(None)
                progress.close()
            outcome.update(marginal_calls=marginal_calls, elapsed_seconds=time.monotonic()-started,
                           saved_null_start_parameters=null_p, tau_max_fixture=settings['tau_max_fixture'],
                           interpretation='exact synthetic exception diagnostic; original draw/statistical law untouched')
            save(directory / 'summary.json', outcome)
            outcomes.append(dict(draw=b, exception=outcome['exception']))
        for p, h in {**guards, **initial}.items():
            require(digest(p) == h, 'immutable source/input changed during diagnostic: ' + str(p))
        save(out / 'summary.json', dict(status='both_exact_draws_diagnosed', draws=outcomes,
             original_role_mutated=False, model_changed=False, calibration_claim=False))
        (out / 'COMPLETE').write_text('diagnostic complete\n')
    except BaseException as exc:
        sys.setprofile(None)
        (out / 'failure.json').write_text(json.dumps(dict(phase=phase, type=type(exc).__name__, message=str(exc),
             traceback=traceback.format_exc(), partial_outputs_retained=True), indent=2) + '\n')
        raise


if __name__ == '__main__':
    main()
