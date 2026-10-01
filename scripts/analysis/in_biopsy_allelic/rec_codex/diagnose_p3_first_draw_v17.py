"""Reproduce only supplied synthetic P3 draw 0; retain the actual exception."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time
import traceback


ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
ROLE = REC / 'bfix_p3_21998519'
ARCHIVE = ROLE / 'executed_sources'
FIXTURE = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
MODULES = ('bfix_v1', 'design_v1', 'likelihood_v1', 'tests_v1')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true', help='hash/syntax inspection only; no arrays or fits')
    args = parser.parse_args()
    receipt = json.loads((ROLE / 'receipt.json').read_text())
    if receipt['spec'] != 'v1.7' or receipt['role'] != 'p3' or receipt['B'] != 50:
        raise ValueError('wrong source role receipt')
    if Path(receipt['fixture']).resolve() != FIXTURE.resolve():
        raise ValueError('wrong fixed synthetic fixture')
    for name, expected in receipt['source_sha256'].items():
        if digest(ARCHIVE / name) != expected:
            raise ValueError('archived source hash differs: ' + name)
    for name, expected in receipt['fixture_sha256'].items():
        if digest(FIXTURE / name) != expected:
            raise ValueError('fixture hash differs: ' + name)
    if digest(ROLE / 'executed_role.py') != receipt['executed_role_sha256']:
        raise ValueError('archived role source hash differs')
    own_source = Path(__file__).read_bytes()
    compile(own_source, str(Path(__file__)), 'exec')
    diagnostic_receipt = {
        'source_role': str(ROLE), 'archive': str(ARCHIVE), 'draw': 0,
        'role_receipt_sha256': digest(ROLE / 'receipt.json'),
        'observed_null_sha256': digest(ROLE / 'observed_null.json'),
        'diagnostic_sha256': hashlib.sha256(own_source).hexdigest(),
        'source_sha256': receipt['source_sha256'],
        'fixture_sha256': receipt['fixture_sha256'],
        'fit_kwargs': {'section4_from_null': True},
        'note': 'Synthetic exception diagnostic; original P3 draw remains unchanged.'}
    if args.inspect:
        print(json.dumps(diagnostic_receipt, indent=2))
        return
    job = os.environ['SLURM_JOB_ID']
    if not job.isdigit():
        raise ValueError('numeric SLURM_JOB_ID required')
    output = Path(os.environ['CODEX_REC_OUTPUT']) / ('p3_first_draw_diagnostic_' + job)
    output.mkdir(exist_ok=False)
    (output / 'executed_diagnostic.py').write_bytes(own_source)
    (output / 'receipt.json').write_text(json.dumps(diagnostic_receipt, indent=2) + '\n')
    started = time.monotonic()
    progress = (output / 'progress.jsonl').open('x')
    phase = 'imports'
    fit = None
    marginal_calls = 0
    latest_marginal = None

    def serial(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(type(value).__name__)

    def event(name, **values):
        progress.write(json.dumps({'event': name, 'elapsed_seconds': time.monotonic()-started,
                                   **values}, default=serial) + '\n')
        progress.flush()
        print(name, values.get('phase', ''), flush=True)

    def save(name, value):
        (output / (name + '.json')).write_text(json.dumps(value, default=serial, indent=2) + '\n')

    try:
        if any(name in sys.modules for name in MODULES):
            raise RuntimeError('numerical module imported before archived-source isolation')
        sys.path.insert(0, str(ARCHIVE))
        import numpy as np
        import scipy
        from bfix_v1 import load, table
        from likelihood_v1 import Model, conditional_draw
        for name in MODULES:
            module_path = Path(sys.modules[name].__file__).resolve()
            if module_path.parent != ARCHIVE.resolve():
                raise RuntimeError('numerical import escaped archived source: ' + name)
            if digest(module_path) != receipt['source_sha256'][name + '.py']:
                raise RuntimeError('archived source changed during import: ' + name)
        event('imports_verified', paths={name: sys.modules[name].__file__ for name in MODULES},
              python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)
        phase = 'load_synthetic_design'
        event('phase', phase=phase)
        design, settings = load()
        model = design.model()
        if settings['B'] != 50:
            raise ValueError('fixture B changed')
        observed = json.loads((ROLE / 'observed_null.json').read_text())
        if digest(ROLE / 'observed_null.json') != diagnostic_receipt['observed_null_sha256']:
            raise RuntimeError('saved observed null changed')
        if not observed['converged']:
            raise ValueError('saved observed null did not converge')
        null_p = np.asarray(observed['parameters'])
        ordered = [(r['individual_id'], r['gene_id']) for r in table('rows.tsv')]
        if ordered != [(r.individual, settings['gene_order'][r.gene]) for r in model.rows]:
            raise ValueError('row ordering differs from original P3 runner')
        uniforms = np.load(FIXTURE / 'draws/p3_uniforms.npy')
        if uniforms.shape != (50, len(model.rows), 2):
            raise ValueError('uniform shape differs from original P3 runner')
        phase = 'conditional_draw_generation'
        event('phase', phase=phase, rows=len(model.rows))
        rows = []
        for index, (row, (u1, u2)) in enumerate(zip(model.rows, uniforms[0])):
            event('row_start', index=index, individual=row.individual, gene=row.gene)
            eta, r, _, _ = model.predictors(row, null_p)
            rows.append(conditional_draw(row, eta, r, u1, u2))
            event('row_complete', index=index)
        phase = 'Model_rows_construction'
        event('phase', phase=phase, rows=len(rows))
        drawn_model = Model(rows)
        phase = 'fit_marginal'
        event('phase', phase=phase)
        methods = {getattr(Model, name).__code__: name for name in
                   ('fit', 'fit_single', 'fit_marginal', 'fit_marginal_single', 'marginal')}

        def profile(frame, action, argument):
            nonlocal marginal_calls, latest_marginal
            if action != 'call' or frame.f_code not in methods:
                return
            method = methods[frame.f_code]
            if method == 'marginal':
                marginal_calls += 1
                latest_marginal = {'call': marginal_calls, 'v': float(frame.f_locals['v']),
                                   'parameters': np.asarray(frame.f_locals['p']).copy()}
                return
            payload = {'method': method}
            if 'start' in frame.f_locals and frame.f_locals['start'] is not None:
                payload['start'] = np.asarray(frame.f_locals['start']).copy()
            if 'fixed_v' in frame.f_locals:
                payload['fixed_v'] = frame.f_locals['fixed_v']
            event('method_entry', **payload)

        sys.setprofile(profile)
        try:
            fit = drawn_model.fit_marginal(settings['tau_max_fixture'], null_p,
                                          section4_from_null=True)
        finally:
            sys.setprofile(None)
        save('fit', fit)
        okay = fit['converged'] and fit['null']['converged']
        lr = 0. if fit['v'] == 0. else max(0., 2*(fit['loglik']-fit['null']['loglik']))
        save('summary', {'draw': 0, 'exception': False, 'converged': bool(okay),
                         'LR': lr, 'marginal_calls': marginal_calls,
                         'elapsed_seconds': time.monotonic()-started,
                         'note': 'Diagnostic only; original P3 draw unchanged.'})
        event('complete', phase='finished', converged=bool(okay))
    except Exception as error:
        sys.setprofile(None)
        full_traceback = traceback.format_exc()
        (output / 'traceback.txt').write_text(full_traceback)
        frames = []
        trace = error.__traceback__
        while trace is not None:
            frame = trace.tb_frame
            if Path(frame.f_code.co_filename).resolve().parent == ARCHIVE.resolve():
                values = {}
                for name in ('gene', 'v', 'mode', 'curvature', 'nodes', 'with_hessian'):
                    value = frame.f_locals.get(name)
                    if value is not None and isinstance(value, (bool, int, float, np.generic)):
                        values[name] = value
                frames.append({'method': frame.f_code.co_name, 'line': trace.tb_lineno,
                               'locals': values})
            trace = trace.tb_next
        # The exception itself, not a guessed numerical explanation, is retained.
        record = {'draw': 0, 'exception': True, 'phase': phase,
                  'exception_type': type(error).__name__, 'message': str(error),
                  'traceback': full_traceback, 'marginal_calls': marginal_calls,
                  'latest_marginal': latest_marginal, 'archived_exception_frames': frames,
                  'elapsed_seconds': time.monotonic()-started,
                  'note': 'Diagnostic only; original P3 draw unchanged.'}
        if 'np' in locals():
            save('summary', record)
        else:
            (output / 'summary.json').write_text(json.dumps(record, indent=2) + '\n')
        print(full_traceback, file=sys.stderr, flush=True)
        raise
    finally:
        progress.close()


if __name__ == '__main__':
    main()
