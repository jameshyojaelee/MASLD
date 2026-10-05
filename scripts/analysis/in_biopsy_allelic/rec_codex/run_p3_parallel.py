"""Execute the original P3 draw body in independent, source-pinned processes."""
import argparse
import ast
import concurrent.futures
import hashlib
import json
import math
import multiprocessing
import os
from pathlib import Path
import sys
import time


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def initialize(campaign):
    global context, uniforms, output
    root = Path(campaign)
    manifest = json.loads((root / 'campaign.json').read_text())
    for name, expected in manifest['files'].items():
        if digest(Path(name)) != expected:
            raise ValueError('Pinned input changed: ' + name)
    sys.path.insert(0, str(root / 'executed_sources'))
    import numpy as np
    from bfix_v1 import load, table, FIXTURE
    from likelihood_v1 import Model, conditional_draw
    design, settings = load()
    model = design.model()
    assert settings['B'] == 50
    ordered = [(r['individual_id'], r['gene_id']) for r in table('rows.tsv')]
    assert ordered == [(r.individual, settings['gene_order'][r.gene]) for r in model.rows]
    uniforms = np.load(FIXTURE / 'draws/p3_uniforms.npy')
    assert uniforms.shape == (50, len(model.rows), 2)
    null = json.loads((root / 'baseline/observed_null.json').read_text())
    marginal = json.loads((root / 'baseline/observed_marginal.json').read_text())
    assert null['converged'] and marginal['converged']
    null_p = np.asarray(null['parameters'])
    observed_lr = 0. if marginal['v'] == 0. else max(0., 2 * (marginal['loglik'] - null['loglik']))
    original = ast.parse((root / 'executed_sources/check_bfix_tests_v13.py').read_text())
    loops = [node for node in original.body if isinstance(node, ast.For)
             and isinstance(node.target, ast.Tuple)
             and [x.id for x in node.target.elts] == ['b', 'draw']]
    assert len(loops) == 1
    function = ast.FunctionDef(name='original_draw',
        args=ast.arguments(posonlyargs=[], args=[ast.arg(arg='b'), ast.arg(arg='draw')],
                           kwonlyargs=[], kw_defaults=[], defaults=[]),
        body=loops[0].body, decorator_list=[])
    module = ast.fix_missing_locations(ast.Module(body=[function], type_ignores=[]))
    assert ast.dump(module.body[0].body[0], include_attributes=False) == ast.dump(loops[0].body[0], include_attributes=False)
    context = dict(np=np, Model=Model, conditional_draw=conditional_draw, model=model,
                   settings=settings, null_p=null_p, observed_lr=observed_lr)
    exec(compile(module, 'unchanged_original_p3_draw', 'exec'), context)
    output = root / 'draws'


def run_draw(b):
    started = time.monotonic()
    def serial(value):
        import numpy as np
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(type(value).__name__)
    def save(name, value):
        assert name == f'p3_{b:03d}'
        temporary = output / f'.{name}.{os.getpid()}.tmp'
        with temporary.open('x') as handle:
            handle.write(json.dumps(value, default=serial, indent=2) + '\n')
        destination = output / (name + '.json')
        # A hard link publishes a complete record without replacing any result.
        os.link(temporary, destination)
        temporary.unlink()
    context.update(save=save, p3_failed=[], p3_lr=[], p3_exceed=[])
    context['original_draw'](b, uniforms[b])
    timing = {'draw': b, 'elapsed_seconds': time.monotonic() - started,
              'job_id': os.environ['SLURM_JOB_ID']}
    with (output / f'p3_{b:03d}_timing.json').open('x') as handle:
        json.dump(timing, handle, indent=2)
    return timing


def collect(root):
    manifest = json.loads((root / 'campaign.json').read_text())
    for name, expected in manifest['files'].items():
        assert digest(Path(name)) == expected, name
    replay = json.loads((root / 'draws/p3_019.json').read_text())
    original = json.loads((root / 'baseline/p3_019.json').read_text())
    # Tolerance fixed before worker execution; no selecting among duplicate fits.
    checks = {'LR': abs(replay['LR'] - original['LR']) <= 1e-6,
              'exceed': replay['exceed'] == original['exceed']}
    for key in ['v', 'loglik']:
        checks[key] = abs(replay['fit'][key] - original['fit'][key]) <= 1e-6
    checks['converged'] = replay['fit']['converged'] == original['fit']['converged']
    checks['null_converged'] = replay['fit']['null']['converged'] == original['fit']['null']['converged']
    checks['kept_start'] = replay['fit']['kept_start'] == original['fit']['kept_start']
    def agree(a, b):
        if isinstance(a, list) or isinstance(b, list):
            return isinstance(a, list) and isinstance(b, list) and len(a) == len(b) and all(agree(x, y) for x, y in zip(a, b))
        if a is None or b is None:
            return a is b
        if math.isnan(a) or math.isnan(b):
            return math.isnan(a) and math.isnan(b)
        return abs(a - b) <= 1e-6
    checks['parameters'] = agree(replay['fit']['parameters'], original['fit']['parameters'])
    checks['null_parameters'] = agree(replay['fit']['null']['parameters'], original['fit']['null']['parameters'])
    checks['null_loglik'] = agree(replay['fit']['null']['loglik'], original['fit']['null']['loglik'])
    assert all(checks.values()), checks
    records = [json.loads((root / ('baseline' if b < 20 else 'draws') / f'p3_{b:03d}.json').read_text())
               for b in range(50)]
    assert [x['b'] for x in records] == list(range(50))
    failed = [x['b'] for x in records if x['fit'] is None or not
              (x['fit']['converged'] and x['fit']['null']['converged'])]
    exceed = [bool(x['exceed']) for x in records]
    p = (1 + sum(exceed)) / 51
    result = {'B': 50, 'p': p, 'mc_se': (p * (1 - p) / 50) ** .5,
              'failed': failed, 'draw_LR': [x['LR'] for x in records],
              'exceedance': exceed, 'replay_checks': checks,
              'assignment': 'original 0..19; new fixed workers 20..49',
              'note': 'Synthetic agreement fixture, not calibration.'}
    with (root / 'parallel_p3_result.json').open('x') as handle:
        json.dump(result, handle, indent=2)
    print(json.dumps(result), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign', type=Path)
    parser.add_argument('--draws', nargs='+', type=int)
    parser.add_argument('--collect', action='store_true')
    args = parser.parse_args()
    if args.collect:
        collect(args.campaign)
        return
    assert args.draws and len(set(args.draws)) == len(args.draws)
    manifest = json.loads((args.campaign / 'campaign.json').read_text())
    assert args.draws in manifest['worker_groups']
    workers = len(args.draws)
    assert workers <= int(os.environ['SLURM_CPUS_PER_TASK']) <= 16
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers,
            mp_context=multiprocessing.get_context('spawn'),
            initializer=initialize, initargs=(str(args.campaign),)) as pool:
        futures = [pool.submit(run_draw, b) for b in args.draws]
        for future in concurrent.futures.as_completed(futures):
            print(json.dumps(future.result()), flush=True)


if __name__ == '__main__':
    main()
