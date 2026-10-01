"""Run the unchanged supplied B-FIX P2 or P3 block in its own process.

This only isolates execution. It does not generate draws, change statistical
rules, read planted classes, or produce a scientific calibration result.
"""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
FIXTURE = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
RUNNER = 'check_bfix_tests_v13.py'
MODULES = ('bfix_v1.py', 'design_v1.py', 'likelihood_v1.py', 'tests_v1.py')
INPUTS = ('globals.json', 'individuals.tsv', 'rows.tsv', 'cohorts.tsv',
          'draws/lambda_draws.tsv', 'draws/p2_resample.tsv', 'draws/p3_uniforms.npy')
PRIOR = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/bfix_fits_21984981/null.json'
OBSERVED = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/bfix_tests_21985999'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def assignment(node, name):
    return isinstance(node, ast.Assign) and any(
        isinstance(target, ast.Name) and target.id == name for target in node.targets)


def source_block(source, nodes):
    # Source segments honor columns, including several statements on one line.
    text = ''.join(ast.get_source_segment(source, node) + '\n' for node in nodes)
    original = ast.dump(ast.Module(body=nodes, type_ignores=[]), include_attributes=False)
    if ast.dump(ast.parse(text), include_attributes=False) != original:
        raise ValueError('isolated statement syntax differs from original')
    return text


def isolate(source, role):
    """Keep original statement text and ordering; omit the other two roles."""
    body = ast.parse(source).body
    positions = {}
    for name in ('se', 'source_names', 'ordered'):
        matches = [i for i, node in enumerate(body) if assignment(node, name)]
        if len(matches) != 1:
            raise ValueError('original runner boundary changed: ' + name)
        positions[name] = matches[0]
    if not positions['se'] < positions['source_names'] < positions['ordered']:
        raise ValueError('original role ordering changed')
    final = body[-1]
    if not (isinstance(final, ast.Expr) and isinstance(final.value, ast.Call)
            and isinstance(final.value.func, ast.Name) and final.value.func.id == 'save'
            and isinstance(final.value.args[0], ast.Constant)
            and final.value.args[0].value == 'result'):
        raise ValueError('original runner final summary changed')
    prefix = body[:positions['se']]
    output_nodes = [node for node in prefix if assignment(node, 'output')]
    mkdir_nodes = [node for node in prefix if isinstance(node, ast.Expr)
                   and isinstance(node.value, ast.Call)
                   and isinstance(node.value.func, ast.Attribute)
                   and isinstance(node.value.func.value, ast.Name)
                   and node.value.func.value.id == 'output'
                   and node.value.func.attr == 'mkdir']
    if len(output_nodes) != 1 or len(mkdir_nodes) != 1:
        raise ValueError('original output setup changed')
    omitted = output_nodes + mkdir_nodes
    prefix = [node for node in prefix if node not in omitted]
    begin = positions['source_names'] if role == 'p2' else positions['ordered']
    end = positions['ordered'] if role == 'p2' else len(body) - 1
    selected = body[begin:end]
    # Draw/statistic expressions and statement order are unchanged, checked by AST.
    return source_block(source, prefix), source_block(source, selected), {
        'selected_lines': [selected[0].lineno, selected[-1].end_lineno],
        'prefix_statement_count': len(prefix), 'role_statement_count': len(selected)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--role', required=True, choices=('p2', 'p3'))
    parser.add_argument('--inspect', action='store_true', help='static inspection only; no fits/output directory')
    args = parser.parse_args()
    source_files = {name: (HERE / name).read_bytes() for name in (RUNNER,) + MODULES}
    source_files[Path(__file__).name] = Path(__file__).read_bytes()
    source = source_files[RUNNER].decode()
    prefix, selected, structure = isolate(source, args.role)
    published = {}
    for line in (FIXTURE / 'SHA256SUMS').read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        published[name.removeprefix('./')] = expected
    inputs = {name: digest((FIXTURE / name).read_bytes()) for name in INPUTS}
    for name, actual in inputs.items():
        if actual != published[name]:
            raise ValueError('synthetic input checksum changed: ' + name)
    settings = json.loads((FIXTURE / 'globals.json').read_text())
    if settings['B'] != 50:
        raise ValueError('exactly fifty supplied draws required')
    reference_hashes = {name: digest((OBSERVED / (name + '.json')).read_bytes())
                        for name in ('observed_null', 'observed_marginal')}
    receipt = {'role': args.role, 'spec': 'v1.7', 'B': 50,
               'fixture': str(FIXTURE), 'fixture_sha256': inputs,
               'prior_null': str(PRIOR), 'prior_null_sha256': digest(PRIOR.read_bytes()),
               'observed_reference': str(OBSERVED), 'observed_reference_sha256': reference_hashes,
               'source_sha256': {name: digest(data) for name, data in source_files.items()},
               **structure, 'note': 'Synthetic agreement fixture, not calibration.'}
    job = 'STATIC' if args.inspect else os.environ['SLURM_JOB_ID']
    if not args.inspect and not job.isdigit():
        raise ValueError('numeric SLURM_JOB_ID required')
    base = Path(os.environ.get('CODEX_REC_OUTPUT', str(PRIOR.parent.parent)))
    output = base / ('bfix_' + args.role + '_' + job)
    snapshot = output / 'executed_sources'
    comparison = ("reference = Path(" + repr(str(OBSERVED)) + ")\n"
                  "comparison = {}\n"
                  "for name, fitted in [('observed_null', null), ('observed_marginal', marginal)]:\n"
                  "    previous = json.loads((reference / (name + '.json')).read_text())\n"
                  "    comparison[name] = {'loglik_difference':abs(fitted['loglik']-previous['loglik']),\n"
                  "                        'kappa_difference':abs(fitted['parameters'][2*model.G]-previous['parameters'][2*model.G])}\n"
                  "    if name == 'observed_marginal':\n"
                  "        comparison[name]['v_difference'] = abs(fitted['v']-previous['v'])\n"
                  "save('observed_comparison', comparison)\n"
                  "if any(value > 1e-6 for differences in comparison.values() for value in differences.values()):\n"
                  "    raise ValueError('observed fits differ from live-original baseline by >1e-6')\n")
    code = ('output = ' + repr(str(output)) + '\nfrom pathlib import Path\n'
            'output = Path(output)\n' + prefix + '\n' + comparison + selected + '\n'
            "save('result', {'fixture':str(FIXTURE),'spec':'v1.7','B':B,"
            "'role':" + repr(args.role) + "," + repr(args.role.upper()) + ':' + args.role + ','
            "'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,"
            "'elapsed_seconds':time.monotonic()-start,"
            "'note':'Synthetic agreement fixture, not calibration.'})\n")
    compile(code, str(output / 'executed_role.py'), 'exec')
    receipt['executed_role_sha256'] = digest(code.encode())
    if args.inspect:
        print(json.dumps(receipt, indent=2))
        return
    output.mkdir(exist_ok=False)
    snapshot.mkdir()
    for name, data in source_files.items():
        (snapshot / name).write_bytes(data)
    (output / 'executed_role.py').write_text(code)
    (output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    for name, data in source_files.items():
        if (HERE / name).read_bytes() != data:
            raise RuntimeError('source changed during snapshot: ' + name)
    for name, actual in inputs.items():
        if digest((FIXTURE / name).read_bytes()) != actual:
            raise RuntimeError('synthetic input changed during snapshot: ' + name)
    if digest(PRIOR.read_bytes()) != receipt['prior_null_sha256']:
        raise RuntimeError('prior null changed during snapshot')
    for name, actual in reference_hashes.items():
        if digest((OBSERVED / (name + '.json')).read_bytes()) != actual:
            raise RuntimeError('observed reference changed during snapshot: ' + name)
    # All local numerical modules import from this archived snapshot.
    if any(Path(name).stem in sys.modules for name in MODULES):
        raise RuntimeError('numerical modules were imported before snapshot isolation')
    sys.path.insert(0, str(snapshot))
    exec(compile(code, str(output / 'executed_role.py'), 'exec'), {'__name__': '__main__'})


if __name__ == '__main__':
    main()
