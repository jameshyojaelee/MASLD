#!/usr/bin/env python3
"""Same-process diagnosis of the first saved GL rule; no fits or acceptance."""
import argparse
import ast
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import resource
import sys
import time
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
HERE = ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
MODULES = ('likelihood_v1.py', 'design_v1.py', 'bfix_v1.py', 'tests_v1.py',
           'repair_p3_bb_mode_anchor.py', 'repair_p3_bb_mode_anchor_cache_v2.py')
POINT = REC/'repair_p3_observed_failure_anchor_22042612'
CACHE_POINT = REC/'repair_p3_bb_mode_anchor_cache_check_22044584'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pins', required=True, type=Path)
    parser.add_argument('--pins-sha', required=True)
    parser.add_argument('--self-sha', required=True)
    parser.add_argument('--inspect', action='store_true')
    args = parser.parse_args()
    require(sha(Path(__file__)) == args.self_sha and sha(args.pins) == args.pins_sha,
            'Reviewed diagnostic or manifest changed')
    pins = {Path(p): h for p, h in json.loads(args.pins.read_text()).items()}
    for path, digest in pins.items():
        require(path.is_absolute() and path.resolve().is_relative_to(ROOT) and sha(path) == digest,
                'Reviewed input changed: '+str(path))
        if path.suffix == '.py':
            ast.parse(path.read_text())
    require(all(HERE/n in pins for n in MODULES), 'Executable dependency missing')
    for directory, names in [(POINT, ('protocol.json', 'GL_candidate.json')),
                             (CACHE_POINT, ('GL_candidate.json', 'failure.json'))]:
        require(all(directory/n in pins for n in names), 'Completed result pin missing')
    protocol = dict(rule=dict(z_bound=8., grid_spacing=.25, order=16, nodes=1024),
        saved_gene=2, rows=24, fixture_seed=20260930, compare_full_and_gene_only_predictors=True,
        same_process_scalar_and_cache=True, raw_byte_comparison=True,
        maximum_saved_row_discrepancies=10, output_cap_bytes=32*1024**2,
        no_fit=True, no_accuracy_acceptance=True, no_calibration=True,
        source_and_input_sha256={str(p): h for p, h in pins.items()})
    if args.inspect:
        print(json.dumps(protocol))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and
            not os.environ.get('SLURM_ARRAY_TASK_ID'), 'Individual compute allocation required')
    output = REC/('repair_p3_cache_rule_diagnosis_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(exist_ok=False)
    archive = output/'executed_sources'
    archive.mkdir()
    for path in (Path(__file__), *(HERE/n for n in MODULES)):
        target = archive/path.name
        target.write_bytes(path.read_bytes())
        require(sha(target) == (args.self_sha if path == Path(__file__) else pins[path]),
                'Executable archive differs')
    began = time.perf_counter()
    try:
        import numpy as np
        import scipy

        def save(name, value):
            def serialize(item):
                if isinstance(item, np.ndarray):
                    return item.tolist()
                if isinstance(item, np.generic):
                    return item.item()
                raise TypeError(type(item).__name__)
            body = json.dumps(value, default=serialize, indent=2, allow_nan=False)+'\n'
            require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file())+len(body.encode())
                    < 32*1024**2-65536, 'Output capacity exceeded')
            with (output/(name+'.json')).open('x') as handle:
                handle.write(body)
                handle.flush()
                os.fsync(handle.fileno())

        protocol['environment'] = dict(python=platform.python_version(), numpy=np.__version__,
                                       scipy=scipy.__version__, host=platform.node(), machine=platform.machine())
        save('protocol', protocol)
        sys.path.insert(0, str(archive))
        require(not any(n[:-3] in sys.modules for n in MODULES[:4]), 'Original modules already imported')
        loader = importlib.import_module('bfix_v1')
        importlib.import_module('tests_v1')
        for name in MODULES[:4]:
            require(Path(sys.modules[name[:-3]].__file__).resolve().parent == archive.resolve(),
                    'Original import escaped archive')
        private = load('rule_diagnosis_original', archive/'likelihood_v1.py')
        anchor = load('rule_diagnosis_anchor', archive/'repair_p3_bb_mode_anchor.py')
        reuse = load('rule_diagnosis_cache', archive/'repair_p3_bb_mode_anchor_cache_v2.py')
        cache = reuse.OnePointBBTerms(anchor.stable_bb_terms)
        design, settings = loader.load()
        full = private.Model(design.model().rows, gene_count=8)
        rows = tuple(row for row in full.rows if row.gene == 2)
        gene_only = private.Model(rows, gene_count=8)
        point = json.loads((POINT/'protocol.json').read_text())
        p, v = np.asarray(point['fixed_p']), point['fixed_v']
        count_sha = hashlib.sha256(np.asarray([row.a for row in full.rows], dtype='<i8').tobytes()).hexdigest()
        require(len(rows) == 24 and len(full.rows) == 173 and full.size == gene_only.size == 28
                and settings['seed'] == 20260930 and count_sha == point['count_vector_sha256']
                and protocol['environment']['numpy'] == point['environment']['numpy']
                and protocol['environment']['scipy'] == point['environment']['scipy'], 'Fixture/environment differs')

        def identical(a, b):
            aa, bb = np.asarray(a), np.asarray(b)
            return aa.dtype == bb.dtype and aa.shape == bb.shape and aa.tobytes() == bb.tobytes()

        nodes, weights = np.polynomial.legendre.leggauss(16)
        edges = np.linspace(-8., 8., 65)
        totals = np.zeros((2, 30))
        corrections = np.zeros_like(totals)
        shifts = np.full(2, -np.inf)
        differences = []
        predictor_mismatches = conditional_mismatches = aggregate_node_mismatches = 0
        ordinal = 0
        for left, right in zip(edges[:-1], edges[1:]):
            for node, weight in zip(nodes, weights):
                z = (left+right)/2+(right-left)*node/2
                u = float(np.sqrt(v)*z)
                ll = np.zeros(2)
                gp = np.zeros((2, 28))
                for row in rows:
                    a = gene_only.predictors(row, p, u, need_hessian=False)
                    b = full.predictors(row, p, u, need_hessian=False)
                    if not all(identical(x, y) for x, y in zip(a[:3], b[:3])):
                        predictor_mismatches += 1
                    eta, dispersion, jac = a[:3]
                    private.bb_terms = anchor.stable_bb_terms
                    scalar = private.conditional_terms(row, eta, dispersion, False)
                    private.bb_terms = cache
                    cached = private.conditional_terms(row, eta, dispersion, False)
                    equal = identical(scalar[0], cached[0]) and identical(scalar[1], cached[1])
                    if not equal:
                        conditional_mismatches += 1
                        if len(differences) < 10:
                            differences.append(dict(node=ordinal, z=z, row=vars(row), eta=eta,
                                dispersion=dispersion, scalar=scalar[:2], cached=cached[:2]))
                    for mode, result in enumerate((scalar, cached)):
                        ll[mode] += result[0]
                        gp[mode] += jac.T@result[1]
                if not identical(ll[0], ll[1]) or not identical(gp[0], gp[1]):
                    aggregate_node_mismatches += 1
                for mode in range(2):
                    exponent = ll[mode]-z*z/2-.5*np.log(2*np.pi)+np.log((right-left)*weight/2)
                    vector = np.r_[1., gp[mode], (z*z-1)/(2*v)]
                    if exponent > shifts[mode]:
                        factor = np.exp(shifts[mode]-exponent)
                        totals[mode] *= factor
                        corrections[mode] *= factor
                        shifts[mode] = exponent
                    increment = np.exp(exponent-shifts[mode])*vector-corrections[mode]
                    updated = totals[mode]+increment
                    corrections[mode] = (updated-totals[mode])-increment
                    totals[mode] = updated
                ordinal += 1
            print(json.dumps(dict(completed_nodes=ordinal, conditional_mismatches=conditional_mismatches)), flush=True)
        results = [dict(logmass=float(shifts[m]+np.log(totals[m, 0])),
                        gradient=totals[m, 1:]/totals[m, 0]) for m in range(2)]
        for path, digest in pins.items():
            require(sha(path) == digest, 'Reviewed source/input changed during diagnostic')
        save('summary', dict(completed=True, nodes=ordinal, conditional_pairs=ordinal*len(rows),
            predictor_mismatches=predictor_mismatches, conditional_mismatches=conditional_mismatches,
            aggregate_node_mismatches=aggregate_node_mismatches, first_discrepancies=differences,
            results=results, rule_logmass_byte_equal=identical(results[0]['logmass'], results[1]['logmass']),
            rule_gradient_byte_equal=identical(results[0]['gradient'], results[1]['gradient']),
            cache_statistics=cache.statistics(), measured_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            elapsed_seconds=time.perf_counter()-began, no_accuracy_acceptance=True, no_fit=True, no_calibration=True))
    except BaseException:
        with (output/'failure.json').open('x') as handle:
            json.dump(dict(traceback=traceback.format_exc(), partials_preserved=True), handle, indent=2)
        raise


if __name__ == '__main__':
    main()
