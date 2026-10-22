#!/usr/bin/env python3
"""Whole original GL score integral with admitted gamma-difference arithmetic."""
import argparse
import ast
from decimal import Decimal, localcontext
import hashlib
import importlib
import importlib.util
import json
import math
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
POINT = REC/'repair_p3_observed_failure_anchor_22042612'
THETA = REC/'repair_p3_theta4_anchored_gl_gene2_22042347'
CASES = REC/'repair_p3_bb_mode_anchor_batch_fast_launch_20261005T172617Z/cases.json'
MODULES = ('likelihood_v1.py', 'design_v1.py', 'bfix_v1.py', 'tests_v1.py')
HELPERS = ('repair_p3_bb_gamma_difference_v3.py','repair_p3_composite_gl.py','repair_p3_integral.py')
ADMISSION = REC/'repair_p3_bb_gamma_difference_v3_check_22076101'
KERNEL_SHA = '461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'
COUNT_SHA = '7009be4adeec6157fa4bbf3f46e48ec24aa398a0bee4d7906c78c231b7d3fbad'
CAP = 32*1024**2


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
    parser.add_argument('--expected-checker-sha', required=True)
    parser.add_argument('--pins', type=Path, required=True)
    parser.add_argument('--expected-pins-sha', required=True)
    parser.add_argument('--inspect', action='store_true')
    args = parser.parse_args()
    require(sha(Path(__file__)) == args.expected_checker_sha and
            sha(args.pins) == args.expected_pins_sha, 'Reviewed checker/manifest differs')
    pins = {Path(p): h for p, h in json.loads(args.pins.read_text()).items()}
    for path, digest in pins.items():
        require(path.is_absolute() and path.resolve().is_relative_to(ROOT.resolve())
                and sha(path) == digest, 'Reviewed source/input differs: '+str(path))
        if path.suffix == '.py':
            ast.parse(path.read_text())
    required = [*(HERE/n for n in (*MODULES, *HELPERS)), CASES,
                *(POINT/n for n in ('protocol.json', 'summary.json', 'GL_candidate.json',
                                   'Decimal_reference_summary.json', 'GL_minus_Decimal.json')),
                THETA/'protocol.json']
    require(all(p in pins for p in required), 'Required source/input pin missing')
    parent = json.loads((POINT/'protocol.json').read_text())
    summary = json.loads((POINT/'summary.json').read_text())
    theta = json.loads((THETA/'protocol.json').read_text())
    packet = json.loads(CASES.read_text())
    require(summary['completed'] is True and summary['GL_candidate_refused'] is False
            and parent['gene'] == 2 and parent['call'] == 12 and parent['target_rows'] == 24
            and parent['fixed_v'] == .0009636462009520182
            and parent['count_vector_sha256'] == COUNT_SHA
            and theta['fixed_v'] == .0004 and len(packet['cases']) == 21,
            'Exact completed synthetic reference family differs')
    require(ADMISSION/'summary.json' in pins and ADMISSION/'protocol.json' in pins,
            'Completed component admission must be pinned')
    admitted = json.loads((ADMISSION/'summary.json').read_text())
    admitted_protocol = json.loads((ADMISSION/'protocol.json').read_text())
    require(admitted_protocol['source_and_input_sha256'][str(HERE/'repair_p3_bb_gamma_difference_v3.py')]
            == pins[HERE/'repair_p3_bb_gamma_difference_v3.py'],
            'Admitted candidate source differs from whole-point candidate')
    require(all(admitted[k] is True for k in ('completed','selected21_criterion_met',
            'independent_rising_criterion_met','conditional_nominal_LP_score_H_criteria_met',
            'native_helper_call_counts_checked')) and admitted['conditional_comparisons'] == 2076,
            'Component admission incomplete')
    protocol = dict(source_and_input_sha256={str(p): h for p, h in pins.items()},
        gene=2,synthetic_participants=60,rows=173,target_rows=24,fixture_seed=20260930,
        count_vector_sha256=COUNT_SHA,fixed_p=parent['fixed_p'],fixed_v=parent['fixed_v'],
        full_GL_with_hessian=False,original_GL_rules_and_caps=True,
        comparison_logmass_absolute_tolerance='1e-9',comparison_gradient_absolute_tolerance='1e-8',
        final_independent_reference=dict(decimal_digits=110,order=80),
        all_six_independent_reference_records_retained=True,old_failures_preserved=True,
        output_cap_bytes=CAP,no_fit=True,no_calibration=True,no_adoption=True)
    if args.inspect:
        print(json.dumps(protocol))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and
            not os.environ.get('SLURM_ARRAY_TASK_ID'), 'Individual compute allocation required')
    output = REC/('repair_p3_gamma_observed_point_GL_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(exist_ok=False)
    archive = output/'executed_sources'
    archive.mkdir()
    for source in (Path(__file__), *(HERE/n for n in (*MODULES, *HELPERS))):
        target = archive/source.name
        target.write_bytes(source.read_bytes())
        require(sha(target) == (args.expected_checker_sha if source == Path(__file__)
                               else pins[source]), 'Executable archive differs')
    began = time.perf_counter()
    try:
        import numpy as np
        import scipy

        def serialize(value):
            if isinstance(value, np.ndarray):
                return value.tolist()
            if isinstance(value, np.generic):
                return value.item()
            raise TypeError(type(value).__name__)

        def save(name, value):
            body = (json.dumps(value, default=serialize, allow_nan=False, indent=2)+'\n').encode()
            require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file())+len(body)
                    < CAP-65536, 'Output capacity exceeded; refusal reserve retained')
            with (output/(name+'.json')).open('xb') as handle:
                handle.write(body)
                handle.flush()
                os.fsync(handle.fileno())
            print(json.dumps(dict(saved=name, elapsed_seconds=time.perf_counter()-began)), flush=True)

        versions = dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)
        require(versions == parent['environment'], 'Numerical environment differs')
        protocol['environment'] = versions
        save('protocol', protocol)
        sys.path.insert(0, str(archive))
        require(not any(n[:-3] in sys.modules for n in MODULES), 'Original modules already imported')
        loader = importlib.import_module('bfix_v1')
        importlib.import_module('tests_v1')
        for name in MODULES:
            module = sys.modules[name[:-3]]
            require(Path(module.__file__).resolve().parent == archive.resolve()
                    and sha(Path(module.__file__)) == pins[HERE/name], 'Original import escaped archive')
        private = load('cache_check_original', archive/'likelihood_v1.py')
        candidate = load('gamma_GL_candidate',archive/'repair_p3_bb_gamma_difference_v3.py')
        gl = load('cache_check_GL', archive/'repair_p3_composite_gl.py')
        design, settings = loader.load()
        model = private.Model(design.model().rows, gene_count=8)
        require(len(model.rows) == 173 and model.size == 28 and model.Z == 7
                and len(design.participants) == 60 and settings['seed'] == 20260930
                and hashlib.sha256(np.asarray([r.a for r in model.rows], dtype='<i8').tobytes()).hexdigest()
                == COUNT_SHA, 'Exact synthetic fixture differs')
        p,v = np.asarray(parent['fixed_p']),parent['fixed_v']
        rows = tuple(row for row in model.rows if row.gene == 2)
        require(len(rows) == 24 and max(row.n for row in model.rows) <= 5000,'Admitted support scope differs')
        calls = 0
        def monitored(*args,**kwargs):
            nonlocal calls
            calls += 1
            return candidate.stable_bb_terms(*args,**kwargs)
        private.bb_terms = monitored
        started = time.perf_counter()
        ll,g,h,report = gl.integrate_gene(model,p,v,rows,private.conditional_terms,
            with_hessian=False,expected_kernel_sha=KERNEL_SHA,kernel_path=archive/'repair_p3_integral.py')
        require(calls > 0 and h is None,'Original GL did not exercise the admitted score helper')
        save('GL_candidate',dict(logmass=ll,gradient=g,Hessian=h,integration=report,
            elapsed_seconds=time.perf_counter()-started,native_BB_calls=calls))
        baseline = json.loads((POINT/'GL_candidate.json').read_text())
        require(len(report['ladder']) == len(baseline['integration']['ladder']) == 9
                and report['likelihood_evaluations'] == baseline['integration']['likelihood_evaluations'] == 51072,
                'Whole original saved refinement family differs')
        comparisons = []
        for current,old in zip(report['ladder'],baseline['integration']['ladder']):
            require(all(current[k] == old[k] for k in ('z_bound','grid_spacing','order','node_evaluations')),
                    'Saved-rule coverage differs')
            delta_g = np.asarray(current['gradient'])-np.asarray(old['gradient'])
            delta_ll = current['logmass']-old['logmass']
            comparisons.append(dict(z_bound=current['z_bound'],grid_spacing=current['grid_spacing'],
                order=current['order'],logmass_difference=delta_ll,gradient_difference=delta_g,
                within_nominal_criteria=abs(delta_ll) <= 1e-9 and bool(np.max(np.abs(delta_g)) <= 1e-8)))
        save('all_saved_rule_comparisons',comparisons)
        reference = json.loads((POINT/'Decimal_reference_summary.json').read_text())
        require(reference['completed'] is True and
                [(r['decimal_digits'],r['order']) for r in reference['records']] ==
                [(d,o) for d in (80,110) for o in (20,40,80)],'Original Decimal reference family differs')
        independent = []
        with localcontext() as context:
            context.prec = 110
            for r in reference['records']:
                dl = Decimal.from_float(float(ll))-Decimal(r['logmass'])
                dg = Decimal.from_float(float(g[-1]))-Decimal(r['variance_score'])
                independent.append(dict(decimal_digits=r['decimal_digits'],order=r['order'],
                    logmass_difference=str(dl),variance_score_difference=str(dg),
                    within_nominal_criteria=abs(dl) <= Decimal('1e-9') and abs(dg) <= Decimal('1e-8')))
        save('all_Decimal_comparisons',independent)
        final = independent[-1]
        for path,digest in pins.items():
            require(sha(path) == digest,'Source/input changed during whole-point integral')
        passed = all(r['within_nominal_criteria'] for r in comparisons) and final['within_nominal_criteria']
        save('summary',dict(completed=True,whole_saved_score_integral_nominal_criteria_met=passed,
            original_rules_and_evaluation_count_unchanged=True,native_BB_calls=calls,
            final_prespecified_Decimal_110_order80_comparison=final,
            no_full_independent_Hessian=True,no_model_fit=True,no_calibration=True,no_adoption=True,
            elapsed_seconds=time.perf_counter()-began,
            measured_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        require(passed,'Whole saved-point nominal score-integral criteria failed')
    except BaseException as exc:
        failure = dict(traceback=traceback.format_exc(),partials_preserved=True,no_fit=True,
            no_calibration=True,no_adoption=True,
            refused_refinement_reports=getattr(exc,'integration_reports',None),
            refused_likelihood_evaluations=getattr(exc,'likelihood_evaluations',None))
        if 'save' in locals():
            save('failure',failure)
        else:
            with (output/'failure.json').open('x') as handle:
                json.dump(failure,handle,indent=2)
        raise


if __name__ == '__main__':
    main()
