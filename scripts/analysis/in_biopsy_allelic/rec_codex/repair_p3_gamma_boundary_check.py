#!/usr/bin/env python3
"""Full synthetic v=0 parent identity and arithmetic-boundary consistency."""
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
HELPERS = ('repair_p3_bb_gamma_difference_v3.py','repair_p3_composite_gl.py',
           'repair_p3_integral.py','repair_p3_gamma_parallel_model.py')
GL_POINT = REC/'repair_p3_gamma_observed_point_GL_22076126'
NULL_POINT = REC/'repair_p3_parallel_observed_fit_22042114'
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
    require(all(p in pins for p in (GL_POINT/'summary.json',NULL_POINT/'protocol.json',
                                   NULL_POINT/'observed_null.json')),
            'Completed GL admission and saved native null point must be pinned')
    gl_admitted = json.loads((GL_POINT/'summary.json').read_text())
    null_identity = json.loads((NULL_POINT/'protocol.json').read_text())
    require(gl_admitted['completed'] is True and gl_admitted['whole_saved_score_integral_nominal_criteria_met'] is True
            and null_identity['fixture_seed'] == 20260930 and null_identity['participants'] == 60
            and null_identity['rows'] == 173 and null_identity['observed_count_vector_sha256'] == COUNT_SHA,
            'Whole-point admission or original null fixture differs')
    protocol = dict(source_and_input_sha256={str(p): h for p, h in pins.items()},
        synthetic_participants=60,rows=173,G=8,Z=7,P=28,fixture_seed=20260930,
        count_vector_sha256=COUNT_SHA,
        points='saved converged native null, original cold_start(), saved call12; each kappa shift 0,-0.1,+0.1 clipped to [-0.9,0.9]',
        point_count=9,zero_variance_parent_byte_identity_required=True,
        boundary_LL_absolute_tolerance='1e-8',boundary_score_atol='1e-8',boundary_score_rtol='1e-6',
        boundary_nuisance_Hessian_absolute_tolerance='1e-5',right_score_restart_cutoff='1e-6',
        original_cold_start_only=True,no_observed_model_fit=True,no_positive_worker_evaluation=True,
        no_calibration=True,no_adoption=True,output_cap_bytes=CAP)
    if args.inspect:
        print(json.dumps(protocol))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and
            not os.environ.get('SLURM_ARRAY_TASK_ID'), 'Individual compute allocation required')
    output = REC/('repair_p3_gamma_boundary_check_'+os.environ['SLURM_JOB_ID'])
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
        kernel = load('boundary_reviewed_kernel',archive/'repair_p3_integral.py')
        model_module = load('boundary_prepared_model',archive/'repair_p3_gamma_parallel_model.py')
        Parent = model_module.gamma_parallel_model_class(private,kernel,gl,candidate,
                                                         kernel_path=archive/'repair_p3_integral.py')
        subject = Parent(model.rows,gene_count=8)
        require(all(getattr(Parent,name) is getattr(private.Model,name) for name in model_module.INHERITED),
                'Inherited fit/design function identity changed')
        original_bb = private.bb_terms
        counterpart = load('boundary_isolated_gamma_original',archive/'likelihood_v1.py')
        counterpart.bb_terms = candidate.stable_bb_terms
        gamma_model = counterpart.Model(model.rows,gene_count=8)
        # Test detached snapshots without altering any biological/synthetic row.
        from types import SimpleNamespace
        toy = SimpleNamespace(value=np.asarray([1.,2.]))
        signature = model_module.row_signatures((toy,),np)
        toy.value[0] = 3.
        require(model_module.row_signatures((toy,),np) != signature,'In-place array change undetected')
        toy.value = np.asarray([1.,2.])
        require(model_module.row_signatures((toy,),np) == signature,'Equal array values spuriously differ')
        save('snapshot_invariants',dict(in_place_mutation_detected=True,equal_rebinding_supported=True))
        null = json.loads((NULL_POINT/'observed_null.json').read_text())['fit']
        require(null['converged'] is True,'Original native null is not converged')
        starts = [('native_null',np.asarray(null['parameters'],dtype=float)),
                  ('native_cold',model.cold_start()),
                  ('saved_call12',np.asarray(parent['fixed_p'],dtype=float))]
        require(int(os.environ.get('SLURM_CPUS_PER_TASK','0')) == 1,
                'Boundary preflight check requires the declared one-CPU job')
        try:
            subject.marginal(starts[0][1],.001)
        except model_module.ModelIntegrityError:
            pass
        else:
            raise ValueError('One-CPU positive call did not stop as an integrity error')
        require(subject.last_integration_reports == [] and subject.last_evaluated_point is None,
                'Preflight refusal retained stale integration state')
        save('preflight_integrity',dict(one_CPU_positive_refusal_is_RuntimeError=True,
            not_ValueError=True,no_worker_started=True,no_stale_integration_state=True))
        results = []
        def byte_equal(a,b):
            return all(np.asarray(x).dtype == np.asarray(y).dtype and np.asarray(x).shape == np.asarray(y).shape
                       and np.asarray(x).tobytes() == np.asarray(y).tobytes() for x,y in zip(a,b))
        for label,start in starts:
            for shift in (0.,-.1,.1):
                p = start.copy()
                p[16] = np.clip(p[16]+shift,-.9,.9)
                native = model.marginal(p,0.,with_hessian=True)
                parent_answer = subject.marginal(p,0.,with_hessian=True)
                gamma = gamma_model.marginal(p,0.,with_hessian=True)
                native_score_only = model.marginal(p,0.,with_hessian=False)
                parent_score_only = subject.marginal(p,0.,with_hessian=False)
                same = byte_equal(native,parent_answer) and byte_equal(native_score_only,parent_score_only)
                delta_ll = gamma[0]-native[0]
                delta_g = gamma[1]-native[1]
                delta_h = gamma[2][:-1,:-1]-native[2][:-1,:-1]
                require(np.isfinite(delta_ll) and np.isfinite(delta_g).all() and np.isfinite(delta_h).all(),
                        'Nonfinite boundary comparison outside native unavailable v-Hessian')
                restart_flags = [bool(x[1][-1] > 1e-6) for x in (native,gamma)]
                passed = (same and abs(delta_ll) <= 1e-8 and
                    np.allclose(gamma[1],native[1],atol=1e-8,rtol=1e-6) and
                    np.max(np.abs(delta_h)) <= 1e-5 and restart_flags[0] == restart_flags[1])
                record = dict(label=label,kappa_shift=shift,parameters=p,
                    original_parent_byte_identity=bool(same),native_loglik=native[0],gamma_loglik=gamma[0],
                    gamma_minus_native_LL=delta_ll,native_score=native[1],gamma_score=gamma[1],
                    gamma_minus_native_score=delta_g,native_nuisance_Hessian=native[2][:-1,:-1],
                    gamma_nuisance_Hessian=gamma[2][:-1,:-1],gamma_minus_native_nuisance_Hessian=delta_h,
                    right_variance_scores=[float(native[1][-1]),float(gamma[1][-1])],
                    above_original_restart_cutoff=restart_flags,passed=bool(passed),
                    variance_Hessian_unavailable_by_original_definition=True)
                results.append(record)
                save(label+'_shift_'+str(shift).replace('.','p'),record)
        require(private.bb_terms is original_bb,'Isolated boundary comparison modified parent BB')
        require(hashlib.sha256(np.asarray([r.a for r in model.rows],dtype='<i8').tobytes()).hexdigest() == COUNT_SHA,
                'Original count vector changed during boundary comparison')
        for path,digest in pins.items():
            require(sha(path) == digest,'Source/input changed during boundary comparison')
        passed = all(r['passed'] for r in results)
        save('summary',dict(completed=True,records=len(results),all_boundary_nominal_criteria_met=passed,
            all_parent_zero_variance_calls_byte_identical=all(r['original_parent_byte_identity'] for r in results),
            original_parent_BB_unchanged=True,original_cold_start_only=True,no_observed_model_fit=True,
            no_positive_worker_evaluation=True,no_full_independent_Hessian=True,
            elapsed_seconds=time.perf_counter()-began,
            measured_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            no_calibration=True,no_adoption=True))
        require(passed,'Boundary arithmetic consistency criteria failed')
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
