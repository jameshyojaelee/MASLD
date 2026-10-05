#!/usr/bin/env python3
"""Fixed synthetic composite-GL diagnostics; all scalar/FD references retained."""
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
KERNEL = HERE / 'repair_p3_composite_gl_integral.py'
KERNEL_SHA = '58de736b9b2d458367d5726e2582d45d7aa5e68136a223335239a5cbbd7663e9'
INSTRUMENTED = HERE / 'repair_p3_integral_instrumented.py'
INSTRUMENTED_SHA = '03a5e8b504ac426b16557c75f7d416d90d17c299cf474b6598a980dd072fa347'
SCALARS = {'positive18':REC/'repair_p3_scalar_prefix_positive18_22039119',
           'failed0gene2':REC/'repair_p3_scalar_prefix_failed0gene2_22039118'}
BASE_KERNEL = HERE / 'repair_p3_integral.py'
BASE_KERNEL_SHA = '461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'
ANCESTOR_SHA = '89434fd7e32b06cb93d9d31f18b646eaa646736ec329ab48cc00fa4409561f72'
ANCESTOR = REC / 'repair_p3_launch_20261005T121740Z/repair_p3_fixed_point_check.py'
CAP = 32 * 1024**2
CASES = {'failed0gene2':0,'positive18':18}
G = 2
STEPS = (1e-3, 3e-4, 1e-4)
FIXED_AXES = (2,10,16,17,24,25,28)
GUARDS = {
    PROTOCOL: '8e7ddf8b0a7bb18afd9692b53c3560f4258fdc373c4b822aee36aa94395eabe9',
    CAMPAIGN / 'campaign.json': 'bb3e82833476f67b0b11d2388e2dc00358eb987367b7a6a666d11938c6baca94',
    KERNEL: KERNEL_SHA,
    INSTRUMENTED: INSTRUMENTED_SHA,
    SCALARS['positive18'] / 'independent_baseline.json': '0184daa0afae8e392272868530d31ad064089c3de06483d6cc8b87b76b109257',
    SCALARS['positive18'] / 'summary.json': 'fba2f38821c59e85ea4973e5ebe8cfe2509bc35d9910fcb7d9638078d6adb43d',
    SCALARS['positive18'] / 'protocol.json': '26f061f8e481a1b9d1e11f2a5dfc0984c82010d3b1dd9c8f270eb1418d701001',
    SCALARS['positive18'] / 'finite_differences.json': 'ebce73a49bbca829ff182e4cb0fc710c18fea794edbe8f8e476a7d145c91ac62',
    SCALARS['failed0gene2'] / 'independent_baseline.json': 'b6f8422966c6db5b10ccbe498a30e1e899f92ec16c66736f71a6c08a7e685d1a',
    SCALARS['failed0gene2'] / 'summary.json': '74f375cf748940e4bf91137e55ab6e2e8d564b7fedb48032c19d60c19e688d53',
    SCALARS['failed0gene2'] / 'protocol.json': 'b27a85f6c11342406e07e99aacafa70839d09fcd870e0449dac40e59aab97a32',
    SCALARS['failed0gene2'] / 'finite_differences.json': '2ae8d822a8e456ad6e247ffaf91bcc21f5a90e06b63c66c8d2846a0cda96451a',
    BASE_KERNEL: BASE_KERNEL_SHA,
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
    ast.parse(KERNEL.read_text())
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
    parser.add_argument('--expected-kernel-sha')
    args = parser.parse_args()
    SCALAR = SCALARS[args.case]
    protocol, campaign, guards = metadata()
    if args.inspect:
        print(json.dumps(dict(case=args.case, draw=CASES[args.case], fixed_gene=G,
              numerical_execution=False, model_fitting=False, calibration=False,
              diagnostic_only=True,z_bound=15.,orders=(16,32,64),refined_z_bound=17.,
              refined_grid_spacing=.125,accepted_candidate_integral=False,
              scalar_integrals_repeated=False,finite_differences_repeated=False,
              kernel_sha256=sha(KERNEL), base_kernel_sha256=BASE_KERNEL_SHA,
              ancestor_sha256=ANCESTOR_SHA, metadata_and_source_checks=True)))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'), 'Individual compute job required')
    require(args.expected_helper_sha == sha(Path(__file__)) and args.expected_kernel_sha == sha(KERNEL), 'Reviewed instrumented source drift')
    output = REC / ('repair_p3_composite_gl_' + args.case + '_' + os.environ['SLURM_JOB_ID'])
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
        own = {Path(__file__):args.expected_helper_sha, KERNEL:args.expected_kernel_sha}
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
        loader = importlib.import_module('bfix_v1')
        original = importlib.import_module('likelihood_v1')
        importlib.import_module('tests_v1')  # explicit before the archived-module assertion
        for name in ('likelihood_v1','bfix_v1','design_v1','tests_v1'):
            require(Path(sys.modules[name].__file__).resolve().parent == archive.resolve() and
                    sha(Path(sys.modules[name].__file__)) == campaign['files'][str(CAMPAIGN/'executed_sources'/ (name+'.py'))],
                    'Original import escaped guarded archive')
        kernel = load_module('repair_p3_archived_instrumented_kernel', archive / KERNEL.name)
        require(Path(kernel.__file__).resolve().parent == archive.resolve(), 'Instrumented import escaped archive')
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
        require(v>0, 'Candidate-only diagnostic requires saved positive variance')
        # Only this archived private module is modified; no candidate marginal/fit callable.
        original.bb_terms = kernel.stable_bb_terms
        candidate = original.Model(rows,gene_count=model.G)
        scalar_protocol=json.loads((SCALAR/'protocol.json').read_text())
        scalar_summary=json.loads((SCALAR/'summary.json').read_text())
        scalar_baseline=json.loads((SCALAR/'independent_baseline.json').read_text())
        scalar_fd=json.loads((SCALAR/'finite_differences.json').read_text())
        require(len(scalar_fd)==21 and {(r['axis'],r['step']) for r in scalar_fd}==
                {(a,h) for a in FIXED_AXES for h in STEPS}, 'Full fixed scalar FD family differs')
        require(scalar_summary['scalar_fixed_point_check_completed'] and scalar_summary['scalar_fd_record_count']==21 and
                scalar_summary['case']==args.case and scalar_summary['draw']==draw, 'Complete scalar reference case differs')
        require(scalar_protocol['fixed_p']==p.tolist() and scalar_protocol['fixed_v']==v and
                scalar_protocol['count_vector_sha256']==count_sha, 'Scalar reference point or row counts differ')
        save('protocol',dict(case=args.case,draw=draw,fixed_gene=G,fixed_gene_rows=24,fixed_p=p,fixed_v=v,
            count_vector_sha256=count_sha,environment=versions,
            source_guards={str(k):d for k,d in guards.items()},own_source_hashes={str(k):d for k,d in own.items()},
            instrumented_kernel_settings=kernel.Settings().__dict__,base_kernel_sha256=BASE_KERNEL_SHA,
            numerical_targets_and_caps_unchanged=True,z_bound=15.,orders=(16,32,64),refined_z_bound=17.,
            refined_grid_spacing=.125,accepted_candidate_integral=False,uncertified_estimates=True,
            unlocated_narrow_peak_may_be_missed=True,
            no_model_fit=True,no_calibration=True,no_scalar_integral_repeat=True,no_finite_difference_repeat=True,
            stochastic_draw_source='Original supplied uniform row only; no new Monte Carlo seed.',
            output_cap_bytes=CAP,error_estimates_not_certified=True))
        phase='composite_gl_diagnostic'
        report_count=0
        def observer(record):
            nonlocal report_count
            report_count += 1
            save('run_event_%04d' % report_count,dict(event_number=report_count,
                elapsed_seconds=time.monotonic()-began,**record))
            event('instrumentation_saved',event_number=report_count,
                  integration_event=record['event'],run_number=record.get('run_number'),
                  unique_evaluations=record['unique_evaluations'])
        returned=False
        comparison_count=0
        completed_order_count=0
        def compare_reference(result):
            nonlocal comparison_count,completed_order_count
            ll,g,h,report=result
            rows=[]
            for reference in scalar_fd:
                axis=reference['axis']
                difference=float(g[axis]-reference['scalar_reference_gradient'])
                tolerance=1e-8+1e-6*abs(float(g[axis]))
                rows.append(dict(axis=axis,step=reference['step'],method=reference['method'],
                    scalar_reference_gradient=reference['scalar_reference_gradient'],moment_gradient=float(g[axis]),
                    moment_minus_scalar_fd=difference,binding_atol=1e-8,binding_rtol=1e-6,
                    binding_tolerance_at_moment_gradient=tolerance,
                    raw_difference_within_binding_tolerance_descriptive_only=abs(difference)<=tolerance,
                    scalar_estimated_integration_plus_tail_FD_error=reference['estimated_integration_plus_tail_FD_error'],
                    FD_truncation_and_arithmetic_error_unbounded=True,no_derivative_accuracy_acceptance=True))
            comparison_count+=len(rows)
            completed_order_count+=1
            save('scalar_comparison_%02d' % completed_order_count,dict(order=report['order'],z_bound=report['z_bound'],
                grid_spacing=report['grid_spacing'],scalar_reference_source_sha256=GUARDS[SCALAR/'independent_baseline.json'],
                scalar_reference_logmass=scalar_baseline['tighter_expanded']['logmass'],
                moment_minus_scalar_logmass=ll-scalar_baseline['tighter_expanded']['logmass'],
                scalar_reference_error_estimate=scalar_baseline['tighter_expanded']['estimated_quad_plus_mass_tail_log_error'],
                finite_difference_comparisons=rows,comparison_is_descriptive_not_accuracy_acceptance=True))
        # Emit reference comparisons as each order finishes, retaining partials
        # even if a later order/refinement reaches the unchanged evaluation cap.
        base_observer=observer
        def observer(record):
            base_observer(record)
            if record['event']=='order_result':
                compare_reference((record['logmass'],record['gradient'],record['hessian'],record['report']))
        try:
            results,report=kernel.integrate_gene(candidate,p,v,gene_rows,original.conditional_terms,
                                                with_hessian=True,observer=observer)
            save('composite_result',dict(diagnostic_only=True,accepted_candidate_integral=False,
                results=[dict(logmass=r[0],gradient=r[1],hessian=r[2],integration=r[3]) for r in results],
                refinement_report=report,complete_scalar_baseline=scalar_baseline))
            require(completed_order_count==4 and comparison_count==84,'Incomplete order/reference comparison family')
            returned=True
            event('composite_gl_diagnostic_returned')
        except (ValueError,FloatingPointError,np.linalg.LinAlgError) as error:
            save('composite_refusal',dict(completed=False,exception_type=type(error).__name__,
                traceback=traceback.format_exc(),instrumentation_events_saved=report_count,
                completed_orders=completed_order_count,scalar_fd_comparisons_saved=comparison_count,
                not_a_replacement_failed_draw=True))
            event('composite_gl_diagnostic_refused',reason=str(error))
        phase='final_source_hashes'
        for path,expected in {**guards,**own}.items():
            require(sha(path)==expected,'Source/fixture changed during instrumented check')
        save('summary',dict(composite_gl_diagnostic_completed=True,case=args.case,draw=draw,fixed_gene=G,
            diagnostic_integrals_returned=returned,instrumentation_events_saved=report_count,
            numerical_targets_and_caps_unchanged=True,z_bound=15.,orders=(16,32,64),refined_z_bound=17.,
            completed_orders=completed_order_count,scalar_fd_comparisons_saved=comparison_count,
            accepted_candidate_integral=False,full_domain_grid_refinement_completed=returned,
            uncertified_estimates=True,unlocated_narrow_peak_may_be_missed=True,
            no_model_fits=True,no_statistical_calibration=True,no_scalar_integrals_repeated=True,
            no_finite_differences_repeated=True,derivative_accuracy_pass_not_asserted=True,
            full_model_stationarity_not_assessed=True,original_records_preserved=True,
            error_estimates_not_certified=True,elapsed_seconds=time.monotonic()-began,output_cap_bytes=CAP,
            saved_bytes_before_summary=sum(q.stat().st_size for q in output.rglob('*') if q.is_file())))
        event('completed')
        require(sum(q.stat().st_size for q in output.rglob('*') if q.is_file())<=CAP,'Final inclusive output cap exceeded')
    except BaseException:
        with (output/'failure.json').open('x') as handle:
            json.dump(dict(completed=False,phase=phase,traceback=traceback.format_exc(),partials_preserved=True,
                           no_fit_or_calibration=True),handle,indent=2)
            handle.write('\n')
        raise


if __name__=='__main__':
    main()



