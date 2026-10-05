#!/usr/bin/env python3
"""Guarded complete observed fit of the original60-person synthetic fixture only.

No draws, calibration, real/private inputs, model claim, or release change.
All numerical imports occur after source admission and the compute guard.
"""
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

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
CAMPAIGN = REC/'p3_parallel_20261003T130323Z'
FIXTURE = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
HERE = Path(__file__).resolve().parent
MODEL = HERE/'repair_p3_batched_model.py'
KERNEL = REC/'repair_p3_launch_20261005T121740Z/repair_p3_integral.py'
GL = HERE/'repair_p3_batched_gl.py'
CONDITIONAL = HERE/'repair_p3_batched_conditional.py'
V2 = HERE/'repair_p3_prepared_conditional_v2.py'
DYNAMIC = HERE/'repair_p3_dynamic_conditional.py'
ORIGINAL_RUNNER = REC/'repair_p3_batched_observed_fit_launch_20261005T152152Z/repair_p3_batched_observed_fit.py'
MODEL_SHA = '2f64fdaffb44e5d78fd327ca6caf8672242ff53e17e53ad89b7e249bd77c0c3c'
FULL_CHECKS = {
    REC/'repair_p3_dynamic_full_failed0gene2_22041698/protocol.json':'67b2592aaff947963ebee81397f32d7827f2c8bf6152efbc5029f26225607aa1',
    REC/'repair_p3_dynamic_full_failed0gene2_22041698/summary.json':'fb4fe6af4f1570982b23efbf26d4440a20c6cc98d4147e50ff7e1130a4ada668',
    REC/'repair_p3_dynamic_full_failed0gene2_22041698/dynamic_full_zero_law_comparison.json':'c69514434cedc793756186a069cfc4c73cfbe2481c97ba2e23e90007fb20e3d1',
    REC/'repair_p3_dynamic_full_positive18_22041699/protocol.json':'eccb2f9094488b8b9617a89d644aeeb57fab054dc0d75358cc545eb25db44d0a',
    REC/'repair_p3_dynamic_full_positive18_22041699/summary.json':'7efb0c5ec298c92c9e7dd3e65361eccf7213b7e63e26060243be5cdbd58cce89',
    REC/'repair_p3_dynamic_full_positive18_22041699/dynamic_full_zero_law_comparison.json':'e0801d8e8141da0dac263905a038ccdc6480a03f809494e80c3f61f06442187e',
}
WARM = REC/'bfix_fits_21984981/null.json'
CAP = 32*1024**2
RESERVE = 8*1024**2
PROGRESS_SOFT_CAP = 8*1024**2
GUARDS = {
    CAMPAIGN/'campaign.json':'bb3e82833476f67b0b11d2388e2dc00358eb987367b7a6a666d11938c6baca94',
    CAMPAIGN/'baseline/observed_null.json':'a617d331fd3fb93576e5bd8588581aaab8f77946936536401905490438660181',
    CAMPAIGN/'baseline/observed_marginal.json':'f2192cc54878147416e212990622d88e3866f2d2cbd478c221aa36c03691cb3d',
    WARM:'8f24cb37b0fa75a604a05066c354cd7f6d9663f0f08f038762a5c42f0de30c5a',
    KERNEL:'461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3',
    GL:'8bffaf8057ace5a2593861f879068f7e4a3e3f1ef72d2fecc4199b94c8372c52',
    CONDITIONAL:'184b3ec43d205290e0c4e2d8df0237a9df770af38116d62762a8b1f9811f3d76',
    V2:'289e93b651411042093a2b987717d87cdb5ed156ab7c736e1ac5cf820670fd4b',
    DYNAMIC:'57ff8c707f17f5cea2a1d7fd8fd76e2e5b57d666cf0512f17e01d8fc66790cea',
    ORIGINAL_RUNNER:'065fd1b0056edadac3842b89e5650e129d652f18a1fa251f75dc41985f4f293d',
    **FULL_CHECKS,
}


def require(ok,message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(1<<20),b''):
            h.update(block)
    return h.hexdigest()


def guard(values):
    for path,digest in values.items():
        require(sha(path)==digest,'Pinned source/input changed: '+str(path))


def metadata(args):
    require(args.expected_model_sha==MODEL_SHA,'Require unchanged reviewed batched model source')
    own = {Path(__file__).resolve():args.expected_runner_sha,MODEL:args.expected_model_sha}
    guard({**GUARDS,**own})
    campaign = json.loads((CAMPAIGN/'campaign.json').read_text())
    originals = {}
    for name in ('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py','check_bfix_tests_v13.py'):
        path = CAMPAIGN/'executed_sources'/name
        originals[path] = campaign['files'][str(path)]
    guard(originals)
    for path in (*own,KERNEL,GL,CONDITIONAL,V2,DYNAMIC,ORIGINAL_RUNNER,*originals):
        ast.parse(path.read_text())
    for directory in {path.parent for path in FULL_CHECKS}:
        summary = json.loads((directory/'summary.json').read_text())
        report = json.loads((directory/'dynamic_full_zero_law_comparison.json').read_text())
        require(summary['dynamic_full_zero_law_comparison_completed'] and summary['no_model_fits']
                and report['comparisons']==14 and report['parameter_points']==7
                and report['registered_rows']==173 and report['genes']==8 and report['exception_count']==0,
                'Pinned full synthetic comparison incomplete or scope changed')
        flags = ('all_LL_within_existing_allowance','all_gradients_within_existing_allowance',
                 'counts_unchanged','row_signatures_unchanged','seed_parameters_unchanged',
                 'final_point_exact_saved_seed','original_callables_and_inherited_method_identities_unchanged',
                 'all_required_refusals_observed','valid_context_after_refusal')
        require(all(report[name] is True for name in flags)
                and all(all(row['exact_equal']) for row in report['full_zero_law_comparisons']),
                'Pinned full scalar arithmetic or guard comparison did not pass')
    return campaign,{**GUARDS,**own,**originals}


def load_module(name,path):
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-runner-sha',required=True)
    parser.add_argument('--expected-model-sha',required=True)
    parser.add_argument('--batch-size',type=int,choices=(32,64),required=True)
    parser.add_argument('--inspect',action='store_true')
    args = parser.parse_args()
    campaign,guards = metadata(args)
    if args.inspect:
        print(json.dumps(dict(source_metadata_checks=True,numerical_execution=False,
            model_fit_run=False,prepared_synthetic_participants=60,prepared_rows=173,prepared_genes=8,
            requested_batch_size=args.batch_size,maximum_model_rows=256,
            maximum_prepared_count_axis_entries_per_gene=100000,
            inherited_v17_fitter_methods=True,explicit_positive_variance_GL_revision=True,
            original_boundary_law_inherited=True,biological_model_validated=False,
            prepared_dynamic_scalar_cache=True,passing_full_comparison_jobs=[22041698,22041699],
            original_observed_runner_sha256=GUARDS[ORIGINAL_RUNNER],
            statistical_calibration=False,independent_participant_validation=False,
            future558_participant_scaling_checked=False)))
        return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    for name in ('globals.json','individuals.tsv','rows.tsv','cohorts.tsv'):
        path = FIXTURE/name
        require(path.resolve().is_relative_to(FIXTURE.resolve()),'Fixture input escaped synthetic root')
        guards[path] = campaign['files'][str(path)]
    guard(guards)
    output = REC/('repair_p3_cached_observed_fit_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(mode=0o750,exist_ok=False)
    archive = output/'executed_sources';archive.mkdir()
    archive_guards = {}
    for path,digest in guards.items():
        if path.suffix=='.py':
            target = archive/path.name
            require(not target.exists(),'Duplicate source basename')
            target.write_bytes(path.read_bytes())
            require(sha(target)==digest,'Archived source copy differs')
            archive_guards[target] = digest
    for path,digest in FULL_CHECKS.items():
        directory = archive/'passing_full_checks'/path.parent.name
        directory.mkdir(parents=True,exist_ok=True)
        target = directory/path.name
        with target.open('xb') as handle:
            handle.write(path.read_bytes())
        require(sha(target)==digest,'Archived passing full-check receipt differs')
        archive_guards[target] = digest
    began,phase = time.monotonic(),'numerical_imports'
    events_seen,events_saved,events_omitted = 0,0,0
    np = None

    def safe(value):
        if np is not None and isinstance(value,np.ndarray):
            return safe(value.tolist())
        if np is not None and isinstance(value,np.generic):
            return safe(value.item())
        if isinstance(value,float):
            return value if value==value and abs(value)!=float('inf') else None
        if isinstance(value,dict):
            return {str(k):safe(v) for k,v in value.items()}
        if isinstance(value,(list,tuple)):
            return [safe(v) for v in value]
        return value

    def bytes_saved():
        return sum(p.stat().st_size for p in output.rglob('*') if p.is_file())

    def save(name,value,terminal=False):
        body = (json.dumps(safe(value),indent=2,allow_nan=False)+'\n').encode()
        require(bytes_saved()+len(body)<=CAP-(0 if terminal else RESERVE),'Observed-fit saved-output cap exceeded')
        with (output/(name+'.json')).open('xb') as handle:
            handle.write(body);handle.flush();os.fsync(handle.fileno())

    def event(fields):
        nonlocal events_seen,events_saved,events_omitted
        events_seen += 1
        fields = dict(phase=phase,elapsed_seconds=time.monotonic()-began,**fields)
        body = (json.dumps(safe(fields),allow_nan=False)+'\n').encode()
        if fields.get('event')=='marginal_evaluation' and bytes_saved()+len(body)>PROGRESS_SOFT_CAP:
            events_omitted += 1
            return
        require(bytes_saved()+len(body)<=CAP-RESERVE,'Terminal/refusal event cannot fit retained reserve')
        with (output/'progress.jsonl').open('ab') as handle:
            handle.write(body);handle.flush()
        events_saved += 1
        if fields.get('event')!='marginal_evaluation':
            print(json.dumps(safe(fields),allow_nan=False),flush=True)

    def coordinate_names(model,settings):
        return ([g+':alpha' for g in settings['gene_order']]+[g+':linear_r' for g in settings['gene_order']]
                +['kappa']+['delta:'+c for c in settings['design_columns']]
                +['omega_S','rho_S','rho_b','rho_d'])

    def diagnostic(model,fit,settings,marginal=False,null=None):
        p = np.asarray(fit['parameters'],float)
        if marginal:
            point = np.r_[p,float(fit['v'])]
            ll,g,h = model.marginal(p,float(fit['v']),with_hessian=True)
            error = model.last_total_error_allowance
            integration = model.last_integration_reports
            tolerance = 1e-10
            bounds = [b if flag else (null['parameters'][j],null['parameters'][j])
                      for j,(b,flag) in enumerate(zip(model.bounds(),model.informative()))]+[(0,None)]
            names = coordinate_names(model,settings)+['v']
        else:
            point = p
            ll,g,h = model.evaluate(p,need_hessian=True)
            error,integration,tolerance = None,[],1e-8
            bounds = [b if flag else (p[j],p[j]) for j,(b,flag) in
                      enumerate(zip(model.bounds(),model.informative()))]
            names = coordinate_names(model,settings)
        require(len(names)==len(point)==len(g) and h.shape==(len(point),len(point))
                and np.isfinite(point).all() and np.isfinite(g).all() and np.isfinite(ll),
                'Final diagnostic shape/finiteness differs')
        active = np.array([(lo is not None and lo==hi)
                           or (lo is not None and point[j]<=lo+tolerance and g[j]<0)
                           or (hi is not None and point[j]>=hi-tolerance and g[j]>0)
                           for j,(lo,hi) in enumerate(bounds)])
        projected = g.copy();projected[active] = 0.
        maximum = float(np.max(np.abs(projected),initial=0))
        boundary = marginal and fit['v']==0.
        nuisance_hessian = h[:-1,:-1] if boundary else h
        require(np.isfinite(nuisance_hessian).all(),'Nonfinite Hessian outside prescribed variance-boundary entries')
        coordinates = [dict(index=j,name=names[j],value=float(point[j]),score=float(g[j]),
            projected_score=float(projected[j]),active_by_original_rule=bool(active[j]),
            lower=bounds[j][0],upper=bounds[j][1],projection_bound_tolerance=tolerance)
            for j in range(len(point))]
        return dict(loglik=float(ll),gradient=g,hessian=h,hessian_finite_mask=np.isfinite(h),
            returned_fit_loglik=fit['loglik'],loglik_recomputation_difference=float(ll-fit['loglik']),
            returned_fit_gradient=fit['gradient'],
            gradient_recomputation_difference=g-np.asarray(fit['gradient']),
            recomputed_gradient_within_existing_fixed_vector_allowance=
                np.allclose(g,fit['gradient'],atol=1e-8,rtol=1e-6),
            coordinates=coordinates,max_projected_gradient=maximum,
            original_stationarity_cutoff=1e-6,projected_gradient_strictly_below_cutoff=maximum<1e-6,
            native_fit_converged=bool(fit['converged']),variance_boundary=boundary,
            native_flag_not_overridden_by_diagnostic=True,
            native_boundary_null_snap_may_use_inherited_null_projection_tolerance=True,
            right_variance_score=float(g[-1]) if marginal else None,
            right_variance_score_at_or_below_original_restart_cutoff=(float(g[-1])<=1e-6 if boundary else None),
            original_null_snap_condition=(bool(null['converged'] and g[-1]<=1e-6) if boundary else None),
            hessian_variance_entries_unavailable_by_inherited_boundary_method=boundary,
            hessian_symmetry_max_absolute_difference=float(np.max(np.abs(nuisance_hessian-nuisance_hessian.T))),
            raw_hessian_information_eigenvalues=np.linalg.eigvalsh(-(nuisance_hessian+nuisance_hessian.T)/2),
            numerical_error_estimates=error,integration_reports=integration,
            null_or_zero_boundary_floating_error_not_certified=error is None,
            stationarity_not_optimization_global_maximum_or_model_validation=True)

    model = None
    scalar_cache = None
    scalar_setup_seconds = None
    def scalar_cache_receipt():
        return None if scalar_cache is None else dict(setup_seconds=scalar_setup_seconds,
            memory_estimate=scalar_cache.memory_estimate,dispersion_cache_ledger=scalar_cache.ledger(),
            persistent_eta_cache=False,seed_used_for_static_initialization_only=True,
            scalar_signature='conditional_terms(row,eta,linear_r,need_hessian=True)',
            exact_latest_finite_float_r_invalidates_each_row=True)
    try:
        sys.path.insert(0,str(archive))
        require(not any(name in sys.modules for name in ('likelihood_v1','bfix_v1','design_v1','tests_v1')),
                'Original numerical modules already imported')
        import numpy as np
        import scipy
        loader = importlib.import_module('bfix_v1')
        original = importlib.import_module('likelihood_v1')
        importlib.import_module('tests_v1')
        for name in ('likelihood_v1','bfix_v1','design_v1','tests_v1'):
            require(Path(sys.modules[name].__file__).resolve().parent==archive.resolve(),
                    'Original import escaped guarded archive')
        environment = dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        require(environment==dict(python='3.10.18',numpy='2.2.6',scipy='1.15.2'),'Original numerical environment changed')
        design,settings = loader.load()
        base = design.model()
        require(len(design.participants)==60 and len(base.rows)==173 and base.G==8
                and base.Z==7 and base.size==28,'Exact observed synthetic fixture shape changed')
        original_bb = original.bb_terms
        original_conditional = original.conditional_terms
        private = load_module('_repair_p3_batched_observed_private_original',archive/'likelihood_v1.py')
        kernel = load_module('_repair_p3_batched_observed_kernel',archive/KERNEL.name)
        gl = load_module('_repair_p3_batched_observed_GL',archive/GL.name)
        conditional = load_module('_repair_p3_batched_observed_conditional',archive/CONDITIONAL.name)
        model_source = load_module('_repair_p3_batched_observed_model',archive/MODEL.name)
        Candidate = model_source.batched_model_class(private,kernel,gl,conditional,
            kernel_path=archive/KERNEL.name,batch_size=args.batch_size,event_sink=event)
        model = Candidate(base.rows,gene_count=base.G)
        require(original.bb_terms is original_bb,'Private candidate leaked into unmodified reference')
        inherited = {name:getattr(Candidate,name) is getattr(private.Model,name) for name in model_source.INHERITED}
        require(all(inherited.values()),'Original fitter/design method identity changed')
        warm_record = json.loads(WARM.read_text())['fit']
        require(warm_record['converged'],'Original observed warm-start fit is not converged')
        warm = np.asarray(warm_record['parameters'],float)
        phase = 'scalar_conditional_preparation'
        dynamic = load_module('_repair_p3_cached_observed_dynamic',archive/DYNAMIC.name)
        setup_start = time.perf_counter()
        scalar_cache = dynamic.DynamicConditional(model,warm,model.rows)
        scalar_setup_seconds = time.perf_counter()-setup_start
        private_bb = private.bb_terms
        private_original_conditional = dynamic.install_private(scalar_cache)
        private_shim = private.conditional_terms
        require(private is not original and scalar_cache.original is private
                and private_shim.__self__ is scalar_cache
                and private_shim.__func__ is dynamic.DynamicConditional.terms
                and private_original_conditional is not private_shim,
                'Scalar conditional escaped intended private model')
        require(original.bb_terms is original_bb and original.conditional_terms is original_conditional
                and all(getattr(Candidate,name) is getattr(private.Model,name) for name in model_source.INHERITED),
                'Scalar install changed unmodified reference or inherited methods')
        old_null = json.loads((CAMPAIGN/'baseline/observed_null.json').read_text())
        old_marginal = json.loads((CAMPAIGN/'baseline/observed_marginal.json').read_text())
        counts = hashlib.sha256(np.asarray([row.a for row in model.rows],dtype='<i8').tobytes()).hexdigest()
        save('protocol',dict(source_guards={str(p):h for p,h in guards.items()},environment=environment,
            participants=60,participant_relative_clusters=len(set(design.clusters.values())),
            rows=173,genes=8,parameters=28,variance_coordinate=28,observed_count_vector_sha256=counts,
            fixture=str(FIXTURE),fixture_seed=settings['seed'],warm_start_source=str(WARM),
            original_observed_runner_sha256=GUARDS[ORIGINAL_RUNNER],
            passing_full_comparison_receipts={str(p):h for p,h in FULL_CHECKS.items()},
            scalar_conditional_source_sha256=GUARDS[DYNAMIC],scalar_cache_initial=scalar_cache_receipt(),
            private_conditional_is_expected_dynamic_shim=True,
            unmodified_reference_conditional_identity_preserved=True,
            observed_recipe='model.fit(original prior observed null parameters); model.fit_marginal(tau_max_fixture,new null parameters,null_fit=new null)',
            tau_max_fixture=settings['tau_max_fixture'],inherited_methods=inherited,requested_batch_size=args.batch_size,
            variance_bounds='v>=0; tau_max_fixture is original cold-start scale, not an imposed upper variance bound',
            original_zero_boundary_law=True,positive_variance_approximation_revision=True,
            original20node_approximation_agreement_not_asserted=True,
            maximum_model_rows=256,maximum_prepared_count_axis_entries_per_gene=100000,
            saved_output_cap_bytes=CAP,only_synthetic_observed_fit=True,
            progress_soft_cap_bytes=PROGRESS_SOFT_CAP,terminal_reserved_bytes=RESERVE,
            no_simulation_draws_or_calibration=True,no_biological_private_inputs=True,
            generating_truth_not_used_as_start=True,future558_participant_scaling_checked=False,
            independent_participant_validation=False,model_validated=False))
        phase = 'zero_law_baseline_comparison'
        require(private.conditional_terms is private_shim and private.bb_terms is private_bb,
                'Private scalar conditional or stable BB changed before baseline')
        old_p = np.asarray(old_null['parameters'],float)
        reference = original.Model(base.rows,gene_count=base.G)
        ref_ll,ref_g,ref_h = reference.evaluate(old_p)
        new_ll,new_g,new_h = model.evaluate(old_p)
        save('original_zero_law_baseline',dict(reference_loglik=ref_ll,candidate_loglik=new_ll,
            loglik_difference=new_ll-ref_ll,gradient_difference=new_g-ref_g,hessian_difference=new_h-ref_h,
            parameters=old_p,existing_fixed_vector_LL_allowance=1e-6,
            existing_fixed_vector_gradient_atol=1e-8,existing_fixed_vector_gradient_rtol=1e-6,
            LL_within_original_allowance=abs(new_ll-ref_ll)<=1e-6,
            gradient_within_original_allowance=np.allclose(new_g,ref_g,atol=1e-8,rtol=1e-6),
            same_zero_variance_probability_law_not_bit_exact_arithmetic=True))
        phase = 'observed_null_fit'
        event(dict(event='fit_started',component='null'))
        start = time.perf_counter()
        null = model.fit(warm)
        save('observed_null',dict(fit=null,seconds=time.perf_counter()-start))
        null_diagnostic = diagnostic(model,null,settings)
        save('observed_null_diagnostic',null_diagnostic)
        require(null['converged'],'Observed null fit nonconverged; original stop retained')
        phase = 'observed_marginal_fit'
        event(dict(event='fit_started',component='marginal'))
        start = time.perf_counter()
        fit = model.fit_marginal(settings['tau_max_fixture'],np.asarray(null['parameters']),null_fit=null)
        save('observed_marginal',dict(fit=fit,seconds=time.perf_counter()-start))
        phase = 'complete_selected_final_diagnostics'
        final = diagnostic(model,fit,settings,marginal=True,null=null)
        save('observed_marginal_diagnostic',final)
        zero_comparison = None
        if fit['v']==0.:
            fixed = np.asarray(fit['parameters'])
            ref_ll,ref_g = reference.marginal(fixed,0.)
            new_ll,new_g = model.marginal(fixed,0.)
            zero_comparison = dict(parameters=fixed,reference_loglik=ref_ll,candidate_loglik=new_ll,
                same_parameter_loglik_difference=new_ll-ref_ll,full_gradient_difference=new_g-ref_g,
                right_boundary_variance_score_reference=ref_g[-1],right_boundary_variance_score_candidate=new_g[-1],
                original_observed_marginal_loglik=old_marginal['loglik'],
                candidate_minus_original_observed_fit_loglik=new_ll-old_marginal['loglik'],
                comparison_of_distinct_fit_points_not_fixed_vector_accuracy=True,
                LL_within_original_fixed_vector_allowance=abs(new_ll-ref_ll)<=1e-6,
                gradient_within_original_fixed_vector_allowance=np.allclose(new_g,ref_g,atol=1e-8,rtol=1e-6))
            save('selected_zero_boundary_comparison',zero_comparison)
        guard(guards);guard(archive_guards)
        require(original.bb_terms is original_bb and original.conditional_terms is original_conditional,
                'Unmodified reference changed during fitting')
        require(private.conditional_terms is private_shim and private.bb_terms is private_bb
                and all(getattr(Candidate,name) is getattr(private.Model,name) for name in model_source.INHERITED),
                'Private scalar conditional or inherited method changed during fitting')
        save('scalar_conditional_cache_final',scalar_cache_receipt())
        event(dict(event='fit_records_completed',native_null_converged=bool(null['converged']),
                   native_marginal_converged=bool(fit['converged'])))
        save('summary',dict(status='complete_synthetic_observed_fit_records',participants=60,rows=173,genes=8,
            native_null_converged=bool(null['converged']),native_marginal_converged=bool(fit['converged']),
            final_projected_gradient=final['max_projected_gradient'],
            final_projected_stationarity=final['projected_gradient_strictly_below_cutoff'],
            final_variance=fit['v'],right_variance_boundary_score=final['right_variance_score'],
            start_logliks=fit['start_logliks'],start_converged=fit['start_converged'],kept_start=fit['kept_start'],
            limitation='Inherited returns retain warm/cold LL and convergence flags, not both start parameter paths',
            original_fitting_methods_unchanged=True,original_reference_module_unchanged=True,inputs_unchanged=True,
            private_conditional_is_expected_dynamic_shim=True,scalar_conditional_cache=scalar_cache_receipt(),
            marginal_evaluation_calls=model.marginal_call_count,events_seen=events_seen,
            events_saved=events_saved,events_omitted_due_to_cap=events_omitted,
            total_seconds=time.monotonic()-began,no_calibration_or_draw_expansion=True,
            no_statistical_tests_or_new_thresholds=True,biological_model_validated=False,
            independent_participant_validation=False,actual558_participant_scaling_checked=False,
            numerical_error_bounds_certified=False,
            saved_output_bytes_before_summary=bytes_saved(),
            output_sha256={str(p.relative_to(output)):sha(p) for p in output.rglob('*') if p.is_file()}),terminal=True)
    except BaseException as exc:
        failure = dict(phase=phase,error_type=type(exc).__name__,message=str(exc),
            traceback=traceback.format_exc(),elapsed_seconds=time.monotonic()-began,outputs_preserved=True,
            events_seen=events_seen,events_saved=events_saved,events_omitted_due_to_cap=events_omitted,
            model_validated=False)
        failure['scalar_conditional_cache'] = scalar_cache_receipt()
        if model is not None:
            failure['last_integration_failure'] = getattr(model,'last_integration_failure',None)
            if failure['last_integration_failure'] is None:
                point = getattr(model,'last_evaluated_point',{})
                failure['last_integration_failure'] = dict(call=point.get('call'),gene=None,
                    variance=point.get('variance'),parameters=point.get('parameters'),
                    completed_gene_reports=getattr(model,'last_integration_reports',[]),
                    refused_gene_ladders=None,error_type=type(exc).__name__,message=str(exc))
        full_body = (json.dumps(safe(failure),indent=2,allow_nan=False)+'\n').encode()
        if bytes_saved()+len(full_body)>CAP:
            detail = failure.get('last_integration_failure')
            failure['numeric_failure_details_omitted_due_to_output_cap'] = True
            failure['full_failure_payload_sha256'] = hashlib.sha256(full_body).hexdigest()
            if detail is not None:
                failure['last_integration_failure'] = {
                    key:detail.get(key) for key in ('call','gene','variance','parameters','error_type','message')}
                failure['last_integration_failure']['completed_gene_report_count'] = len(
                    detail.get('completed_gene_reports') or [])
                failure['last_integration_failure']['refused_gene_ladder_count'] = len(
                    detail.get('refused_gene_ladders') or [])
            failure['omission_limit'] = 'Detailed numerical ladders may be absent; terminal refusal and counts retained; no successful-fit claim'
        else:
            failure['numeric_failure_details_omitted_due_to_output_cap'] = False
        save('failure',failure,terminal=True)
        raise


if __name__=='__main__':
    main()
