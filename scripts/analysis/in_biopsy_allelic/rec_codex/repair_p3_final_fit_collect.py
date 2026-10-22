#!/usr/bin/env python3
"""Compare saved sequential/parallel synthetic fits using bounded JSON only.

No progress trajectories, molecular arrays, numerical imports or scheduler calls.
Missing/failed terminal results are retained; accuracy is never autoaccepted.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path

REC=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z')
READ_CAP=64*1024**2
OUTPUT_CAP=32*1024**2
JOBS={
 'sequential':dict(job='22040873',directory='repair_p3_batched_observed_fit_22040873',
    launch='repair_p3_batched_observed_fit_launch_20261005T152152Z',script='observed_fit.sbatch',
    script_sha='9e0f764a8b3b16fe147a398a6964645cdcd7b6cd9d49c3a0e61b81498426cbd7',
    runner='repair_p3_batched_observed_fit.py',runner_sha='065fd1b0056edadac3842b89e5650e129d652f18a1fa251f75dc41985f4f293d',
    protocol_sha='bd38c3d80bf62fd32f8a642874cd31586cb6df21407ade9d08ea47ef4ab28e4c'),
 'parallel':dict(job='22042114',directory='repair_p3_parallel_observed_fit_22042114',
    launch='repair_p3_parallel_observed_fit_launch_20261005T161858Z',script='parallel_fit.sbatch',
    script_sha='f7a182bbca2f13416e71374a471da59b38737da7ad2eb9b3444ef52504997197',
    runner='repair_p3_parallel_observed_fit.py',runner_sha='ded982c1114e9ba4241e37472cc05a2ff4682f254f267595357b0c7419019966',
    protocol_sha='05c69cdb3aa293a7b1150bfc16c0a28487d76d731307236341e4e2c0de5b0604'),
}
IDENTITY=('environment','participants','participant_relative_clusters','rows','genes','parameters','variance_coordinate',
 'observed_count_vector_sha256','fixture','fixture_seed','warm_start_source','observed_recipe','tau_max_fixture',
 'inherited_methods','requested_batch_size','variance_bounds','original_zero_boundary_law','positive_variance_approximation_revision')
RESULT_NAMES=('original_zero_law_baseline.json','observed_null.json','observed_null_diagnostic.json',
 'observed_marginal.json','observed_marginal_diagnostic.json','selected_zero_boundary_comparison.json',
 'scalar_conditional_cache_final.json','parallel_gene_backend_final.json','summary.json','failure.json')


def require(ok,message):
    if not ok:raise ValueError(message)


def finite(value):
    return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)


class Reader:
    def __init__(self):self.bytes_read=0

    def bytes(self,path,expected=None):
        require(path.resolve().is_relative_to(REC.resolve()) and not path.is_symlink(),'Input escaped synthetic REC')
        size=path.stat().st_size
        require(size<=8*1024**2 and self.bytes_read+size<=READ_CAP,'Small receipt/source input cap reached')
        body=path.read_bytes();self.bytes_read+=len(body)
        require(len(body)==size and self.bytes_read<=READ_CAP,'Receipt changed size during reading')
        digest=hashlib.sha256(body).hexdigest()
        require(expected is None or digest==expected,'Pinned source/protocol differs: '+str(path))
        return body,dict(sha256=digest,bytes=len(body))

    def read(self,path,expected=None):
        require(path.suffix=='.json','Only JSON result files allowed')
        body,receipt=self.bytes(path,expected)
        def invalid(value):raise ValueError('Nonfinite JSON constant '+value)
        value=json.loads(body,parse_constant=invalid)
        require(isinstance(value,dict),'JSON receipt object required')
        return value,receipt


def metadata(reader):
    jobs={}
    for label,spec in JOBS.items():
        launch=REC/spec['launch'];directory=REC/spec['directory']
        _,script_receipt=reader.bytes(launch/spec['script'],spec['script_sha'])
        _,runner_receipt=reader.bytes(launch/spec['runner'],spec['runner_sha'])
        protocol,receipt=reader.read(directory/'protocol.json',spec['protocol_sha'])
        require(protocol['participants']==60 and protocol['participant_relative_clusters']==58 and protocol['rows']==173 and
            protocol['genes']==8 and protocol['parameters']==28 and protocol['variance_coordinate']==28 and
            protocol['observed_count_vector_sha256']=='7009be4adeec6157fa4bbf3f46e48ec24aa398a0bee4d7906c78c231b7d3fbad' and
            protocol['requested_batch_size']==32 and all(protocol['inherited_methods'].values()),'Frozen observed fit design differs')
        require(protocol['environment']==dict(python='3.10.18',numpy='2.2.6',scipy='1.15.2') and
            protocol['only_synthetic_observed_fit'] and protocol['no_simulation_draws_or_calibration'] and
            protocol['no_biological_private_inputs'],'Fit environment/scientific scope differs')
        require(protocol['source_guards'][str(launch/spec['runner'])]==spec['runner_sha'],'Executed runner declaration differs')
        jobs[label]=dict({**spec,'directory':str(directory)},protocol=protocol,protocol_receipt=receipt,
            runner_receipt=runner_receipt,launcher_receipt=script_receipt)
    identities={k:dict(sequential=jobs['sequential']['protocol'][k],parallel=jobs['parallel']['protocol'][k],
        equal=jobs['sequential']['protocol'][k]==jobs['parallel']['protocol'][k]) for k in IDENTITY}
    common={}
    for label,job in jobs.items():
        normalized={Path(p).name:d for p,d in job['protocol']['source_guards'].items()}
        common[label]=normalized
    source_comparisons=[dict(name=n,sequential_sha256=common['sequential'][n],parallel_sha256=common['parallel'][n],
        equal=common['sequential'][n]==common['parallel'][n]) for n in sorted(common['sequential'].keys()&common['parallel'].keys())]
    return jobs,identities,source_comparisons


def diagnostic_checks(fit,diagnostic,marginal):
    errors=[];checks={}
    if fit is None:return dict(available=False,reason='fit_missing'),errors
    try:
        p=fit['parameters'];g=fit['gradient'];width=29 if marginal else 28
        require(len(p)==28 and all(finite(x) for x in p) and len(g)==width and all(finite(x) for x in g) and
            finite(fit['loglik']) and isinstance(fit['converged'],bool),'Final fit shape/finiteness differs')
        if marginal:require(finite(fit['v']) and fit['v']>=0,'Invalid saved variance')
        starts=fit.get('start_logliks',{});start_flags=fit.get('start_converged',{})
        expected_kept=None
        if set(starts)=={'warm','cold'} and all(finite(x) for x in starts.values()):
            expected_kept='cold' if starts['cold']-starts['warm']>=1e-6 else 'warm'
        checks.update(available=True,native_converged=fit['converged'],kept_start=fit.get('kept_start'),
            start_logliks=starts,start_converged=start_flags,two_initial_start_records_present=expected_kept is not None,
            expected_initial_start_from_saved_LL=expected_kept,
            boundary_restart_label=fit.get('kept_start')=='boundary_restart',
            final_kept_label_matches_initial_rule=(fit.get('kept_start')==expected_kept if expected_kept else None),
            both_start_parameter_trajectories_not_saved=True)
        if diagnostic is None:
            checks['diagnostic_available']=False
            return checks,errors
        coords=diagnostic['coordinates'];dg=diagnostic['gradient']
        require(len(coords)==width and len(dg)==width and all(finite(x) for x in dg),'Diagnostic full coordinate shape differs')
        point=p+([fit['v']] if marginal else [])
        projected=[]
        for index,coord in enumerate(coords):
            require(coord['index']==index and coord['value']==point[index] and coord['score']==dg[index],
                'Diagnostic coordinate join differs')
            low,high,tol=coord['lower'],coord['upper'],coord['projection_bound_tolerance']
            active=((low is not None and low==high) or (low is not None and point[index]<=low+tol and dg[index]<0) or
                (high is not None and point[index]>=high-tol and dg[index]>0))
            score=0. if active else dg[index]
            require(coord['active_by_original_rule'] is active and coord['projected_score']==score,'Saved KKT projection differs')
            projected.append(score)
        maximum=max(abs(x) for x in projected)
        require(diagnostic['max_projected_gradient']==maximum and diagnostic['original_stationarity_cutoff']==1e-6 and
            diagnostic['projected_gradient_strictly_below_cutoff'] is (maximum<1e-6),'Stationarity summary differs')
        delta=diagnostic['loglik']-fit['loglik']
        require(diagnostic['loglik_recomputation_difference']==delta and diagnostic['native_fit_converged'] is fit['converged'],
            'Diagnostic fit identity differs')
        boundary=marginal and fit['v']==0
        if boundary:require(diagnostic['right_variance_score']==dg[-1] and
            diagnostic['right_variance_score_at_or_below_original_restart_cutoff'] is (dg[-1]<=1e-6),'Boundary score condition differs')
        differences=[a-b for a,b in zip(dg,g)]
        require(diagnostic['gradient_recomputation_difference']==differences,'Recomputed gradient differences differ')
        checks.update(diagnostic_available=True,recomputed_loglik_difference=delta,recomputed_full_gradient_differences=differences,
            max_projected_gradient_rederived=maximum,projected_strictly_below_original_cutoff=maximum<1e-6,
            variance_boundary=boundary,right_variance_score=dg[-1] if marginal else None,
            boundary_right_score_at_or_below_original_cutoff=dg[-1]<=1e-6 if boundary else None,
            native_flag_not_overridden=True,diagnostic_same_candidate_not_independent_reference=True)
    except Exception as error:errors.append(dict(reason='saved_fit_or_diagnostic_invalid',error=str(error)))
    return checks,errors


def collect_job(reader,job):
    directory=Path(job['directory']);documents={};hashes={};missing=[];errors=[]
    optional_names={'selected_zero_boundary_comparison.json','scalar_conditional_cache_final.json','parallel_gene_backend_final.json','failure.json'}
    for name in RESULT_NAMES:
        path=directory/name
        if not path.exists():
            if name not in optional_names:missing.append(name)
            continue
        try:
            documents[name],hashes[name]=reader.read(path)
        except Exception as error:errors.append(dict(file=name,error=str(error)))
    fits={name:documents.get('observed_'+name+'.json',{}).get('fit') for name in ('null','marginal')}
    checks={}
    for name,fit in fits.items():
        checks[name],issues=diagnostic_checks(fit,documents.get('observed_'+name+'_diagnostic.json'),name=='marginal')
        errors.extend(dict(component=name,**issue) for issue in issues)
    summary=documents.get('summary.json');failure=documents.get('failure.json')
    if summary is not None:
        try:
            require(summary['status']=='complete_synthetic_observed_fit_records' and summary['participants']==60 and
                summary['rows']==173 and summary['genes']==8 and summary['original_fitting_methods_unchanged'] and
                summary['original_reference_module_unchanged'] and summary['inputs_unchanged'] and
                summary['biological_model_validated'] is False,'Final summary scope differs')
            require(summary['native_null_converged'] is fits['null']['converged'] and
                summary['native_marginal_converged'] is fits['marginal']['converged'] and
                summary['kept_start']==fits['marginal']['kept_start'] and summary['start_logliks']==fits['marginal']['start_logliks'] and
                summary['start_converged']==fits['marginal']['start_converged'] and summary['final_variance']==fits['marginal']['v'],
                'Final summary/fit identity differs')
            for name,receipt in hashes.items():
                if name!='summary.json' and name in summary['output_sha256']:
                    require(summary['output_sha256'][name]==receipt['sha256'],'Final saved result hash differs')
        except Exception as error:errors.append(dict(reason='summary_invalid',error=str(error)))
    state=('failure_receipt' if failure is not None else 'final_summary_saved' if summary is not None else
        'marginal_saved_summary_missing' if fits['marginal'] is not None else 'null_saved_marginal_missing' if fits['null'] is not None else
        'fits_missing')
    return dict(**job,documents=documents,result_hashes=hashes,missing_required_receipts=missing,observed_receipt_state=state,
        saved_fit_checks=checks,structural_or_arithmetic_errors=errors,final_summary_present=summary is not None,
        failure_receipt_present=failure is not None,scheduler_state_not_checked=True,accuracy_accepted=False)


def compare_component(jobs,name):
    a=jobs['sequential']['documents'].get('observed_'+name+'.json',{}).get('fit')
    b=jobs['parallel']['documents'].get('observed_'+name+'.json',{}).get('fit')
    if a is None or b is None:return dict(comparison_available=False,missing=[label for label,f in (('sequential',a),('parallel',b)) if f is None])
    result=dict(comparison_available=True,component=name,raw_sequential_fit=a,raw_parallel_fit=b,
        convergence_flags=dict(sequential=a['converged'],parallel=b['converged']),
        start_metadata_equal={k:a.get(k)==b.get(k) for k in ('kept_start','start_logliks','start_converged')},
        comparisons_not_independent_probability_accuracy=True,full_trajectory_agreement_unverified=True)
    try:
        require(len(a['parameters'])==len(b['parameters'])==28,'Parameter comparison shape differs')
        result['parameter_differences_parallel_minus_sequential']=[y-x for x,y in zip(a['parameters'],b['parameters'])]
        result['loglik_difference_parallel_minus_sequential']=b['loglik']-a['loglik']
        result['loglik_within_existing_absolute_1e_minus6_allowance']=abs(b['loglik']-a['loglik'])<=1e-6
        result['kappa_difference_parallel_minus_sequential']=b['parameters'][16]-a['parameters'][16]
        result['kappa_within_existing_absolute_1e_minus4_allowance']=abs(result['kappa_difference_parallel_minus_sequential'])<=1e-4
        require(len(a['gradient'])==len(b['gradient'])==(29 if name=='marginal' else 28),'Gradient comparison shape differs')
        result['gradient_comparisons']=[dict(index=i,sequential=x,parallel=y,difference=y-x,
            existing_fixed_vector_tolerance=1e-8+1e-6*abs(x),within_nominal_tolerance=abs(y-x)<=1e-8+1e-6*abs(x))
            for i,(x,y) in enumerate(zip(a['gradient'],b['gradient']))]
        result['identical_parameter_vectors']=a['parameters']==b['parameters']
        if name=='marginal':
            require(finite(a['v']) and finite(b['v']) and a['v']>=0 and b['v']>=0,'Invalid variance comparison')
            result['variance_difference_parallel_minus_sequential']=b['v']-a['v']
            result['tau_difference_parallel_minus_sequential']=math.sqrt(b['v'])-math.sqrt(a['v'])
            result['same_fixed_parameter_point']=result['identical_parameter_vectors'] and a['v']==b['v']
            nulls=[jobs[label]['documents'].get('observed_null.json',{}).get('fit') for label in ('sequential','parallel')]
            if all(nulls):
                lrs=[2*(max(f['loglik'],n['loglik'])-n['loglik']) for f,n in zip((a,b),nulls)]
                result['saved_fit_LRs_rederived']=dict(sequential=lrs[0],parallel=lrs[1])
                result['LR_difference_parallel_minus_sequential']=lrs[1]-lrs[0]
                result['existing_tau_or_LR_agreement_rule']=('tau_absolute1e-3' if max(lrs)>1 else 'LR_absolute1e-6')
                result['tau_or_LR_nominal_agreement']=(abs(result['tau_difference_parallel_minus_sequential'])<=1e-3
                    if max(lrs)>1 else abs(lrs[1]-lrs[0])<=1e-6)
        else:result['same_fixed_parameter_point']=result['identical_parameter_vectors']
        result['score_comparison_at_distinct_fit_points_not_fixed_vector_accuracy']=not result['same_fixed_parameter_point']
    except Exception as error:result['comparison_error']=str(error)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--inspect',action='store_true')
    args=parser.parse_args();reader=Reader();specs,identities,sources=metadata(reader)
    if args.inspect:
        print(json.dumps(dict(job_ids=[j['job'] for j in specs.values()],protocol_and_runner_pins_verified=True,
            common_scientific_identity_equal=all(r['equal'] for r in identities.values()),
            common_source_hashes_equal=all(r['equal'] for r in sources),numerical_imports=False,
            result_collection_not_run=True,progress_not_read=True,input_cap=READ_CAP,output_cap=OUTPUT_CAP)))
        return
    jobs={label:collect_job(reader,spec) for label,spec in specs.items()}
    comparisons={name:compare_component(jobs,name) for name in ('null','marginal')}
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    result=dict(collected_utc=stamp,collector_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),jobs=jobs,
        scientific_identity_comparisons=identities,common_source_comparisons=sources,final_fit_comparisons=comparisons,
        both_final_summaries_present=all(j['final_summary_present'] for j in jobs.values()),
        failed_jobs_with_saved_receipts=[label for label,j in jobs.items() if j['failure_receipt_present']],
        missing_final_result_jobs=[label for label,j in jobs.items() if not j['final_summary_present']],
        progress_files_not_read=True,full_trajectory_agreement_unverified=True,both_start_parameter_paths_not_saved=True,
        positive_variance_integral_is_approximation_revision=True,frozen20node_equivalence_not_claimed=True,
        final_score_recomputation_same_candidate_not_independent_verification=True,
        numerical_uncertainties_not_certified=True,accuracy_accepted=False,statistical_calibration=False,biological_validation=False,
        scheduler_state_not_checked=True,input_bytes=reader.bytes_read,input_cap=READ_CAP,output_cap=OUTPUT_CAP)
    body=json.dumps(result,indent=2,allow_nan=False)+'\n'
    require(len(body.encode())<=OUTPUT_CAP,'Comparison output cap reached')
    output=REC/('repair_p3_final_fit_collection_'+stamp);output.mkdir(mode=0o750)
    partial=output/'collection.json.partial'
    with partial.open('x') as handle:handle.write(body)
    partial.rename(output/'collection.json')
    print(json.dumps(dict(output=str(output),both_final_summaries_present=result['both_final_summaries_present'],
        failed_jobs=result['failed_jobs_with_saved_receipts'],missing_final_jobs=result['missing_final_result_jobs'],accuracy_accepted=False)))


if __name__=='__main__':main()
