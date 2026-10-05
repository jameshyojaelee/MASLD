#!/usr/bin/env python3
"""JSON-only collection of the 14 fixed complex-axis jobs; never accepts accuracy."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
LAUNCH = REC/'repair_p3_complex_launch_20261005T141156Z'
AXES = (2,10,16,17,24,25,28)
STEPS = (1e-10,1e-15,1e-20)
BOUNDS = (24.,26.)
FD_STEPS = (1e-3,3e-4,1e-4)
CAP = 16*1024**2
JSON_CAP = 4*1024**2
SUBMISSIONS = {
 '22039439':'d3f6ecb93df288351bce982b0a320ac378c4e8b094ebc064bf214bcefe06a069',
 '22039440':'ef4221362789f65fe3db2390695b376b59c1fcdc26c16e9c2376b952a767c320',
 '22039441':'b4df2207630a78c509f4d266f39054c7c47f1213a8184e41a223d27ebd47fca0',
 '22039442':'d323ecc405205e4a61b0f05b5d7757ed555f5e791dd1cbd62125b83215e190c6',
 '22039443':'7425d2c15f8bf33d5a72d93dba1a703f0c3ab75f1fab589941bbbea24caf8f05',
 '22039444':'6716515420417962f208e2d022ea6dc153d34fff2ced850e4434b7b0a8b51f17',
 '22039445':'003dbb04f8566cc8cf871b0ace633c4684bec62be5620a62fae20ae26f79b779',
 '22039446':'ba0310b64f0ad59cc67cb6e43dc260b4d81f4c4ff7cb0e2b22fdf98f0e9c631a',
 '22039447':'2024c19a1c7f7445cc4d9d92ddd07e5e2e7f2535890722a5dad2ece7a1f4b97b',
 '22039448':'0e6d9569fe1ac77b1fae66d77521e61680f90a9e336c1479e9712153ad9150ff',
 '22039449':'8f65ade6f386ce2b5108bca9dfa466e041f185b5cc256770f0c2ec3a96737225',
 '22039450':'d74ac5bf8d5985c0b4799e32e93b01e498e449eefa08e9bdd9aadbcf15529f81',
 '22039451':'17dca64ecd9ff67941be9b66e528479eddf4cec5ae4994828153c8805a4903a7',
 '22039452':'c3496cc6baa1b7f75fa79cec95825986af40a4b0ae8d71fde31474e12683a895',
}
CASE_DIRS = {
 'failed0gene2':('repair_p3_scalar_prefix_failed0gene2_22039118','repair_p3_fixed_point_failed0gene2_22039200',0),
 'positive18':('repair_p3_scalar_prefix_positive18_22039119','repair_p3_fixed_point_positive18_22039201',18),
}
CASE_PINS = {
 'failed0gene2': {
  'scalar_protocol':'b27a85f6c11342406e07e99aacafa70839d09fcd870e0449dac40e59aab97a32',
  'scalar_summary':'74f375cf748940e4bf91137e55ab6e2e8d564b7fedb48032c19d60c19e688d53',
  'scalar_fd':'2ae8d822a8e456ad6e247ffaf91bcc21f5a90e06b63c66c8d2846a0cda96451a',
  'composite_protocol':'7b5d793be5485a25e5fc7f5b3d97c4ddf30a1a497d697728cd64468ed61f8d40',
  'composite_summary':'5ad1ea571c71525e6ae1eb18bf625baf1e6779934d64a149e0d16ca5f3137a0b',
  'composite_baseline':'edec063451383f9a1386619186c2089d3ef8e6ee34d854e80a2894fdd28e0f9a'},
 'positive18': {
  'scalar_protocol':'26f061f8e481a1b9d1e11f2a5dfc0984c82010d3b1dd9c8f270eb1418d701001',
  'scalar_summary':'fba2f38821c59e85ea4973e5ebe8cfe2509bc35d9910fcb7d9638078d6adb43d',
  'scalar_fd':'ebce73a49bbca829ff182e4cb0fc710c18fea794edbe8f8e476a7d145c91ac62',
  'composite_protocol':'ed6c190424fa29501e1377b5b47b12762917e53f46338c3c3f1a408c3f9e4774',
  'composite_summary':'1daa3275c8d524857c90369b7686fc2234229df765e41f139a4b232c9e323613',
  'composite_baseline':'4a576370be2f33ac28bf23052fe283b4a5f614bc30e983bd6a84224da857c722'},
}


def require(ok,message):
    if not ok:
        raise ValueError(message)


def finite_number(value):
    require(isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value),'Nonfinite numeric receipt field')
    return value


def read_json(path,expected=None):
    require(path.suffix=='.json' and path.stat().st_size<=JSON_CAP,'Only bounded JSON receipts permitted')
    body=path.read_bytes()
    digest=hashlib.sha256(body).hexdigest()
    require(expected is None or digest==expected,'Pinned JSON receipt differs: '+str(path))
    def invalid_constant(value):
        raise ValueError('Nonfinite JSON constant: '+value)
    return json.loads(body,parse_constant=invalid_constant),digest


def metadata():
    pins,_=read_json(LAUNCH/'source_pins.json','155e40b826a57db669eabf19367486103ed0bdebbd2552c55a583d72cd20cc8a')
    scripts={r['path']:r for r in pins['scripts']}
    require(len(scripts)==14 and {(r['case'],r['axis']) for r in scripts.values()}==
            {(c,a) for c in CASE_DIRS for a in AXES},'Frozen 14-script identity family differs')
    jobs=[]
    for job,expected in SUBMISSIONS.items():
        receipt,digest=read_json(LAUNCH/('submission_'+job+'.json'),expected)
        require(str(receipt['job_id'])==job and receipt['script'] in scripts,'Submission identity differs')
        spec=scripts[receipt['script']]
        require(receipt['source_pins']==pins['source_sha256'] and receipt['no_fit'] and receipt['no_calibration'],
                'Submission source/scientific scope differs')
        require(receipt['cpus']==1 and receipt['memory']=='4G' and receipt['time']=='90:00:00','Submission resource contract differs')
        require(Path(receipt['output']).name==f"{spec['case']}_axis{spec['axis']:02d}_{job}.out",'Submission log identity differs')
        jobs.append(dict(job_id=job,case=spec['case'],axis=spec['axis'],submission=receipt,submission_sha256=digest,
                         result_directory=str(REC/('repair_p3_complex_'+spec['case']+'_'+job))))
    require(len({(r['case'],r['axis']) for r in jobs})==14,'Duplicate submitted case/axis')
    return pins,jobs


def references():
    result={}
    for case,(scalar,composite,draw) in CASE_DIRS.items():
        paths={'scalar_protocol':REC/scalar/'protocol.json','scalar_summary':REC/scalar/'summary.json',
               'scalar_fd':REC/scalar/'finite_differences.json','composite_protocol':REC/composite/'protocol.json',
               'composite_summary':REC/composite/'summary.json','composite_baseline':REC/composite/'candidate_baseline.json'}
        loaded={key:read_json(path,CASE_PINS[case][key])[0] for key,path in paths.items()}
        require(loaded['scalar_summary']['scalar_fixed_point_check_completed'] and
                loaded['composite_summary']['candidate_integral_returned'],'Fixed baseline not complete')
        require(len(loaded['scalar_fd'])==21 and Counter((r['axis'],r['step']) for r in loaded['scalar_fd'])==
                Counter((a,h) for a in AXES for h in FD_STEPS),'Original full FD family differs')
        identity=('case','draw','fixed_gene','fixed_gene_rows','fixed_p','fixed_v','count_vector_sha256','environment')
        require(all(loaded['scalar_protocol'][k]==loaded['composite_protocol'][k] for k in identity),'Scalar/composite fixed point differs')
        require(loaded['scalar_protocol']['case']==case and loaded['scalar_protocol']['draw']==draw,'Wrong baseline case')
        result[case]=dict(protocol=loaded['scalar_protocol'],original_FD_records=loaded['scalar_fd'],
            composite_gradient=loaded['composite_baseline']['gradient'],
            composite_uncertified_gradient_errors=loaded['composite_baseline']['integration']['empirical_refinement_plus_tail_allowances'][1],
            source_paths={k:str(p) for k,p in paths.items()},source_sha256=CASE_PINS[case])
    return result


def collect_job(job,pins,reference):
    path=Path(job['result_directory'])
    expected={(job['axis'],h,b) for h in STEPS for b in BOUNDS}
    errors=[]
    receipts={}
    def optional(name):
        p=path/name
        if not p.exists():
            return None
        try:
            value,digest=read_json(p)
            receipts[name]=dict(sha256=digest,record=value)
            require(isinstance(value,list) and all(isinstance(r,dict) for r in value) if name=='all_derivatives.json'
                    else isinstance(value,dict),'Unexpected receipt JSON schema')
            return value
        except Exception as error:
            errors.append(dict(file=name,error=str(error)))
            return None
    protocol=optional('protocol.json')
    summary=optional('summary.json')
    phase_failure=optional('failure.json')
    if phase_failure is not None:
        errors.append(dict(reason='phase_failure_record_present'))
    mass_records={bound:optional('real_mass_%02d.json' % bound) for bound in BOUNDS}
    optional('real_mass_refinement.json')
    for p in sorted(path.glob('real_mass_refusal_*.json')):
        optional(p.name)
    if protocol is None:
        errors.append(dict(reason='missing_protocol'))
    else:
        try:
            require(protocol['requested_axes']==[job['axis']] and tuple(protocol['all_required_axes'])==AXES,
                    'Protocol axis family differs')
            require(tuple(protocol['complex_steps'])==STEPS and tuple(protocol['bounds'])==BOUNDS,'Protocol step/domain family differs')
            for key in ('case','draw','fixed_gene','fixed_gene_rows','fixed_p','fixed_v','count_vector_sha256','environment'):
                require(protocol[key]==reference['protocol'][key],'Complex/reference point differs: '+key)
            require({Path(k).name:v for k,v in protocol['own_source_hashes'].items()}==
                    {k:pins['source_sha256'][k] for k in ('repair_p3_complex_check.py','repair_p3_complex_reference.py')},
                    'Executed helper/reference source declarations differ')
        except Exception as error:
            errors.append(dict(reason='protocol_refusal',error=str(error)))
    records=[]
    for p in sorted(path.glob('derivative_*.json')):
        record=optional(p.name)
        if record is not None:
            records.append(record)
    aggregate=optional('all_derivatives.json')
    identity=lambda r:(r.get('axis'),r.get('step'),r.get('bound'))
    counts=Counter(identity(r) for r in records)
    if set(counts)!=expected or any(n!=1 for n in counts.values()):
        errors.append(dict(reason='missing_duplicate_or_unexpected_derivative',
            missing=[list(k) for k in sorted(expected-set(counts))],unexpected=[list(k) for k in set(counts)-expected],
            duplicate_keys=[list(k) for k,n in counts.items() if n!=1]))
    if aggregate is None:
        errors.append(dict(reason='missing_final_derivative_aggregate'))
    elif Counter(identity(r) for r in aggregate)!=Counter({k:1 for k in expected}) or any(r not in records for r in aggregate):
        errors.append(dict(reason='aggregate_identity_or_record_disagreement',aggregate_retained=True))
    comparisons=[]
    completed=[]
    for record in records:
        if record.get('completed') is not True:
            continue  # Raw refusal remains in records and receipts, never omitted.
        try:
            axis,step,bound=identity(record)
            require((axis,step,bound) in expected,'Unexpected successful record identity')
            result=record['result']
            require(result['axis']==axis and result['step']==step,'Nested result identity differs')
            require(mass_records[bound] is not None and
                    mass_records[bound]['mass']['geometry']['bound']==bound and
                    mass_records[bound]['mass']['value']==result['real_mass'],'Real denominator identity differs')
            g=finite_number(result['gradient'])
            uncertainty=finite_number(result['estimated_combined_gradient_uncertainty'])
            require(uncertainty>=0,'Negative complex uncertainty')
            cg=finite_number(reference['composite_gradient'][axis])
            ce=finite_number(reference['composite_uncertified_gradient_errors'][axis])
            require(ce>=0,'Negative composite uncertainty')
            difference=cg-g
            tolerance=1e-8+1e-6*abs(g)
            originals=[r for r in reference['original_FD_records'] if r['axis']==axis]
            comparisons.append(dict(axis=axis,step=step,bound=bound,composite_gradient=cg,complex_gradient=g,
                composite_minus_complex=difference,binding_tolerance=tolerance,
                nominal_tolerance_within_descriptive_only=abs(difference)<=tolerance,
                summed_uncertified_gradient_error=ce+uncertainty,
                uncertainty_overlap_descriptive_only=abs(difference)<=ce+uncertainty,
                complex_error_components={k:result[k] for k in ('estimated_numerator_denominator_gradient_error',
                    'relative_derivative_tail_allowance','finite_domain_variance_endpoint_allowance')},
                all_original_FD_comparisons=[dict(original_record=r,complex_minus_scalar_FD=g-r['scalar_reference_gradient'],
                    original_FD_truncation_and_arithmetic_unbounded=True) for r in originals],
                accuracy_not_accepted=True))
            completed.append(record)
        except Exception as error:
            errors.append(dict(reason='successful_record_invalid',identity=list(identity(record)),error=str(error)))
    pairs=[]
    for i,left in enumerate(completed):
        for right in completed[i+1:]:
            if left['bound']!=right['bound'] and left['step']!=right['step']:
                continue
            a,b=left['result'],right['result']
            difference=b['gradient']-a['gradient']
            uncertainty=a['estimated_combined_gradient_uncertainty']+b['estimated_combined_gradient_uncertainty']
            quad_only=a['estimated_numerator_denominator_gradient_error']+b['estimated_numerator_denominator_gradient_error']
            pairs.append(dict(first_identity=list(identity(left)),second_identity=list(identity(right)),
                gradient_difference=difference,summed_uncertified_gradient_error=uncertainty,
                summed_quad_only_estimated_error=quad_only,
                difference_outside_uncertified_estimated_error=abs(difference)>uncertainty,
                difference_outside_quad_only_estimated_error=abs(difference)>quad_only,
                estimated_errors_not_certified=True,agreement_alone_not_accuracy_acceptance=True))
    if summary is None:
        errors.append(dict(reason='missing_final_summary'))
    else:
        try:
            require(summary['complex_reference_records_completed'] and summary['case']==job['case'], 'Final summary identity differs')
            require(summary['requested_axes']==[job['axis']] and tuple(summary['complex_steps'])==STEPS and
                    tuple(summary['bounds'])==BOUNDS,'Final summary requested family differs')
            require(summary['requested_record_count']==6 and summary['saved_record_count']==len(records)==6,
                    'Final summary record count differs')
            require(summary['successful_record_count']==sum(r.get('completed') is True for r in records),
                    'Final summary successful count differs')
            require(summary['derivative_accuracy_accepted'] is False,'Checker falsely claims accepted accuracy')
        except Exception as error:
            errors.append(dict(reason='summary_refusal',error=str(error)))
    numerical_complete=not errors and len(completed)==6 and summary is not None and summary.get('requested_family_numerically_completed') is True
    return dict(**job,receipts=receipts,raw_derivative_records=records,structural_errors=errors,
        phase_failure_retained=phase_failure is not None,refusal_records=[r for r in records if r.get('completed') is not True],
        family_metadata_complete=not errors,successful_derivative_count=len(completed),
        all_requested_complex_derivatives_completed=numerical_complete,
        composite_comparisons=comparisons,step_and_domain_comparisons=pairs,
        all_original_FD_records_retained_at_case_level=True,scheduler_state_not_checked=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect',action='store_true')
    parser.add_argument('--output-directory',type=Path)
    args=parser.parse_args()
    pins,jobs=metadata()
    if args.inspect:
        print(json.dumps(dict(submission_mapping_complete=True,individual_jobs=len(jobs),cases=list(CASE_DIRS),
            axes=AXES,steps=STEPS,bounds=BOUNDS,json_only=True,numerical_imports=False,results_collected=False)))
        return
    refs=references()
    results=[collect_job(job,pins,refs[job['case']]) for job in jobs]
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    output=args.output_directory or REC/('repair_p3_complex_collection_'+stamp)
    require(output.parent.resolve()==REC.resolve(),'Collection output must be a fresh owned REC directory')
    output.mkdir(mode=0o750)
    collection=dict(collected_utc=stamp,collector_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        submitted_jobs=list(SUBMISSIONS),source_pins=pins,case_references=refs,jobs=results,
        family_metadata_complete=all(r['family_metadata_complete'] for r in results),
        all_requested_complex_derivatives_completed=all(r['all_requested_complex_derivatives_completed'] for r in results),
        successful_derivative_count=sum(r['successful_derivative_count'] for r in results),expected_derivative_count=84,
        original_FD_record_count=sum(len(r['original_FD_records']) for r in refs.values()),
        phase_failures_retained=sum(r['phase_failure_retained'] for r in results),
        derivative_refusals_retained=sum(len(r['refusal_records']) for r in results),
        nominal_tolerance_outside_count=sum(not c['nominal_tolerance_within_descriptive_only'] for r in results for c in r['composite_comparisons']),
        uncertainty_overlap_outside_count=sum(not c['uncertainty_overlap_descriptive_only'] for r in results for c in r['composite_comparisons']),
        step_domain_outside_estimated_error_count=sum(c['difference_outside_uncertified_estimated_error'] for r in results for c in r['step_and_domain_comparisons']),
        accuracy_accepted=False,model_fit=False,calibration=False,adoption=False,scheduler_state_not_checked=True,
        estimates_not_certified=True,nominal_agreement_does_not_bound_arithmetic_or_complex_truncation=True,
        output_cap_bytes=CAP)
    body=json.dumps(collection,indent=2,allow_nan=False)+'\n'
    require(len(body.encode())<=CAP,'Collection JSON output cap exceeded')
    partial=output/'collection.json.partial'
    with partial.open('x') as handle:
        handle.write(body)
    partial.rename(output/'collection.json')
    print(json.dumps(dict(output=str(output),family_metadata_complete=collection['family_metadata_complete'],
        successful_derivatives=collection['successful_derivative_count'],expected_derivatives=84,
        accuracy_accepted=False,refusals_preserved=True)))


if __name__=='__main__':
    main()
