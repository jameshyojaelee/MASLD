#!/usr/bin/env python3
"""Bounded JSON-only Hessian-family collection; no accuracy acceptance or jobs."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
LAUNCH=REC/'repair_p3_independent_hessian_launch_20261005T144646Z'
AXES=(2,10,16,17,24,25,28)
REAL_STEPS=(1e-3,3e-4,1e-4,3e-5,1e-5)
COMPLEX_STEPS=(1e-10,1e-15,1e-20)
BOUNDS=(24.,26.)
FD_STEPS=(1e-3,3e-4,1e-4)
READ_CAP=64*1024**2
OUTPUT_CAP=64*1024**2
SUBMISSIONS={
 '22040175':'989e6f537e2c4168705346663a1ce44607d3fac947a17b2ae1881caad6a61fcf',
 '22040176':'4ee22993023916c8c9b6eccf9944825556f60154354bae694a01cc1effde1f2b',
 '22040177':'f096513d5bbd4fbcefb5eb57217f73b363b4a8cda3474dba3cb4e6d6f5192a44',
 '22040178':'bf79e98d3b981c8097ac0c97a51143fea204985d7db6ea2cbf9e0349655edfa0',
 '22040179':'fedcfafeee57ba8772e98a4bc368e204ce6d4d38fcaf01545439be94359d0deb',
 '22040180':'1f56c01f036f4928bee2ef6ac57af3426cd054cbaa1d952d8445918e3581ddf5',
 '22040181':'0bdd0d6ee7b0c326b6241fba1c85c8653f22ee5ee05a92147e112d37c9208a14',
 '22040182':'e104eb8a523c41c071548004b587c3aea63417d0c9030094a2adf3fa4254320f',
 '22040183':'feade0ca3b8c1a720b6ce4e7acc3b44da592b95deeb8554396ad3dcc167132b0',
 '22040184':'f4bffd66d00771c5f3d65d291178c2844a3bb0a89af16d69c3c792aee0827ee2',
 '22040185':'808d8b38aefdc582609a4861d50bc41538f69d00dd644604af54f3f5edad4688',
 '22040186':'076cba2f0aa591bdf6eb06f791cfe0b7733bc4a9764a18fa3404aef938fc2f2f',
 '22040187':'316762b0d2cdbf87111f1bf617481c7812f79451e95b725ec620d367ef209cce',
 '22040188':'e08612c2c7f528bcad9e04317bdc242457854f0035a101b06be7f3d4b3403fe0',
}
REFERENCES={
 'failed0gene2':dict(draw=0,scalar='repair_p3_scalar_prefix_failed0gene2_22039118',
   protocol_sha='b27a85f6c11342406e07e99aacafa70839d09fcd870e0449dac40e59aab97a32',
   fd_sha='2ae8d822a8e456ad6e247ffaf91bcc21f5a90e06b63c66c8d2846a0cda96451a'),
 'positive18':dict(draw=18,scalar='repair_p3_scalar_prefix_positive18_22039119',
   protocol_sha='26f061f8e481a1b9d1e11f2a5dfc0984c82010d3b1dd9c8f270eb1418d701001',
   fd_sha='ebce73a49bbca829ff182e4cb0fc710c18fea794edbe8f8e476a7d145c91ac62'),
}


def require(ok,message):
    if not ok:
        raise ValueError(message)


def number(value):
    require(isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value),'Invalid numeric receipt field')
    return value


class Reader:
    def __init__(self):
        self.bytes_read=0

    def read(self,path,expected=None):
        require(path.suffix=='.json' and path.stat().st_size<=4*1024**2,'Only small JSON receipts allowed')
        body=path.read_bytes()
        self.bytes_read+=len(body)
        require(self.bytes_read<=READ_CAP,'Explicit JSON input byte cap reached')
        digest=hashlib.sha256(body).hexdigest()
        require(expected is None or expected==digest,'Pinned JSON differs: '+str(path))
        def invalid(value):
            raise ValueError('Nonfinite JSON constant: '+value)
        return json.loads(body,parse_constant=invalid),digest


def submission_metadata(reader):
    pins,_=reader.read(LAUNCH/'source_pins.json','18ade218da26f1e18db84226bb09a6b6c117cf793bf16954048f784873312b7a')
    normalized={Path(k).name:v for k,v in pins.items()}
    require(normalized=={'repair_p3_independent_hessian_check.py':'423442e1123593385fffa085e42cad3978b61b53f3e8e0596c2441fe50c7c8c3',
                        'repair_p3_complex_reference.py':'f722eb76a2fb288f4f4e697dcd149d714b89221ec791ba6a43db8cfebb2f7ec3'},
            'Frozen checker/reference declarations differ')
    jobs=[]
    for job,digest in SUBMISSIONS.items():
        receipt,_=reader.read(LAUNCH/('submission_'+job+'.json'),digest)
        match=re.fullmatch(r'(failed0gene2|positive18)_axis(\d{2})\.sbatch',Path(receipt['script']).name)
        require(match is not None and Path(receipt['script']).parent==LAUNCH,'Submission script identity differs')
        case,column=match.group(1),int(match.group(2))
        require(receipt['job_id']==job and receipt['source_pins']==pins and receipt['independent_seven_score_axes'],
                'Submission source/scope differs')
        require(tuple(receipt['real_steps'])==REAL_STEPS and receipt['no_fit'] and receipt['no_calibration'],
                'Submission scientific scope differs')
        require(receipt['cpus']==1 and receipt['memory']=='4G' and receipt['time']=='90:00:00','Resource contract differs')
        jobs.append(dict(job_id=job,case=case,column=column,submission=receipt,submission_sha256=digest,
                         directory=str(REC/f'repair_p3_independent_hessian_{case}_{job}')))
    require(Counter((j['case'],j['column']) for j in jobs)==Counter((c,a) for c in REFERENCES for a in AXES),
            'Submitted full two-case/seven-column family differs')
    return normalized,jobs


def identity(row):
    return tuple(row.get(k) for k in ('column','score_axis','real_step','complex_step','bound'))


def collect_job(reader,job,pins,reference):
    directory=Path(job['directory'])
    errors=[]
    files={}
    def optional(name,kind=dict):
        path=directory/name
        if not path.exists():
            return None
        try:
            value,digest=reader.read(path)
            files[name]=dict(sha256=digest,bytes=path.stat().st_size)
            require(isinstance(value,kind),'Unexpected JSON schema')
            return value
        except Exception as error:
            errors.append(dict(file=name,error=str(error)))
            return None
    protocol=optional('protocol.json')
    summary=optional('summary.json')
    failure=optional('failure.json')
    if failure is not None:
        errors.append(dict(reason='phase_failure_present'))
    expected={(job['column'],a,d,h,b) for a in AXES for d in REAL_STEPS for h in COMPLEX_STEPS for b in BOUNDS}
    records=optional('independent_hessian_before_candidate_comparison.json',list)
    if records is None:
        records=[]
    require(all(isinstance(r,dict) for r in records),'Malformed Hessian record list')
    counts=Counter(identity(r) for r in records)
    missing=expected-set(counts)
    duplicates=[list(k) for k,n in counts.items() if n!=1]
    unexpected=[list(k) for k in set(counts)-expected]
    if missing or duplicates or unexpected:
        errors.append(dict(reason='incomplete_duplicate_or_unexpected_210_family',missing=[list(k) for k in sorted(missing)],
                           duplicates=duplicates,unexpected=unexpected))
    original_fd=optional('original_scalar_FD_records.json',list)
    if original_fd!=reference['original_FD_records']:
        errors.append(dict(reason='original21FD_not_fully_retained'))
    if protocol is None:
        errors.append(dict(reason='missing_protocol'))
    else:
        try:
            require(protocol['requested_columns']==[job['column']] and tuple(protocol['score_axes'])==AXES,
                    'Requested column/score-axis scope differs')
            require(tuple(protocol['real_steps'])==REAL_STEPS and tuple(protocol['complex_steps'])==COMPLEX_STEPS and
                    tuple(protocol['bounds'])==BOUNDS,'Step/domain family differs')
            for key in ('case','draw','fixed_gene','fixed_gene_rows','fixed_p','fixed_v','count_vector_sha256','environment'):
                require(protocol[key]==reference['protocol'][key],'Fixed reference identity differs: '+key)
            require({Path(k).name:v for k,v in protocol['own_source_hashes'].items()}==pins,'Executed source declarations differ')
            require(protocol['per_realpoint_scoreaxis_evaluation_cap']==100000 and
                    protocol['aggregate_successful_evaluation_cap']==protocol['distinct_realpoints']*7*100000,
                    'Evaluation cap declaration differs')
        except Exception as error:
            errors.append(dict(reason='protocol_invalid',error=str(error)))
    point_metadata=[]
    score_metadata=[]
    point_refusals=[]
    physical_calls=0
    for path in sorted(directory.glob('point_*_metadata.json')):
        value=optional(path.name)
        if value is None:
            continue
        try:
            count=number(value['shared_geometry_mass_successful_probability_evaluations'])
            require(0<=count<=100000,'Shared realpoint cap invalid')
            physical_calls+=count
            point_metadata.append(dict(point_id=value['point']['point_id'],vector_sha256=value['point']['vector_sha256'],
                shared_calls=count,shared_counter_including_refusal=value['shared_counter_including_refused_attempt']))
            point_refusals.extend(value['refusals'])
        except Exception as error:
            errors.append(dict(file=path.name,error=str(error)))
    for path in sorted(directory.glob('point_*_score_*.json')):
        value=optional(path.name)
        if value is None:
            continue
        try:
            require(value['score_axis'] in AXES and value['maximum_axis_budget_evaluations']==100000,'Score-axis cap differs')
            count=number(value['derivative_successful_probability_evaluations'])
            shared=number(value['shared_baseline_budget_evaluations'])
            require(0<=count and 0<=shared and count+shared<=100000,'Score-axis successful budget exceeded')
            require(len(value['records'])==6 and Counter((r['complex_step'],r['bound']) for r in value['records'])==
                    Counter((h,b) for h in COMPLEX_STEPS for b in BOUNDS),'Realpoint score-axis step family differs')
            physical_calls+=count
            score_metadata.append(dict(point_id=value['point_id'],score_axis=value['score_axis'],
                derivative_calls=count,shared_scope_calls=shared,budget_counter_including_refusal=value['axis_budget_counter_including_refused_attempt']))
            point_refusals.extend(r for r in value['records'] if r.get('completed') is not True)
        except Exception as error:
            errors.append(dict(file=path.name,error=str(error)))
    successful=[r for r in records if r.get('completed') is True]
    for row in successful:
        try:
            number(row['independent_hessian'])
            require(number(row['abs_stencil_weighted_uncertified_error'])>=0,'Negative stencil uncertainty')
        except Exception as error:
            errors.append(dict(identity=list(identity(row)),error=str(error)))
    saved_comparisons=optional('candidate_comparisons.json',list)
    if saved_comparisons is None:
        saved_comparisons=[]
    if Counter(identity(r) for r in saved_comparisons)!=Counter(identity(r) for r in successful):
        errors.append(dict(reason='candidate_comparison_family_differs_from_successful_reference'))
    refinements=optional('real_complex_step_domain_refinement.json',list)
    if refinements is None:
        refinements=[]
    if Counter((r['column'],r['score_axis']) for r in refinements)!=Counter((job['column'],a) for a in AXES):
        errors.append(dict(reason='missing_step_domain_refinement_axes'))
    if summary is None:
        errors.append(dict(reason='missing_final_summary'))
    else:
        try:
            require(summary['independent_hessian_records_completed'] and summary['case']==job['case'] and
                    summary['requested_columns']==[job['column']] and tuple(summary['score_axes'])==AXES,'Summary identity differs')
            require(tuple(summary['real_steps'])==REAL_STEPS and tuple(summary['complex_steps'])==COMPLEX_STEPS and
                    tuple(summary['bounds'])==BOUNDS,'Summary step/domain family differs')
            require(summary['expected_hessian_record_count']==summary['saved_hessian_record_count']==len(records)==210,
                    'Summary210-record family differs')
            require(summary['successful_hessian_record_count']==len(successful),'Successful-record summary differs')
            require(summary['actual_successful_probability_evaluations']==physical_calls,'Physical evaluation sum differs')
            require(summary['aggregate_successful_probability_evaluation_cap']==protocol['aggregate_successful_evaluation_cap'] and
                    physical_calls<=summary['aggregate_successful_probability_evaluation_cap'],'Aggregate evaluation cap exceeded/differs')
            roster=protocol['point_roster']
            require(Counter(r['point_id'] for r in point_metadata)==Counter(p['point_id'] for p in roster), 'Point metadata roster incomplete')
            require(Counter((r['point_id'],r['score_axis']) for r in score_metadata)==
                    Counter((p['point_id'],a) for p in roster for a in AXES),'Realpoint score-axis roster incomplete')
            require(summary['hessian_accuracy_accepted'] is False and summary['no_model_fits'] and summary['no_statistical_calibration'],
                    'Unexpected scientific acceptance claim')
        except Exception as error:
            errors.append(dict(reason='summary_invalid',error=str(error)))
    numerical_complete=(not errors and not point_refusals and not failure and len(successful)==210 and
        summary is not None and summary.get('requested_columns_numerically_complete') is True)
    return dict(**job,receipt_file_hashes=files,protocol=protocol,summary=summary,phase_failure=failure,
        structural_errors=errors,missing_record_identities=[list(k) for k in sorted(missing)],
        all210_records=records,reference_refusal_records=[r for r in records if r.get('completed') is not True],
        point_gradient_refusals=point_refusals,point_metadata=point_metadata,score_axis_evaluation_metadata=score_metadata,
        physical_probability_evaluations_rederived=physical_calls,
        candidate_comparisons=saved_comparisons,full_step_domain_refinement_comparisons=refinements,
        record_family_metadata_complete=not errors,numerical_family_complete=numerical_complete,
        successful_record_count=len(successful),accuracy_accepted=False,scheduler_state_not_checked=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect',action='store_true')
    parser.add_argument('--output-directory',type=Path)
    args=parser.parse_args()
    reader=Reader()
    pins,jobs=submission_metadata(reader)
    if args.inspect:
        print(json.dumps(dict(json_only=True,numerical_imports=False,results_collected=False,submitted_jobs=len(jobs),
            source_pins=pins,score_axes=AXES,columns=AXES,real_steps=REAL_STEPS,complex_steps=COMPLEX_STEPS,bounds=BOUNDS,
            expected_records_per_column=210,expected_two_case_records=2940,accuracy_accepted=False)))
        return
    refs={}
    for case,spec in REFERENCES.items():
        protocol,_=reader.read(REC/spec['scalar']/'protocol.json',spec['protocol_sha'])
        fd,_=reader.read(REC/spec['scalar']/'finite_differences.json',spec['fd_sha'])
        require(len(fd)==21 and Counter((r['axis'],r['step']) for r in fd)==Counter((a,h) for a in AXES for h in FD_STEPS),
                'Original complete scalar FD family differs')
        require(protocol['case']==case and protocol['draw']==spec['draw'],'Original case differs')
        refs[case]=dict(protocol=protocol,original_FD_records=fd,source_metadata=spec)
    results=[collect_job(reader,job,pins,refs[job['case']]) for job in jobs]
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    output=args.output_directory or REC/('repair_p3_independent_hessian_collection_'+stamp)
    require(output.parent.resolve()==REC.resolve(),'Output must be a fresh owned REC directory')
    output.mkdir(mode=0o750)
    collection=dict(collected_utc=stamp,collector_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_pins=pins,case_references=refs,jobs=results,
        expected_records=2940,saved_records=sum(len(r['all210_records']) for r in results),
        successful_records=sum(r['successful_record_count'] for r in results),
        family_metadata_complete=all(r['record_family_metadata_complete'] for r in results),
        full_two_case7x7_numerically_complete=all(r['numerical_family_complete'] for r in results),
        reference_refusals_retained=sum(len(r['reference_refusal_records']) for r in results),
        point_gradient_refusals_retained=sum(len(r['point_gradient_refusals']) for r in results),
        phase_failures_retained=sum(r['phase_failure'] is not None for r in results),original_FD_records_retained=42,
        source_scope='Only selected7x7; all29moment entries outside this block remain independently unverified.',
        nominal_hessian_budget=1e-5,accuracy_accepted=False,model_fit=False,calibration=False,adoption=False,
        estimated_errors_not_certified=True,real_complex_stencil_truncation_arithmetic_unbounded=True,
        scheduler_state_not_checked=True,json_bytes_read=reader.bytes_read,json_read_cap=READ_CAP,output_cap=OUTPUT_CAP)
    body=json.dumps(collection,indent=2,allow_nan=False)+'\n'
    require(len(body.encode())<=OUTPUT_CAP,'Collection output byte cap exceeded')
    partial=output/'collection.json.partial'
    with partial.open('x') as handle:
        handle.write(body)
    partial.rename(output/'collection.json')
    print(json.dumps(dict(output=str(output),saved_records=collection['saved_records'],expected_records=2940,
        successful_records=collection['successful_records'],family_metadata_complete=collection['family_metadata_complete'],
        numerical_family_complete=collection['full_two_case7x7_numerically_complete'],accuracy_accepted=False)))


if __name__=='__main__':
    main()
