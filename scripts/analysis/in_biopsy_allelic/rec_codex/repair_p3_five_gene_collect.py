#!/usr/bin/env python3
"""Bounded JSON collection of all fifteen prescribed synthetic mass checks.

Retain incomplete points, numerical refusals and frozen summary flags separately.
No array reads, numerical imports, scheduler access, fits or accuracy acceptance.
"""
import argparse
from collections import Counter
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path

REC=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z')
LAUNCH=REC/'repair_p3_five_gene_accuracy_launch_20261005T143525Z'
JOBS=tuple(range(22039932,22039945))+(22039946,22039947)
COUNTS=(5,25,75,200,500)
MULTIPLES=(.25,1.,4.)
READ_CAP=64*1024**2
OUTPUT_CAP=32*1024**2
SOURCES={
 'repair_p3_five_gene_accuracy.py':'4aec076799a157551d86c1c2fc3da8901e436f8779cca190643c3398fd0f0b21',
 'repair_p3_integral.py':'461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3',
 'repair_p3_composite_gl.py':'45b4423811fded58712bd5e85666f9ffc2b6972abbd9bd71964a584ec6d3c59c',
 'check_quadrature_v1.py':'40a51fd0090d797b631c897f1c35c046b750e040a780432f0f614960b9a17735',
 'likelihood_v1.py':'41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6',
}
SUBMISSION_SHA='7f9d9b51de8d3b4a44acbb3f7cc5776d7d8ef1e37115fad64c1476454d690694'
SOURCE_RECEIPT_SHA='56efa9ae67f85c921409f04db63b7844b1e066f4b2e635031da9094448da226a'


def require(ok,message):
    if not ok:raise ValueError(message)


def finite(value):
    return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)


class Reader:
    def __init__(self):self.bytes_read=0

    def bytes(self,path,expected=None,cap=4*1024**2):
        require(path.resolve().is_relative_to(REC.resolve()) and not path.is_symlink(),'Input escaped synthetic REC')
        size=path.stat().st_size
        require(size<=cap and self.bytes_read+size<=READ_CAP,'Explicit input byte cap reached')
        body=path.read_bytes();self.bytes_read+=len(body)
        require(len(body)==size and self.bytes_read<=READ_CAP,'Input grew while collecting')
        digest=hashlib.sha256(body).hexdigest()
        require(expected is None or expected==digest,'Pinned receipt/source differs: '+str(path))
        return body,digest

    def json(self,path,expected=None):
        require(path.suffix=='.json','Only JSON result receipts allowed')
        body,digest=self.bytes(path,expected)
        def invalid(value):raise ValueError('Nonfinite JSON constant '+value)
        return json.loads(body,parse_constant=invalid),dict(sha256=digest,bytes=len(body))


def metadata(reader):
    sources,_=reader.json(LAUNCH/'source_hashes.json',SOURCE_RECEIPT_SHA)
    require(sources==SOURCES,'Frozen source declaration differs')
    for name,digest in SOURCES.items():reader.bytes(LAUNCH/name,digest)
    body,_=reader.bytes(LAUNCH/'submissions.jsonl',SUBMISSION_SHA)
    submissions=[json.loads(line) for line in body.decode().splitlines()]
    require([int(s['job_id']) for s in submissions]==list(JOBS),'Submitted fifteen-job roster differs')
    points=[(g,m) for g in range(5) for m in MULTIPLES]
    for s,(g,m) in zip(submissions,points):
        require(Path(s['script'])==LAUNCH/f'gene{g}_tau{m}.sbatch','Submitted Cartesian point differs')
        require(s['cpu']==1 and s['memory']=='4G' and s['walltime']=='90:00:00','Resource receipt differs')
    return submissions,points


def point_checks(record):
    """Recompute only arithmetic from saved scalar likelihood receipts."""
    checks={};errors=[];refusals=[]
    for name in ('candidate','original_GH20','original_GH40'):
        entry=record.get(name,{})
        if entry.get('refused') is True:refusals.append(dict(method=name,raw_refusal=entry))
    if 'scalar_refusal' in record:refusals.append(dict(method='independent_scalar',raw_refusal=record['scalar_refusal']))
    scalar=record.get('scalar',{})
    first,expanded=scalar.get('first',{}),scalar.get('expanded',{})
    scalar_available=(finite(first.get('logmass')) and finite(expanded.get('logmass')) and
        finite(expanded.get('estimated_logmass_error')) and finite(first.get('relative_prior_mass_tail_allowance')))
    scalar_ok=False
    if scalar_available:
        refinement=abs(expanded['logmass']-first['logmass'])+expanded['estimated_logmass_error']
        scalar_ok=('scalar_refusal' not in record and refinement<=1e-9 and
            first['relative_prior_mass_tail_allowance']<=1e-10)
        checks['scalar']=dict(saved_first=first,saved_expanded=expanded,
            rederived_refinement_plus_estimate=refinement,refinement_budget=1e-9,
            first_tail_budget=1e-10,within_recorded_scalar_requirements=scalar_ok,
            uncertainty_not_certified=True)
    else:checks['scalar']=dict(available=False,refused='scalar_refusal' in record)
    for name in ('candidate','original_GH20','original_GH40'):
        entry=record.get(name,{})
        result=dict(available=finite(entry.get('logmass')),refused=entry.get('refused') is True)
        if result['available'] and scalar_available:
            difference=entry['logmass']-expanded['logmass'];limit=1e-8*abs(expanded['logmass'])
            saved_delta=record.get('candidate_scalar_difference') if name=='candidate' else entry.get('scalar_difference')
            saved_flag=record.get('candidate_scalar_relative_LL_pass') if name=='candidate' else entry.get('scalar_relative_LL_pass')
            raw_pass=abs(difference)<=limit
            result.update(saved_difference=saved_delta,saved_point_flag=saved_flag,rederived_difference=difference,
                formal_relative_loglikelihood_limit=limit,rederived_raw_point_flag=raw_pass,
                qualified_point_comparison=scalar_ok and raw_pass,
                scalar_plus_candidate_error_estimate_not_certified=True)
            if scalar_ok and (saved_delta!=difference or saved_flag is not raw_pass):
                errors.append(dict(method=name,reason='saved_point_arithmetic_or_flag_differs'))
        checks[name]=result
    a,b=record.get('original_GH20',{}),record.get('original_GH40',{})
    if finite(a.get('logmass')) and finite(b.get('logmass')):
        difference=a['logmass']-b['logmass'];passed=abs(difference)<=1e-8
        checks['GH20_GH40']=dict(saved_difference=record.get('original_GH20_GH40_difference'),
            saved_point_flag=record.get('original_GH20_GH40_absolute_pass'),rederived_difference=difference,
            absolute_limit=1e-8,rederived_raw_point_flag=passed)
        if record.get('original_GH20_GH40_difference')!=difference or record.get('original_GH20_GH40_absolute_pass') is not passed:
            errors.append(dict(method='GH20_GH40',reason='saved_point_arithmetic_or_flag_differs'))
    else:checks['GH20_GH40']=dict(available=False)
    return checks,errors,refusals


def collect_job(reader,submission,point):
    gene,multiple=point;directory=REC/('repair_p3_five_gene_accuracy_'+submission['job_id'])
    files={};missing=[];errors=[]
    def optional(name):
        path=directory/name
        if not path.exists():missing.append(name);return None
        try:
            value,receipt=reader.json(path);files[name]=receipt
            require(isinstance(value,dict),'JSON object required')
            return value
        except Exception as error:
            errors.append(dict(file=name,error=str(error)));return None
    protocol=optional('protocol.json');fixture=optional('fixture.json');summary=optional('summary.json')
    name=f'gene{gene}_tau{multiple}.json';record=optional(name)
    phase_failure=optional('failure.json') if (directory/'failure.json').exists() else None
    if protocol is not None:
        try:
            require(protocol['seed']==20260923 and protocol['rows_per_gene']==list(COUNTS) and protocol['tau_max']==.05 and
                protocol['tau_multiples']==list(MULTIPLES) and protocol['required_points']==15 and
                protocol['requested_points']==[[multiple,gene]] and protocol['requested_point_count']==1 and
                protocol['combined_family_all15_remains_required'] and protocol['no_model_fits'] and protocol['no_calibration'],
                'Protocol scientific point/family differs')
            require({Path(k).name:v for k,v in protocol['source_sha256'].items()}==SOURCES,'Executed source pins differ')
        except Exception as error:errors.append(dict(reason='protocol_invalid',error=str(error)))
    if fixture is not None:
        try:
            require(fixture['rows_per_gene']==list(COUNTS) and len(fixture['drawn_counts'])==sum(COUNTS) and
                len(fixture['parameters'])==16 and all(finite(v) for v in fixture['parameters']) and
                fixture['parameters'][:5]==[-1.2,-.7,.4,.9,1.3] and fixture['parameters'][10:12]==[.08,.04],
                'Frozen synthetic fixture identity differs')
        except Exception as error:errors.append(dict(reason='fixture_invalid',error=str(error)))
    checks={};refusals=[];summary_disagreements=[]
    if record is not None:
        try:
            require(record['gene']==gene and record['rows']==COUNTS[gene] and record['tau']==multiple*.05 and
                record['variance']==(multiple*.05)**2,'Saved point identity differs')
            checks,arithmetic_errors,refusals=point_checks(record);errors.extend(arithmetic_errors)
        except Exception as error:errors.append(dict(reason='point_record_invalid',error=str(error)))
    if summary is not None:
        try:
            require(summary['requested_point_count']==1 and summary['all_requested_points_recorded'] and
                summary['combined_family_all15_remains_required'] and summary['all15points_recorded'] is False and
                summary['candidate_scalar_relative_LL_pass_all15'] is False and summary['original_GH20_GH40_pass_all15'] is False and
                summary['original_GH20_scalar_pass_all15'] is False and summary['independent_derivative_accuracy_verified'] is False and
                summary['no_model_fits'] and summary['no_calibration'],'Shard summary scope differs')
            raw_flag=record.get('candidate_scalar_relative_LL_pass') if record is not None else None
            saved=summary['candidate_scalar_relative_LL_pass_requested_points']
            if raw_flag is not None and saved is not raw_flag:
                summary_disagreements.append(dict(saved_requested_summary_flag=saved,saved_raw_point_flag=raw_flag,
                    consistent_with_frozen_numpy_bool_identity_bug=saved is False and raw_flag is True,
                    raw_summary_preserved_not_corrected=True))
        except Exception as error:errors.append(dict(reason='summary_invalid',error=str(error)))
    state=('phase_failure' if phase_failure is not None else 'point_missing' if record is None else
        'point_saved_summary_missing' if summary is None else 'point_and_summary_saved')
    return dict(job_id=submission['job_id'],gene=gene,rows=COUNTS[gene],tau_multiple=multiple,directory=str(directory),
        observed_receipt_state=state,raw_protocol=protocol,raw_summary=summary,raw_point_record=record,
        file_hashes=files,missing_receipts=missing,structural_or_arithmetic_errors=errors,phase_failure=phase_failure,
        numerical_refusals=refusals,rederived_point_comparisons=checks,
        requested_summary_disagreements=summary_disagreements,scheduler_state_not_checked=True,
        point_record_present=record is not None,final_summary_present=summary is not None,accuracy_accepted=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect',action='store_true')
    args=parser.parse_args();reader=Reader();submissions,points=metadata(reader)
    if args.inspect:
        print(json.dumps(dict(source_and_submission_pins_verified=True,required_points=15,
            cartesian_points=points,source_pins=SOURCES,input_cap=READ_CAP,output_cap=OUTPUT_CAP,
            numerical_imports=False,results_collected=False,accuracy_accepted=False)))
        return
    results=[collect_job(reader,s,p) for s,p in zip(submissions,points)]
    fixture_hashes={r['file_hashes']['fixture.json']['sha256'] for r in results if 'fixture.json' in r['file_hashes']}
    environments=[r['raw_protocol']['environment'] for r in results if r['raw_protocol'] is not None and 'environment' in r['raw_protocol']]
    cross_shard_errors=[]
    if len(fixture_hashes)>1:cross_shard_errors.append('Different synthetic fixture bytes across shards')
    if environments and any(e!=environments[0] for e in environments):cross_shard_errors.append('Different numerical environments across shards')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    result=dict(collected_utc=stamp,collector_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        submitted_jobs=list(JOBS),required_cartesian_points=points,required_points=15,results=results,
        saved_point_records=sum(r['point_record_present'] for r in results),final_summaries=sum(r['final_summary_present'] for r in results),
        observed_receipt_states=dict(Counter(r['observed_receipt_state'] for r in results)),
        numerical_refusals_retained=sum(len(r['numerical_refusals']) for r in results),
        summary_boolean_disagreements_retained=sum(len(r['requested_summary_disagreements']) for r in results),
        cross_shard_errors=cross_shard_errors,source_pins=SOURCES,shared_fixture_sha256=list(fixture_hashes),
        environment=environments[0] if environments else None,
        frozen_summary_bug='Executed4aec uses NumPy comparison booleans followed by identity with True; retain raw false requested-summary flags alongside serialized point flags.',
        relative_loglikelihood_comparison_not_relative_probability_accuracy=True,
        uncertainty_estimates_not_certified=True,independent_derivative_accuracy_verified=False,
        complete_fit_accuracy_verified=False,accuracy_accepted=False,statistical_calibration=False,
        scheduler_state_not_checked=True,json_and_source_bytes_read=reader.bytes_read,input_cap=READ_CAP,output_cap=OUTPUT_CAP)
    body=json.dumps(result,indent=2,allow_nan=False)+'\n'
    require(len(body.encode())<=OUTPUT_CAP,'Collector output byte cap exceeded')
    output=REC/('repair_p3_five_gene_collection_'+stamp);output.mkdir(mode=0o750)
    partial=output/'collection.json.partial'
    with partial.open('x') as handle:handle.write(body)
    partial.rename(output/'collection.json')
    print(json.dumps(dict(output=str(output),saved_point_records=result['saved_point_records'],required_points=15,
        final_summaries=result['final_summaries'],numerical_refusals=result['numerical_refusals_retained'],accuracy_accepted=False)))


if __name__=='__main__':main()
