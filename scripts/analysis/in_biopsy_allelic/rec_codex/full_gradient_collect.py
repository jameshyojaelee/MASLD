#!/usr/bin/env python3
"""JSON-only complete/partial full-gradient collection; no accuracy acceptance."""
import argparse
from collections import Counter
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import sys

REC=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z')
LAUNCH=REC/'repair_p3_full_gradient_launch_20261005T155545Z'
BASE=REC/'repair_p3_full_gradient_baseline_22041766'
HELPER_SHA='24e77619aef4da42e4c7d9e280de1a41e2a48764bf53b3b25e6dd1b99ce08f48'
BASE_SHA='855670c22753501710d4dc2d02da14841116e03132155ef0606046e5de2f1624'
STEPS=(1e-10,1e-15,1e-20)
BOUNDS=(24.,26.)
READ_CAP=64*1024**2
OUTPUT_CAP=32*1024**2
SUBMISSION_SHAS=(
 'e8a3d181312dfa3242e8721590631ab75bd72914427feadbc2593cef8db9d31c',
 '0dfba6dda66925e3c132f155f5ff9c898a2758a08885c89d1441de4936a6a123',
 '1951c6f5e3d05ecee19fa1fae8b17229c81e8f1aba4d55e8b4638adca25d353a',
 'caf2a06408707992a0e78c5fc4171ffacdabcea032792adc6905bc975110e655',
 '238502ae824548911ceb36152e109be382c02b8982cda11ea1c5e9706e2f0368',
 'd2d065b6c70732e6d996414faf12fd6477d849df832aaf7bf557596d849475c3',
 'fcd6e7c4091dfe3d9a502aa7bb763b37117b82eef71689f2cffb725b9724e1e8',
 '3106350d928cab66d76b357d858393fcdecee693676c60f273061f9b78afa129',
 '9e3c0f9a1947cf198fd674ab8db49c29d3e3a3d523dbf704dd12191a3ebfefa6',
 '7197348dc76dbb9c0a85d49df016fa3b1901852ec240e1f212c265d2ad8bea4f',
 'a81898618b61e1af4bfab54e3bc155ff73efecf1875374cefa5c1a029da6e1c0',
 '48819cca9b58792aa52ad19f7e6d015b7a0f081a7a75d24705a52e79f79620b4',
 '7de7fbd0960ae16a56e060dfcf3a9d8d79bd29fb3d8c89e9ad6bdafc762e1b22',
 '784fc78b6d90f78034afa86ec10b71eff10b85a6bf8c45989e105eb1b2a106b1',
 '00b35c91e45cd364099041a6e2480dd52b30776bbe2cda61b3769b45043b81e1',
 '872bbd77746967247fb30b56c25c8a794623d50b03c2f6b2be8396fb3b4d8d52',
 '0df78d6888e52a50e8cd96040cb76fd68ea0f99c61d97a41b26916ca21fe5599',
 '12d5be382ca81c6fd2ecda66de6ce3b4e6bf97dfd4bce6435b23795a1c5b69c3',
 '9ae704165bd2a3de31d7fd8a2344e28473d04d33b3cac7d91a4b1905e79856b7',
 'ab3c5e894bc3b640b223bd8fe579f3c485f927aa9c7864cb690e186f31454fa7',
 '70cd8327e392ad8feaa69e93fae318e7ba5b0c1c9d22742b14cf4e4ed875e22b',
 'c95a6cb86d8a70a8f19747ee162febcefab8e573306cbb878370e4105c9fdaf1',
 'a13e256946ebf4fb0bd2aa1582c60612e47a63c175d45110af447b84f3ca2db4',
 '8fc406d8bdfd772bc6ba96d928a1b8801cc33e46db64293a98c10b2ed7fc9c3c',
 'b19939dc984a62a241493fbe0425a4364ca444f75d9f33bbdcf33d13a41a0db0',
 '51dbd15ba6b9bd2cc15aa06e08ec718c7c55f38b49195a434431bfc029578a52',
 'cd1f7281d058d0d4697f476d7e25081664e6645bff42e995ee53805a8cafe16c',
 'cfe3d468e948354b897a24886aae0e631a451f89a81f2295f3efacf8e7eb10de',
 '2565804d904abe6188dcfb00b1e48cdf09c5488c8331d13511da9f8989a0347d',
)


def require(ok,message):
    if not ok:raise ValueError(message)


def finite(x):return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)


def close_arithmetic(a,b,scale):
    return finite(a) and finite(b) and abs(a-b)<=32*sys.float_info.epsilon*max(1.,scale)


class Reader:
    def __init__(self):self.bytes_read=0

    def read(self,path,expected=None):
        require(path.resolve().is_relative_to(REC.resolve()) and not path.is_symlink(),'Input escaped synthetic REC')
        size=path.stat().st_size
        require(size<=4*1024**2 and self.bytes_read+size<=READ_CAP,'Explicit JSON input cap reached')
        body=path.read_bytes();self.bytes_read+=len(body)
        digest=hashlib.sha256(body).hexdigest()
        require(len(body)==size and self.bytes_read<=READ_CAP,'Receipt grew during collection')
        require(expected is None or digest==expected,'Pinned receipt differs: '+str(path))
        def invalid(value):raise ValueError('Nonfinite JSON constant '+value)
        return json.loads(body,parse_constant=invalid),dict(sha256=digest,bytes=len(body))


def selected_genes(axis):return [axis%8] if axis<16 else list(range(8))


def key(row,pergene=False):
    return tuple(row.get(k) for k in (('gene','axis','step','bound') if pergene else ('axis','step','bound')))


def metadata(reader):
    baseline,receipt=reader.read(BASE/'shared_baseline.json',BASE_SHA)
    protocol,proto_receipt=reader.read(BASE/'protocol.json','4dea1a3ce1e85a9ecd866921748192337a04b191958325f69e5085ab5acfd4b7')
    summary,summary_receipt=reader.read(BASE/'summary.json','bec335139565733e8ca85e66c91e072d8aed6a3df9626843458cd0c761cd0808')
    require(baseline['helper_sha256']==HELPER_SHA and baseline['reference_mass_complete'] and baseline['candidate_complete'] and
        baseline['candidate_gene_count']==8 and not baseline['mass_refusals'] and not baseline['candidate_refusals'],
        'Completed full scalar/candidate baseline required')
    require(summary['mode']=='baseline' and summary['numerical_records_completed'] and summary['accuracy_accepted'] is False,
        'Baseline summary scope differs')
    require(all(protocol[k]==v for k,v in baseline['identity'].items()) and protocol['source_pins']==baseline['source_pins'],
        'Baseline fixedpoint/source identity differs')
    submissions=[]
    for axis,digest in enumerate(SUBMISSION_SHAS):
        job=str(22041901+axis)
        submission,_=reader.read(LAUNCH/f'submission_axis{axis:02d}_{job}.json',digest)
        require(submission['axis']==axis and submission['job_id']==job and
            Path(submission['script'])==LAUNCH/f'axis{axis:02d}.sbatch' and submission['helper_sha256']==HELPER_SHA and
            submission['baseline_sha256']==BASE_SHA and submission['cpus']==1 and submission['memory']=='4G' and
            submission['time']=='90:00:00' and submission['no_fit'] and submission['no_calibration'],
            'Exact submission/axis/source/resource family differs')
        submissions.append(submission)
    return baseline,dict(protocol=protocol,summary=summary,shared_baseline_receipt=receipt,
        protocol_receipt=proto_receipt,summary_receipt=summary_receipt),submissions


def collect_axis(reader,submission,baseline):
    axis=submission['axis'];directory=REC/f"repair_p3_full_gradient_axis{axis:02d}_{submission['job_id']}"
    hashes={};missing=[];errors=[]
    def optional(path):
        if not path.exists():missing.append(str(path));return None
        try:
            value,digest=reader.read(path);hashes[str(path)]=digest
            require(isinstance(value,dict),'Object receipt required')
            return value
        except Exception as error:errors.append(dict(path=str(path),error=str(error)));return None
    protocol=optional(directory/'protocol.json');summary=optional(directory/'summary.json')
    final=optional(directory/'axis_results.json')
    failure_receipts=[]
    # Preserve only existing small JSON failure/accounting receipts, never infer
    # scheduler failure from absent output or call the scheduler here.
    paths={p for p in directory.glob('*.json') if any(s in p.name.lower() for s in ('failure','scheduler','sacct'))}
    paths.update(p for p in LAUNCH.glob('*'+submission['job_id']+'*.json')
        if any(s in p.name.lower() for s in ('failure','scheduler','sacct')))
    for path in sorted(paths):
        value=optional(path)
        if value is not None:failure_receipts.append(dict(path=str(path),receipt=value))
    if protocol is not None:
        try:
            require(all(protocol[k]==v for k,v in baseline['identity'].items()) and
                protocol['source_pins']==baseline['source_pins'] and protocol['helper_sha256']==HELPER_SHA and
                protocol['mode']=='axis' and protocol['axis']==axis and protocol['steps']==list(STEPS) and
                protocol['domains']==[[24.,.5],[26.,.25]] and protocol['per_gene_axis_cap']==100000 and
                protocol['required_combined_axes']==list(range(29)) and protocol['accuracy_accepted'] is False,
                'Axis fixedpoint/source/step/cap declaration differs')
        except Exception as error:errors.append(dict(reason='protocol_invalid',error=str(error)))
    genes=selected_genes(axis)
    expected={(g,axis,h,b) for g in genes for h in STEPS for b in BOUNDS}
    standalone=[]
    for gene in genes:
        for step in STEPS:
            for bound in BOUNDS:
                path=directory/f'derivative_gene{gene:02d}_{step}_B{bound:02.0f}.json'
                row=optional(path)
                if row is not None:
                    standalone.append(row)
                    if key(row,True)!=(gene,axis,step,bound):errors.append(dict(path=str(path),reason='Derivative identity differs'))
    standalone_counts=Counter(key(r,True) for r in standalone)
    if any(n!=1 for n in standalone_counts.values()) or set(standalone_counts)-expected:
        errors.append(dict(reason='Duplicate/unexpected per-gene records'))
    sums=[];pairs=[];counter_receipts=[];reference_refusals=[r for r in standalone if r.get('completed') is not True]
    rederived=[];physical_new=0
    if final is not None:
        try:
            require(final['axis']==axis and final['baseline_sha256']==BASE_SHA and final['per_gene_axis_cap']==100000 and
                final['aggregate_axis_cap']==len(genes)*100000 and final['accuracy_accepted'] is False,
                'Final axis source/cap differs')
            require(Counter(key(r,True) for r in final['records'])==Counter(expected) and
                Counter(key(r,True) for r in standalone)==Counter(expected),'Final/standalone per-gene family incomplete')
            require(sorted(final['records'],key=lambda r:key(r,True))==sorted(standalone,key=lambda r:key(r,True)),
                'Final records differ from standalone receipts')
            sums=final['full_sum_results'];pairs=final['all_step_domain_comparisons'];counter_receipts=final['probability_counters']
            require(Counter(key(r) for r in sums)==Counter((axis,h,b) for h in STEPS for b in BOUNDS),
                'Six full-sum step/domain family differs')
            for row in sums:
                require(row['required_genes']==genes,'Full-sum contributing genes differ')
                inputs=[r for r in standalone if r['step']==row['step'] and r['bound']==row['bound']]
                successful=all(r.get('completed') is True for r in inputs)
                require(row['completed'] is successful,'Full-sum completion differs from reference inputs')
                if not successful:continue
                gradient=sum(r['result']['gradient'] for r in inputs)
                uncertainty=sum(r['result']['estimated_combined_gradient_uncertainty'] for r in inputs)
                candidate=baseline['candidate_total']['gradient'][axis]
                difference=gradient-candidate;limit=1e-8+1e-6*abs(gradient)
                combined=uncertainty+baseline['candidate_total']['gradient_uncertified_error'][axis]
                require(row['gradient']==gradient and row['summed_uncertified_error']==uncertainty and
                    row['candidate_gradient']==candidate and row['independent_minus_candidate']==difference and
                    row['formal_gradient_tolerance']==limit and row['combined_uncertified_error']==combined and
                    row['nominal_tolerance_within'] is (abs(difference)<=limit) and
                    row['uncertainty_overlap'] is (abs(difference)<=combined),'Full-sum saved arithmetic differs')
                if axis<16:require(row['other_gene_contributions_exact_zero']==[g for g in range(8) if g not in genes] and
                    row['other_genes_zero_by_disjoint_gene_parameters'],'Disjoint-gene zero mapping differs')
                rederived.append(dict(axis=axis,step=row['step'],bound=row['bound'],gradient=gradient,
                    difference=difference,formal_tolerance=limit,nominal_tolerance_within=abs(difference)<=limit,
                    summed_uncertified_error=combined,uncertainty_overlap=abs(difference)<=combined,accuracy_accepted=False))
            completed=[r for r in sums if r['completed']]
            expected_pairs={(left['step'],left['bound'],right['step'],right['bound'])
                for i,left in enumerate(completed) for right in completed[i+1:]}
            actual_pairs=Counter((p['first']['step'],p['first']['bound'],p['second']['step'],p['second']['bound']) for p in pairs)
            require(actual_pairs==Counter(expected_pairs),'All saved pair refinements not retained')
            for pair in pairs:
                first=next(r for r in completed if (r['step'],r['bound'])==(pair['first']['step'],pair['first']['bound']))
                second=next(r for r in completed if (r['step'],r['bound'])==(pair['second']['step'],pair['second']['bound']))
                difference=second['gradient']-first['gradient'];error=first['summed_uncertified_error']+second['summed_uncertified_error']
                require(pair['difference']==difference and pair['summed_uncertified_error']==error and
                    pair['disagreement_outside_estimate'] is (abs(difference)>error),'Pair refinement arithmetic differs')
            require(Counter(c['gene'] for c in counter_receipts)==Counter(genes),'Per-gene counter family differs')
            for counter in counter_receipts:
                shared=baseline['masses'][str(counter['gene'])]['shared_geometry_mass_counter']
                require(counter['shared_counter']==shared and counter['successful_new_derivative_calls']==
                    max(0,min(counter['total_counter_including_refused_attempt'],100000)-shared),'Evaluation scope accounting differs')
                require(0<=shared<=100000,'Shared counter cap exceeded')
                physical_new+=counter['successful_new_derivative_calls']
            require(final['independent_total_scalar_masses_reused']==baseline['independent_total_scalar_masses'] and
                final['likelihood_comparisons_reused']==baseline['likelihood_comparisons'], 'Full scalar baseline reuse differs')
            require(final['refusals']==reference_refusals,'Final reference refusals differ')
        except Exception as error:errors.append(dict(reason='final_axis_invalid',error=str(error)))
    if summary is not None:
        try:
            require(summary['mode']=='axis' and summary['axis']==axis and summary['required_combined_axes']==list(range(29)) and
                summary['accuracy_accepted'] is False and summary['all29_accuracy_accepted'] is False and
                summary['no_fits'] and summary['no_calibration'] and summary['source_hashes_checked_before_after'],
                'Summary source/scientific scope differs')
            require(final is not None and summary['numerical_records_completed'] is
                (not reference_refusals and len(sums)==6 and all(r['completed'] for r in sums)),'Summary completion differs')
        except Exception as error:errors.append(dict(reason='summary_invalid',error=str(error)))
    return dict(submission=submission,axis=axis,directory=str(directory),protocol=protocol,summary=summary,
        raw_axis_results=final,standalone_per_gene_records=standalone,full_sum_records=sums,
        all_step_domain_pairs=pairs,rederived_comparisons=rederived,source_receipt_hashes=hashes,
        missing_receipts=missing,missing_per_gene_identities=[list(k) for k in sorted(expected-set(standalone_counts))],
        structural_or_arithmetic_errors=errors,reference_refusals=reference_refusals,
        failure_and_scheduler_receipts=failure_receipts,probability_counter_receipts=counter_receipts,
        physical_new_derivative_calls_rederived=physical_new,
        full_record_family_present=len(standalone)==len(expected) and len(sums)==6 and summary is not None and not errors,
        accuracy_accepted=False,scheduler_state_not_checked=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--inspect',action='store_true')
    args=parser.parse_args();reader=Reader();baseline,base_metadata,submissions=metadata(reader)
    if args.inspect:
        print(json.dumps(dict(submitted_axes=29,expected_full_sum_records=174,expected_per_gene_records=720,
            expected_pairs=435,helper_sha256=HELPER_SHA,baseline_sha256=BASE_SHA,
            input_cap=READ_CAP,output_cap=OUTPUT_CAP,results_collected=False,numerical_imports=False,accuracy_accepted=False)))
        return
    candidate_receipts=[];base_errors=[]
    for gene in range(8):
        value,receipt=reader.read(BASE/f'candidate_gene{gene:02d}.json')
        require(value['gene']==gene and len(value['gradient'])==29,'Baseline candidate gene shape differs')
        candidate_receipts.append(dict(record=value,receipt=receipt))
    for axis in range(29):
        values=[c['record']['gradient'][axis] for c in candidate_receipts]
        actual=baseline['candidate_total']['gradient'][axis]
        if not close_arithmetic(sum(values),actual,sum(abs(v) for v in values)):
            base_errors.append(dict(axis=axis,reason='Candidate full score sum differs',saved=actual,rederived=sum(values)))
    require(sum(c['record']['logmass'] for c in candidate_receipts)==baseline['candidate_total']['logmass'],
        'Candidate full likelihood sum differs')
    axes=[collect_axis(reader,s,baseline) for s in submissions]
    sum_records=[r for a in axes for r in a['full_sum_records']]
    comparisons=[r for a in axes for r in a['rederived_comparisons']]
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    collection=dict(collected_utc=stamp,collector_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        baseline=baseline,baseline_metadata=base_metadata,baseline_candidate_receipts=candidate_receipts,
        baseline_arithmetic_errors=base_errors,axes=axes,expected_axes=29,expected_full_sum_records=174,
        expected_per_gene_records=720,expected_step_domain_pairs=435,
        saved_per_gene_records=sum(len(a['standalone_per_gene_records']) for a in axes),
        saved_full_sum_records=len(sum_records),successful_full_sum_records=sum(r.get('completed') is True for r in sum_records),
        nominal_tolerance_outside=sum(r['nominal_tolerance_within'] is False for r in comparisons),
        uncertainty_overlap_outside=sum(r['uncertainty_overlap'] is False for r in comparisons),
        saved_step_domain_pairs=sum(len(a['all_step_domain_pairs']) for a in axes),
        pairs_outside_uncertified_estimates=sum(p['disagreement_outside_estimate'] for a in axes for p in a['all_step_domain_pairs']),
        reference_refusals_retained=sum(len(a['reference_refusals']) for a in axes),
        failure_and_scheduler_receipts_retained=sum(len(a['failure_and_scheduler_receipts']) for a in axes),
        full_record_family_present=all(a['full_record_family_present'] for a in axes) and not base_errors,
        shared_baseline_probability_calls_counted_once=baseline['shared_probability_evaluations'],
        new_derivative_probability_calls_from_final_receipts=sum(a['physical_new_derivative_calls_rederived'] for a in axes),
        counters_for_axes_without_final_receipts_not_invented=True,
        fixed_point_scope='Saved positive18 full173-row point; no kept-fit or convergence validation.',
        all_errors_uncertified=True,complex_step_truncation_and_arithmetic_unbounded=True,
        accuracy_accepted=False,full_fit_accuracy_verified=False,calibration=False,adoption=False,
        scheduler_state_not_checked=True,json_bytes_read=reader.bytes_read,input_cap=READ_CAP,output_cap=OUTPUT_CAP)
    body=json.dumps(collection,indent=2,allow_nan=False)+'\n'
    require(len(body.encode())<=OUTPUT_CAP,'Collection output cap exceeded')
    output=REC/('repair_p3_full_gradient_collection_'+stamp);output.mkdir(mode=0o750)
    partial=output/'collection.json.partial'
    with partial.open('x') as handle:handle.write(body)
    partial.rename(output/'collection.json')
    print(json.dumps(dict(output=str(output),saved_sum_records=len(sum_records),expected_sum_records=174,
        saved_per_gene_records=collection['saved_per_gene_records'],expected_per_gene_records=720,
        reference_refusals=collection['reference_refusals_retained'],accuracy_accepted=False)))


if __name__=='__main__':main()
