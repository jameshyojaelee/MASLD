#!/usr/bin/env python3
"""JSON-only descriptive supplement; no model evaluation/new acceptance rule."""
import hashlib
import json
import math
from pathlib import Path

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
INPUTS={
 'repair_p3_batched_gl_failed0gene2_22040257/candidate_baseline.json':'b2207a47757ed50d63d59802bc6d81566c23529337d6e5465cee6c35554bf6c5',
 'repair_p3_batched_gl_positive18_22040259/candidate_baseline.json':'8e84dcabf557dd71824936f9a581ed8a1405068330602359f6431bb3038425bf',
 'repair_p3_local_gl_failed0gene2_22041594/protocol.json':'7c7b7b696b4efc92d7194a1e5a64ce0ecd957ba45c8ee80fc559851c159136a2',
 'repair_p3_local_gl_failed0gene2_22041594/summary.json':'af820d623fcb2fe92b19c4ef5010f1dedafbb9e3191f04ec65a9d6bd86c63348',
 'repair_p3_local_gl_failed0gene2_22041594/local_gl_comparison.json':'c587ca41b27a9cf623b2fd7de9eb342d6700e90373fc90d05fd9249f3d50293b',
 'repair_p3_local_gl_positive18_22041595/protocol.json':'76d7d884fdc8307375dcf23e1d7ac3d0e6053a5d9627791d8b919f46b83241af',
 'repair_p3_local_gl_positive18_22041595/summary.json':'af0e24522951f837e59756e9c38b24f459312ff3f0600c485982f7cc3df05ec6',
 'repair_p3_local_gl_positive18_22041595/local_gl_comparison.json':'82215c38468332263d33ac11f2e30217ffe0f8c1e959d799524b7c5a8f03305a'}


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def flat(a):return [z for x in a for z in flat(x)] if isinstance(a,list) else [a]
def maxima(a):return max(flat(a))
def within(errors,settings,fraction):
    return (errors[0]<=settings['logmass_tolerance']*fraction
        and maxima(errors[1])<=settings['gradient_tolerance']*fraction
        and maxima(errors[2])<=settings['hessian_tolerance']*fraction)
def decisions(ladder,settings):
    tail=(ladder['relative_mass_tail_allowance']<=settings['logmass_tolerance']/8
        and maxima(ladder['relative_gradient_tail_allowance'])<=settings['gradient_tolerance']/8
        and maxima(ladder['relative_hessian_tail_allowance'])<=settings['hessian_tolerance']/8)
    errors=ladder.get('order_refinement_32_64_allowances')
    return dict(tail_within_existing_eighth_budget=tail,
        order_within_existing_quarter_budget=(within(errors,settings,.25) if errors is not None else None))
def projection(array,indices,dimension):
    return ([array[j] for j in indices] if dimension==1 else [[array[i][j] for j in indices] for i in indices])
def scales(a,b):
    pairs=list(zip(flat(a),flat(b)))
    differences=[abs(x-y) for x,y in pairs]
    i=max(range(len(pairs)),key=lambda k:differences[k])
    x,y=pairs[i]
    return dict(max_absolute_difference=differences[i],
        old_at_max_absolute_difference=y,new_at_max_absolute_difference=x,
        old_maximum_absolute_magnitude=max(abs(y) for x,y in pairs),
        new_maximum_absolute_magnitude=max(abs(x) for x,y in pairs),
        maximum_scale_normalized_difference=max(abs(x-y)/max(1.,abs(x),abs(y)) for x,y in pairs),
        descriptive_scale='max(1,abs(old),abs(new)); no acceptance criterion')


def main():
    for name,digest in INPUTS.items():
        if sha(REC/name)!=digest:raise ValueError('Completed synthetic JSON changed: '+name)
    records=[]
    for case,job,oldjob in [('failed0gene2',22041594,22040257),('positive18',22041595,22040259)]:
        directory=REC/('repair_p3_local_gl_'+case+'_'+str(job))
        new=json.loads((directory/'local_gl_comparison.json').read_text())
        settings=json.loads((directory/'protocol.json').read_text())['candidate_settings']
        if (settings['logmass_tolerance'],settings['gradient_tolerance'],settings['hessian_tolerance'])!=(1e-9,1e-8,1e-5):
            raise ValueError('Existing numerical budgets differ')
        old=json.loads((REC/('repair_p3_batched_gl_'+case+'_'+str(oldjob))/'candidate_baseline.json').read_text())['integration']
        indices=new['records'][0]['augmented_indices']
        for record in new['records']:
            actual=record['full_local_integration_report']
            ladders=[]
            for index,(a,b) in enumerate(zip(actual['ladder'],old['ladder'])):
                grids=('z_bound','grid_spacing','panels','order','node_evaluations','cumulative_evaluations')
                olddecision,newdecision=decisions(b,settings),decisions(a,settings)
                tail={key:scales(a[key],projection(b[key],indices,dim)) for key,dim in
                    [('relative_gradient_tail_allowance',1),('relative_hessian_tail_allowance',2)]}
                ladders.append(dict(ordinal=index,old_existing_decisions=olddecision,new_existing_decisions=newdecision,
                    existing_decisions_equal=olddecision==newdecision,grids_and_counts_equal=all(a[k]==b[k] for k in grids),
                    z_bound=a['z_bound'],order=a['order'],cumulative_evaluations=a['cumulative_evaluations'],
                    tail_allowance_scale_descriptions=tail,
                    original_absolute_checker_pass=record['ladder_comparisons'][index]['within_unchanged_allowances']))
            final=dict(old_within_existing_full_budget=within(old['empirical_refinement_plus_tail_allowances'],settings,1.),
                new_within_existing_full_budget=within(actual['empirical_refinement_plus_tail_allowances'],settings,1.),
                old_error_maxima=[old['empirical_refinement_plus_tail_allowances'][0],
                    maxima(old['empirical_refinement_plus_tail_allowances'][1]),maxima(old['empirical_refinement_plus_tail_allowances'][2])],
                new_error_maxima=[actual['empirical_refinement_plus_tail_allowances'][0],
                    maxima(actual['empirical_refinement_plus_tail_allowances'][1]),maxima(actual['empirical_refinement_plus_tail_allowances'][2])])
            records.append(dict(case=case,job=job,shape=record['shape'],ladder_count=len(ladders),ladder_records=ladders,
                all_saved_decisions_equal=all(x['existing_decisions_equal'] for x in ladders),
                all_grids_and_counts_equal=all(x['grids_and_counts_equal'] for x in ladders),final_existing_decisions=final,
                original_all_comparisons_flag=record['all_comparisons_within_unchanged_allowances'],
                primary_saved_differences=record['primary_differences'],integration_seconds=record['integration_seconds'],
                original_callables_unchanged=record['original_callables_unchanged']))
    for name,digest in INPUTS.items():
        if sha(REC/name)!=digest:raise ValueError('JSON input changed during collection')
    result=dict(method='Read saved JSON; no evaluator, arrays, fitting or quadrature',
        source_sha256=sha(Path(__file__)),input_sha256=INPUTS,records=records,
        new_accuracy_acceptance_rule=False,original_failure_flags_preserved=True,
        no_fullscale_validation=True,no_independent_probability_reference=True,no_calibration=True)
    output=REC/'repair_p3_local_gl_json_supplement_22041594_22041595.json'
    with output.open('x') as handle:json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(str(output),sha(output))
    for r in records:
        print(r['case'],r['shape'],'same_decisions',r['all_saved_decisions_equal'],
              'same_grids',r['all_grids_and_counts_equal'],'final',r['final_existing_decisions'])


if __name__=='__main__':main()
