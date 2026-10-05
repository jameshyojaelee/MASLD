#!/usr/bin/env python3
"""Stdlib reconciliation of one four-point synthetic mixture-policy experiment."""
import argparse
import hashlib
import json
import math
from pathlib import Path

REC = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z')
LAUNCH = REC/'repair_p3_mixture_gh_policies_launch_20261005T151726Z'
KEYS = [('failed0gene2','scale_precision'),('failed0gene2','scale_sqrt'),
        ('positive18','scale_precision'),('positive18','scale_sqrt')]
SOURCES = {
    'repair_p3_fixed_point_check.py':'97687ec12978a66a4dc0ded0299a5ea0b326696716a8a1b896a5ddb59ef7e027',
    'repair_p3_integral.py':'461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3',
    'repair_p3_mixture_gh.py':'68f8311749795c145ce836bc56876255ba1ad967af04df2342a19eb121cae24b',
    'repair_p3_reference_fast.py':'2eaaf549739f64322d394f0674aadba98f0f2ae6bedda1df56db8e9a8a3949fd',
    'codex_stable_bb_probe.py':'723781113c17061b9a28341a91149bbd4ab71f195121a67f88d8cbc710a21488',
}
PARENTS = {
    'failed0gene2':(22039200,{
        'protocol.json':'7b5d793be5485a25e5fc7f5b3d97c4ddf30a1a497d697728cd64468ed61f8d40',
        'candidate_baseline.json':'edec063451383f9a1386619186c2089d3ef8e6ee34d854e80a2894fdd28e0f9a',
        'independent_baseline.json':'f621bcc6067bc0577f7cb35ef9c82197c35c49f1fc4e53a7d2b321cfd238fab0',
        'summary.json':'5ad1ea571c71525e6ae1eb18bf625baf1e6779934d64a149e0d16ca5f3137a0b'}),
    'positive18':(22039201,{
        'protocol.json':'ed6c190424fa29501e1377b5b47b12762917e53f46338c3c3f1a408c3f9e4774',
        'candidate_baseline.json':'4a576370be2f33ac28bf23052fe283b4a5f614bc30e983bd6a84224da857c722',
        'independent_baseline.json':'3c9342dd88fc449ef91458e534116a6c721a94f384356b6ea33d7880fec199b2',
        'summary.json':'1daa3275c8d524857c90369b7686fc2234229df765e41f139a4b232c9e323613'}),
}
FIELDS = ('case','draw','fixed_gene','fixed_gene_rows','fixed_p','fixed_v','count_vector_sha256',
          'environment','candidate_settings','unchanged_outer_methods','finite_difference_steps',
          'required_combined_derivative_family_axes','reference_arithmetic')


def require(ok,message):
    if not ok:
        raise ValueError(message)


def read(path):
    require(path.resolve().is_relative_to(REC.resolve()) and not path.is_symlink(),'Input escaped synthetic REC')
    require(path.stat().st_size <= 16*1024**2,'Bounded result file cap exceeded')
    body = path.read_bytes()
    return json.loads(body),hashlib.sha256(body).hexdigest()


def metadata():
    pins,_ = read(LAUNCH/'source_hashes.json')
    require(pins == SOURCES,'Frozen source manifest differs')
    for name,digest in SOURCES.items():
        require(hashlib.sha256((LAUNCH/name).read_bytes()).hexdigest() == digest,'Frozen mathematical source differs')
    entries = [json.loads(line) for line in (LAUNCH/'submissions.jsonl').read_text().splitlines()]
    require([int(e['job_id']) for e in entries] == list(range(22040801,22040805)),'Exact four submitted jobs required')
    for entry,(case,policy) in zip(entries,KEYS):
        require(Path(entry['script']) == LAUNCH/f'{case}_{policy}.sbatch','Point/policy ownership differs')
    return entries


def vec(value):
    return isinstance(value,list) and len(value) == 29 and all(isinstance(x,(float,int)) and math.isfinite(x) for x in value)


def matrix(value):
    return isinstance(value,list) and len(value) == 29 and all(vec(row) for row in value)


def max_vector_difference(a,b):
    return max(abs(x-y) for x,y in zip(a,b))


def reconcile(order,policy,settings):
    components = order['component_contributions']
    require(len(components) == order['component_count'] and len(components) > 0 and
            [c['index'] for c in components] == list(range(len(components))), 'Component ownership/count differs')
    require(order['prior_policy'] == policy and order['complete_mixture_denominator_used'] is True and
            order['prior_normalization_included_exactly_once'] is True, 'Mixture probability identity differs')
    require(vec(order['gradient']) and matrix(order['hessian']) and vec(order['centered_moment_anchor']) and
            vec(order['posterior_centered_score_mean']), 'Global moment dimensions/values invalid')
    for c in components:
        require(math.isfinite(c['log_mass_contribution']) and
                math.isfinite(c['posterior_normalized_mass_contribution']) and c['posterior_normalized_mass_contribution'] > 0 and
                c['proposal_weight'] > 0 and math.isfinite(c['proposal_weight']) and
                vec(c['posterior_weighted_score_contribution']) and vec(c['posterior_weighted_centered_score_contribution']) and
                matrix(c['posterior_weighted_conditional_H_contribution']) and
                matrix(c['posterior_weighted_centered_score_outer_contribution']), 'Component numerical values invalid')
    anchor = order['centered_moment_anchor']
    mass = math.fsum(c['posterior_normalized_mass_contribution'] for c in components)
    maximum = max(c['log_mass_contribution'] for c in components)
    ll = maximum+math.log(math.fsum(math.exp(c['log_mass_contribution']-maximum) for c in components))
    scores = [math.fsum(c['posterior_weighted_score_contribution'][j] for c in components) for j in range(29)]
    centered = [math.fsum(c['posterior_weighted_centered_score_contribution'][j] for c in components) for j in range(29)]
    hessian = [[math.fsum(c['posterior_weighted_conditional_H_contribution'][i][j]+
                            c['posterior_weighted_centered_score_outer_contribution'][i][j] for c in components)
                -centered[i]*centered[j] for j in range(29)] for i in range(29)]
    anchor_relation = max(abs(c['posterior_weighted_score_contribution'][j]-
        (c['posterior_weighted_centered_score_contribution'][j]+anchor[j]*c['posterior_normalized_mass_contribution']))
        for c in components for j in range(29))
    peaks = len(components)-1
    epsilon = order['prior_policy_precision_epsilon']
    if peaks == 0:
        expected_weights = [1.]
    else:
        require(components[0]['proposal']['kind'] == 'prior' and components[0]['proposal']['mean_z'] == 0 and
                components[0]['proposal']['sd_z'] == 1,'Prior proposal differs')
        scale = min(1.,min(c['proposal']['sd_z'] for c in components[1:]))
        require(scale > 0 and math.isfinite(scale),'Invalid peak scale')
        if policy == 'scale_sqrt':
            require(epsilon is None,'Unexpected precision epsilon for square-root policy')
            prior = min(1/(peaks+1),scale*math.sqrt(settings['logmass_tolerance']))
        else:
            require(isinstance(epsilon,(float,int)) and math.isfinite(epsilon) and
                    0 < epsilon <= .01*settings['logmass_tolerance'],'Invalid saved precision epsilon')
            prior = min(1/(peaks+1),scale*epsilon)
        expected_weights = [prior]+[(1-prior)/peaks]*peaks
    weights = [c['proposal_weight'] for c in components]
    require(weights == order['proposal_weights'],'Component and global proposal weights differ')
    return dict(order=order['order'],z_bound=order['z_bound'],grid_spacing=order['grid_spacing'],
        component_mass_fraction_sum=mass,mass_fraction_sum_minus_one=mass-1,
        component_logsumexp_minus_global_logmass=ll-order['logmass'],
        weighted_score_max_reconciliation_residual=max_vector_difference(scores,order['gradient']),
        centered_score_max_reconciliation_residual=max_vector_difference(centered,order['posterior_centered_score_mean']),
        anchor_plus_centered_mean_max_residual=max_vector_difference([a+b for a,b in zip(anchor,centered)],order['gradient']),
        per_component_anchor_relation_max_residual=anchor_relation,
        hessian_max_reconciliation_residual=max(max_vector_difference(a,b) for a,b in zip(hessian,order['hessian'])),
        proposal_recipe_max_residual=max_vector_difference(weights,expected_weights),
        proposal_weight_sum_minus_one=math.fsum(weights)-1,
        precision_epsilon_from_saved_metadata_not_independently_rederived=policy=='scale_precision',
        residuals_are_arithmetic_diagnostics_not_accuracy_certificates=True)


def collect(entries):
    missing = [dict(job_id=e['job_id'],case=c,policy=p) for e,(c,p) in zip(entries,KEYS)
               if not (REC/f"repair_p3_fixed_point_{c}_{e['job_id']}"/'summary.json').exists()]
    if missing:
        return dict(ready=False,missing_completed_summaries=missing,no_output_written=True)
    parents,hashes = {},{}
    for case,(job,pins) in PARENTS.items():
        directory = REC/f'repair_p3_fixed_point_{case}_{job}'
        parents[case] = {}
        for name,digest in pins.items():
            value,actual = read(directory/name)
            require(actual == digest,'Original GL baseline bytes changed')
            parents[case][name] = value
            hashes[str(directory/name)] = actual
    results = []
    for entry,(case,policy) in zip(entries,KEYS):
        directory = REC/f"repair_p3_fixed_point_{case}_{entry['job_id']}"
        protocol,summary,independent = [read(directory/name)[0] for name in ('protocol.json','summary.json','independent_baseline.json')]
        parent = parents[case]
        for key in FIELDS:
            require(protocol[key] == parent['protocol.json'][key],'Fixed point differs: '+key)
        require(protocol['mixture_gh'] is True and protocol['composite_gl'] is False and
                protocol['mixture_prior_policy'] == policy and protocol['baseline_only'] is True and
                protocol['baseline_reuse'] is False,'Policy/implementation identity differs')
        pins = {Path(path).name:digest for section in ('source_guards','own_source_hashes')
                for path,digest in protocol[section].items() if Path(path).name in SOURCES}
        require(pins == SOURCES,'Shard mathematical source pins differ')
        for name in ('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
                     'repair_p3_integral.py','repair_p3_reference_fast.py','codex_stable_bb_probe.py'):
            current = {digest for section in ('source_guards','own_source_hashes')
                       for path,digest in protocol[section].items() if Path(path).name == name}
            previous = {digest for section in ('source_guards','own_source_hashes')
                        for path,digest in parent['protocol.json'][section].items() if Path(path).name == name}
            require(len(current) == 1 and current == previous,'Shared mathematical source differs: '+name)
        require(all(math.isfinite(independent[key]['logmass']) for key in ('first','tighter_expanded')),
                'Independent scalar baseline nonfinite')
        require(summary['fixed_point_check_completed'] is True and summary['baseline_only'] is True and
                summary['executed_derivative_axes'] == [] and summary['derivative_checks_completed'] is False and
                summary['no_model_fits'] is True and summary['no_statistical_calibration'] is True,'Baseline summary differs')
        record = dict(job_id=entry['job_id'],case=case,policy=policy,raw_summary=summary,independent_baseline=independent)
        candidate_path = directory/'candidate_baseline.json'
        if candidate_path.exists():
            candidate,digest = read(candidate_path)
            hashes[str(candidate_path)] = digest
            require(summary['candidate_integral_returned'] is True and vec(candidate['gradient']) and
                    matrix(candidate['hessian']) and math.isfinite(candidate['logmass']), 'Accepted MGH values invalid')
            report = candidate['integration']
            record['candidate_baseline'] = candidate
            record['component_reconciliations'] = [reconcile(order,policy,protocol['candidate_settings']) for order in report['ladder']]
            gl = parent['candidate_baseline.json']
            scalar = independent['tighter_expanded']['logmass']
            record['comparison'] = dict(MGH_minus_GL_logmass=candidate['logmass']-gl['logmass'],
                MGH_minus_own_scalar_logmass=candidate['logmass']-scalar,
                MGH_minus_original_scalar_logmass=candidate['logmass']-parent['independent_baseline.json']['tighter_expanded']['logmass'],
                gradient_max_difference_vs_GL=max_vector_difference(candidate['gradient'],gl['gradient']),
                hessian_max_difference_vs_GL=max(max_vector_difference(a,b) for a,b in zip(candidate['hessian'],gl['hessian'])),
                empirical_refinement_plus_tail_allowances=report['empirical_refinement_plus_tail_allowances'],
                differences_and_allowances_not_independent_derivative_accuracy=True)
        else:
            refusal,digest = read(directory/'candidate_refusal.json')
            hashes[str(directory/'candidate_refusal.json')] = digest
            require(summary['candidate_integral_returned'] is False,'Refusal and summary disagree')
            record['candidate_refusal'] = refusal
            record['component_reconciliations'] = [
                reconcile(order,policy,protocol['candidate_settings'])
                for order in refusal.get('integration_reports',[]) or []
                if isinstance(order,dict) and 'component_contributions' in order]
        for name in ('protocol.json','summary.json','independent_baseline.json'):
            hashes[str(directory/name)] = read(directory/name)[1]
        stderr = Path(entry['stderr'])
        require(stderr.parent == LAUNCH and stderr.stat().st_size <= 1024**2,'Unexpected shard stderr')
        record['stderr'] = stderr.read_text(errors='replace')
        results.append(record)
    return dict(ready=True,execution_complete=True,completed_points=4,required_points=4,results=results,input_hashes=hashes,
        all_four_candidate_integrals_returned=all('candidate_baseline' in r for r in results),
        numerical_policy_selection_not_automatic=True,independent_derivative_accuracy_verified=False,
        full_model_fit_verified=False,statistical_calibration_verified=False,
        empirical_error_allowances_not_certified=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect',action='store_true')
    args = parser.parse_args()
    entries = metadata()
    if args.inspect:
        print(json.dumps(dict(metadata_verified=True,required_points=4,point_policies=KEYS,
                              numerical_array_reads=False,output_written=False)))
        return
    result = collect(entries)
    if not result['ready']:
        print(json.dumps(result))
        return
    body = json.dumps(result,indent=2,allow_nan=False)+'\n'
    require(len(body.encode()) <= 32*1024**2,'Collection output cap exceeded')
    target = LAUNCH/'mixture_policy_collection.json'
    with target.open('x') as handle:
        handle.write(body)
    print(json.dumps(dict(output=str(target),execution_complete=True,independent_derivative_accuracy_verified=False)))


if __name__ == '__main__':
    main()
