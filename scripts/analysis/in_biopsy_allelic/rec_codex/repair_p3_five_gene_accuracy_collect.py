#!/usr/bin/env python3
"""Collect one frozen fifteen-point synthetic likelihood accuracy experiment."""
import argparse
import hashlib
import json
import math
from pathlib import Path

REC = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z')
LAUNCH = REC/'repair_p3_five_gene_accuracy_launch_20261005T143525Z'
JOBS = list(range(22039932,22039945))+[22039946,22039947]
MULTIPLES = (.25,1.,4.)
COUNTS = (5,25,75,200,500)
SOURCES = {
    'repair_p3_five_gene_accuracy.py':'4aec076799a157551d86c1c2fc3da8901e436f8779cca190643c3398fd0f0b21',
    'repair_p3_integral.py':'461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3',
    'repair_p3_composite_gl.py':'45b4423811fded58712bd5e85666f9ffc2b6972abbd9bd71964a584ec6d3c59c',
    'check_quadrature_v1.py':'40a51fd0090d797b631c897f1c35c046b750e040a780432f0f614960b9a17735',
    'likelihood_v1.py':'41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6',
}


def require(ok,message):
    if not ok:
        raise ValueError(message)


def read(path):
    require(path.resolve().is_relative_to(REC.resolve()) and not path.is_symlink(),'Input escaped synthetic REC')
    require(path.stat().st_size <= 4*1024**2,'Small result file cap exceeded')
    body = path.read_bytes()
    return json.loads(body),hashlib.sha256(body).hexdigest()


def metadata():
    source,_ = read(LAUNCH/'source_hashes.json')
    require(source == SOURCES,'Frozen source manifest differs')
    for name,digest in SOURCES.items():
        require(hashlib.sha256((LAUNCH/name).read_bytes()).hexdigest() == digest,'Frozen mathematical source differs')
    submissions = [json.loads(line) for line in (LAUNCH/'submissions.jsonl').read_text().splitlines()]
    require([int(e['job_id']) for e in submissions] == JOBS,'Exact fifteen submitted jobs required')
    keys = [(gene,multiple) for gene in range(5) for multiple in MULTIPLES]
    for entry,(gene,multiple) in zip(submissions,keys):
        require(Path(entry['script']) == LAUNCH/f'gene{gene}_tau{multiple}.sbatch','Cartesian point ownership differs')
    return submissions,keys


def collect(submissions,keys):
    missing = [dict(job_id=e['job_id'],gene=g,tau_multiple=m)
               for e,(g,m) in zip(submissions,keys)
               if not (REC/f"repair_p3_five_gene_accuracy_{e['job_id']}"/'summary.json').exists()]
    if missing:
        return dict(ready=False,missing_completed_summaries=missing,no_output_written=True)
    results,hashes = [],{}
    fixture_sha = environment = None
    for entry,(gene,multiple) in zip(submissions,keys):
        directory = REC/f"repair_p3_five_gene_accuracy_{entry['job_id']}"
        documents = {}
        for name in ('protocol.json','fixture.json','summary.json',f'gene{gene}_tau{multiple}.json'):
            documents[name],digest = read(directory/name)
            hashes[str(directory/name)] = digest
        protocol,fixture,summary = (documents[name] for name in ('protocol.json','fixture.json','summary.json'))
        require(protocol['seed'] == 20260923 and protocol['rows_per_gene'] == list(COUNTS) and
                protocol['tau_max'] == .05 and protocol['tau_multiples'] == list(MULTIPLES) and
                protocol['required_points'] == 15 and protocol['requested_points'] == [[multiple,gene]] and
                protocol['requested_point_count'] == 1 and protocol['combined_family_all15_remains_required'] is True and
                protocol['no_model_fits'] is True and protocol['no_calibration'] is True,'Shard protocol differs')
        pins = {Path(path).name:digest for path,digest in protocol['source_sha256'].items()}
        require(pins == SOURCES,'Shard mathematical source pins differ')
        require(summary['all15points_recorded'] is False and summary['all_requested_points_recorded'] is True and
                summary['requested_point_count'] == 1 and summary['combined_family_all15_remains_required'] is True and
                summary['candidate_scalar_relative_LL_pass_all15'] is False and
                summary['original_GH20_GH40_pass_all15'] is False and summary['original_GH20_scalar_pass_all15'] is False and
                summary['independent_derivative_accuracy_verified'] is False and summary['no_model_fits'] is True and
                summary['no_calibration'] is True,'Shard subset/family summary differs')
        current_fixture_sha = hashes[str(directory/'fixture.json')]
        if fixture_sha is None:
            fixture_sha,environment = current_fixture_sha,protocol['environment']
        require(current_fixture_sha == fixture_sha and protocol['environment'] == environment,
                'Fixture parameters/drawn counts/environment differ across shards')
        p = fixture['parameters']
        require(len(p) == 16 and all(math.isfinite(x) for x in p) and
                p[:5] == [-1.2,-.7,.4,.9,1.3] and p[10] == .08 and p[11] == .04 and
                fixture['rows_per_gene'] == list(COUNTS) and len(fixture['drawn_counts']) == sum(COUNTS),
                'Prespecified full synthetic fixture differs')
        record = documents[f'gene{gene}_tau{multiple}.json']
        require(record['gene'] == gene and record['rows'] == COUNTS[gene] and
                record['tau'] == multiple*.05 and record['variance'] == (multiple*.05)**2,
                'Saved point differs')
        arithmetic_errors = []
        summary_boolean_type_mismatches = []
        checks = {}
        scalar = record.get('scalar',{}).get('expanded',{})
        scalar_ok = ('scalar_refusal' not in record and 'logmass' in scalar and math.isfinite(scalar['logmass']))
        if scalar_ok:
            first = record['scalar']['first']
            scalar_refinement = abs(scalar['logmass']-first['logmass'])+scalar['estimated_logmass_error']
            scalar_ok = (math.isfinite(scalar_refinement) and scalar_refinement <= 1e-9 and
                         first['relative_prior_mass_tail_allowance'] <= 1e-10)
            checks['scalar_refinement_rederived'] = scalar_refinement
        for name in ('candidate','original_GH20','original_GH40'):
            value = record.get(name,{})
            valid = not value.get('refused',False) and 'logmass' in value and math.isfinite(value['logmass'])
            passed = False
            if scalar_ok and valid:
                delta = value['logmass']-scalar['logmass']
                limit = 1e-8*abs(scalar['logmass'])
                passed = abs(delta) <= limit
                checks[name] = dict(scalar_difference=delta,likelihood_relative_absolute_limit=limit,
                                   likelihood_relative_pass=passed)
                saved_delta = record.get('candidate_scalar_difference') if name == 'candidate' else value.get('scalar_difference')
                saved_pass = record.get('candidate_scalar_relative_LL_pass') if name == 'candidate' else value.get('scalar_relative_LL_pass')
                if saved_delta != delta or saved_pass is not passed:
                    arithmetic_errors.append(name+' saved scalar arithmetic differs')
            else:
                checks[name] = dict(likelihood_relative_pass=False,missing_or_refused_comparison=True)
            if name == 'candidate' and summary['candidate_scalar_relative_LL_pass_requested_points'] is not passed:
                known_numpy_boolean_pattern = (passed is True and
                    record.get('candidate_scalar_relative_LL_pass') is True and
                    summary['candidate_scalar_relative_LL_pass_requested_points'] is False)
                summary_boolean_type_mismatches.append(dict(
                    saved_requested_summary=summary['candidate_scalar_relative_LL_pass_requested_points'],
                    saved_point_criterion=record.get('candidate_scalar_relative_LL_pass'),
                    rederived_criterion=passed,
                    consistent_with_frozen_numpy_boolean_identity_bug=known_numpy_boolean_pattern))
                if not known_numpy_boolean_pattern:
                    arithmetic_errors.append('Unexpected requested-point summary mismatch')
        gh20,gh40 = record['original_GH20'],record['original_GH40']
        gh_pass = False
        if all(not x.get('refused',False) and 'logmass' in x and math.isfinite(x['logmass']) for x in (gh20,gh40)):
            delta = gh20['logmass']-gh40['logmass']
            gh_pass = abs(delta) <= 1e-8
            checks['original_GH20_GH40'] = dict(difference=delta,absolute_limit=1e-8,absolute_pass=gh_pass)
            if record.get('original_GH20_GH40_difference') != delta or record.get('original_GH20_GH40_absolute_pass') is not gh_pass:
                arithmetic_errors.append('Saved GH20/GH40 arithmetic differs')
        else:
            checks['original_GH20_GH40'] = dict(absolute_pass=False,missing_or_refused_comparison=True)
        stderr = Path(entry['stderr'])
        require(stderr.parent == LAUNCH and stderr.stat().st_size <= 1024**2,'Unexpected shard stderr')
        results.append(dict(job_id=entry['job_id'],gene=gene,tau_multiple=multiple,record=record,
                            rederived_checks=checks,arithmetic_errors=arithmetic_errors,
                            raw_summary=summary,
                            summary_boolean_type_mismatches=summary_boolean_type_mismatches,
                            stderr=stderr.read_text(errors='replace')))
    consistent = all(not r['arithmetic_errors'] for r in results)
    return dict(ready=True,execution_complete=True,completed_points=15,required_points=15,
        shared_fixture_sha256=fixture_sha,environment=environment,results=results,input_hashes=hashes,
        arithmetic_consistent=consistent,
        requested_summary_consistent=all(not r['summary_boolean_type_mismatches'] for r in results),
        frozen_requested_summary_boolean_bug='Frozen checker uses numpy.bool_ identity with True; preserved raw summaries may incorrectly report a false requested-point pass.',
        candidate_scalar_likelihood_accuracy_all15=consistent and all(r['rederived_checks']['candidate']['likelihood_relative_pass'] for r in results),
        original_GH20_scalar_likelihood_accuracy_all15=consistent and all(r['rederived_checks']['original_GH20']['likelihood_relative_pass'] for r in results),
        original_GH40_scalar_likelihood_accuracy_all15=consistent and all(r['rederived_checks']['original_GH40']['likelihood_relative_pass'] for r in results),
        original_GH20_GH40_absolute_agreement_all15=consistent and all(r['rederived_checks']['original_GH20_GH40']['absolute_pass'] for r in results),
        criterion='abs(candidate_LL-reference_LL) <= 1e-8*abs(reference_LL); GH20/GH40 absolute difference <=1e-8',
        relative_probability_accuracy_not_claimed=True,independent_derivative_accuracy_verified=False,
        full_model_fit_verified=False,statistical_calibration_verified=False,
        empirical_quadrature_and_floating_error_not_certified=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect',action='store_true')
    args = parser.parse_args()
    submissions,keys = metadata()
    if args.inspect:
        print(json.dumps(dict(metadata_verified=True,required_points=15,cartesian_keys=keys,
                              numerical_array_reads=False,output_written=False)))
        return
    result = collect(submissions,keys)
    if not result['ready']:
        print(json.dumps(result))
        return
    body = json.dumps(result,indent=2,allow_nan=False)+'\n'
    require(len(body.encode()) <= 32*1024**2,'Collection output cap exceeded')
    target = LAUNCH/'five_gene_accuracy_collection.json'
    with target.open('x') as handle:
        handle.write(body)
    print(json.dumps(dict(output=str(target),execution_complete=True,
                          candidate_likelihood_accuracy=result['candidate_scalar_likelihood_accuracy_all15'])))


if __name__ == '__main__':
    main()
