#!/usr/bin/env python3
"""Bounded stdlib collector for four synthetic vectors across eight genes."""
import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
LAUNCH = REC/'repair_p3_bfix_vectors_launch_20261005T144634Z'
THETAS = ('theta1_start_like','theta2_truth_v0','theta3_truth_vpos','theta4_perturbed')
SOURCES = {
    'repair_p3_bfix_vectors.py':'08021479792640cac64ba32770498697d02ea4d6769e18b1285f9aab098e6ba5',
    'repair_p3_integral.py':'461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3',
    'repair_p3_composite_gl.py':'45b4423811fded58712bd5e85666f9ffc2b6972abbd9bd71964a584ec6d3c59c',
    'repair_p3_complex_reference.py':'f722eb76a2fb288f4f4e697dcd149d714b89221ec791ba6a43db8cfebb2f7ec3',
    'likelihood_v1.py':'41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6',
    'design_v1.py':'cb3163896f982f39ba94cdbf227a8e7d4c7f85a0a646b26542868e6fa14c9c5e',
    'bfix_v1.py':'324594ce783ba2464921eaee5924065f579060eded6043a48bdf4b7e392bfa56',
    'tests_v1.py':'52cd14d43b649c20347330b64b467dbd81320f10252fe85ce3c538ceef2f092e',
}
ENVIRONMENT = dict(python='3.10.18',numpy='2.2.6',scipy='1.15.2')


def require(ok,message):
    if not ok:
        raise ValueError(message)


def read(path):
    require(not path.is_symlink() and (path.resolve().is_relative_to(REC.resolve()) or
            path.resolve().is_relative_to(FIXTURE.resolve())),'Input escaped synthetic directories')
    require(path.stat().st_size <= 4*1024**2,'Small result file cap exceeded')
    body = path.read_bytes()
    return json.loads(body),hashlib.sha256(body).hexdigest()


def metadata():
    source,_ = read(LAUNCH/'source_hashes.json')
    require(source == SOURCES,'Frozen source manifest differs')
    for name,digest in SOURCES.items():
        require(hashlib.sha256((LAUNCH/name).read_bytes()).hexdigest() == digest,'Frozen mathematical source differs')
    entries = [json.loads(line) for line in (LAUNCH/'submissions.jsonl').read_text().splitlines()]
    keys = [(theta,gene) for theta in THETAS for gene in range(8)]
    require([int(e['job_id']) for e in entries] == list(range(22040192,22040224)), 'Exact thirty-two jobs required')
    for entry,(theta,gene) in zip(entries,keys):
        require(Path(entry['script']) == LAUNCH/f'{theta}_gene{gene}.sbatch','Vector/gene ownership differs')
    return entries,keys


def vector(values,size):
    return isinstance(values,list) and len(values) == size and all(isinstance(x,(int,float)) and math.isfinite(x) for x in values)


def collect(entries,keys):
    missing = [dict(job_id=e['job_id'],theta=t,gene=g) for e,(t,g) in zip(entries,keys)
               if not (REC/f"repair_p3_bfix_vectors_{e['job_id']}"/'summary.json').exists()]
    if missing:
        return dict(ready=False,missing_completed_summaries=missing,no_output_written=True)
    results,hashes,parameters = [],{},{}
    shared_inputs = None
    for theta in THETAS:
        path = FIXTURE/'params'/(theta+'.json')
        data,digest = read(path)
        p = data['alpha']+data['r']+[data['kappa']]+data['delta']+[data[k] for k in ('omega_S','rho_S','rho_b','rho_d')]
        require(vector(p,28) and math.isfinite(data['v']) and data['v'] >= 0,'Invalid prespecified vector')
        parameters[theta] = (p,float(data['v']))
        hashes[str(path)] = digest
    rows_by_gene = {}
    for entry,(theta,gene) in zip(entries,keys):
        directory = REC/f"repair_p3_bfix_vectors_{entry['job_id']}"
        documents = {}
        for name in ('protocol.json','summary.json',f'{theta}_gene{gene}.json'):
            documents[name],digest = read(directory/name)
            hashes[str(directory/name)] = digest
        protocol,summary = documents['protocol.json'],documents['summary.json']
        require(protocol['requested_points'] == [[theta,gene]] and protocol['requested_point_count'] == 1 and
                protocol['required_combined_points'] == 32 and protocol['all32_remain_required'] is True and
                protocol['participants'] == 60 and protocol['genes'] == 8 and protocol['rows'] == 173 and
                protocol['nuisance_coordinates'] == 28 and protocol['positive_variance_coordinates'] == 29 and
                protocol['environment'] == ENVIRONMENT and protocol['no_model_fits'] is True and
                protocol['no_calibration'] is True and protocol['independent_gradient_accuracy_verified'] is False and
                protocol['independent_Hessian_accuracy_verified'] is False,'Shard fixed protocol differs')
        all_pins = protocol['source_and_input_sha256']
        source_pins = {Path(path).name:digest for path,digest in all_pins.items() if Path(path).suffix == '.py'}
        require(source_pins == SOURCES,'Shard mathematical source pins differ')
        if shared_inputs is None:
            shared_inputs = all_pins
        require(all_pins == shared_inputs,'Synthetic fixture/input pins differ across shards')
        for vector_theta in THETAS:
            path = FIXTURE/'params'/(vector_theta+'.json')
            require(all_pins.get(str(path)) == hashes[str(path)],'Prespecified vector parameter bytes differ')
        require(summary['all32points_recorded'] is False and summary['all_requested_points_recorded'] is True and
                summary['requested_point_count'] == 1 and summary['required_combined_points'] == 32 and
                summary['all32_remain_required'] is True and summary['vector_totals'] == {} and
                summary['independent_gradient_accuracy_verified'] is False and
                summary['independent_Hessian_accuracy_verified'] is False and
                summary['frozen20_node_approximation_agreement_claimed'] is False and
                summary['no_model_fits'] is True and summary['no_calibration'] is True,'Shard subset summary differs')
        record = documents[f'{theta}_gene{gene}.json']
        p,v = parameters[theta]
        require(record['theta'] == theta and record['gene'] == gene and record['p'] == p and record['v'] == v and
                record['independent_gradient_accuracy_verified'] is False and
                record['independent_Hessian_accuracy_verified'] is False,'Saved vector point differs')
        require(isinstance(record['rows'],int) and record['rows'] > 0,'Invalid gene row count')
        rows_by_gene.setdefault(gene,record['rows'])
        require(rows_by_gene[gene] == record['rows'],'Gene rows changed across vectors')
        candidate = record.get('candidate')
        independent = record.get('independent_scalar')
        arithmetic_errors = []
        complete = (candidate is not None and independent is not None and
                    'candidate_refusal' not in record and 'independent_scalar_refusal' not in record)
        if candidate is not None:
            size = 28 if v == 0 else 29
            require(math.isfinite(candidate['logmass']) and vector(candidate['gradient'],29) and
                    len(candidate['hessian']) == size and all(vector(row,size) for row in candidate['hessian']),
                    'Accepted candidate numerical values/dimensions invalid')
            if v == 0:
                require(candidate['integration'].get('exact_variance_zero_branch') is True and
                        candidate['integration'].get('positive_variance_Hessian_not_applicable') is True,
                        'Boundary nuisance Hessian scope differs')
        delta = per_gene_pass = None
        if independent is not None:
            require(math.isfinite(independent['logmass']),'Independent scalar likelihood nonfinite')
        if complete:
            delta = candidate['logmass']-independent['logmass']
            per_gene_pass = abs(delta) <= 1e-6/8
            if record.get('candidate_scalar_logmass_difference') != delta or record.get('candidate_scalar_per_gene_budget_pass') is not per_gene_pass:
                arithmetic_errors.append('Saved per-gene scalar arithmetic differs')
        stderr = Path(entry['stderr'])
        require(stderr.parent == LAUNCH and stderr.stat().st_size <= 1024**2,'Unexpected shard stderr')
        results.append(dict(job_id=entry['job_id'],theta=theta,gene=gene,record=record,
                            numerical_pair_available=complete,rederived_scalar_difference=delta,
                            per_gene_budget_descriptive_pass=per_gene_pass,arithmetic_errors=arithmetic_errors,
                            stderr=stderr.read_text(errors='replace')))
    require(sum(rows_by_gene.values()) == 173,'Combined gene rows differ from original fixture')
    totals = {}
    for theta in THETAS:
        points = [r for r in results if r['theta'] == theta]
        valid = len(points) == 8 and all(r['numerical_pair_available'] and not r['arithmetic_errors'] for r in points)
        total = dict(all_eight_numerical_pairs_available=valid,likelihood_only_pass=False)
        if valid:
            v = parameters[theta][1]
            size = 28 if v == 0 else 29
            ll = math.fsum(r['record']['candidate']['logmass'] for r in points)
            scalar = math.fsum(r['record']['independent_scalar']['logmass'] for r in points)
            total.update(candidate_logmass=ll,independent_scalar_logmass=scalar,
                gradient=[math.fsum(r['record']['candidate']['gradient'][j] for r in points) for j in range(29)],
                hessian=[[math.fsum(r['record']['candidate']['hessian'][i][j] for r in points) for j in range(size)] for i in range(size)],
                hessian_coordinate_scope='nuisance28_only' if v == 0 else 'all29_positive_variance',
                scalar_logmass_difference=ll-scalar,likelihood_only_pass=abs(ll-scalar) <= 1e-6,
                per_gene_budgets_descriptive_all_pass=all(r['per_gene_budget_descriptive_pass'] for r in points))
        totals[theta] = total
    return dict(ready=True,execution_complete=True,completed_points=32,required_points=32,
        shared_source_and_input_pins=shared_inputs,environment=ENVIRONMENT,input_hashes=hashes,
        results=results,vector_totals=totals,
        all_four_vector_scalar_likelihood_pass=all(t['likelihood_only_pass'] for t in totals.values()),
        combined_likelihood_absolute_tolerance=1e-6,per_gene_budget_descriptive=1e-6/8,
        independent_gradient_accuracy_verified=False,independent_Hessian_accuracy_verified=False,
        full_model_fit_verified=False,statistical_calibration_verified=False,
        frozen20_node_approximation_agreement_claimed=False,
        empirical_quadrature_and_floating_error_bounds_not_certified=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect',action='store_true')
    args = parser.parse_args()
    entries,keys = metadata()
    if args.inspect:
        print(json.dumps(dict(metadata_verified=True,required_points=32,cartesian_keys=keys,
                              numerical_array_reads=False,output_written=False)))
        return
    result = collect(entries,keys)
    if not result['ready']:
        print(json.dumps(result))
        return
    body = json.dumps(result,indent=2,allow_nan=False)+'\n'
    require(len(body.encode()) <= 32*1024**2,'Collection output cap exceeded')
    target = LAUNCH/'bfix_vectors_collection.json'
    with target.open('x') as handle:
        handle.write(body)
    print(json.dumps(dict(output=str(target),execution_complete=True,
                          likelihood_only_pass=result['all_four_vector_scalar_likelihood_pass'])))


if __name__ == '__main__':
    main()
