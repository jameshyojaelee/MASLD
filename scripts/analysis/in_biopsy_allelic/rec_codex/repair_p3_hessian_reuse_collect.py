#!/usr/bin/env python3
"""Bounded stdlib collector for one frozen synthetic derivative campaign."""
import argparse
import hashlib
import json
import math
from pathlib import Path

REC = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z')
LAUNCH = REC/'repair_p3_hessian_reuse_launch_20261005T142226Z'
AXES = (2,10,16,17,24,25,28)
STEPS = (.001,.0003,.0001)
SOURCES = {
    'repair_p3_fixed_point_check.py':'9ce9ae5ab4fd72312ffdd7e1f88cfd5917bbd25da2ab93e1366575d9bd200dcf',
    'repair_p3_integral.py':'461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3',
    'repair_p3_baseline_reuse.py':'9de3fc38f954f5b101d1e03a057ce5501943910b3c32fbccbb6dc9ff23a96327',
    'repair_p3_composite_gl.py':'45b4423811fded58712bd5e85666f9ffc2b6972abbd9bd71964a584ec6d3c59c',
    'repair_p3_reference_fast.py':'2eaaf549739f64322d394f0674aadba98f0f2ae6bedda1df56db8e9a8a3949fd',
    'codex_stable_bb_probe.py':'723781113c17061b9a28341a91149bbd4ab71f195121a67f88d8cbc710a21488',
}
FIELDS = ('case','draw','fixed_gene','fixed_gene_rows','fixed_p','fixed_v','count_vector_sha256',
          'environment','candidate_settings','unchanged_outer_methods','finite_difference_steps',
          'required_combined_derivative_family_axes','reference_arithmetic','composite_gl')


def require(ok,message):
    if not ok:
        raise ValueError(message)


def read(path):
    require(path.resolve().is_relative_to(REC.resolve()) and not path.is_symlink(), 'Input escaped synthetic REC')
    require(path.stat().st_size <= 2*1024**2, 'Small result file cap exceeded')
    body = path.read_bytes()
    return json.loads(body),hashlib.sha256(body).hexdigest()


def finite(value):
    if isinstance(value,dict):
        return all(finite(v) for v in value.values())
    if isinstance(value,list):
        return all(finite(v) for v in value)
    return math.isfinite(value) if isinstance(value,(int,float)) else True


def metadata():
    sources,_ = read(LAUNCH/'source_hashes.json')
    require(sources == SOURCES, 'Frozen source manifest differs')
    for name,digest in SOURCES.items():
        require(hashlib.sha256((LAUNCH/name).read_bytes()).hexdigest() == digest, 'Frozen mathematical source differs')
    entries = [json.loads(line) for line in (LAUNCH/'submissions.jsonl').read_text().splitlines()]
    require(len(entries) == 14 and [int(x['job_id']) for x in entries] == list(range(22039656,22039670)),
            'Exact fourteen submitted jobs required')
    for ordinal,entry in enumerate(entries):
        case = 'failed0gene2' if ordinal < 7 else 'positive18'
        axis = AXES[ordinal%7]
        require(Path(entry['script']) == LAUNCH/f'{case}_axis{axis:02d}.sbatch','Shard ownership differs')
    return entries


def collect(entries):
    missing = []
    for ordinal,entry in enumerate(entries):
        case = 'failed0gene2' if ordinal < 7 else 'positive18'
        directory = REC/f"repair_p3_fixed_point_{case}_{entry['job_id']}"
        if not (directory/'summary.json').exists():
            missing.append(dict(job_id=entry['job_id'],case=case,axis=AXES[ordinal%7],
                                failure_file_present=(directory/'failure.json').exists()))
    if missing:
        return dict(ready=False,missing_completed_summaries=missing,no_output_written=True)
    parent_protocols,parent_pins = {},{}
    hashes = {}
    for case,parent_job in (('failed0gene2',22039200),('positive18',22039201)):
        pins,pins_sha = read(LAUNCH/f'{case}_parent_hashes.json')
        require(set(pins) == {'protocol.json','candidate_baseline.json','independent_baseline.json','summary.json'},
                'Four parent hashes required')
        parent = REC/f'repair_p3_fixed_point_{case}_{parent_job}'
        documents = {}
        for name,digest in pins.items():
            documents[name],actual = read(parent/name)
            require(actual == digest,'Parent result bytes changed')
            hashes[str(parent/name)] = actual
        require(documents['summary.json']['baseline_only'] is True and
                documents['summary.json']['fixed_point_check_completed'] is True and
                documents['summary.json']['candidate_integral_returned'] is True,
                'Parent baseline computation incomplete')
        parent_protocols[case],parent_pins[case] = documents['protocol.json'],pins
        hashes[str(LAUNCH/f'{case}_parent_hashes.json')] = pins_sha
    shards,refusals = [],[]
    for ordinal,entry in enumerate(entries):
        case = 'failed0gene2' if ordinal < 7 else 'positive18'
        axis = AXES[ordinal%7]
        directory = REC/f"repair_p3_fixed_point_{case}_{entry['job_id']}"
        protocol,summary,records = [read(directory/name)[0] for name in ('protocol.json','summary.json','finite_differences.json')]
        for key in FIELDS:
            require(protocol[key] == parent_protocols[case][key],'Shard fixed input differs: '+key)
        require(protocol.get('baseline_reuse') is True and protocol.get('reused_parent_file_hashes') == parent_pins[case] and
                protocol['requested_derivative_axes'] == [axis], 'Shard parent/axis differs')
        pins = {Path(path).name:digest for section in ('source_guards','own_source_hashes')
                for path,digest in protocol[section].items() if Path(path).name in SOURCES}
        require(pins == SOURCES,'Shard mathematical source pins differ')
        require(summary['fixed_point_check_completed'] is True and summary['candidate_integral_returned'] is True and
                summary['no_model_fits'] is True and summary['no_statistical_calibration'] is True and
                summary['baseline_only'] is False and summary['derivative_checks_completed'] is False and
                summary['executed_derivative_axes'] == [axis] and
                summary['finite_difference_steps'] == list(STEPS) and
                summary['required_combined_derivative_family_axes'] == list(AXES),
                'Shard execution/family summary differs')
        require(len(records) == 3 and [r['step'] for r in records] == list(STEPS) and
                all(r['axis'] == axis for r in records) and finite(records), 'Three finite original-step records required')
        diagnostics = []
        for record in records:
            suffix = str(record['step']).replace('.','_')
            separately,record_sha = read(directory/f'finite_difference_{axis:02d}_{suffix}.json')
            require(separately == record,'Individual and assembled finite differences disagree')
            hashes[str(directory/f'finite_difference_{axis:02d}_{suffix}.json')] = record_sha
            if record['step'] in STEPS[1:]:
                if 'hessian_score_difference_refused' in record:
                    refusals.append(dict(job_id=entry['job_id'],axis=axis,step=record['step'],
                                         refusal=record['hessian_score_difference_refused']))
                else:
                    vectors = [record[name] for name in ('hessian_moment_column','score_FD_hessian_column',
                               'Hessian_column_difference','estimated_score_integration_FD_error')]
                    require(all(len(v) == 29 for v in vectors),'Complete 29-coordinate Hessian columns required')
                    reconstructed = [b-a for a,b in zip(vectors[0],vectors[1])]
                    diagnostics.append(dict(step=record['step'],
                        max_absolute_column_difference=max(abs(x) for x in vectors[2]),
                        max_saved_difference_arithmetic_residual=max(abs(a-b) for a,b in zip(reconstructed,vectors[2])),
                        max_estimated_integration_FD_error=max(vectors[3]),
                        truncation_and_arithmetic_error_not_bounded=True))
        for name in ('protocol.json','summary.json','finite_differences.json'):
            hashes[str(directory/name)] = read(directory/name)[1]
        stderr_path = Path(entry['stderr'])
        require(stderr_path.parent == LAUNCH and stderr_path.stat().st_size <= 2*1024**2,
                'Unexpected or excessive shard stderr')
        stderr_body = stderr_path.read_bytes()
        hashes[str(stderr_path)] = hashlib.sha256(stderr_body).hexdigest()
        shards.append(dict(job_id=entry['job_id'],case=case,axis=axis,records=records,
                           hessian_column_diagnostics=diagnostics,
                           stderr=stderr_body.decode('utf-8',errors='replace')))
    return dict(ready=True,execution_complete=True,completed_shards=14,required_shards=14,
        scalar_gradient_records=42,hessian_columns_complete=not refusals,
        hessian_column_refusals=refusals,shards=shards,input_hashes=hashes,
        accuracy_verified=False,no_automatic_accuracy_pass=True,
        finite_difference_truncation_and_arithmetic_unbounded=True,
        empirical_integral_allowances_not_certified=True,no_model_fit=True,no_calibration=True,
        full_gene_accuracy_and_optimizer_checks_still_required=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect',action='store_true')
    args = parser.parse_args()
    entries = metadata()
    if args.inspect:
        print(json.dumps(dict(metadata_verified=True,required_shards=14,numerical_array_reads=False,
                              output_written=False,fixed_campaign=str(LAUNCH))))
        return
    result = collect(entries)
    if not result['ready']:
        print(json.dumps(result))
        return
    output = LAUNCH/'hessian_reuse_collection.json'
    body = json.dumps(result,indent=2,allow_nan=False)+'\n'
    require(len(body.encode()) <= 4*1024**2,'Collection output cap exceeded')
    with output.open('x') as handle:
        handle.write(body)
    print(json.dumps(dict(output=str(output),execution_complete=True,accuracy_verified=False)))


if __name__ == '__main__':
    main()
