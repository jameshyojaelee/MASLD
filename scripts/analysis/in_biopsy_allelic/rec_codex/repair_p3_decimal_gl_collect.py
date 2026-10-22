#!/usr/bin/env python3
"""Collect every declared Decimal GL chunk before forming mixed curvature."""
import argparse
from decimal import Decimal,localcontext
import hashlib
import json
from pathlib import Path
import subprocess
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
RULES = ((24,'0.25',16),(26,'0.125',16),(26,'0.125',32))
CASES = {'failed0gene2':22039200,'positive18':22039201}


def require(ok,message):
    if not ok:
        raise ValueError(message)


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--launch',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    require(args.launch.resolve().is_relative_to(REC) and args.out.resolve().is_relative_to(REC),
            'Collection must remain in the existing synthetic execution tree')
    args.out.mkdir(exist_ok=False)
    try:
        jobs = read(args.launch/'jobs.json')
        expected = {(case,rule,chunk) for case in CASES for rule in range(3) for chunk in range(16)}
        require(len(jobs) == 96 and {(j['case'],j['rule'],j['chunk']) for j in jobs} == expected
                and len({j['job_id'] for j in jobs}) == 96,'Declared 96-job family differs')
        pins = read(args.launch/'pins.json')
        for path,digest in pins.items():
            require(sha(Path(path)) == digest,'Source/input differs: '+path)
        receipt = read(args.launch/'source_receipt.json')
        for path,digest in receipt.items():
            require(sha(Path(path)) == digest,'Frozen source differs: '+path)
        status = subprocess.run(['sacct','-j',','.join(str(j['job_id']) for j in jobs),
            '--noheader','--parsable2','--format=JobIDRaw,State,ExitCode'],
            capture_output=True,text=True,check=True).stdout
        (args.out/'sacct.txt').write_text(status)
        states = {}
        for line in status.splitlines():
            parts = line.split('|')
            if len(parts) >= 3 and parts[0].isdigit():
                states[parts[0]] = parts[1:3]
        missing,refusals,loaded = [],[],{}
        evidence = {}
        for job in jobs:
            key = (job['case'],job['rule'],job['chunk'])
            out = Path(job['output_directory'])
            require(out.resolve().is_relative_to(REC),'Chunk output escaped execution tree')
            state = states.get(str(job['job_id']),['MISSING','MISSING'])
            if state != ['COMPLETED','0:0']:
                (refusals if state[0] in ('FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY') else missing).append(
                    dict(job=job,state=state))
                continue
            require((out/'summary.json').is_file(),'Completed chunk lacks summary')
            summary,protocol = read(out/'summary.json'),read(out/'protocol.json')
            require(summary['completed'] is True and summary['records'] == 2
                    and summary['raw_moments_only'] is True and summary['no_chunk_covariance'] is True,
                    'Chunk summary does not prove complete raw moments')
            require(protocol['case'] == job['case'] and protocol['rule_index'] == job['rule']
                    and protocol['chunk'] == job['chunk'] and protocol['decimal_digits'] == [80,110]
                    and protocol['fixed_log_rescaling'] == '0'
                    and protocol['source_and_input_sha256'] == pins,'Chunk protocol differs')
            parent = read(REC/f"repair_p3_fixed_point_{job['case']}_{CASES[job['case']]}/protocol.json")
            require(all(protocol[k] == parent[k] for k in ('fixed_p','fixed_v','count_vector_sha256')),
                    'Chunk changed original synthetic point/counts')
            records = []
            for digits in (80,110):
                path = out/f'moments_{digits}.json'
                value = read(path)
                bound,step,order = RULES[job['rule']]
                count = int(Decimal(2*bound)/Decimal(step))
                first,last = count*job['chunk']//16,count*(job['chunk']+1)//16
                require((value['digits'],value['rule_index'],value['chunk'],value['bound'],
                         value['spacing'],value['order'],value['total_panels'],value['first_panel'],
                         value['last_panel_exclusive'],value['node_calls']) ==
                        (digits,job['rule'],job['chunk'],bound,step,order,count,first,last,(last-first)*order),
                        'Chunk panel coverage or node count differs')
                require(len(value['moments']) == 6 and len(value['GL_polynomial_absolute_errors']) == 2*order
                        and all(Decimal(x).is_finite() for x in value['moments'])
                        and Decimal(value['moments'][0]) > 0
                        and max(map(Decimal,value['GL_polynomial_absolute_errors'])) <=
                            Decimal(10)**(-digits+13),'Invalid moments or polynomial rule')
                records.append(value)
                evidence[str(path)] = sha(path)
            loaded[key] = records
        report = dict(complete_family=len(loaded) == 96,completed_chunks=len(loaded),
            total_chunks=96,missing_or_running=missing,refusals=refusals,
            input_result_sha256=evidence,no_full_Hessian_accuracy=True,
            no_certified_quadrature_error=True,no_fit=True,no_calibration=True)
        if len(loaded) == 96:
            families = {}
            with localcontext() as context:
                context.prec = 110
                for case in CASES:
                    candidate_path = REC/f'repair_p3_fixed_point_{case}_{CASES[case]}/candidate_baseline.json'
                    require(str(candidate_path) in pins,'Candidate comparison not pinned')
                    candidate = read(candidate_path)
                    values = []
                    for rule in range(3):
                        for precision_index,digits in enumerate((80,110)):
                            chunks = [loaded[(case,rule,k)][precision_index] for k in range(16)]
                            coverage = [j for chunk in chunks for j in range(
                                chunk['first_panel'],chunk['last_panel_exclusive'])]
                            require(coverage == list(range(chunks[0]['total_panels'])),
                                    'Whole-domain panels omitted or duplicated')
                            raw = [sum((Decimal(chunk['moments'][j]) for chunk in chunks),Decimal(0))
                                   for j in range(6)]
                            curvature = [raw[4]/raw[0]-raw[1]*raw[2]/raw[0]**2,
                                         raw[5]/raw[0]-raw[1]*raw[3]/raw[0]**2]
                            differences = [Decimal.from_float(float(candidate['hessian'][10][j]))-x
                                           for j,x in zip((24,28),curvature)]
                            tail = chunks[0]['analytic_absolute_tail_bounds_not_interval_evaluated']
                            require(all(c['analytic_absolute_tail_bounds_not_interval_evaluated'] == tail
                                        for c in chunks),'Whole-line tail bounds changed between chunks')
                            values.append(dict(rule=rule,digits=digits,raw_moments=list(map(str,raw)),
                                mixed_curvature_r_omega_and_r_variance=list(map(str,curvature)),
                                candidate_minus_reference=list(map(str,differences)),
                                candidate_within_nominal_1e_minus5=all(abs(x) <= Decimal('1e-5') for x in differences),
                                raw_absolute_tail_bounds_not_interval_evaluated=tail))
                    stability = [[str(abs(Decimal(a)-Decimal(b))) for a,b in zip(
                        values[2*r]['mixed_curvature_r_omega_and_r_variance'],
                        values[2*r+1]['mixed_curvature_r_omega_and_r_variance'])] for r in range(3)]
                    refinement = [[str(abs(Decimal(a)-Decimal(b))) for a,b in zip(
                        values[2*r+1]['mixed_curvature_r_omega_and_r_variance'],
                        values[2*(r+1)+1]['mixed_curvature_r_omega_and_r_variance'])] for r in range(2)]
                    families[case] = dict(all_rule_and_precision_results=values,
                        precision_differences=stability,refinement_differences=refinement,
                        precision_within_1e_minus30=all(Decimal(x) <= Decimal('1e-30') for pair in stability for x in pair),
                        refinement_within_1e_minus8=all(Decimal(x) <= Decimal('1e-8') for pair in refinement for x in pair),
                        candidate_within_nominal_1e_minus5_for_entire_family=all(
                            value['candidate_within_nominal_1e_minus5'] for value in values))
            report['families'] = families
        for path,digest in pins.items():
            require(sha(Path(path)) == digest,'Source/input changed during collection')
        with (args.out/'collection.json').open('x') as handle:
            json.dump(report,handle,indent=2,allow_nan=False)
            handle.write('\n')
        print(json.dumps(dict(complete_family=report['complete_family'],completed_chunks=len(loaded),
                             missing_or_running=len(missing),refusals=len(refusals))))
    except BaseException:
        (args.out/'failure.txt').write_text(traceback.format_exc())
        raise


if __name__ == '__main__':
    main()
