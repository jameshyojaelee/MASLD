#!/usr/bin/env python3
"""Run the nine saved-rule reuse comparisons in eight compute subprocesses."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

REC = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z')
POINT = REC/'repair_p3_observed_failure_anchor_22042612'
NODE_COUNTS = (1024,2048,4096,1920,3840,7680,4352,8704,17408)


def require(ok,message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checker',type=Path,required=True)
    parser.add_argument('--checker-sha',required=True)
    parser.add_argument('--pins',type=Path,required=True)
    parser.add_argument('--pins-sha',required=True)
    parser.add_argument('--self-sha',required=True)
    parser.add_argument('--inspect',action='store_true')
    args = parser.parse_args()
    require(sha(Path(__file__)) == args.self_sha and sha(args.checker) == args.checker_sha
            and sha(args.pins) == args.pins_sha,'Reviewed orchestration source differs')
    pins = json.loads(args.pins.read_text())
    for path,digest in pins.items():
        require(sha(Path(path)) == digest,'Reviewed source/input differs: '+path)
    command = [sys.executable,str(args.checker),'--pins',str(args.pins),
               '--pins-sha',args.pins_sha,'--self-sha',args.checker_sha]
    if args.inspect:
        for rule in range(9):
            subprocess.run(command+['--rule',str(rule),'--inspect'],check=True,capture_output=True,text=True)
        print(json.dumps(dict(static_checks=True,rule_nodes=NODE_COUNTS,workers=8,no_fit=True)))
        return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID')
            and int(os.environ.get('SLURM_CPUS_PER_TASK','0')) >= 8 and len(os.sched_getaffinity(0)) >= 8,
            'Individual eight-CPU allocation required')
    require(all(os.environ.get(key) == '1' for key in
                ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
            'Each numerical subprocess must use one numerical thread')
    require(int(os.environ.get('SLURM_MEM_PER_NODE','0')) >= 8192,'Eight GiB allocation required')
    output = REC/('repair_p3_cache_rule_family_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(exist_ok=False)
    began = time.perf_counter()
    def save(name,value):
        body = json.dumps(value,indent=2,allow_nan=False)+'\n'
        require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file())+len(body.encode())
                < 32*1024**2,'Combined nine-rule job output cap exceeded')
        with (output/(name+'.json')).open('x') as handle:
            handle.write(body)
            handle.flush()
            os.fsync(handle.fileno())
    save('protocol',dict(workers=8,rule_nodes=NODE_COUNTS,total_nodes=sum(NODE_COUNTS),
        source_sha256={str(Path(__file__)):args.self_sha,str(args.checker):args.checker_sha},
        source_and_input_sha256=pins,no_fit=True,no_calibration=True,no_accuracy_acceptance=True,
        original_GL_callable_not_reexecuted=True,complete_original_saved_rule_family=True))
    def run(rule):
        with (output/f'rule_{rule}.out').open('x') as stdout, (output/f'rule_{rule}.err').open('x') as stderr:
            result = subprocess.run(command+['--rule',str(rule)],stdout=stdout,stderr=stderr)
        return dict(rule=rule,exit_code=result.returncode,elapsed_seconds=time.perf_counter()-began)
    results = []
    try:
        # Longest comparisons start first; each subprocess constructs its own
        # process/thread-owned one-point cache and archives its own sources.
        with ThreadPoolExecutor(max_workers=8) as pool:
            futures = [pool.submit(run,rule) for rule in sorted(range(9),key=lambda r:-NODE_COUNTS[r])]
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                save('terminal_rule_'+str(result['rule']),result)
                print(json.dumps(dict(completed_rules=len(results),total_rules=9,**result)),flush=True)
        baseline = json.loads((POINT/'GL_candidate.json').read_text())
        records = []
        for result in sorted(results,key=lambda r:r['rule']):
            rule = result['rule']
            path = output/f'rule_{rule}/summary.json'
            record = dict(result,summary_available=path.is_file())
            if path.is_file():
                summary = json.loads(path.read_text())
                passed = (summary['completed'] is True and summary['rule_index'] == rule
                    and summary['nodes'] == NODE_COUNTS[rule]
                    and summary['conditional_pairs'] == 24*NODE_COUNTS[rule]
                    and all(summary[k] == 0 for k in ('predictor_mismatches','conditional_mismatches',
                                                      'aggregate_node_mismatches'))
                    and summary['rule_logmass_byte_equal'] is True
                    and summary['rule_gradient_byte_equal'] is True)
                record.update(summary_sha256=sha(path),same_process_byte_checks_passed=passed)
                scalar = summary['results'][0]
                old = baseline['integration']['ladder'][rule]
                difference = [a-b for a,b in zip(scalar['gradient'],old['gradient'])]
                record.update(cross_job_logmass_difference=scalar['logmass']-old['logmass'],
                    cross_job_gradient_difference=difference,
                    cross_job_gradient_maximum_absolute_difference=max(map(abs,difference)))
            records.append(record)
        for path,digest in pins.items():
            require(sha(Path(path)) == digest,'Reviewed source/input changed during comparison')
        require(sha(Path(__file__)) == args.self_sha and sha(args.checker) == args.checker_sha,
                'Reviewed orchestration source changed')
        passed = len(records) == 9 and all(r['exit_code'] == 0 and
                    r.get('same_process_byte_checks_passed') is True for r in records)
        save('summary',dict(completed=True,records=records,all_nine_rules_same_process_byte_equal=passed,
            total_nodes=sum(NODE_COUNTS),conditional_pairs=24*sum(NODE_COUNTS),
            elapsed_seconds=time.perf_counter()-began,no_fit=True,no_calibration=True,
            no_accuracy_acceptance=True,original_GL_callable_not_reexecuted=True,
            old_cross_job_failure_preserved=True))
        require(passed,'Complete nine-rule scalar/reuse comparison did not pass')
    except BaseException:
        save('failure',dict(traceback=traceback.format_exc(),partials_preserved=True,
            completed_subprocesses=results,no_fit=True,no_calibration=True))
        raise


if __name__ == '__main__':
    main()
