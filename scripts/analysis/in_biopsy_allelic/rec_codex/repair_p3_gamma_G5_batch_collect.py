#!/usr/bin/env python3
"""Collect fifteen separate G5 synthetic mass checks; no fit or calibration."""
import argparse
import ast
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
LAUNCH=REC/'repair_p3_gamma_G5_batch_family_launch_20261007T072543Z'
RUNNER_SHA='64b130f6ba6d569c97dda84b063b7e2c0b0b2db7dbe34bc28a399261733c935c'
PINS_SHA='ea630bba40d99d6aaf773a59d3b900f91fc8103fe2fff232edda4df7db507e36'
COUNTS=(5,25,75,200,500)


def require(ok,message):
    if not ok:raise ValueError(message)


def finite(x):return type(x) in (int,float) and math.isfinite(x)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    require(args.out.resolve().parent==REC.resolve(),'New direct child of REC required')
    args.out.mkdir(exist_ok=False);digests={};bytes_read=0
    def body(p):
        nonlocal bytes_read
        require(p.resolve().is_relative_to(ROOT.resolve()) and not p.is_symlink(),'Input escaped repository')
        data=p.read_bytes();bytes_read+=len(data)
        require(len(data)<=32*1024**2 and bytes_read<=256*1024**2,'Read allowance exceeded')
        digests[str(p)]=hashlib.sha256(data).hexdigest();return data
    def read(p):
        def invalid(x):raise ValueError('Nonfinite JSON constant: '+x)
        return json.loads(body(p),parse_constant=invalid)
    jobs=read(LAUNCH/'jobs.json');pins=read(LAUNCH/'pins.json')
    require(digests[str(LAUNCH/'pins.json')]==PINS_SHA,'Input pins differ')
    require([j['jobid'] for j in jobs]==[str(n) for n in range(22083938,22083953)]
        and [(j['gene'],j['tau_multiple']) for j in jobs]==[(g,m) for m in (.25,1.,4.) for g in range(5)],
        'Fixed fifteen-point roster differs')
    runner=body(LAUNCH/'repair_p3_gamma_G5_batch_family_check.py')
    require(hashlib.sha256(runner).hexdigest()==RUNNER_SHA,'Runner differs')
    tree=ast.parse(runner);main_node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
    sources=next(ast.literal_eval(n.value) for n in main_node.body if isinstance(n,ast.Assign)
        and any(isinstance(t,ast.Name) and t.id=='sources' for t in n.targets))
    require(len(sources)==12 and len(set(sources))==12,'Helper family differs')
    for path,h in pins.items():require(hashlib.sha256(body(Path(path))).hexdigest()==h,'Pinned input changed')
    collection=read(REC/'repair_p3_five_gene_accuracy_launch_20261005T143525Z/five_gene_accuracy_collection.json')
    require(collection['execution_complete'],'Independent scalar references incomplete')
    refs={(x['gene'],x['tau_multiple']):x['record'] for x in collection['results']}
    fixture=read(REC/('repair_p3_five_gene_accuracy_'+collection['results'][0]['job_id'])/'fixture.json')
    count_sha=hashlib.sha256(struct.pack('<805q',*fixture['drawn_counts'])).hexdigest()
    settings_class=next(n for n in ast.parse(body(HERE/'repair_p3_integral.py')).body
        if isinstance(n,ast.ClassDef) and n.name=='Settings')
    settings={n.target.id:ast.literal_eval(n.value) for n in settings_class.body if isinstance(n,ast.AnnAssign)}
    status=subprocess.check_output(['sacct','-j',','.join(j['jobid'] for j in jobs),'-n','-P',
        '--format=JobIDRaw,State,ExitCode,Elapsed,MaxRSS'],text=True)
    (args.out/'sacct.txt').write_text(status)
    states={p[0]:p[1:] for line in status.splitlines() if (p:=line.split('|'))[0].isdigit()}
    records=[];issues=[]
    for job in jobs:
        gene=job['gene'];multiple=job['tau_multiple'];state=states.get(job['jobid'],['MISSING'])
        record=dict(jobid=job['jobid'],gene=gene,tau_multiple=multiple,state=state,completed=False);records.append(record)
        if state[:2]!=['COMPLETED','0:0']:
            if state[0] not in ('RUNNING','PENDING','CONFIGURING','COMPLETING','SUSPENDED'):
                issues.append(dict(jobid=job['jobid'],state=state))
            continue
        d=Path(job['output_directory']);require(d==REC/('repair_p3_gamma_G5_batch_family_check_'+job['jobid']),
            'Output identity differs')
        protocol=read(d/'protocol.json');summary=read(d/'summary.json')
        answer=read(d/f'gene{gene}_tau{.05*multiple}.json');raw=read(d/'GL_returned_raw.json')
        require(protocol['source_and_input_sha256']==pins and protocol['seed']==20260923
            and protocol['requested_gene']==gene and protocol['requested_tau_multiple']==multiple
            and protocol['rows_per_gene']==list(COUNTS) and protocol['full_drawn_count_vector_sha256']==count_sha
            and protocol['executed_GL_settings']==settings and protocol['separate_G5_fixture_context'] is True
            and protocol['requested_node_batch_size']==64,'Executed fixture/settings differ')
        for name in (*sources,'repair_p3_gamma_G5_batch_family_check.py'):
            require(hashlib.sha256(body(d/'executed_sources'/name)).hexdigest()==
                (RUNNER_SHA if name=='repair_p3_gamma_G5_batch_family_check.py' else pins[str(HERE/name)]),
                'Archived source differs')
        require(summary['completed'] is True and summary['points']==1 and summary['required_family_points']==15
            and summary['requested_gene']==gene and summary['requested_tau_multiple']==multiple,'Summary scope differs')
        require(answer['gene']==gene and answer['tau']==.05*multiple and answer['variance']==(.05*multiple)**2
            and answer['scalar']==refs[gene,multiple]['scalar'],'Point/reference differs')
        native=read(d/'native_GH.json');require(native=={k:v for k,v in answer.items() if k!='GL'},'Native output changed')
        scalar=answer['scalar']['expanded']['logmass'];require(finite(scalar),'Invalid scalar reference')
        require([o['order'] for o in answer['orders']]==[20,40],'Native order family differs')
        native_ok=True
        for o in answer['orders']:
            if not o['returned']:native_ok=False;continue
            delta=o['logmass']-scalar;ok=abs(delta)<=1e-8*abs(scalar)
            require(finite(o['logmass']) and len(o['gradient'])==17 and all(finite(x) for x in o['gradient'])
                and delta==o['scalar_difference'] and ok==o['scalar_relative_LL_pass'],'Native scalar comparison differs')
            native_ok=native_ok and ok
        if all(o['returned'] for o in answer['orders']):
            delta=answer['orders'][0]['logmass']-answer['orders'][1]['logmass'];ok=abs(delta)<=1e-8
            require(delta==answer['GH20_GH40_difference'] and ok==answer['GH20_GH40_absolute_pass'],'Native refinement differs')
            native_ok=native_ok and ok
        gl=answer['GL'];require(gl['returned'] is True and summary['GL_returned'] is True,'GL did not return')
        require(raw=={k:gl[k] for k in ('logmass','gradient','hessian','report')},'Raw GL output differs')
        g=gl['gradient'];h=gl['hessian'];r=gl['report'];width=9
        require(finite(gl['logmass']) and len(g)==17 and all(finite(x) for x in g)
            and len(h)==17 and all(len(row)==17 and all(finite(x) for x in row) for row in h),'GL arrays invalid')
        active=[gene,5+gene,10,11,12,13,14,15,16]
        require(r['augmented_global_indices']==active and r['full_parameter_dimension']==16
            and r['local_augmented_dimension']==width and r['requested_node_batch_size']==64
            and r['estimated_GL_workspace_bytes']<=256*1024**2,'GL scope differs')
        require(r['refined']['logmass']==gl['logmass'] and r['refined']['gradient']==[g[i] for i in active]
            and r['refined']['hessian']==[[h[i][j] for j in active] for i in active],'Local/full embedding differs')
        inactive=[i for i in range(17) if i not in active]
        require(all(g[i]==0 and all(x==0 for x in h[i]) and all(row[i]==0 for row in h) for i in inactive),
            'Inactive derivatives differ from exact zero')
        e=r['empirical_refinement_plus_tail_allowances'];ref=r['refined']
        first=r['first']
        require(r['logmass_refinement_difference']==abs(first['logmass']-ref['logmass'])
            and r['gradient_refinement_difference']==[abs(a-b) for a,b in zip(first['gradient'],ref['gradient'])]
            and r['hessian_refinement_difference']==[[abs(a-b) for a,b in zip(ra,rb)]
                for ra,rb in zip(first['hessian'],ref['hessian'])],
            'First/final integration refinement differs from raw arrays')
        require(len(e)==3 and len(e[1])==width and len(e[2])==width and all(len(row)==width for row in e[2]),
            'Empirical allowance dimensions differ')
        rederived=[r['logmass_refinement_difference']+ref['logmass_error_allowance'],
            [a+b for a,b in zip(r['gradient_refinement_difference'],ref['gradient_error_allowance'])],
            [[a+b for a,b in zip(ra,rb)] for ra,rb in zip(r['hessian_refinement_difference'],ref['hessian_error_allowance'])]]
        require(e==rederived,'Empirical refinement arithmetic differs')
        flat=[e[0],*e[1],*(x for row in e[2] for x in row)]
        require(all(finite(x) and x>=0 for x in flat) and e[0]<=settings['logmass_tolerance']
            and max(e[1])<=settings['gradient_tolerance']
            and max(x for row in e[2] for x in row)<=settings['hessian_tolerance'],'Unresolved empirical allowance')
        ledger=gl['conditional_chunk_ledger'];require(ledger and {x['row_ordinal'] for x in ledger}==set(range(COUNTS[gene])),
            'Conditional row coverage differs')
        for x in ledger:
            require(type(x['calls']) is int and x['calls']>0 and x['need_hessian'] is True
                and x['all_requested_nodes_returned'] is True and 0<x['requested_nodes']<=64
                and x['effective_chunk_limit']==x['requested_nodes']
                and x['chunk_lengths']==[x['requested_nodes']],'Shallow-count chunk coverage differs')
        nodes=sum(x['requested_nodes']*x['calls'] for x in ledger)
        require(nodes==COUNTS[gene]*r['likelihood_evaluations'],'Conditional node count differs')
        require(all(sum(x['requested_nodes']*x['calls'] for x in ledger if x['row_ordinal']==i)
            ==r['likelihood_evaluations'] for i in range(COUNTS[gene])),
            'Per-row conditional node coverage differs')
        batches=sum(x['calls']*len(x['chunk_lengths']) for x in ledger)
        require(type(gl['Gamma_BB_batch_calls']) is int and 2*batches<=gl['Gamma_BB_batch_calls']<=3*batches,
            'Gamma batch counter outside support call bounds')
        probe=read(d/'conditional_probe.json')
        require(probe['rows']==COUNTS[gene] and probe['comparisons']==6*COUNTS[gene]
            and probe['maxima']['logmass']<=1e-11 and probe['maxima']['gradient']<=1e-11
            and probe['maxima']['hessian']<=1e-10,'Conditional component comparison differs')
        dl=gl['logmass']-scalar;gl_ok=abs(dl)<=1e-8*abs(scalar);passed=native_ok and gl_ok
        require(dl==gl['scalar_difference'] and gl_ok==gl['scalar_relative_LL_pass']
            and passed==summary['selected_point_likelihood_criteria_passed'],'Likelihood criterion differs')
        record.update(completed=True,point_likelihood_criteria_passed=passed,GL_scalar_difference=dl,
            GL_seconds=gl['seconds'],quadrature_nodes=r['likelihood_evaluations'],
            selected_orders=[r['selected_first_order'],r['selected_refined_order']],
            conditional_nodes=nodes,Gamma_BB_batch_calls=gl['Gamma_BB_batch_calls'],
            Gamma_support_calls_per_batch=gl['Gamma_BB_batch_calls']/batches,
            empirical_allowance_maxima=[e[0],max(e[1]),max(x for row in e[2] for x in row)],
            conditional_probe_maxima=probe['maxima'])
    complete=all(x['completed'] for x in records)
    report=dict(complete_family=complete,completed_points=sum(x['completed'] for x in records),required_points=15,
        all15_likelihood_criteria_passed=all(x['point_likelihood_criteria_passed'] for x in records) if complete else None,
        records=records,issues=issues,scheduler_states=dict(Counter(x['state'][0] for x in records)),
        collector_sha256=hashlib.sha256(body(Path(__file__).resolve())).hexdigest(),input_sha256=digests,
        error_allowances_not_certified=True,no_independent_full_derivative_validation=True,
        conditional_comparisons_not_all_integration_nodes=True,no_chunk_splitting_test=True,
        no_fits=True,no_calibration=True,no_human_data=True)
    (args.out/'collection.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(completed=report['completed_points'],required=15,complete=complete,issues=len(issues),
        all15_likelihood_criteria_passed=report['all15_likelihood_criteria_passed'])))
    if issues or (complete and not report['all15_likelihood_criteria_passed']):raise SystemExit(1)


if __name__=='__main__':main()
