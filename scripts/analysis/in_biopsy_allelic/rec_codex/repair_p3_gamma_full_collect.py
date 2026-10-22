#!/usr/bin/env python3
"""Read completed fixed-point jobs; never infer fit or independent accuracy."""
import argparse
import ast
import hashlib
import json
import math
from pathlib import Path
import subprocess


def require(ok,message):
    if not ok:
        raise RuntimeError(message)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def flat(x):
    if isinstance(x,list):
        for v in x:
            yield from flat(v)
    else:
        yield x


def difference(a,b):
    x,y=list(flat(a)),list(flat(b))
    require(len(x)==len(y) and all(math.isfinite(v) for v in x+y),'Invalid finite numerical comparison')
    return max((abs(v-w) for v,w in zip(x,y)),default=0.)


def verify_settings(pins,checker):
    """Derive criteria from the pinned executable expressions, not prose."""
    paths={Path(p).name:Path(p) for p in pins
           if Path(p).parent.name=='rec_codex'}
    kernel=ast.parse(paths['repair_p3_integral.py'].read_text())
    settings=next(n for n in kernel.body if isinstance(n,ast.ClassDef) and n.name=='Settings')
    defaults={n.target.id:ast.literal_eval(n.value) for n in settings.body
              if isinstance(n,ast.AnnAssign)}
    expected_wrapper='replace(settings,logmass_tolerance=min(settings.logmass_tolerance,1e-8/max(1,present)),gradient_tolerance=min(settings.gradient_tolerance,1e-7/max(1,present)),hessian_tolerance=min(settings.hessian_tolerance,1e-5/max(1,present)))'
    expected_serial='replace(kernel.Settings(),logmass_tolerance=1e-9,gradient_tolerance=1e-8,hessian_tolerance=1e-5/8)'
    for source,expression in ((paths['repair_p3_gamma_parallel_model.py'],expected_wrapper),
                              (checker,expected_serial)):
        actual=[n for n in ast.walk(ast.parse(source.read_text())) if isinstance(n,ast.Call)
                and isinstance(n.func,ast.Name) and n.func.id=='replace']
        require(len(actual)==1 and ast.dump(actual[0])==ast.dump(ast.parse(expression,mode='eval').body),
                'Executable numerical settings expression changed')
    parallel=dict(defaults)
    parallel.update(logmass_tolerance=min(defaults['logmass_tolerance'],1e-8/8),
        gradient_tolerance=min(defaults['gradient_tolerance'],1e-7/8),
        hessian_tolerance=min(defaults['hessian_tolerance'],1e-5/8))
    serial=dict(defaults)
    serial.update(logmass_tolerance=1e-9,gradient_tolerance=1e-8,hessian_tolerance=1e-5/8)
    require(parallel==serial,'Executable sequential/parallel criteria differ')
    return parallel


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--launch',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    jobs=json.loads((args.launch/'jobs.json').read_text())
    require(len(jobs)==2 and {j['mode'] for j in jobs}=={'parallel','sequential'}
            and len({j['job_id'] for j in jobs})==2,'Expected two independent jobs')
    checker=args.launch/'repair_p3_gamma_full_point.py'
    manifest=args.launch/'pins.json'
    require(all(sha(checker)==j['checker_sha256'] and sha(manifest)==j['manifest_sha256']
                for j in jobs),'Submitted executable/manifest changed')
    pins=json.loads(manifest.read_text())
    for path,digest in pins.items():
        require(sha(Path(path))==digest,'Pinned settings/input source changed')
    settings=verify_settings(pins,checker)
    receipt=subprocess.run(['sacct','-j',','.join(j['job_id'] for j in jobs),
        '--format=JobIDRaw,State,ExitCode','-n','-P'],text=True,capture_output=True,check=True).stdout
    states={}
    for line in receipt.splitlines():
        fields=line.split('|')
        if len(fields)>=3 and fields[0].isdigit():
            states[fields[0]]=fields[1:3]
    args.out.mkdir(exist_ok=False)
    (args.out/'sacct.txt').write_text(receipt)
    result=dict(complete_pair=False,states=states,no_fit=True,no_calibration=True,
                no_independent_full_Hessian=True,not_a_certified_numerical_error_bound=True,
                executable_settings_AST_checked=True,per_gene_settings=settings)
    def save():
        (args.out/'collection.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print(json.dumps(result))
    if any(states.get(j['job_id'])!=['COMPLETED','0:0'] for j in jobs):
        result['waiting_or_failed_job_ids']=[j['job_id'] for j in jobs
            if states.get(j['job_id'])!=['COMPLETED','0:0']]
        save();return
    protocols={};answers={};records={};inputs={}
    for j in jobs:
        mode=j['mode'];out=Path(j['output_directory'])
        require(out.name=='repair_p3_gamma_full_'+mode+'_'+j['job_id'],
                'Output directory is not tied to the declared job')
        require(sha(args.launch/'repair_p3_gamma_full_point.py')==j['checker_sha256']
                and sha(args.launch/'pins.json')==j['manifest_sha256'],'Submitted executable/manifest changed')
        def read(p):
            inputs[str(p)]=sha(p)
            return json.loads(p.read_text())
        protocol=read(out/'protocol.json');summary=read(out/'summary.json');answer=read(out/'answer.json')
        require(protocol['mode']==mode and summary['mode']==mode and summary['completed'] is True
                and summary['all_fixed_inputs_unchanged'] is True
                and summary['returned']==answer['returned'],'Incomplete actual result')
        require(sha(out/'executed_sources/repair_p3_gamma_full_point.py')==j['checker_sha256'],
                'Executed checker differs')
        for path,digest in protocol['source_and_input_sha256'].items():
            require(sha(Path(path))==digest,'Source/input changed: '+path)
            archived=out/'executed_sources'/Path(path).name
            if Path(path).suffix=='.py' and archived.exists():
                require(sha(archived)==digest,'Archived source differs')
        for name in ('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
                     'repair_p3_integral.py','repair_p3_composite_gl.py',
                     'repair_p3_bb_gamma_difference_v3.py','repair_p3_gamma_parallel_model.py'):
            require((out/'executed_sources'/name).is_file(),'Required executed source absent: '+name)
        require(not answer.get('integrity_failure',False),'Recorded evaluator integrity failure')
        if mode=='parallel':
            rows=read(out/'workers.json')
        else:
            rows=[read(out/('gene'+str(g)+'.json')) for g in range(8)]
        require(len(rows)==8 and [r['gene'] for r in rows]==list(range(8)),
                'Expected every gene in ascending order')
        protocols[mode]=protocol;answers[mode]=answer;records[mode]=rows
    a,b=protocols['sequential'],protocols['parallel']
    require({k:v for k,v in a.items() if k!='mode'}=={k:v for k,v in b.items() if k!='mode'},
            'Different fixed point, numerical criteria or sources')
    result.update(complete_pair=True,input_result_sha256=inputs,
        both_returned=all(x['returned'] for x in answers.values()),
        per_gene_dispositions=[dict(gene=g,sequential=records['sequential'][g]['success'],
                                   parallel=records['parallel'][g]['success']) for g in range(8)])
    if not result['both_returned']:
        for r in records['parallel']:
            if not r['success']:
                require(r.get('numerical_integration_refusal') is True
                        and r.get('value_error') is True
                        and isinstance(r.get('refused_gene_ladders'),list)
                        and isinstance(r.get('refused_gene_evaluations'),int),
                        'Failed worker is not a documented numerical integration refusal')
        result.update(nominal_comparison_passed=False,all_refusals_retained=True)
        save();return
    comparisons=[]
    for g,(x,y) in enumerate(zip(records['sequential'],records['parallel'])):
        require(x['success'] and y['success'],'Successful total with failed gene')
        require(all(len(r['gradient'])==29 and len(r['hessian'])==29
                    and all(len(row)==29 for row in r['hessian']) for r in (x,y)),
                'Expected 29 scores and a 29 by 29 Hessian')
        dg=difference(x['gradient'],y['gradient']);dh=difference(x['hessian'],y['hessian'])
        dl=abs(x['logmass']-y['logmass'])
        inactive=[j for j in range(16) if j not in (g,8+g)]
        zeros=all(r['gradient'][j]==0 and all(v==0 for v in r['hessian'][j])
                  and all(row[j]==0 for row in r['hessian']) for r in (x,y) for j in inactive)
        require(zeros,'Cross-gene derivative layout/integrity defect')
        comparisons.append(dict(gene=g,LL_abs=dl,score_max_abs=dg,Hessian_max_abs=dh,
            reports_exact=x['report']==y['report'],cross_gene_alpha_r_exact_zero=zeros,
            passed=dl<=a['comparison_LL_atol'] and dg<=a['comparison_score_atol']
                and dh<=a['comparison_Hessian_atol'] and zeros))
    x,y=answers['sequential'],answers['parallel']
    for mode,answer in answers.items():
        rows=records[mode]
        require(len(answer['gradient'])==29 and len(answer['hessian'])==29
                and all(len(row)==29 for row in answer['hessian']),'Invalid total dimensions')
        require(answer['logmass']==sum(r['logmass'] for r in rows)
                and answer['gradient']==[sum(r['gradient'][j] for r in rows) for j in range(29)]
                and answer['hessian']==[[sum(r['hessian'][i][j] for r in rows) for j in range(29)]
                                       for i in range(29)],'Reported total differs from gene sums')
        require(answer['gene_terms']=={str(r['gene']):r['logmass'] for r in rows},
                'Gene contributions differ from retained records')
    total=dict(LL_abs=abs(x['logmass']-y['logmass']),score_max_abs=difference(x['gradient'],y['gradient']),
               Hessian_max_abs=difference(x['hessian'],y['hessian']))
    for answer in (x,y):
        for key,limit in (('logmass',a['global_LL_allowance']),
                          ('gradient',a['global_score_allowance']),
                          ('hessian',a['global_Hessian_allowance'])):
            values=list(flat(answer['error_allowances'][key]))
            require(all(math.isfinite(v) and 0<=v<=limit for v in values),'Accumulated allowance exceeded')
    result.update(per_gene_comparisons=comparisons,total_differences=total,
        nominal_comparison_passed=all(c['passed'] for c in comparisons)
            and total['LL_abs']<=a['comparison_LL_atol'] and total['score_max_abs']<=a['comparison_score_atol']
            and total['Hessian_max_abs']<=a['comparison_Hessian_atol'])
    save()


if __name__=='__main__':
    main()
