#!/usr/bin/env python3
"""Bounded profile of existing synthetic G5 conditional evaluations only."""
import argparse
import ast
import cProfile
from dataclasses import replace
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pstats
import sys
import time

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
BASE=REC/'repair_p3_gamma_G5_batch_family_check_22083952'
FIXTURE=REC/'repair_p3_five_gene_accuracy_22039932/fixture.json'


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pins',type=Path,required=True);ap.add_argument('--expected-pins-sha',required=True)
    ap.add_argument('--expected-runner-sha',required=True);ap.add_argument('--inspect',action='store_true');a=ap.parse_args()
    require(sha(Path(__file__))==a.expected_runner_sha and sha(a.pins)==a.expected_pins_sha,'Runner/pins differ')
    pins=json.loads(a.pins.read_text())
    for p,h in pins.items():
        path=Path(p);require(path.resolve().is_relative_to(ROOT.resolve()) and sha(path)==h,'Input differs')
        if path.suffix=='.py':ast.parse(path.read_text())
    require(str(FIXTURE) in pins and str(BASE/'summary.json') in pins,'Missing completed fixture evidence')
    require(json.loads((BASE/'summary.json').read_text())['selected_point_likelihood_criteria_passed'],
        'Completed selected G5 point required')
    if a.inspect:print(json.dumps(dict(rows=500,node_batches=3,nodes_per_batch=64,paired_prefix=True,no_integral=True,no_fit=True)));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
        'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
        'Single numerical thread required')
    out=REC/('profile_G5_conditional_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p,h in pins.items():
        p=Path(p)
        if p.suffix=='.py':
            q=archive/p.name;q.write_bytes(p.read_bytes());require(sha(q)==h,'Archive differs')
    def load(name,filename):
        spec=importlib.util.spec_from_file_location(name,archive/filename);m=importlib.util.module_from_spec(spec)
        sys.modules[name]=m;spec.loader.exec_module(m);return m
    import numpy as np
    original=load('profile_original','likelihood_v1.py')
    prepared=load('profile_prepared','repair_p3_gamma_batched_conditional_g5.py')
    fixture=json.loads(FIXTURE.read_text());require(fixture['rows_per_gene']==[5,25,75,200,500],'Fixture differs')
    rows=[]
    for gene,count in enumerate(fixture['rows_per_gene']):
        for j in range(count):
            rows.append(original.Row(gene,f'{j:04d}','fixture',20+j%17,10,.6,.3,.1,.01,.02,.05,
                (-1)**(gene+j),float(j%3-1),(float(j%5-2)/10,),float(j%7-3)/10,float(j%5-2)/10,10,2,10))
    require(len(fixture['drawn_counts'])==len(rows)==805,'Saved counts differ')
    # The original fixture records counts in Model's individual/gene order.
    rows=[replace(row,a=count) for row,count in zip(original.Model(rows).rows,fixture['drawn_counts'])]
    rows=tuple(r for r in rows if r.gene==4);model=original.Model(rows,gene_count=5)
    p=np.asarray(fixture['parameters']);p.setflags(write=False)
    Context=prepared.context_class(archive);context=Context(model,p,rows)
    inputs=[]
    for nodes in (np.linspace(-.3,-.1,64),np.linspace(.1,.3,64),np.linspace(-20,20,64)):
        for row in rows:
            predictors=[model.predictors(row,p,float(u),need_hessian=False) for u in nodes]
            inputs.append((row,np.asarray([x[0] for x in predictors]),predictors[0][1]))
    profiler=cProfile.Profile();began=time.perf_counter();profiler.enable()
    sums=np.zeros(3);references=[]
    for row,etas,r in inputs:
        answer=context.callback(row,etas,r,need_hessian=True)
        require(all(np.isfinite(x).all() for x in answer),'Nonfinite profile output')
        sums+=[float(np.sum(x)) for x in answer]
        references.append(answer)
    profiler.disable();seconds=time.perf_counter()-began;profiler.dump_stats(str(out/'profile.pstats'))
    stats=pstats.Stats(profiler);records=[]
    for (file,line,name),(primitive,calls,self_time,cumulative,callers) in stats.stats.items():
        records.append(dict(file=file,line=line,name=name,primitive_calls=primitive,calls=calls,
            self_seconds=self_time,cumulative_seconds=cumulative))
    records.sort(key=lambda r:r['cumulative_seconds'],reverse=True)
    prefix=load('profile_prefix','repair_p3_batched_conditional_g5.py').BatchedConditional(model,p,rows)
    maxima=np.zeros(3);prefix_profiler=cProfile.Profile();began=time.perf_counter();prefix_profiler.enable()
    actuals=[]
    for (row,etas,r),reference in zip(inputs,references):
        answer=prefix.terms_batch(row,etas,r,p,need_hessian=True)
        require(all(np.isfinite(x).all() for x in answer),'Nonfinite prefix output')
        maxima=np.maximum(maxima,[float(np.max(np.abs(x-y))) for x,y in zip(answer,reference)])
        actuals.append(answer)
    prefix_profiler.disable();prefix_seconds=time.perf_counter()-began
    prefix_profiler.dump_stats(str(out/'prefix_profile.pstats'))
    # Small synthetic arrays retained for independent arithmetic inspection.
    with (out/'outputs.npz').open('xb') as f:
        np.savez_compressed(f,**{label+'_'+str(j):np.stack([x[j] for x in values])
            for label,values in (('Gamma',references),('prefix',actuals)) for j in range(3)})
        f.flush();os.fsync(f.fileno())
    require(maxima[0]<=1e-10 and maxima[1]<=1e-10 and maxima[2]<=1e-10,'Prefix/Gamma conditional differs')
    for path,h in pins.items():require(sha(Path(path))==h,'Input changed during profile')
    report=dict(completed=True,node=os.environ.get('SLURMD_NODENAME'),seconds=seconds,
        rows=500,node_batches=3,nodes_per_batch=64,conditional_calls=len(inputs),
        output_sums=sums.tolist(),Gamma_BB_calls=context.gamma_BB_calls,
        prefix_seconds=prefix_seconds,prefix_Gamma_max_absolute_differences=maxima.tolist(),
        paired_same_job=True,profiling_ratio_not_end_to_end_speed=True,
        records=records,source_and_input_sha256=pins,
        no_integral=True,no_fit=True,no_calibration=True,no_human_data=True,
        timing_includes_profiler_overhead=True)
    with (out/'summary.json').open('x') as f:json.dump(report,f,indent=2);f.flush();os.fsync(f.fileno())
    print(json.dumps(dict(completed=True,seconds=seconds,conditional_calls=len(inputs))))


if __name__=='__main__':main()
