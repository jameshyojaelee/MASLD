"""Paired P2 execution check on synthetic participants, counts and error draws.

Two draws exercise the driver only. This is not a size, power, or coverage
estimate and supplies no Gate A information about MASLD cohorts.
"""
import json
import os
import platform
import time
from pathlib import Path

import numpy as np
import scipy

from likelihood_v1 import Row,conditional_draw
from design_v1 import Participant,Design
from lambda_v1 import PairedReliability
from tests_v1 import p2,p2_tau


seed=20260923
rng=np.random.default_rng(seed)
people=[Participant(f'P{i:03d}','synthetic',i%3,
                    (float(np.log(.01)+float(i%5-2)/10),.3+float(i%7-3)/10),
                    float(.01*np.exp(float(i%5-2)/10)),35.,4)
        for i in range(36)]
rows=[Row(g,p.individual,p.cohort,30+i%11,15,.6,.3,.1,p.e,.02,.05,(-1)**(g+i),
          0.,(),0.,0.,10,2,10) for g in range(4) for i,p in enumerate(people)]
clusters={people[0].individual:'fixture_family',people[1].individual:'fixture_family'}
design=Design(people,rows,clusters)
model=design.model();planted=model.initial();planted[:4]=[-1.2,-.6,.7,1.3];planted[8]=.08
drawn=[]
for row in model.rows:
    eta,r,_,_=model.predictors(row,planted)
    drawn.append(conditional_draw(row,eta,r,*rng.random(2)))
design=Design(people,drawn,clusters)
y5_cycle=(0,2,3,1,2,4)
target_y5={p.individual:y5_cycle[i%6] for i,p in enumerate(people)}
assert all(int(np.searchsorted((1.5,2.5),target_y5[p.individual]))==p.stage for p in people)
reliability=PairedReliability({'fixture_source_1':'S3','fixture_source_2':'Y5'},
                            {'fixture_source_1':[.4,.45],'fixture_source_2':[.2,.25]},
                            target_y5,{'synthetic':1.},observed_r_squared=design.stage_r_squared())
traces=[]
def paired(sampled,b):
    values=reliability(sampled,b)
    traces.append({'b':b,'individuals':[p.individual for p in sampled.participants],
                   'clusters':len(set(sampled.clusters.values())),
                   'stage_mean':float(sampled.stage.mean()),
                   'centered_stage_sum':float(sampled.stage_centered.sum()),
                   'resampled_r_squared_sensitivity':sampled.stage_r_squared(),
                   'fixed_observed_r_squared':reliability.observed_r_squared,'lambda':values.tolist()})
    assert abs(sampled.stage_centered.sum())<1e-10
    assert np.isfinite(values).all() and np.all(values>.05),values
    return values
paired.source_specs=reliability.source_specs
start=time.monotonic()
output=Path(os.environ['CODEX_REC_OUTPUT'])/('p2_v13_'+os.environ['SLURM_JOB_ID'])
output.mkdir(exist_ok=False)
events=[]
def serial(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)
def progress(event,payload):
    name=f'{len(events):03d}_{event}.json'
    (output/name).write_text(json.dumps(payload,default=serial,indent=2)+'\n')
    events.append(name)
    print(event,'saved',name,flush=True)
result=p2(design,paired,draws=2,progress=progress)
out={'spec':'v1.3','seed':seed,'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
     'participants':36,'relative_clusters':35,'genes':4,'rows':144,
     'draw_count':2,'label_errors':'synthetic fixed draws; not estimates from public pairs',
     'elapsed_seconds':time.monotonic()-start,'result':result,'pairing':traces}
(output/'result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out),flush=True)
assert result['computed'],result
assert not result['failed'] and not result['lambda_failed'],result
assert len(traces)==len(result['kappa_draws'])==2
assert all('signed_root' in source for source in result['tau']['sources']),result['tau']
assert 0<=result['p']<=1
# v1.2 invalid reliability must make the tau component non-rejecting.
invalid=p2_tau(design.model(),{'converged':True},result['tau_max'],[[.04]])
assert invalid['p_tau']==1.,invalid
print('v1.2 invalid-lambda non-rejection confirmed.',flush=True)
