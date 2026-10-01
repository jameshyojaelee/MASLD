"""Diagnose the single nonconverged synthetic P3 draw without extending B."""
import json
import os
from pathlib import Path

import numpy as np
from likelihood_v1 import Row,Model,conditional_draw
from design_v1 import Participant,Design
from tests_v1 import rng_for


output=Path(os.environ['CODEX_REC_OUTPUT'])
saved=json.loads((output/'fit_v1_smoke.json').read_text())
rng=np.random.default_rng(20260923)
people=[Participant(f'P{i:03d}','synthetic',i%3,(float(i%5-2)/10,float(i%7-3)/10),.01,35.,4) for i in range(36)]
rows=[Row(g,p.individual,p.cohort,30+i%11,15,.6,.3,.1,.01,.02,.05,(-1)**(g+i),
          0.,(),0.,0.,10,2,10) for g in range(4) for i,p in enumerate(people)]
model=Design(people,rows).model();planted=model.initial();planted[:4]=[-1.2,-.6,.7,1.3];planted[8]=.08
drawn=[]
for row in model.rows:
    eta,r,_,_=model.predictors(row,planted)
    drawn.append(conditional_draw(row,eta,r,*rng.random(2)))
model=Design(people,drawn).model();observed=np.asarray(saved['null_parameters'])
np.testing.assert_allclose(model.evaluate(observed)[0],saved['null_loglik'],atol=1e-10,rtol=0)
simulated=[]
for row,(u1,u2) in zip(model.rows,rng_for(3,0).random((len(model.rows),2))):
    eta,r,_,_=model.predictors(row,observed)
    simulated.append(conditional_draw(row,eta,r,u1,u2))
fitted=Model(simulated).fit_marginal(saved['tau_max'],observed)
def details(fit,marginal=False):
    p=np.r_[fit['parameters'],fit['v']] if marginal else fit['parameters']
    bounds=model.bounds()+[(0,None)] if marginal else model.bounds()
    gradient=fit['gradient'].copy()
    active=[]
    for j,(lo,hi) in enumerate(bounds):
        if (lo is not None and p[j]<=lo+1e-10 and gradient[j]<0) or (hi is not None and p[j]>=hi-1e-10 and gradient[j]>0):
            gradient[j]=0;active.append(j)
    return {'converged':fit['converged'],'loglik':fit['loglik'],'parameters':p.tolist(),
            'gradient':fit['gradient'].tolist(),'projected_gradient':gradient.tolist(),
            'max_projected_gradient':float(np.abs(gradient).max()),'active':active}
result={'null':details(fitted['null']),'marginal':details(fitted,True)}
(output/os.environ.get('CODEX_REC_DIAGNOSTIC_FILE','failed_p3_draw0_diagnostics.json')).write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2),flush=True)
