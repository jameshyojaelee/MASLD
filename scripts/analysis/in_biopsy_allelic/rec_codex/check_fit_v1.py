"""Synthetic fit and small end-to-end inference smoke; no calibration claim."""
import json
import os
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from likelihood_v1 import Row,Model,conditional_draw
from design_v1 import Participant,Design
from tests_v1 import p1_design,p3_design,ordinal_stage_probabilities,rng_for


seed=20260923
rng=np.random.default_rng(seed)
people=[Participant(f'P{i:03d}','synthetic',i%3,(float(i%5-2)/10,float(i%7-3)/10),.01,35.,4) for i in range(36)]
rows=[Row(g,p.individual,p.cohort,30+i%11,15,.6,.3,.1,.01,.02,.05,(-1)**(g+i),
          0.,(),0.,0.,10,2,10) for g in range(4) for i,p in enumerate(people)]
design=Design(people,rows);model=design.model();planted=model.initial();planted[:4]=[-1.2,-.6,.7,1.3]
planted[8]=.08
drawn=[]
for row in model.rows:
    eta,r,_,_=model.predictors(row,planted)
    drawn.append(conditional_draw(row,eta,r,*rng.random(2)))
design=Design(people,drawn);model=design.model()
start=time.monotonic();null=model.fit();null_time=time.monotonic()-start
summary={'seed':seed,'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
         'participants':len(people),'genes':4,'rows':len(rows),'null_fit_seconds':null_time,
         'null_converged':null['converged'],'null_loglik':null['loglik'],
         'null_parameters':null['parameters'].tolist()}
print(json.dumps(summary),flush=True)
assert null['converged'],summary
tau_max=.125*float(np.median(np.abs(null['parameters'][:4])))
start=time.monotonic();marginal=model.fit_marginal(tau_max,null['parameters'])
summary.update({'tau_max':tau_max,'marginal_fit_seconds':time.monotonic()-start,
                'marginal_converged':marginal['converged'],'marginal_loglik':marginal['loglik'],
                'tau_hat':float(np.sqrt(marginal['v']))})
print(json.dumps(summary),flush=True)
assert marginal['converged'],summary
assert marginal['loglik']>=null['loglik']-1e-8
# Small supplied stage/uniform fixtures exercise resampling and failure handling.
# Draw counts here test execution only; scientific P1/P3 use B=1000.
stages=[rng_for(1,b).permutation(design.stage.astype(int)) for b in range(2)]
summary['p1_smoke']=p1_design(design,stages)
summary['p3_smoke']=p3_design(design,tau_max,draws=2)
assert summary['p1_smoke']['computed'],summary['p1_smoke']
assert summary['p3_smoke']['computed'],summary['p3_smoke']
summary['stage_law']=design.stage_law()['converged']
assert summary['stage_law']
output=Path(os.environ['CODEX_REC_OUTPUT'])/'fit_v1_smoke.json'
output.write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary),flush=True)
