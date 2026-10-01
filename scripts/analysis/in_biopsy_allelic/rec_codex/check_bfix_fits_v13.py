"""Observed shared-fixture fits with each completed fit saved separately."""
import json
import os
import platform
import time
from pathlib import Path
import numpy as np
import scipy
from bfix_v1 import FIXTURE,load


def serial(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)


design,settings=load();model=design.model()
output=Path(os.environ['CODEX_REC_OUTPUT'])/('bfix_fits_'+os.environ['SLURM_JOB_ID'])
output.mkdir(exist_ok=False)
inventory={'fixture':str(FIXTURE),'spec':'v1.3','participants':60,'rows':173,'genes':8,
           'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__}
start=time.monotonic();null=model.fit()
null_out={**inventory,'elapsed_seconds':time.monotonic()-start,'fit':null}
(output/'null.json').write_text(json.dumps(null_out,default=serial,indent=2)+'\n')
print('section3',null['converged'],null['loglik'],null['parameters'][2*model.G],flush=True)
assert null['converged'],null['gradient']
start=time.monotonic();marginal=model.fit_marginal(settings['tau_max_fixture'],null['parameters'])
lr=max(0.,2*(marginal['loglik']-null['loglik']))
marginal_out={**inventory,'elapsed_seconds':time.monotonic()-start,'fit':marginal,'LR':lr}
(output/'marginal.json').write_text(json.dumps(marginal_out,default=serial,indent=2)+'\n')
print('section4',marginal['converged'],marginal['loglik'],marginal['v'],'LR',lr,flush=True)
assert marginal['converged'] and marginal['null']['converged'],marginal['gradient']
