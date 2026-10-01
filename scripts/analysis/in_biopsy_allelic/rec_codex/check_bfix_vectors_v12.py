"""Fixed-vector shared-fixture results; synthetic counts only, no MLE claims."""
import hashlib
import json
import os
import platform
import time
from pathlib import Path
import numpy as np
import scipy
from bfix_v1 import FIXTURE,load,parameters


for line in (FIXTURE/'SHA256SUMS').read_text().splitlines():
    expected,name=line.split(maxsplit=1)
    assert hashlib.sha256((FIXTURE/name).read_bytes()).hexdigest()==expected,name
design,settings=load();model=design.model()
start=time.monotonic();results={}
for path in sorted((FIXTURE/'params').glob('theta*.json')):
    p,v=parameters(path)
    if v==0:
        ll,gradient,_=model.evaluate(p)
        pieces={g:model.evaluate(p,[r for r in model.rows if r.gene==g])[0] for g in range(model.G)}
    else:
        pieces={};ll,gradient=model.marginal(p,v,gene_terms=pieces)
        check_ll,check_gradient=model.marginal(p,v,nodes=40)
        np.testing.assert_allclose(ll,check_ll,rtol=1e-8,atol=1e-8)
        np.testing.assert_allclose(gradient,check_gradient,rtol=1e-6,atol=1e-8)
    results[path.stem]={'loglik':ll,'gradient':gradient.tolist(),'v':v,'per_gene_loglik':pieces}
    print(path.stem,ll,flush=True)
out={'fixture':str(FIXTURE),'spec':'v1.2','rows':len(model.rows),'participants':len(design.participants),
     'design_columns':settings['design_columns'],'vectors':results,
     'elapsed_seconds':time.monotonic()-start,'python':platform.python_version(),
     'numpy':np.__version__,'scipy':scipy.__version__,
     'note':'Full three-class numerator retained; planted classes are not used by the likelihood.'}
output=Path(os.environ['CODEX_REC_OUTPUT'])/'bfix_fixed_vectors_v12.json'
output.write_text(json.dumps(out,indent=2)+'\n')
print('Shared fixed-vector results written; MLE and draw fits not run by this script.',flush=True)
