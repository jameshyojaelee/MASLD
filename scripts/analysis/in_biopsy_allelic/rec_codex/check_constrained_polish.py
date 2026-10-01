"""Check the recovered failed draw and invariants of constrained polishing."""
import json
import os
from pathlib import Path
import numpy as np
from bfix_v1 import load,table


root=Path(os.environ['CODEX_REC_OUTPUT'])
design,settings=load();draw=table('draws/p1_Sstar.tsv')[2]
sample=design.with_stage([int(draw[p.individual]) for p in design.participants])
model=sample.model();old=json.loads((root/'bfix_tests_21985999/p1_002.json').read_text())
reference=json.loads((root/'failed_p1_21991584/summary.json').read_text())
parameters,history=model.constrained_polish(old['fit']['parameters'],model.bounds())
ll,gradient,hessian=model.evaluate(parameters)
projected=gradient.copy()
for j,(lo,hi) in enumerate(model.bounds()):
    assert (lo is None or parameters[j]>=lo-1e-10) and (hi is None or parameters[j]<=hi+1e-10)
    if (lo is not None and parameters[j]<=lo+1e-8 and gradient[j]<0) or (hi is not None and parameters[j]>=hi-1e-8 and gradient[j]>0):projected[j]=0
assert np.max(np.abs(projected))<1e-6
np.testing.assert_allclose(ll,reference['loglik'],atol=1e-6,rtol=0)
np.testing.assert_allclose(parameters[2*model.G],reference['parameters'][2*model.G],atol=1e-4,rtol=0)
assert all(b['loglik']>=a['loglik'] for a,b in zip(history,history[1:]))
se=model.sandwich(parameters,sample.clusters);assert np.isfinite(se) and se>0
# An already converged draw must remain at its saved optimum.
first=table('draws/p1_Sstar.tsv')[0]
first_sample=design.with_stage([int(first[p.individual]) for p in design.participants])
first_model=first_sample.model();first_fit=json.loads((root/'bfix_tests_21985999/p1_000.json').read_text())['fit']
same,same_history=first_model.constrained_polish(first_fit['parameters'],first_model.bounds())
np.testing.assert_array_equal(same,first_fit['parameters'])
assert len(same_history)==1
# Explicitly fixed coordinates must remain fixed during the fallback.
fixed_bounds=model.bounds();fixed_bounds[0]=(parameters[0],parameters[0])
fixed,_=model.constrained_polish(parameters,fixed_bounds)
assert fixed[0]==parameters[0]
out=root/('constrained_polish_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
report={'draw':2,'loglik':ll,'kappa':float(parameters[2*model.G]),
        'max_projected_gradient':float(np.max(np.abs(projected))),
        'sandwich_se':se,'z':float(parameters[2*model.G]/se),
        'history':history,'checks':['original bounds','original convergence threshold',
                                   'independent diagnostic agreement','likelihood monotonicity',
                                   'already converged fit unchanged','fixed coordinates unchanged'],
        'limit':'Polish verification; full warm/cold draw refit and B-MODEL comparison remain pending.'}
(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
