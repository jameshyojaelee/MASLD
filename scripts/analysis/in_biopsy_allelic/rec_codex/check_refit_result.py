"""Re-evaluate a saved supplied-draw refit independently of optimization."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import numpy as np
from bfix_v1 import load,table


root=Path(os.environ['CODEX_REC_OUTPUT'])
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source',type=Path,default=root/'refit_failed_p1_21992790/summary.json')
source=parser.parse_args().source
assert source.resolve().is_relative_to(root.resolve())
result=json.loads(source.read_text());fit=result['fit']
assert result['source_sha256']==hashlib.sha256(Path('likelihood_v1.py').read_bytes()).hexdigest()
assert fit['converged'] and all(fit['start_converged'].values())
assert fit['kept_start']=='warm' and fit['start_logliks']['cold']-fit['start_logliks']['warm']<1e-6
draw_index=int(result['draw'])
design,settings=load();draw=table('draws/p1_Sstar.tsv')[draw_index]
assert int(draw['draw'])==draw_index
sample=design.with_stage([int(draw[p.individual]) for p in design.participants]);model=sample.model()
p=np.asarray(fit['parameters']);ll,gradient,hessian=model.evaluate(p)
np.testing.assert_allclose(ll,fit['loglik'],atol=1e-8,rtol=0)
np.testing.assert_allclose(gradient,fit['gradient'],atol=1e-8,rtol=0)
projected=gradient.copy()
for j,(lo,hi) in enumerate(model.bounds()):
    assert (lo is None or p[j]>=lo-1e-10) and (hi is None or p[j]<=hi+1e-10)
    if (lo is not None and p[j]<=lo+1e-8 and gradient[j]<0) or (hi is not None and p[j]>=hi-1e-8 and gradient[j]>0):projected[j]=0
assert np.max(np.abs(projected))<1e-6
se=model.sandwich(p,sample.clusters);assert np.isfinite(se) and se>0
np.testing.assert_allclose(se,result['sandwich_se'],atol=1e-10,rtol=1e-8)
z=float(p[2*model.G]/se)
np.testing.assert_allclose(z,result['z'],atol=1e-8,rtol=1e-8)
assert result['exceed']==bool(abs(z)>=abs(result['observed_z'])-1e-6)
report={'saved_result':str(source),'saved_result_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'draw':draw_index,'loglik':ll,'kappa':float(p[2*model.G]),'max_projected_gradient':float(np.max(np.abs(projected))),
        'sandwich_se':float(se),'z':z,'exceed':result['exceed'],
        'start_converged':fit['start_converged'],'start_logliks':fit['start_logliks'],'kept_start':fit['kept_start'],
        'checks':'Source hash, both starts, bounds, convergence, likelihood, gradient, sandwich and exceedance rechecked.',
        'saved_result_note':result['note'],
        'limit':'Synthetic draw recovery only; blind comparison and full 50-draw tests remain pending.'}
out=root/('verified_refit_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
