"""Repeat a failed supplied P1 draw using the prescribed starts and bounds."""
import argparse
import hashlib
import json
import os
import time
from pathlib import Path
import numpy as np
from bfix_v1 import load,table


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--draw',type=int,default=2)
args=parser.parse_args()
root=Path(os.environ['CODEX_REC_OUTPUT'])
draws=table('draws/p1_Sstar.tsv')
assert 0<=args.draw<len(draws)
design,settings=load();draw=draws[args.draw]
assert int(draw['draw'])==args.draw
sample=design.with_stage([int(draw[p.individual]) for p in design.participants])
model=sample.model()
null=json.loads((root/'bfix_tests_21985999/observed_null.json').read_text())
assert null['converged']
out=root/('refit_failed_p1_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
start=time.monotonic()
fit=model.fit(np.asarray(null['parameters']))
def serial(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)
se=model.sandwich(fit['parameters'],sample.clusters) if fit['converged'] else float('nan')
observed_se=design.model().sandwich(null['parameters'],design.clusters)
observed_z=null['parameters'][2*model.G]/observed_se
z=fit['parameters'][2*model.G]/se if np.isfinite(se) and se>0 else float('nan')
okay=bool(fit['converged'] and np.isfinite(z))
result={'draw':args.draw,'fit':fit,'sandwich_se':se,'z':z,'observed_z':observed_z,
        'exceed':bool(not okay or abs(z)>=abs(observed_z)-1e-6),
        'elapsed_seconds':time.monotonic()-start,
        'source_sha256':hashlib.sha256(Path('likelihood_v1.py').read_bytes()).hexdigest(),
        'original_file':f'bfix_tests_21985999/p1_{args.draw:03d}.json',
        'note':'Original retained. Corrected numerical implementation uses the prescribed warm/cold starts; no B-MODEL comparison or full-test result claimed.'}
(out/'summary.json').write_text(json.dumps(result,default=serial,indent=2)+'\n')
print(json.dumps({'draw':args.draw,'converged':okay,'loglik':fit['loglik'],
                  'kappa':float(fit['parameters'][2*model.G]),'kept_start':fit['kept_start'],
                  'start_converged':fit['start_converged'],'z':z,'exceed':result['exceed']},default=serial),flush=True)
assert okay
