"""Check optional first-derivative evaluations against the full RNA likelihood."""
import json
import os
import time
from dataclasses import replace
from pathlib import Path
import numpy as np
from bfix_v1 import FIXTURE, load, parameters
from likelihood_v1 import conditional_terms
from test_likelihood_v1 import fixture


design, settings = load()
model = design.model()
root = Path(os.environ['CODEX_REC_OUTPUT'])
out = root/('score_path_'+os.environ['SLURM_JOB_ID'])
out.mkdir(exist_ok=False)
maximum_ll = 0.
maximum_score = 0.
for path in sorted((FIXTURE/'params').glob('theta*.json')):
    p, v = parameters(path)
    for u in (0., -.3, .3):
        full = model.evaluate(p,u=u)
        lean = model.evaluate(p,u=u,need_hessian=False)
        maximum_ll = max(maximum_ll, abs(full[0]-lean[0]))
        maximum_score = max(maximum_score, float(np.max(np.abs(full[1]-lean[1]))))
        np.testing.assert_allclose(full[0], lean[0], atol=1e-12, rtol=0)
        np.testing.assert_allclose(full[1], lean[1], atol=1e-12, rtol=0)
        assert lean[2] is None
for base in fixture():
    for eta in (-5.,0.,5.):
        for r in (-15.,-3.,5.):
            for row in (base, replace(base,minor=0,percent=0), replace(base,orientation=-base.orientation)):
                full = conditional_terms(row,eta,r)
                lean = conditional_terms(row,eta,r,False)
                np.testing.assert_allclose(full[0],lean[0],atol=1e-12,rtol=0)
                np.testing.assert_allclose(full[1],lean[1],atol=1e-12,rtol=0)
# Same parameters, same process, warmed caches, five repeated evaluations each.
p = np.array(json.loads((root/'bfix_tests_21985999/observed_null.json').read_text())['parameters'])
timing = {}
for mode in (True,False):
    model.evaluate(p,need_hessian=mode)
    start = time.monotonic()
    for repeat in range(5):model.evaluate(p,need_hessian=mode)
    timing['full' if mode else 'score_only'] = (time.monotonic()-start)/5
arithmetic = {'fixed_vector_max_ll_difference':maximum_ll,
              'fixed_vector_max_score_difference':maximum_score,
              'seconds_per_evaluation':timing,
              'evaluation_speed_ratio':timing['full']/timing['score_only']}
(out/'arithmetic.json').write_text(json.dumps(arithmetic,indent=2)+'\n')
print(json.dumps({'arithmetic_complete':arithmetic}),flush=True)
# Freeze the pre-change fit as a start and verify the fit remains within the
# independent agreement tolerance, including its final Hessian and convergence.
fit = model.fit_single(p)
reference = json.loads((root/'bfix_tests_21985999/observed_null.json').read_text())
assert fit['converged']
np.testing.assert_allclose(fit['loglik'],reference['loglik'],atol=1e-6,rtol=0)
np.testing.assert_allclose(fit['parameters'][2*model.G],p[2*model.G],atol=1e-4,rtol=0)
summary = {**arithmetic,
           'observed_warm_fit_loglik':fit['loglik'],'observed_warm_fit_kappa':float(fit['parameters'][2*model.G]),
           'observed_warm_fit_converged':fit['converged'],
           'limits':'Synthetic arithmetic and one warm-start fit only; no full-fit speed or real-scale calibration claim.'}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary),flush=True)
