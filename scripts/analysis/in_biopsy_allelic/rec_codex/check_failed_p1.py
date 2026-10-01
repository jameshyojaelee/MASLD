"""Diagnose and polish shared P1 draw 2 without replacing its original result."""
import json
import os
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from bfix_v1 import load, table


design, settings = load()
stages = table('draws/p1_Sstar.tsv')[2]
assert int(stages['draw']) == 2
sample = design.with_stage([int(stages[p.individual]) for p in design.participants])
model = sample.model()
base = Path(os.environ['CODEX_REC_OUTPUT'])
original = json.loads((base/'bfix_tests_21985999/p1_002.json').read_text())
parameters = np.asarray(original['fit']['parameters'],float)
out = base/('failed_p1_'+os.environ['SLURM_JOB_ID'])
out.mkdir(exist_ok=False)
ll, gradient, hessian = model.evaluate(parameters)
np.testing.assert_allclose(ll,original['fit']['loglik'],atol=1e-8,rtol=0)
np.testing.assert_allclose(gradient,original['fit']['gradient'],atol=1e-8,rtol=0)
checks = []
for coordinate in (4,12,16,20,24):
    shift = np.zeros(model.size)
    shift[coordinate] = 1e-4 if coordinate == 12 else 1e-5
    numeric = (model.evaluate(parameters+shift,need_hessian=False)[0]
               - model.evaluate(parameters-shift,need_hessian=False)[0])/(2*shift[coordinate])
    checks.append({'coordinate':coordinate,'analytic':float(gradient[coordinate]),'numeric':float(numeric)})
(out/'gradient_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps({'gradient_checks':checks}),flush=True)
bounds = model.bounds()
history = []
for iteration in range(100):
    ll, gradient, hessian = model.evaluate(parameters)
    active = np.array([(lo is not None and parameters[j] <= lo+1e-8 and gradient[j] < 0)
                       or (hi is not None and parameters[j] >= hi-1e-8 and gradient[j] > 0)
                       for j,(lo,hi) in enumerate(bounds)])
    projected = gradient.copy()
    projected[active] = 0
    maximum = float(np.max(np.abs(projected)))
    history.append({'iteration':iteration,'loglik':ll,'max_projected_gradient':maximum})
    if maximum < 1e-6:
        break
    free = np.flatnonzero(~active)
    curvature = -hessian[np.ix_(free,free)]
    values, vectors = np.linalg.eigh((curvature+curvature.T)/2)
    # Positive curvature regularization only proposes a direction. Acceptance
    # is judged on the unchanged exact likelihood and original parameter bounds.
    floor = max(1e-9, float(np.max(np.abs(values)))*1e-10)
    curvature = (vectors*np.maximum(values,floor))@vectors.T
    score = gradient[free]
    limits = [(None if bounds[j][0] is None else bounds[j][0]-parameters[j],
               None if bounds[j][1] is None else bounds[j][1]-parameters[j]) for j in free]
    def quadratic(step):
        return .5*step@curvature@step-score@step, curvature@step-score
    direction_fit = minimize(quadratic,np.zeros(len(free)),jac=True,method='SLSQP',bounds=limits,
                             options={'ftol':1e-14,'maxiter':1000})
    step = np.zeros(model.size)
    step[free] = direction_fit.x
    improved = False
    for damping in 2.**-np.arange(24):
        trial = parameters+damping*step
        trial_ll = model.evaluate(trial,need_hessian=False)[0]
        if trial_ll > ll:
            parameters = trial
            improved = True
            break
    if not improved:
        break
ll, gradient, hessian = model.evaluate(parameters)
active = np.array([(lo is not None and parameters[j] <= lo+1e-8 and gradient[j] < 0)
                   or (hi is not None and parameters[j] >= hi-1e-8 and gradient[j] > 0)
                   for j,(lo,hi) in enumerate(bounds)])
projected = gradient.copy();projected[active] = 0
converged = bool(np.max(np.abs(projected)) < 1e-6)
se = model.sandwich(parameters,sample.clusters) if converged else float('nan')
result = {'original_file':'bfix_tests_21985999/p1_002.json','draw':2,
          'loglik':ll,'original_loglik':original['fit']['loglik'],'parameters':parameters.tolist(),
          'gradient':gradient.tolist(),'max_projected_gradient':float(np.max(np.abs(projected))),
          'converged':converged,'sandwich_se':se,
          'z':float(parameters[2*model.G]/se) if np.isfinite(se) else None,
          'history':history,'gradient_checks':checks,
          'limit':'Numerical diagnostic only. Original draw preserved; no B-MODEL agreement verdict or change to convergence threshold.'}
(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('loglik','max_projected_gradient','converged','z')}),flush=True)
