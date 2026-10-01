"""Scientific invariants of shared-nuisance expected-information projection."""
import json
import os
from pathlib import Path
import numpy as np
from information_projection import project_information


scores = {'c1': np.array([[1.,2.,0.], [2.,0.,1.], [-1.,1.,1.]]),
          'c2': np.array([[3.,1.,2.], [1.,-1.,0.], [2.,2.,-1.]])}
matrices = {name: x.T@x for name, x in scores.items()}
result = project_information(matrices, 0)
# Independently fit the target score against nuisance scores at row level.
stacked = np.vstack(list(scores.values()))
beta = np.linalg.lstsq(stacked[:,1:], stacked[:,0], rcond=None)[0]
residual = stacked[:,0]-stacked[:,1:]@beta
np.testing.assert_allclose(result['information'], residual@residual, atol=1e-12)
for name, x in scores.items():
    r = x[:,0]-x[:,1:]@beta
    np.testing.assert_allclose(result['cohort_information'][name], r@r, atol=1e-12)
    np.testing.assert_allclose(result['cohort_local_response'][name], r@x[:,0], atol=1e-12)
# Nuisance basis changes, large unit changes, and duplicate nuisance columns
# must preserve the target's efficient information and its cohort attribution.
for transform in (np.array([[1.,2.],[-.5,1.]]), np.diag([1e-4,1e4])):
    changed = {name: np.column_stack((x[:,0], x[:,1:]@transform)) for name,x in scores.items()}
    got = project_information({name:x.T@x for name,x in changed.items()}, 0)
    np.testing.assert_allclose(got['information'], result['information'], atol=1e-10)
    np.testing.assert_allclose(list(got['cohort_information'].values()), list(result['cohort_information'].values()), atol=1e-10)
duplicate = {name: np.column_stack((x, x[:,1]*2)) for name,x in scores.items()}
got = project_information({name:x.T@x for name,x in duplicate.items()}, 0)
np.testing.assert_allclose(got['information'], result['information'], atol=1e-10)
assert got['nuisance_rank'] == 2
# Complete confounding of stage and reference bias must report no information.
confounded = {name: np.column_stack((x[:,0],x[:,0],x[:,1:])) for name,x in scores.items()}
got = project_information({name:x.T@x for name,x in confounded.items()}, 0)
assert not got['identifiable'] and got['information'] < 1e-20
# Splitting a cohort leaves total information unchanged and contributions additive.
split = dict(matrices);original = split.pop('c1');split['c1a'] = .3*original;split['c1b'] = .7*original
got = project_information(split, 0)
np.testing.assert_allclose(got['information'], result['information'], atol=1e-12)
np.testing.assert_allclose(got['cohort_information']['c1a']+got['cohort_information']['c1b'], result['cohort_information']['c1'], atol=1e-12)
# A simple independent analytic counterexample: shared nuisance coefficient
# beta = (1*2 + 3*1)/(2**2 + 1**2) = 1. Residuals (-1,2) give information 5,
# cohort variance contributions (1,4), response contributions (-1,6).
example_scores = {'c1': np.array([[1.,2.]]), 'c2': np.array([[3.,1.]])}
example = project_information({name:x.T@x for name,x in example_scores.items()}, 0)
np.testing.assert_allclose(list(example['cohort_information'].values()), [1.,4.], atol=1e-12)
np.testing.assert_allclose(list(example['cohort_local_response'].values()), [-1.,6.], atol=1e-12)
example['efficient_score_coefficients'] = example['efficient_score_coefficients'].tolist()
out = Path(os.environ['CODEX_REC_OUTPUT'])/('information_projection_'+os.environ['SLURM_JOB_ID'])
out.mkdir(exist_ok=False)
report = {**result, 'efficient_score_coefficients': result['efficient_score_coefficients'].tolist(),
          'checks': ['independent row-level nuisance regression', 'nuisance basis and units',
                     'duplicate nuisance rank', 'complete confounding', 'cohort split additivity'],
          'shared_nuisance_response_example': example,
          'limit': 'Synthetic invariants only; full approved ancestry-weighted matrices still needed.'}
(out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report), flush=True)
