"""Verify identical RNA scores without unused large derivative matrices."""
import hashlib
import json
import os
from pathlib import Path
import time
from unittest.mock import patch

import numpy as np

from bfix_v1 import load
from likelihood_v1 import Model
from test_likelihood_v1 import fixture


class LegacyPredictorAllocation(Model):
    def predictors(self,row,p,u=0,need_hessian=True):
        return super().predictors(row,p,u,need_hessian=True)


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Compute allocation required.')
    out=Path(os.environ['CODEX_REC_OUTPUT'])/('predictor_workspace_'+os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)
    # Positive-variance and boundary integrals must retain identical values/scores.
    small=Model(fixture())
    p=small.initial();p[:2]=[.4,-.5]
    maximum_ll=maximum_gradient=0.
    for variance in (0.,.0025):
        full=small.marginal(p,variance,with_hessian=True)
        lean=small.marginal(p,variance,with_hessian=False)
        np.testing.assert_allclose(full[0],lean[0],atol=1e-12,rtol=0)
        np.testing.assert_allclose(full[1],lean[1],atol=1e-12,rtol=0)
        maximum_ll=max(maximum_ll,abs(full[0]-lean[0]))
        maximum_gradient=max(maximum_gradient,float(np.max(np.abs(full[1]-lean[1]))))
    design,_=load()
    original=design.model()
    # Expand coordinate count only. This remains the 173-row agreement fixture,
    # not an actual 379-gene study or a representative full-fit benchmark.
    large=Model(original.rows,gene_count=379)
    legacy=LegacyPredictorAllocation(original.rows,gene_count=379)
    q=large.initial();q[:8]=[.4,-.5,.7,-.6,.8,-.3,.2,-.9]
    for u in (-.3,0.,.3):
        before=legacy.evaluate(q,u=u,need_hessian=False)
        after=large.evaluate(q,u=u,need_hessian=False)
        np.testing.assert_array_equal(before[0],after[0])
        np.testing.assert_array_equal(before[1],after[1])
    # Check full derivative compatibility on one row at the expanded dimension.
    full=large.evaluate(q,rows=large.rows[:1])
    lean=large.evaluate(q,rows=large.rows[:1],need_hessian=False)
    np.testing.assert_array_equal(full[0],lean[0])
    np.testing.assert_array_equal(full[1],lean[1])
    counts={}
    for label,model in [('legacy_predictor_allocation',legacy),('score_only',large)]:
        zero=np.zeros
        seen=[]
        def counted(shape,*args,**kwargs):
            if isinstance(shape,tuple) and shape==(model.size,model.size):seen.append(shape)
            return zero(shape,*args,**kwargs)
        with patch('likelihood_v1.np.zeros',side_effect=counted):
            model.evaluate(q,need_hessian=False)
        counts[label]=len(seen)
    assert counts['legacy_predictor_allocation']==len(large.rows)
    assert counts['score_only']==0
    timings={}
    # Alternate run order and preserve each timing; no full optimizer speed claim.
    for repeat in range(4):
        order=[('legacy_predictor_allocation',legacy),('score_only',large)]
        if repeat%2:order.reverse()
        for label,model in order:
            start=time.monotonic();model.evaluate(q,need_hessian=False)
            timings.setdefault(label,[]).append(time.monotonic()-start)
    report={'status':'predictor_workspace_checks_passed','coordinate_count':large.size,
            'synthetic_rows':len(large.rows),'unused_dense_allocation_count':counts,
            'bytes_per_avoided_dense_array':large.size**2*8,
            'marginal_max_ll_difference':maximum_ll,'marginal_max_gradient_difference':maximum_gradient,
            'timings_seconds':timings,'python_source_sha256':hashlib.sha256(Path('likelihood_v1.py').read_bytes()).hexdigest(),
            'numpy':np.__version__,'real_outcomes_read':False,
            'limits':'Identical arithmetic and coordinate-allocation benchmark only; no real-scale fit speed, '
                     'size/power/coverage, full-information or scientific validation claim.'}
    (out/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':main()
