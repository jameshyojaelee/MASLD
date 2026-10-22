"""Separate synthetic GL model with checked null score-integral polish.

Run the original constrained solver before the score-integral refinement.
Only constrained_polish changes. Inherited starts, warm/cold ties, bounds,
positive-variance integration and the original convergence cutoff remain.
The empirical directional-integral sign check is not a certified bound.
"""
import hashlib
from pathlib import Path


def model_class(base,original,gamma,polish,*,expected_polish_sha,event_sink=None):
    if hashlib.sha256(Path(polish.__file__).read_bytes()).hexdigest()!=expected_polish_sha:
        raise RuntimeError('Checked polish source differs')
    if base.constrained_polish is not original.Model.constrained_polish:
        raise RuntimeError('Expected original constrained polish')

    class PolishedModel(base):
        def constrained_polish(self,start,bounds,maxiter=100):
            if maxiter!=100:
                raise RuntimeError('Unexpected original polish iteration request')
            baseline,baseline_history=original.Model.constrained_polish(self,start,bounds,maxiter=maxiter)
            p,evidence=polish.polish(self,baseline,bounds,original,gamma,maxiter=20)
            evidence['original_constrained_baseline_history']=baseline_history
            evidence['original_constrained_baseline_parameters']=baseline.tolist()
            evidence['original_constrained_solver_retained']=True
            if not hasattr(self,'null_polish_evidence'):
                self.null_polish_evidence=[]
            self.null_polish_evidence.append(evidence)
            if event_sink:
                event_sink(dict(event='null_score_integral_polish',
                    converged=evidence['converged'],
                    original_max_projected_gradient=evidence['original_max_projected_gradient'],
                    Gamma_max_projected_gradient=evidence['Gamma_max_projected_gradient'],
                    Gamma_evaluations=evidence['Gamma_evaluations']))
            return p,evidence['history']

    for name in ('fit','fit_single','fit_marginal','fit_marginal_single','bounds',
                 'initial','cold_start','informative','predictors','evaluate','sandwich','marginal'):
        if getattr(PolishedModel,name) is not getattr(base,name):
            raise RuntimeError('Unexpected changed method: '+name)
    return PolishedModel
