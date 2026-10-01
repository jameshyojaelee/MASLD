"""Outcome-free population diagnostic of label correction in selected RNA counts.

One gene, constant depths/priors, no covariates, known dispersions. This is
not full Model A calibration, empirical power, or a MASLD outcome analysis.
"""
import json
import os
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.special import logsumexp, ndtr
from likelihood_v1 import EPS, allowed_counts, bb_terms, fixed_bb


def scenario(alpha, n, true_kappa, var_e, order=60):
    nodes, weights = np.polynomial.hermite.hermgauss(order)
    weights /= np.sqrt(np.pi)
    latent = 2 + np.sqrt(2) * nodes  # Synthetic T ~ N(2, 1).
    cumulative = ndtr((np.array([1.5, 2.5])[None, :] - latent[:, None]) / np.sqrt(var_e))
    stage_prob = np.column_stack((cumulative[:, 0], cumulative[:, 1]-cumulative[:, 0], 1-cumulative[:, 1]))
    true_score = stage_prob @ np.arange(3.)
    stage_marginal = weights @ stage_prob
    mean = float(stage_marginal @ np.arange(3.))
    stage_variance = float(stage_marginal @ (np.arange(3.)-mean)**2)
    reliability = float(weights @ (true_score-mean)**2 / stage_variance)
    counts = allowed_counts(n, 10, 1, 10)
    rho = .02
    r = np.log((rho-EPS)/(1-EPS-rho))
    prior = np.array([.4, .35, .25])
    hom = np.stack((fixed_bb(counts, n, .01, .04), fixed_bb(counts, n, .99, .04)))

    def mixture(eta):
        hmass, score, _ = bb_terms(counts, n, eta, r)
        parts = np.vstack((hmass, hom)) + np.log(prior)[:, None]
        mass = logsumexp(parts, axis=0)
        derivative = np.exp(parts[0]-mass) * score[:, 0]
        return mass, derivative

    # Expected frequency of each included (recorded S, RNA allele count).
    joint = np.zeros((3, len(counts)))
    for w, probabilities, score in zip(weights, stage_prob, true_score):
        mass, _ = mixture(alpha * (1+true_kappa*(score-mean)))
        joint += w * probabilities[:, None] * np.exp(mass)[None, :]
    included = float(joint.sum())
    joint /= included

    def objective(parameters):
        a, kappa = parameters
        value = 0.
        gradient = np.zeros(2)
        for stage in range(3):
            centered = stage-mean
            mass, score = mixture(a*(1+kappa*centered))
            normalizer = logsumexp(mass)
            normalized_score = score - np.exp(mass-normalizer) @ score
            value -= joint[stage] @ (mass-normalizer)
            gradient -= (joint[stage] @ normalized_score) * np.array([1+kappa*centered, a*centered])
        return value, gradient

    fits = [minimize(objective, start, jac=True, method='BFGS',
                     options={'gtol': 1e-10, 'maxiter': 500})
            for start in ([alpha, true_kappa*reliability], [alpha, 0.])]
    fit = min(fits, key=lambda item: item.fun)
    residual = float(np.max(np.abs(objective(fit.x)[1])))
    assert residual < 1e-7, (alpha, n, true_kappa, residual)
    corrected = float(fit.x[1]/reliability)
    if true_kappa == 0:
        assert abs(corrected) < 1e-7
    return {'alpha_true': alpha, 'depth': n, 'kappa_true': true_kappa,
            'label_var_e': var_e, 'quadrature_order': order, 'lambda': reliability,
            'inclusion_probability': included, 'alpha_fitted': float(fit.x[0]),
            'kappa_recorded': float(fit.x[1]), 'kappa_corrected': corrected,
            'corrected_minus_true': corrected-true_kappa, 'max_gradient': residual}


out = Path(os.environ['CODEX_REC_OUTPUT']) / ('label_attenuation_'+os.environ['SLURM_JOB_ID'])
out.mkdir(exist_ok=False)
# Public-source planning error variances already independently rederived in C8.
errors = {'GSE193066_S3': .4246578312855484, 'PXD051911_A1_Y5': .2100726140806411}
results = []
for source, var_e in errors.items():
    for alpha in (.3, 1., 3.):
        for n in (20, 80):
            for kappa in (0., -.0001, .0001, -.125, .125):
                result = scenario(alpha, n, kappa, var_e)
                result['source_error'] = source
                results.append(result)
refinements = []
for result in results:
    if result['alpha_true'] == 3 and result['depth'] == 80 and abs(result['kappa_true']) == .125:
        fine = scenario(3., 80, result['kappa_true'], result['label_var_e'], 120)
        difference = abs(fine['kappa_corrected']-result['kappa_corrected'])
        assert difference < 1e-5, difference
        refinements.append({'source': result['source_error'], 'kappa': result['kappa_true'],
                            'corrected_difference_GH60_GH120': difference})
summary = {'results': results, 'quadrature_checks': refinements,
           'numpy': np.__version__, 'scipy': scipy.__version__,
           'design': 'Synthetic T~N(2,1), score E[S|T]; RNA beta-binomial H/R/A mixture; selected counts.',
           'limits': 'One gene, fixed n/priors/errors/dispersion, no z or stage-dependent technical effects. Population bias only; not sampling coverage, power or full Model A calibration.'}
(out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps({'output': str(out), 'scenarios': len(results),
                  'max_abs_corrected_bias_at_margin': max(abs(r['corrected_minus_true']) for r in results if abs(r['kappa_true']) == .125)}), flush=True)
