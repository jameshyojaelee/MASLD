"""Profile expected information and attribute efficient-score variance to cohorts.

Input matrices must be expected information for the same retained gene set,
plant, coordinate order and complete-F centering. No participant data needed.
"""
import numpy as np


def project_information(cohort_matrices, target, tolerance=1e-10):
    names = sorted(cohort_matrices)
    matrices = []
    for name in names:
        matrix = np.asarray(cohort_matrices[name], float)
        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or not np.isfinite(matrix).all():
            raise ValueError('finite square information matrix required: '+name)
        np.testing.assert_allclose(matrix, matrix.T, atol=1e-10, rtol=1e-10)
        matrices.append((matrix+matrix.T)/2)
    if not matrices or len({m.shape for m in matrices}) != 1:
        raise ValueError('nonempty matching cohort matrices required')
    total = np.sum(matrices, axis=0)
    p = len(total)
    if not 0 <= target < p:
        raise ValueError('target coordinate outside matrix')
    nuisance = np.flatnonzero(np.arange(p) != target)
    # Scale each coordinate before deciding numerical rank, so unit changes do
    # not change which estimable nuisance directions are retained.
    diagonal = np.diag(total)
    if np.any(diagonal < -tolerance):
        raise ValueError('negative expected-information diagonal')
    scale = np.sqrt(np.maximum(diagonal, 0))
    active = scale > 0
    if np.any(np.abs(total[~active]) > tolerance):
        raise ValueError('zero-diagonal information coordinate has nonzero coupling')
    active_indices = np.flatnonzero(active)
    if active.any():
        normalized = total[np.ix_(active_indices, active_indices)] / np.outer(scale[active], scale[active])
        if np.linalg.eigvalsh(normalized)[0] < -tolerance:
            raise ValueError('total matrix is not positive semidefinite')
    informative_nuisance = nuisance[active[nuisance]]
    coefficients = np.zeros(p)
    coefficients[target] = 1
    rank = 0
    if active[target] and len(informative_nuisance):
        a = total[np.ix_(informative_nuisance, informative_nuisance)]
        a /= np.outer(scale[informative_nuisance], scale[informative_nuisance])
        eigenvalues, vectors = np.linalg.eigh(a)
        retained = eigenvalues > tolerance*max(1., float(eigenvalues[-1]))
        rank = int(retained.sum())
        cross = total[informative_nuisance, target]/scale[informative_nuisance]
        if np.max(np.abs(vectors[:, ~retained].T @ cross), initial=0) > tolerance*max(1., np.linalg.norm(cross)):
            raise ValueError('target coupling enters a numerically unidentifiable nuisance direction')
        beta = vectors[:, retained] @ ((vectors[:, retained].T @ cross)/eigenvalues[retained])
        coefficients[informative_nuisance] = -beta/scale[informative_nuisance]
    contributions = {}
    response = {}
    for name, matrix in zip(names, matrices):
        if active.any():
            scaled = matrix[np.ix_(active_indices, active_indices)] / np.outer(scale[active], scale[active])
            if np.linalg.eigvalsh(scaled)[0] < -tolerance:
                raise ValueError('cohort matrix is not positive semidefinite: '+name)
        value = float(coefficients @ matrix @ coefficients)
        if value < -tolerance*max(1., diagonal[target]):
            raise ValueError('negative efficient-score contribution: '+name)
        contributions[name] = max(value, 0.)
        # Local response to a cohort-specific target perturbation differs from
        # the nonnegative variance attribution when nuisances are shared.
        response[name] = float(coefficients @ matrix[:, target])
    information = float(sum(contributions.values()))
    residual = total @ coefficients
    np.testing.assert_allclose(residual[nuisance], 0., atol=1e-8*max(1., np.max(np.abs(total))), rtol=0)
    np.testing.assert_allclose(information, residual[target], atol=1e-8*max(1., diagonal[target]), rtol=1e-8)
    np.testing.assert_allclose(sum(response.values()), information, atol=1e-8*max(1., diagonal[target]), rtol=1e-8)
    identifiable = information > tolerance*max(1., diagonal[target])
    return {'information': information, 'identifiable': bool(identifiable),
            'cohort_information': contributions, 'nuisance_rank': rank,
            'cohort_local_response': response,
            'efficient_score_coefficients': coefficients}
