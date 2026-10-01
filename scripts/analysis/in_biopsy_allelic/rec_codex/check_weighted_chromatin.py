#!/usr/bin/env python3
"""Independent primal checks of chosen donor loss and rank-constrained fit."""
import json
import os
import platform
from pathlib import Path

import numpy as np
from weighted_chromatin import donor_weights, fit_predict


def main():
    job = os.environ.get("SLURM_JOB_ID")
    if not job:
        raise RuntimeError("CPU SLURM allocation required")
    rng = np.random.default_rng(20260930)
    x, y, query = rng.normal(size=(13, 5)), rng.normal(size=(13, 7)), rng.normal(size=(4, 5))
    strata = [(0,)] * 5 + [(1,)] * 3 + [(i,) for i in range(2, 7)]
    w = donor_weights(strata, 0.1)
    np.testing.assert_allclose(w @ np.ones(13), np.ones(13), atol=1e-14)
    assert np.linalg.eigvalsh(w).min() > 0
    errors = rng.normal(size=(13, 7))
    pair_loss = sum(np.sum((errors[i] - errors[j]) ** 2) / (sum(key == strata[i] for key in strata) - 1)
                    for i in range(13) for j in range(i + 1, 13) if strata[i] == strata[j])
    np.testing.assert_allclose(np.sum(errors * ((w - np.eye(13)) @ errors)), 0.1 * 13 / 8 * pair_loss, rtol=1e-12)
    mu, sd = x.mean(0), x.std(0)
    a, q = (x - mu) / sd, (query - mu) / sd
    yc = y - y.mean(0)
    fraction = 0.03
    penalty = fraction * np.linalg.eigvalsh(a.T @ w @ a / 5).max() * 5
    normal = a.T @ w @ a + penalty * np.eye(5)
    rhs = a.T @ w @ yc
    beta = np.linalg.solve(normal, rhs)
    expected = q @ beta + y.mean(0)
    observed = fit_predict(x, y, query, w, fraction)
    np.testing.assert_allclose(observed, expected, atol=1e-11)
    # Independent primal whitening of the coefficient quadratic obtains the
    # rank-constrained optimum. This differs from the implementation's dual SVD.
    chol = np.linalg.cholesky(normal)
    whitened = np.linalg.solve(chol, rhs)
    left, singular, right = np.linalg.svd(whitened, full_matrices=False)
    truncated = (left[:, :2] * singular[:2]) @ right[:2]
    beta_rank = np.linalg.solve(chol.T, truncated)
    expected_rank = q @ beta_rank + y.mean(0)
    observed_rank = fit_predict(x, y, query, w, fraction, rank=2)
    np.testing.assert_allclose(observed_rank, expected_rank, atol=1e-11)
    order = rng.permutation(13)
    reordered = fit_predict(x[order], y[order], query, w[np.ix_(order, order)], fraction, rank=2)
    np.testing.assert_allclose(reordered, observed_rank, atol=1e-11)
    shifted = fit_predict(x, y + 3.0, query, w, fraction, rank=2)
    np.testing.assert_allclose(shifted, observed_rank + 3.0, atol=1e-11)
    identity = donor_weights([(i,) for i in range(13)], 0)
    np.testing.assert_array_equal(identity, np.eye(13))
    try:
        donor_weights([(i,) for i in range(13)], 0.1)
    except ValueError:
        pass
    else:
        raise AssertionError("Unidentified positive contrast loss accepted")
    result = {"checks": {"graph_pair_identity": True, "intercept_preserved": True,
                         "positive_loss_weights": True, "primal_ridge_agreement": True,
                         "primal_penalized_rank_optimum": True, "donor_permutation": True,
                         "target_translation": True, "no_partner_rule": True},
              "ridge_max_abs_error": float(np.max(np.abs(observed - expected))),
              "rank2_max_abs_error": float(np.max(np.abs(observed_rank - expected_rank))),
              "seed": 20260930, "python": platform.python_version(), "numpy": np.__version__,
              "status": "synthetic_numerical_checks_only_no_model_accuracy"}
    out = Path(os.environ["CODEX_REC_OUTPUT"]) / f"weighted_chromatin_checks_{job}"
    out.mkdir(exist_ok=False)
    (out / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
