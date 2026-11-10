"""Cache alpha-independent ridge arithmetic within each inner training split."""
import numpy as np
import pandas as pd


def tune_path(core, designs, y, donors, train, inner, grid):
    """Use the original solve and donor-loss arithmetic for every ridge alpha.

    The supplied core defines the original equal-donor weights. Each split has
    its own already selected design; no feature selection or donor split changes.
    """
    losses = {alpha: [] for alpha in grid}
    for k in sorted(set(inner)):
        it, iv = train[inner != k], train[inner == k]
        design = designs[k]
        x, predict = design[it], design[iv]
        target = y[it]
        weights = core.donor_weights(donors[it])
        mean = np.average(x, axis=0, weights=weights) if x.shape[1] else np.empty(0)
        scale = np.sqrt(np.average((x - mean) ** 2, axis=0, weights=weights)) if x.shape[1] else np.empty(0)
        scale = np.maximum(scale, 1e-8)
        z, test = (x - mean) / scale, (predict - mean) / scale
        ym = np.average(target, weights=weights)
        if x.shape[1]:
            gram = z.T @ (weights[:, None] * z)
            rhs = z.T @ (weights * (target - ym))
            identity = np.eye(z.shape[1])
        for alpha in grid:
            if x.shape[1]:
                coefficient = np.linalg.solve(gram + alpha * identity, rhs)
                prediction = test @ coefficient + ym
            else:
                prediction = np.repeat(ym, len(predict))
            values = pd.DataFrame({"donor": donors[iv], "loss": (prediction - y[iv]) ** 2}).groupby("donor").loss.mean()
            losses[alpha].extend(values.tolist())
    return {alpha: float(np.mean(values)) for alpha, values in losses.items()}
