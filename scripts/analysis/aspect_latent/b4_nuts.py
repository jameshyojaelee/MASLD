#!/usr/bin/env python3
"""B4: final reported fit (PRESPEC_B section 4 and Amendment 2): NUTS at the kept K.

All seven development cohorts, primary (ALT-corrected) inputs unless --inputs says
otherwise. 4 chains x 1,000 warmup + 1,000 draws, seed 20260923. The fit is accepted
only with R-hat < 1.01 for every sampled site and no divergent transition; otherwise it
is reported as not converged and nothing downstream reads it.

Writes: posterior draws of the global sites (npz), per-site R-hat / n_eff, divergences,
each label's coefficient vector W @ lambda_j over the 117 programs (posterior mean and
95% interval), and the singular values of W. These do not depend on the rotation.
"""
import argparse
import json
import sys
from pathlib import Path

import numpyro

numpyro.set_host_device_count(4)
import jax  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from numpyro.diagnostics import summary  # noqa: E402
from numpyro.infer import MCMC, NUTS  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit  # noqa: E402
from b3_tests import DEV, check_prespec  # noqa: E402
from model import LABELS, aspect_model  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=2)
    ap.add_argument("--warmup", type=int, default=1000)
    ap.add_argument("--draws", type=int, default=1000)
    a_ = ap.parse_args()
    check_prespec()
    out = Path(a_.out); out.mkdir(parents=True, exist_ok=True)
    d, prog = fit.load_inputs(a_.inputs)
    d = d[d["cohort"].isin(DEV)].reset_index(drop=True)
    a = fit.arrays(d, prog, DEV)
    mcmc = MCMC(NUTS(aspect_model), num_warmup=a_.warmup, num_samples=a_.draws, num_chains=4,
                chain_method="parallel", progress_bar=False)
    mcmc.run(jax.random.PRNGKey(fit.SEED), k=a_.k, extra_fields=("diverging",), **a)
    samples = mcmc.get_samples(group_by_chain=True)
    div = int(np.asarray(mcmc.get_extra_fields()["diverging"]).sum())
    diag = summary({k: v for k, v in samples.items() if k != "z"}, group_by_chain=True)
    rows = []
    for site, st in diag.items():
        rhat, neff = np.atleast_1d(st["r_hat"]), np.atleast_1d(st["n_eff"])
        rows.append({"site": site, "size": int(rhat.size), "max_r_hat": float(np.nanmax(rhat)),
                     "min_n_eff": float(np.nanmin(neff))})
    diag_df = pd.DataFrame(rows).sort_values("max_r_hat", ascending=False)
    diag_df.to_csv(out / "nuts_diagnostics.tsv", sep="\t", index=False)
    max_rhat = float(diag_df["max_r_hat"].max())
    converged = bool(max_rhat < 1.01 and div == 0)

    flat = {k: np.asarray(v).reshape((-1,) + np.asarray(v).shape[2:]) for k, v in samples.items() if k not in ("eps", "z")}
    np.savez_compressed(out / f"nuts_K{a_.k}_global_draws.npz", **flat)
    np.save(out / f"nuts_K{a_.k}_z_mean.npy", np.asarray(samples["z"]).reshape((-1,) + np.asarray(samples["z"]).shape[2:]).mean(0))
    W = flat["W"]  # draws x P x K
    lam = np.asarray(fit._lam({k: v for k, v in flat.items()}, a_.k))  # draws x labels x K
    coef = []
    for j, name in enumerate(LABELS):
        c = np.einsum("dpk,dk->dp", W, lam[:, j, :])
        coef.append(pd.DataFrame({"label": name, "program": prog, "mean": c.mean(0),
                                  "lo95": np.percentile(c, 2.5, 0), "hi95": np.percentile(c, 97.5, 0)}))
    coef = pd.concat(coef)
    coef["excludes_zero"] = (coef["lo95"] > 0) | (coef["hi95"] < 0)
    coef.to_csv(out / f"label_coefficients_K{a_.k}.tsv", sep="\t", index=False)
    sv = np.linalg.svd(W, compute_uv=False)  # draws x K
    between = {}
    for i, li in enumerate(LABELS):
        for j2, lj in enumerate(LABELS):
            if j2 > i:
                ci = np.einsum("dpk,dk->dp", W, lam[:, i, :]); cj = np.einsum("dpk,dk->dp", W, lam[:, j2, :])
                cos = (ci * cj).sum(1) / np.linalg.norm(ci, axis=1) / np.linalg.norm(cj, axis=1)
                between[f"{li}|{lj}"] = [float(np.mean(cos)), float(np.percentile(cos, 2.5)), float(np.percentile(cos, 97.5))]
    summ = {"k": a_.k, "participants": int(len(d)), "divergences": div, "max_r_hat": max_rhat, "converged": converged,
            "singular_values_mean": sv.mean(0).tolist(), "singular_values_95": np.percentile(sv, [2.5, 97.5], 0).tolist(),
            "label_coefficient_cosine_mean_and_95": between,
            "programs_excluding_zero": coef.groupby("label")["excludes_zero"].sum().to_dict()}
    (out / f"nuts_K{a_.k}_summary.json").write_text(json.dumps(summ, indent=2))
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main()
