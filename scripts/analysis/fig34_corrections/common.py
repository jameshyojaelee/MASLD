"""Donor-level statistics for the Figures 3–4 correction candidate."""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
SEED = 20260929
FOCAL = {
    "hotspot_hepatocytes_f05c535ae5bbc0b9": "Stromal ECM",
    "hotspot_hepatocytes_48f39dd4d817a10e": "Ductular injury",
}
HAC = ROOT / "RNA-seq/results/histology_anchored_continuum"
ML = HAC / "molecular_layers/hac-molecular-layers-20260818T173348Z"
AXES = HAC / "candidates/hac-continuum-20260818T024923Z/projection/score_registry.tsv"
PROGRAM = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot"


def read(path):
    return pd.read_csv(path, sep="\t")


def write(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(path)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def bh(values, family):
    values = np.asarray(values, dtype=float)
    result = np.full(len(values), np.nan)
    idx = np.flatnonzero(np.isfinite(values))
    assert len(idx) <= family
    order = idx[np.argsort(values[idx])]
    if len(order):
        result[order] = np.minimum(1, np.minimum.accumulate(
            (values[order] * family / np.arange(1, len(order) + 1))[::-1])[::-1])
    return result


def design(d, stage="stage_group", composition=False):
    parts = [np.ones((len(d), 1)), d[["axis_z"]].to_numpy(float)]
    for col in [stage, "inferred_sex"]:
        parts.append(pd.get_dummies(d[col].astype(str), drop_first=True).to_numpy(float))
    if composition:
        parts.append(d.filter(regex="^ilr_").to_numpy(float))
    return np.column_stack(parts)


def fit(y, x, term=1):
    y, x = np.asarray(y, float), np.asarray(x, float)
    result = dict(n=len(y), estimable=False, beta=np.nan, se=np.nan, p=np.nan,
                  low=np.nan, high=np.nan, hc3_se=np.nan, hc3_p=np.nan,
                  hc3_low=np.nan, hc3_high=np.nan, df=len(y)-x.shape[1],
                  failure_reason="nonfinite_or_singular_design")
    if not np.isfinite(y).all() or not np.isfinite(x).all() or len(y) <= x.shape[1]:
        return result
    if np.linalg.matrix_rank(x) != x.shape[1]:
        return result
    inv = np.linalg.inv(x.T @ x)
    beta = inv @ x.T @ y
    residual = y - x @ beta
    leverage = np.einsum("ij,jk,ik->i", x, inv, x)
    if np.max(leverage) >= 1 - 1e-10:
        result["failure_reason"] = "unit_leverage"
        return result
    df = result["df"]
    se = np.sqrt(np.sum(residual**2) / df * inv[term, term])
    meat = (x * (residual / (1 - leverage))[:, None]).T @ (x * (residual / (1 - leverage))[:, None])
    robust = np.sqrt(max(0, (inv @ meat @ inv)[term, term]))
    b = beta[term]
    critical = stats.t.ppf(.975, df)
    result.update(estimable=True, beta=b, se=se, p=2*stats.t.sf(abs(b/se), df),
                  low=b-critical*se, high=b+critical*se,
                  hc3_se=robust, hc3_p=2*stats.t.sf(abs(b/robust), df),
                  hc3_low=b-critical*robust, hc3_high=b+critical*robust,
                  max_leverage=float(max(leverage)), failure_reason="")
    return result


def meta(d, expected=2, robust=False):
    secol = "hc3_se" if robust else "se"
    z = d[d.estimable & np.isfinite(d[secol]) & (d[secol] > 0)]
    result = dict(n=int(z.n.sum()), n_cohorts=len(z), estimable=False,
                  beta=np.nan, se=np.nan, p=np.nan, low=np.nan, high=np.nan,
                  direction_consistent=False)
    if len(z) != expected or z.dataset.nunique() != expected:
        return result
    w = 1/z[secol]**2
    b = np.sum(w*z.beta)/sum(w)
    se = np.sqrt(1/sum(w))
    result.update(estimable=True, beta=b, se=se, p=2*stats.norm.sf(abs(b/se)),
                  low=b-1.96*se, high=b+1.96*se,
                  direction_consistent=bool((np.sign(z.beta)==np.sign(b)).all()))
    return result


def snapshot(paths, out):
    rows = []
    for p in paths:
        p = Path(p)
        digest = hashlib.sha256()
        with p.open("rb") as handle:
            for block in iter(lambda: handle.read(1024*1024), b""):
                digest.update(block)
        rows.append(dict(path=str(p.resolve()), sha256=digest.hexdigest()))
    write(pd.DataFrame(rows), out)


def plot_setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 6, "axes.labelsize": 6, "xtick.labelsize": 6,
                        "ytick.labelsize": 6, "legend.fontsize": 6,
                        "axes.titlesize": 6, "pdf.fonttype": 42,
                        "axes.spines.top": False, "axes.spines.right": False})
    return plt
