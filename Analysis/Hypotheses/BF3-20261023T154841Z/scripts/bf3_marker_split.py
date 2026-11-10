#!/usr/bin/env python
"""BF3 post hoc biliary versus stromal marker split (DEVIATIONS.md 2026-10-23T16:49Z). Descriptive."""
import json
import sys

import numpy as np
import pandas as pd

BUILD, PRIMARY, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CG0 = f"{ROOT}/Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z/genes.tsv"
SEED, NBOOT = 20261023, 2000
MARKERS = ["KRT7", "KRT19", "EPCAM", "SOX9", "COL1A1", "COL1A2", "COL3A1", "LUM", "DCN"]
DIAG_ORD = {"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2, "Cirrhosis": 3}
TERMS = ["COMP", "WITHIN", "INTER"]
rng = np.random.default_rng(SEED)

z = np.load(f"{BUILD}/donor_type_counts.npz", allow_pickle=True)  # own build output; genes saved as object array
counts, donors, types, genes = z["counts"], z["donors"].astype(str), z["types"].astype(str), z["genes"].astype(str)
dt = pd.read_csv(f"{BUILD}/donors.tsv", sep="\t").set_index("donor").loc[donors]
ok = (dt.n_nuclei >= 500) & dt.diagnosis.isin(list(DIAG_ORD)) & dt.sex.isin(["Male", "Female"])
counts, dt = counts[ok.values], dt[ok]
HEP, FIB, CHOL = [int(np.where(types == t)[0][0]) for t in ["Hepatocytes", "Fibroblasts", "Cholangiocytes"]]

g = pd.read_csv(CG0, sep="\t")
sig = g[g.signature_gene.fillna(False).astype(bool)]
gpos = {s: i for i, s in enumerate(genes)}
U = [gpos[s] for s in sig.gene_name[sig.cont_z > 0] if s in gpos]
M = [gpos[s] for s in MARKERS if s in gpos]

DS = {}
for ds in sorted(dt.dataset.unique()):
    m = (dt.dataset == ds).values
    d = dt[m]
    X = pd.get_dummies(d.diagnosis, drop_first=True, dtype=float)
    X.insert(0, "sex", (d.sex == "Male").astype(float).values)
    X.insert(0, "const", 1.0)
    DS[ds] = dict(n=counts[m], Ndc=counts[m].sum(2), cov=X.values,
                  diag=d.diagnosis.map(DIAG_ORD).values.astype(float), male=(d.sex == "Male").values.astype(float))


def decompose(nn, Ndc):
    """Same as bf3_decompose.py."""
    Nd = Ndc.sum(1)
    w = Ndc / Nd[:, None]
    x = nn.sum(1) / Nd[:, None]
    keep = x.mean(0) >= 1e-6
    nn, x = nn[:, :, keep], x[:, keep]
    xbar = x.mean(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = nn / Ndc[:, :, None]
    has = Ndc > 0
    pbar = np.where(has[:, :, None], p, 0).sum(0) / np.maximum(has.sum(0), 1)[:, None]
    p = np.where(has[:, :, None], p, pbar[None])
    a = (pbar / xbar).mean(1)
    b = ((p - pbar[None]) / xbar).mean(2)
    wbar = w.mean(0)
    S = (x / xbar).mean(1)
    T = {"COMP": (w - wbar) * a, "WITHIN": wbar * b, "INTER": (w - wbar) * b}
    assert np.allclose(S - (wbar * a).sum(), sum(t.sum(1) for t in T.values()), atol=1e-9)
    return S, T, int(keep.sum())


def resid(y, X):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def shares(idx, rows=None, by_type=False):
    num, den, per, ngenes = {k: 0.0 for k in TERMS}, 0.0, {}, {}
    num_t = {k: np.zeros(len(types)) for k in TERMS}
    for ds, D in DS.items():
        r = np.arange(len(D["n"])) if rows is None else rows[ds]
        cov = D["cov"][r]
        cov = cov[:, np.abs(cov).sum(0) > 0]
        S, T, ngenes[ds] = decompose(D["n"][:, :, idx][r], D["Ndc"][r])
        Sr = resid(S, cov)
        dden = (Sr ** 2).sum()
        per[ds] = {k: float((resid(T[k].sum(1), cov) * Sr).sum() / dden) for k in TERMS}
        den += dden
        for k in TERMS:
            num[k] += (resid(T[k].sum(1), cov) * Sr).sum()
            if by_type:
                num_t[k] += np.array([(resid(T[k][:, c], cov) * Sr).sum() for c in range(len(types))])
    out = {"pooled": {k: num[k] / den for k in TERMS}, "per_dataset": per, "n_genes_used": ngenes}
    if by_type:
        out["by_type"] = {k: dict(zip(types, (num_t[k] / den).round(5))) for k in TERMS}
    return out


SLOPE_COLS = ["S"] + TERMS + [f"{k}_{t}" for k in ["COMP", "WITHIN"] for t in types]


def slope_split(idx, rows=None):
    """Diagnosis-ordinal slope with dataset and sex; S, the three terms, and COMP and WITHIN by cell type."""
    Ys, Xs = [], []
    for j, (ds, D) in enumerate(DS.items()):
        r = np.arange(len(D["n"])) if rows is None else rows[ds]
        S, T, _ = decompose(D["n"][:, :, idx][r], D["Ndc"][r])
        Ys.append(np.column_stack([S] + [T[k].sum(1) for k in TERMS] + [T["COMP"], T["WITHIN"]]))
        Xs.append(np.column_stack([D["diag"][r], D["male"][r], np.full(len(r), j)]))
    X = np.vstack(Xs)
    X = np.column_stack([X[:, 0], X[:, 1], (X[:, 2][:, None] == np.arange(len(DS))).astype(float)])
    return dict(zip(SLOPE_COLS, np.linalg.lstsq(X, np.vstack(Ys), rcond=None)[0][0]))


def boot(fun, idx, **kw):
    draws = []
    for _ in range(NBOOT):
        rows = {ds: rng.integers(0, len(D["n"]), len(D["n"])) for ds, D in DS.items()}
        draws.append(fun(idx, rows=rows, **kw))
    return draws


def ci(vals):
    return [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))]


SETS = {"B": [gpos[s] for s in ["KRT7", "KRT19", "EPCAM", "SOX9"]],
        "St": [gpos[s] for s in ["COL1A1", "COL1A2", "COL3A1", "LUM", "DCN"]]}


def stats(rows=None):
    """Shares by type, slope by type, and the two contrasts, for B and St on the same donor draw."""
    out = {}
    for name, idx in SETS.items():
        sh = shares(idx, rows=rows, by_type=True)
        sl = slope_split(idx, rows=rows)
        out[name] = dict(share={k: float(v) for k, v in sh["pooled"].items()},
                         share_by_type={k: {t: float(v) for t, v in d.items()} for k, d in sh["by_type"].items()},
                         slope={k: float(v) for k, v in sl.items()}, n_genes_used=sh["n_genes_used"])
        out[name]["hep_within_share"] = out[name]["share_by_type"]["WITHIN"]["Hepatocytes"]
        out[name]["hep_slope_fraction"] = out[name]["slope"]["WITHIN_Hepatocytes"] / out[name]["slope"]["S"]
    out["contrast_hep_within_share"] = out["B"]["hep_within_share"] - out["St"]["hep_within_share"]
    out["contrast_hep_slope_fraction"] = out["B"]["hep_slope_fraction"] - out["St"]["hep_slope_fraction"]
    return out


pt = stats()
draws = []
for _ in range(NBOOT):
    rows = {ds: rng.integers(0, len(D["n"]), len(D["n"])) for ds, D in DS.items()}
    draws.append(stats(rows))

RES = {"genes": {k: genes[v].tolist() for k, v in SETS.items()}}
for name in SETS:
    x = pt[name]
    RES[name] = dict(
        n_genes_used=x["n_genes_used"],
        share={k: dict(est=v, ci=ci([d[name]["share"][k] for d in draws])) for k, v in x["share"].items()},
        share_by_type={k: {t: dict(est=v, ci=ci([d[name]["share_by_type"][k][t] for d in draws]))
                           for t, v in dd.items() if abs(v) >= 0.01} for k, dd in x["share_by_type"].items()},
        slope={k: dict(est=v, ci=ci([d[name]["slope"][k] for d in draws])) for k, v in x["slope"].items()
               if k in ["S"] + TERMS or abs(v) >= 0.01},
        hep_slope_fraction=dict(est=x["hep_slope_fraction"], ci=ci([d[name]["hep_slope_fraction"] for d in draws])))
for k in ["contrast_hep_within_share", "contrast_hep_slope_fraction"]:
    RES[k] = dict(est=pt[k], ci=ci([d[k] for d in draws]))

amb = {}
for name, idx in SETS.items():
    amb[name] = {}
    for src, cols in {"cholangiocyte_share": [CHOL], "fibroblast_share": [FIB]}.items():
        a_all, b_all, per = [], [], {}
        for ds, D in DS.items():
            S, T, _ = decompose(D["n"][:, :, idx], D["Ndc"])
            w = D["Ndc"] / D["Ndc"].sum(1)[:, None]
            hw = resid(T["WITHIN"][:, HEP], D["cov"])
            pw = resid(w[:, cols].sum(1), D["cov"])
            per[ds] = float(np.corrcoef(hw, pw)[0, 1])
            a_all.append(hw)
            b_all.append(pw)
        amb[name][src] = dict(per_dataset_r=per,
                              pooled_r=float(np.corrcoef(np.concatenate(a_all), np.concatenate(b_all))[0, 1]))
RES["hep_within_vs_compartment_share"] = amb

json.dump(RES, open(f"{OUT}/marker_split.json", "w"), indent=1, default=float)
print(json.dumps(RES, indent=1, default=float))
