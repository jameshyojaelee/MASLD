#!/usr/bin/env python
"""BF3 post hoc ambient-RNA check (DEVIATIONS.md 2026-10-23T16:00Z). Descriptive; cannot change the primary decision."""
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


RES = {}
# 1. positive-control set by cell type
RES["M_by_type"] = shares(M, by_type=True)["by_type"]

# 2. per-gene enrichment by per-UMI rate, pooled over all donors
Ndc = counts.sum(2)
has = Ndc > 0
with np.errstate(invalid="ignore", divide="ignore"):
    p = counts[:, :, U] / Ndc[:, :, None]
pbar = np.where(has[:, :, None], p, 0).sum(0) / has.sum(0)[:, None]
enr = pd.DataFrame(pbar.T, index=genes[U], columns=types)
enr = enr[enr.max(1) > 0]
enr["enriched_type"] = enr[list(types)].idxmax(1)
enr["hep_ratio"] = enr["Hepatocytes"] / enr[list(types)].max(1)
pg = pd.read_csv(f"{PRIMARY}/s5_per_gene.tsv", sep="\t").set_index("gene")
enr = enr.join(pg[["share_COMP", "share_WITHIN", "share_INTER", "within_top_type", "origin_top_type"]])
enr.to_csv(f"{OUT}/gene_enrichment.tsv", sep="\t")
RES["enriched_type_counts"] = enr.enriched_type.value_counts().to_dict()
d5 = enr.dropna(subset=["share_WITHIN"])
RES["s5_by_enriched_type"] = d5.groupby("enriched_type").agg(
    n=("share_WITHIN", "size"), median_share_WITHIN=("share_WITHIN", "median"),
    n_within_over_half=("share_WITHIN", lambda s: int((s > 0.5).sum())),
    n_within_over_half_carried_by_hep=("within_top_type", lambda s: int(
        ((s == "Hepatocytes") & (d5.loc[s.index, "share_WITHIN"] > 0.5)).sum()))).to_dict(orient="index")

# 3. hepatocyte-enriched and other rising genes as separate sets
SETS = {"U": U, "U_hep_enriched": [gpos[s] for s in enr.index[enr.enriched_type == "Hepatocytes"]],
        "U_other": [gpos[s] for s in enr.index[enr.enriched_type != "Hepatocytes"]], "M": M}
for name in ["U_hep_enriched", "U_other"]:
    pt = shares(SETS[name], by_type=True)
    bs = boot(shares, SETS[name])
    pt["pooled_ci"] = {k: ci([b["pooled"][k] for b in bs]) for k in TERMS}
    pt["n_genes"] = len(SETS[name])
    RES[name] = pt

# 4. does the hepatocyte within term track portal-tract abundance?
amb = {}
for name, idx in SETS.items():
    a_all, b_all, per = [], [], {}
    for ds, D in DS.items():
        S, T, _ = decompose(D["n"][:, :, idx], D["Ndc"])
        w = D["Ndc"] / D["Ndc"].sum(1)[:, None]
        hw = resid(T["WITHIN"][:, HEP], D["cov"])
        pw = resid(w[:, FIB] + w[:, CHOL], D["cov"])
        per[ds] = float(np.corrcoef(hw, pw)[0, 1])
        a_all.append(hw)
        b_all.append(pw)
    amb[name] = dict(per_dataset_r=per, pooled_r=float(np.corrcoef(np.concatenate(a_all), np.concatenate(b_all))[0, 1]),
                     n=int(sum(len(x) for x in a_all)))
RES["hep_within_vs_portal_share"] = amb

# 5. diagnosis slope by cell type
for name in ["U", "U_hep_enriched", "U_other"]:
    sl = slope_split(SETS[name])
    bs = pd.DataFrame(boot(slope_split, SETS[name]))
    RES[f"slope_{name}"] = {k: dict(est=float(sl[k]), ci=ci(bs[k])) for k in SLOPE_COLS
                            if k in ["S"] + TERMS or abs(sl[k]) >= 0.01}

json.dump(RES, open(f"{OUT}/posthoc.json", "w"), indent=1, default=float)
print(json.dumps(RES, indent=1, default=float))
