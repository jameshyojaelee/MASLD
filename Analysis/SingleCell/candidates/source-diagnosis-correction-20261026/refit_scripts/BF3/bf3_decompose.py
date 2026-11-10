#!/usr/bin/env python
"""BF3 step 2: composition / within-cell-type decomposition of continuum gene sets (PRESPEC_BF3.md)."""
import json
import sys

import numpy as np
import pandas as pd

BUILD, OUT = sys.argv[1], sys.argv[2]
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CG0 = f"{ROOT}/Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z/genes.tsv"
SEED, NBOOT, NNULL = 20261023, 2000, 1000
MARKERS = ["KRT7", "KRT19", "EPCAM", "SOX9", "COL1A1", "COL1A2", "COL3A1", "LUM", "DCN"]
DIAG_ORD = {"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2, "Cirrhosis": 3}
TERMS = ["COMP", "WITHIN", "INTER"]
rng = np.random.default_rng(SEED)

z = np.load(f"{BUILD}/donor_type_counts.npz", allow_pickle=True)  # own build output; genes saved as object array
counts, donors, types, genes = z["counts"], z["donors"].astype(str), z["types"].astype(str), z["genes"].astype(str)
dt = pd.read_csv(f"{BUILD}/donors.tsv", sep="\t").set_index("donor").loc[donors]
ok = (dt.n_nuclei >= 500) & dt.diagnosis.isin(list(DIAG_ORD)) & dt.sex.isin(["Male", "Female"])
RES = {"excluded_donors": dt.index[~ok].tolist()}
counts, dt = counts[ok.values], dt[ok]
RES["n_donors"] = dt.groupby(["dataset", "diagnosis"]).size().unstack(fill_value=0).to_dict(orient="index")

g = pd.read_csv(CG0, sep="\t")
sig = g[g.signature_gene.fillna(False).astype(bool)]
gpos = {s: i for i, s in enumerate(genes)}
U = [gpos[s] for s in sig.gene_name[sig.cont_z > 0] if s in gpos]
F = [gpos[s] for s in sig.gene_name[sig.cont_z < 0] if s in gpos]
M = [gpos[s] for s in MARKERS if s in gpos]
sig_idx = {gpos[s] for s in sig.gene_name if s in gpos}
RES["n_genes_matched"] = {"U": len(U), "U_signature": int((sig.cont_z > 0).sum()), "F": len(F), "M": len(M)}

DS = {}
for ds in sorted(dt.dataset.unique()):
    m = (dt.dataset == ds).values
    n = counts[m]
    Ndc = n.sum(2)
    d = dt[m]
    X = pd.get_dummies(d.diagnosis, drop_first=True, dtype=float)
    X.insert(0, "sex", (d.sex == "Male").astype(float).values)
    X.insert(0, "const", 1.0)
    DS[ds] = dict(n=n, Ndc=Ndc, cov=X.values, diag=d.diagnosis.map(DIAG_ORD).values.astype(float),
                  male=(d.sex == "Male").values.astype(float), xg=(n.sum(1) / Ndc.sum(1)[:, None]).mean(0))


def decompose(nn, Ndc):
    """Exact split S_d - R = COMP + WITHIN + INTER, each donors x types. nn = donors x types x set genes."""
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
    R = (wbar * a).sum()
    assert np.allclose(S - R, sum(t.sum(1) for t in T.values()), atol=1e-9)
    return S, T, int(keep.sum())


def resid(y, X):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def shares(idx, rows=None, by_type=False, drop_cirrhosis=False):
    """Pooled within-diagnosis variance shares; rows = per-dataset bootstrap index arrays."""
    num, den, per = {k: 0.0 for k in TERMS}, 0.0, {}
    num_t = {k: np.zeros(len(types)) for k in TERMS}
    ngenes = {}
    for ds, D in DS.items():
        r = np.arange(len(D["n"])) if rows is None else rows[ds]
        if drop_cirrhosis:
            r = r[D["diag"][r] != 3]
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


def slope_split(idx, rows=None):
    """Diagnosis-ordinal slope of S with dataset and sex; exact split into term slopes."""
    S_all, T_all, Xs = [], {k: [] for k in TERMS}, []
    for j, (ds, D) in enumerate(DS.items()):
        r = np.arange(len(D["n"])) if rows is None else rows[ds]
        S, T, _ = decompose(D["n"][:, :, idx][r], D["Ndc"][r])
        S_all.append(S)
        for k in TERMS:
            T_all[k].append(T[k].sum(1))
        Xs.append(np.column_stack([D["diag"][r], D["male"][r], np.full(len(r), j)]))
    X = np.vstack(Xs)
    X = np.column_stack([X[:, 0], X[:, 1], (X[:, 2][:, None] == np.arange(len(DS))).astype(float)])
    co = lambda y: np.linalg.lstsq(X, y, rcond=None)[0][0]
    out = {"S": co(np.concatenate(S_all))}
    out.update({k: co(np.concatenate(T_all[k])) for k in TERMS})
    return out


def boot(fun, idx, **kw):
    draws = []
    for _ in range(NBOOT):
        rows = {ds: rng.integers(0, len(D["n"]), len(D["n"])) for ds, D in DS.items()}
        draws.append(fun(idx, rows=rows, **kw))
    return draws


def ci(vals):
    return [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))]


# ---------------------------------------------------------------- primary and positive control
for name, idx in [("U", U), ("M_positive_control", M)]:
    pt = shares(idx, by_type=(name == "U"))
    bs = boot(shares, idx)
    pt["pooled_ci"] = {k: ci([b["pooled"][k] for b in bs]) for k in TERMS}
    RES[name] = pt
rows_out = []
for name, idx in [("U", U), ("M", M), ("F", F)]:
    for ds, D in DS.items():
        S, T, _ = decompose(D["n"][:, :, idx], D["Ndc"])
        df = pd.DataFrame({"set": name, "dataset": ds, "donor": dt.index[(dt.dataset == ds).values], "S": S,
                           **{k: T[k].sum(1) for k in TERMS}})
        for k in ["COMP", "WITHIN"]:
            for c, t in enumerate(types):
                df[f"{k}_{t}"] = T[k][:, c]
        rows_out.append(df)
pd.concat(rows_out).merge(dt[["diagnosis", "sex", "n_nuclei"]], left_on="donor", right_index=True).to_csv(
    f"{OUT}/donor_terms.tsv", sep="\t", index=False)
u = RES["U"]
lo = {k: u["pooled_ci"][k][0] for k in TERMS}
RES["U_decision"] = ("mostly composition" if lo["COMP"] > 0.5 else
                     "mostly within-cell-type" if lo["WITHIN"] > 0.5 else "mixed")
lead = max(TERMS, key=lambda k: u["pooled"][k])
RES["U_consistent_across_datasets"] = bool(all((v[lead] > 0.5) == (u["pooled"][lead] > 0.5)
                                               for v in u["per_dataset"].values()))
RES["positive_control_passed"] = bool(RES["M_positive_control"]["pooled"]["COMP"] > 0.5)

# ---------------------------------------------------------------- matched random sets
xg = np.mean([D["xg"] for D in DS.values()], axis=0)
cand = np.array([i for i in range(len(genes)) if i not in sig_idx and xg[i] >= 1e-6])
Udet = [i for i in U if xg[i] >= 1e-6]
edges = np.quantile(xg[cand], np.linspace(0, 1, 11))
ubin = np.clip(np.searchsorted(edges, xg[Udet], side="right") - 1, 0, 9)
cbin = np.clip(np.searchsorted(edges, xg[cand], side="right") - 1, 0, 9)
need = pd.Series(ubin).value_counts()
null = []
for _ in range(NNULL):
    idx = np.concatenate([rng.choice(cand[cbin == b_], k, replace=False) for b_, k in need.items()])
    null.append(shares(list(idx))["pooled"])
null = pd.DataFrame(null)
null.to_csv(f"{OUT}/null_matched_shares.tsv", sep="\t", index=False)
RES["null_set_size"] = len(Udet)
RES["U_vs_matched"] = {k: dict(null_median=float(null[k].median()), null_ci=ci(null[k]),
                               pct_of_U=float((null[k] < u["pooled"][k]).mean())) for k in TERMS}

# ---------------------------------------------------------------- secondaries
sl = slope_split(U)
bs = boot(slope_split, U)
RES["S2_diagnosis_slope_U"] = {k: dict(est=float(sl[k]), ci=ci([b[k] for b in bs])) for k in ["S"] + TERMS}
RES["S2_diagnosis_slope_U"]["share_of_slope"] = {k: float(sl[k] / sl["S"]) for k in TERMS}
RES["S3_F"] = shares(F)
RES["S6_U_no_cirrhosis"] = shares(U, drop_cirrhosis=True)

Uboth = [i for i in U if all(D["xg"][i] >= 1e-6 for D in DS.values())]
RES["n_U_detected_both"] = len(Uboth)
allc = counts.sum(0)
origin = pd.DataFrame(allc[:, Uboth].T / allc[:, Uboth].sum(0)[:, None], index=genes[Uboth], columns=types)
origin["top_type"] = origin[list(types)].idxmax(1)
origin.to_csv(f"{OUT}/s4_gene_origin.tsv", sep="\t")
RES["S4_top_type_counts"] = origin.top_type.value_counts().to_dict()

pg = []
for i in Uboth:
    s = shares([i], by_type=True)
    wt = s["by_type"]["WITHIN"]
    pg.append(dict(gene=genes[i], **{f"share_{k}": s["pooled"][k] for k in TERMS},
                   within_top_type=max(wt, key=wt.get), within_top_value=max(wt.values()),
                   origin_top_type=origin.loc[genes[i], "top_type"]))
pg = pd.DataFrame(pg).sort_values("share_WITHIN", ascending=False)
pg.to_csv(f"{OUT}/s5_per_gene.tsv", sep="\t", index=False)
RES["S5_per_gene"] = dict(n=len(pg), median_share_COMP=float(pg.share_COMP.median()),
                          median_share_WITHIN=float(pg.share_WITHIN.median()),
                          n_COMP_over_half=int((pg.share_COMP > 0.5).sum()),
                          n_WITHIN_over_half=int((pg.share_WITHIN > 0.5).sum()),
                          within_top_type_counts=pg[pg.share_WITHIN > 0.5].within_top_type.value_counts().to_dict())
json.dump(RES, open(f"{OUT}/results.json", "w"), indent=1, default=float)
print(json.dumps(RES, indent=1, default=float))
