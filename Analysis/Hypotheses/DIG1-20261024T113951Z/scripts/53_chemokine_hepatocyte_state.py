"""DIG1 job 53 (exploratory, post hoc to job 52). What are the CCL20+ (and CXCL10+, SOX9+) hepatocytes?
Job 52: in single livers CCL20 tracks RELB and SOX9 beyond the injury score, but CCL20+ hepatocytes carry
2-5x the UMIs of CCL20- hepatocytes, so raw co-detection rates are depth-confounded.
Here, per donor, anchor+ hepatocytes are compared with anchor- hepatocytes of the same donor matched on
depth (within-donor percentile of log UMI and log genes; 10 nearest neighbours with replacement) and,
second mode, on depth plus the depth-residualised KRT8+SQSTM1 injury score.
Per donor and mode: pseudobulk log2 CPM (prior 1) of anchor+ cells minus the matched set. The same is done
for 20 random anchor genes per anchor (cohort hepatocyte detection 0.67-1.5x the anchor's; listed, score
and zonation genes excluded; seed 20261024); excess = anchor LFC - mean random-anchor LFC.
Cohorts: GSE298719 Flex hepatocytes (anchors CCL20, CXCL10, SOX9) and GSE256398 snRNA hepatocytes
(CCL20, SOX9; job 45L labels). A donor counts for an anchor when >= 10 hepatocytes detect it.
Primary family (Holm over 2), exact two-sided sign-flip over eligible donors, depth-matched mode:
  Z1 Flex CCL20+ zonation excess, (portal - central score) difference vs matched, minus random anchors.
  Z2 snRNA CCL20+ zonation excess, same.
  Portal = CPS1, HAL, ASS1, SDS, PCK1, ALDOB; central = CYP2E1, GLUL, CYP1A2, OAT, CYP3A4; scores are the
  mean within-donor z of log1p CP10k. Direction not fixed in advance (periportal if the cells sit with
  the ductular reaction; pericentral if they are zone-3 injured hepatocytes).
Descriptive: injury-score difference; every gene's LFC and excess with sign-flip p and BH within cohort x
anchor x mode; depth-matched detection rates for listed genes (doublet and spillover markers included);
CCL20 / SOX9 genes concordant across the two cohorts.
Usage: python -I 53_chemokine_hepatocyte_state.py <new_out_dir>
"""
import sys
from itertools import product
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.spatial import cKDTree

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
V5 = ROOT / "Analysis/Hypotheses/V5-20261021T122723Z"
V3PB = ROOT / "Analysis/Hypotheses/V3-20261008T015239Z/pseudobulk"
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
K_MATCH, N_RAND, MIN_POS = 10, 20, 10
SCORE = ["KRT8", "SQSTM1"]
ZON = {"central": ["CYP2E1", "GLUL", "CYP1A2", "OAT", "CYP3A4"], "portal": ["CPS1", "HAL", "ASS1", "SDS", "PCK1", "ALDOB"]}
ANCHORS = {"flex": ["CCL20", "CXCL10", "SOX9"], "snrna": ["CCL20", "SOX9"]}
LISTED = ["CCL20", "CXCL10", "CXCL9", "CXCL11", "CXCL1", "CXCL2", "CXCL8", "CCL2", "SOX9", "RELB", "NFKB2", "NFKBIA", "TNFAIP3", "BIRC3",
          "MAP3K14", "LTBR", "TNFRSF12A", "IL32", "ICAM1", "STAT1", "IRF1", "CD24", "EPCAM", "KRT7", "KRT19", "KRT23", "SPP1", "AKR1B10",
          "LCN2", "SAA1", "SAA2", "CRP", "GDF15", "FGF21", "CCN1", "CCN2", "AMOTL2", "ANKRD1", "HES1", "JAG1", "NOTCH2", "ONECUT1", "HNF1B",
          "PROM1", "ANXA4", "ALB", "APOA1", "HNF4A", "CEBPA", "CYP3A4", "PTPRC", "CD68", "CD3E", "PECAM1", "COL1A1", "HLA-DRA", "CD74"]
PRIMARY = {"flex": "Z1", "snrna": "Z2"}


def signflip(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) < 3:
        return np.nan
    S = np.array(list(product([1, -1], repeat=len(v))), dtype=np.int8)
    return np.mean(np.abs((S * v).mean(1)) >= abs(v.mean()) - 1e-12)


def bh(p):
    p = np.asarray(p, float); out = np.full(len(p), np.nan); ok = np.isfinite(p)
    if ok.any():
        q = p[ok]; o = np.argsort(q); r = q[o] * len(q) / (np.arange(len(q)) + 1)
        out[np.where(ok)[0][o]] = np.minimum(1, np.minimum.accumulate(r[::-1])[::-1])
    return out


def pct(x):
    return (stats.rankdata(x) - 0.5) / len(x)


def cohort(tag, X, gnames, donor, dis_of, seed):
    rng = np.random.default_rng(seed)
    X = sparse.csr_matrix(X); genes = {}
    for i, g in enumerate(gnames):
        genes.setdefault(g, i)
    umi = np.asarray(X.sum(axis=1)).ravel(); ngene = np.asarray((X > 0).sum(axis=1)).ravel()
    det_all = np.asarray((X > 0).mean(axis=0)).ravel()
    zon = {k: [x for x in v if x in genes] for k, v in ZON.items()}
    excl = set(LISTED) | set(SCORE) | set(x for v in zon.values() for x in v)
    gsel = np.where((det_all >= 0.01) | np.isin(gnames, LISTED))[0]
    anchors = {}
    for a in [x for x in ANCHORS[tag] if x in genes]:
        da = det_all[genes[a]]
        cand = [g for g, d in zip(gnames, det_all) if 0.67 * da <= d <= 1.5 * da and g not in excl]
        anchors[a] = list(rng.choice(cand, size=min(N_RAND, len(cand)), replace=False))
    listed = [g for g in LISTED if g in genes]
    side = sorted(set(listed) | set(SCORE) | set(x for v in zon.values() for x in v))
    Ls = np.log1p(X[:, [genes[g] for g in side]].toarray() * (1e4 / umi[:, None])); si = {g: i for i, g in enumerate(side)}
    gene_rows, cell_rows, det_rows = [], [], []
    for d in sorted(np.unique(donor)):
        m = np.where(donor == d)[0]; Xd = X[m]
        S = Ls[m]; sd = S.std(0); sd[sd == 0] = np.nan; Zs = (S - S.mean(0)) / sd
        sc = lambda gs: np.nan_to_num(np.nanmean(Zs[:, [si[g] for g in gs]], 1))
        inj = sc(SCORE)
        B = np.column_stack([np.ones(len(m)), stats.rankdata(np.log(umi[m])), stats.rankdata(np.log(ngene[m]))])
        inj_r = inj - B @ np.linalg.lstsq(B, inj, rcond=None)[0]
        zdiff = sc(zon["portal"]) - sc(zon["central"])
        feats = {"depth": np.column_stack([pct(np.log(umi[m])), pct(np.log(ngene[m]))]),
                 "depth+injury": np.column_stack([pct(np.log(umi[m])), pct(np.log(ngene[m])), pct(inj_r)])}
        Xg = Xd[:, gsel].tocsr()
        Db = (Xd[:, [genes[g] for g in listed]] > 0).toarray()
        tot_d = np.asarray(Xd.sum(axis=1)).ravel()
        for a, rands in anchors.items():
            for anc in [a] + rands:
                pos = np.asarray((Xd[:, genes[anc]] > 0).todense()).ravel()
                if pos.sum() < MIN_POS or (~pos).sum() < 5 * MIN_POS:
                    continue
                P = np.where(pos)[0]; N = np.where(~pos)[0]
                for mode, F in feats.items():
                    tree = cKDTree(F[N]); _, nn = tree.query(F[P], k=K_MATCH)
                    w = np.bincount(N[nn.ravel()], minlength=len(m)).astype(float)
                    sp = np.asarray(Xg[P].sum(axis=0)).ravel(); sn = Xg.T @ w
                    lfc = np.log2(sp / tot_d[P].sum() * 1e6 + 1) - np.log2(sn / (w @ tot_d) * 1e6 + 1)
                    gene_rows.append(dict(anchor=a, anc=anc, mode=mode, donor=d, lfc=lfc))
                    cell_rows.append(dict(cohort=tag, donor=d, disease=dis_of[d], anchor=a, anc=anc, is_anchor=anc == a, mode=mode, n_pos=len(P),
                                          zon_diff=zdiff[P].mean() - (w @ zdiff) / w.sum(), inj_diff=inj_r[P].mean() - (w @ inj_r) / w.sum(),
                                          logumi_diff=np.log(umi[m][P]).mean() - (w @ np.log(umi[m])) / w.sum()))
                    if anc == a:
                        for j, g in enumerate(listed):
                            det_rows.append(dict(cohort=tag, donor=d, disease=dis_of[d], anchor=a, mode=mode, gene=g, n_pos=len(P),
                                                 det_pos=Db[P, j].mean(), det_matched=(w @ Db[:, j]) / w.sum(), det_all_other=Db[N, j].mean()))
        print(tag, d, len(m), flush=True)
    C = pd.DataFrame(cell_rows)
    gsym = np.array(gnames)[gsel]
    out = []
    for (a, mode), grp in pd.DataFrame(gene_rows).groupby(["anchor", "mode"]):
        dons = sorted(grp[grp.anc == a].donor)
        A_ = np.vstack([grp[(grp.anc == a) & (grp.donor == d)].lfc.iloc[0] for d in dons])
        Nl = np.vstack([np.mean(np.vstack(grp[(grp.anc != a) & (grp.donor == d)].lfc.tolist()), 0) if ((grp.anc != a) & (grp.donor == d)).any()
                        else np.full(len(gsel), np.nan) for d in dons])
        Ex = A_ - Nl
        dis = np.array([dis_of[d] for d in dons])
        for grp_name, sel in [("all", np.ones(len(dons), bool)), ("masld", dis == 1)]:
            if sel.sum() < 3:
                continue
            p = np.array([signflip(Ex[sel, j]) for j in range(Ex.shape[1])])
            out.append(pd.DataFrame(dict(cohort=tag, anchor=a, mode=mode, donors=grp_name, gene=gsym, n_donors=int(sel.sum()),
                                         mean_lfc=A_[sel].mean(0), mean_null_lfc=np.nanmean(Nl[sel], 0), mean_excess=np.nanmean(Ex[sel], 0),
                                         n_pos_excess=(Ex[sel] > 0).sum(0), p_signflip=p, q_bh=bh(p), listed=np.isin(gsym, LISTED))))
    return C, pd.concat(out, ignore_index=True), pd.DataFrame(det_rows), anchors


# Flex
A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad")
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str}); dq = dq[dq.included]
A = A[(A.obs.label == "Hepatocytes").to_numpy() & A.obs.gsm.isin(dq.gsm).to_numpy()].copy()
cf, gf, df_, af = cohort("flex", A.X, list(A.var_names), A.obs.gsm.astype(str).to_numpy(), dict(zip(dq.gsm, dq.disease)), 20261024)
del A
# snRNA
H = ad.read_h5ad(RUN / "out_snrna_labels/gse256398_hepatocyte_counts.h5ad")
md = pd.read_csv(V3PB / "donor_metadata.tsv", sep="\t")
cs, gs, ds, as_ = cohort("snrna", H.X, list(H.var_names), H.obs.donor.astype(str).to_numpy(), dict(zip(md.donor.astype(str), md.disease)), 20261025)
del H

C = pd.concat([cf, cs], ignore_index=True); G = pd.concat([gf, gs], ignore_index=True); D = pd.concat([df_, ds], ignore_index=True)
C.to_csv(OUT / "donor_anchor_zonation_injury.tsv", sep="\t", index=False)
G.to_csv(OUT / "gene_lfc_matched.tsv.gz", sep="\t", index=False)
D.to_csv(OUT / "listed_detection_matched.tsv", sep="\t", index=False)
pd.DataFrame([dict(cohort=c, anchor=a, random_anchors=",".join(r)) for c, an in [("flex", af), ("snrna", as_)] for a, r in an.items()]).to_csv(
    OUT / "random_anchors.tsv", sep="\t", index=False)
# zonation and injury summaries
rows = []
for (c, a, mode), o in C.groupby(["cohort", "anchor", "mode"]):
    for var in ["zon_diff", "inj_diff", "logumi_diff"]:
        piv = o.pivot_table(index="donor", columns="is_anchor", values=var, aggfunc="mean")
        piv = piv.dropna(subset=[True]) if True in piv.columns else piv.iloc[0:0]
        if piv.empty:
            continue
        ex = piv[True] - piv.get(False, np.nan)
        for grp, sel in [("all", np.ones(len(piv), bool)), ("masld", np.array([dict(zip(o.donor, o.disease))[x] == 1 for x in piv.index]))]:
            rows.append(dict(cohort=c, anchor=a, mode=mode, var=var, donors=grp, n_donors=int(sel.sum()), mean_anchor=piv[True][sel].mean(),
                             mean_random=piv.get(False, pd.Series(np.nan, index=piv.index))[sel].mean(), mean_excess=ex[sel].mean(),
                             n_pos=int((ex[sel] > 0).sum()), p_signflip=signflip(ex[sel].to_numpy())))
Z = pd.DataFrame(rows); Z.to_csv(OUT / "summary_zonation_injury.tsv", sep="\t", index=False)
prim = Z[(Z.anchor == "CCL20") & (Z["mode"] == "depth") & (Z["var"] == "zon_diff") & (Z.donors == "all")].copy()
prim.insert(0, "id", prim.cohort.map(PRIMARY))
o = np.argsort(prim.p_signflip.to_numpy()); ph = np.maximum.accumulate(prim.p_signflip.to_numpy()[o] * (len(o) - np.arange(len(o))))
prim.loc[prim.index[o], "p_holm"] = np.minimum(1, ph)
prim.to_csv(OUT / "primary_holm.tsv", sep="\t", index=False)
# cross-cohort concordance for shared anchors
X2 = []
for a in ["CCL20", "SOX9"]:
    for mode in ["depth", "depth+injury"]:
        f = G[(G.cohort == "flex") & (G.anchor == a) & (G["mode"] == mode) & (G.donors == "all")].drop_duplicates("gene").set_index("gene")
        s = G[(G.cohort == "snrna") & (G.anchor == a) & (G["mode"] == mode) & (G.donors == "all")].drop_duplicates("gene").set_index("gene")
        j = f[["mean_excess", "n_pos_excess", "n_donors", "p_signflip"]].join(s[["mean_excess", "n_pos_excess", "n_donors", "p_signflip"]],
                                                                             lsuffix="_flex", rsuffix="_snrna", how="inner").reset_index()
        j.insert(0, "mode", mode); j.insert(0, "anchor", a)
        X2.append(j)
X2 = pd.concat(X2, ignore_index=True); X2.to_csv(OUT / "cross_cohort_excess.tsv.gz", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as fh:
    import importlib.metadata as im
    fh.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 400)
print(prim.round(4).to_string()); print(Z.round(3).to_string())
print(G[G.listed & (G.donors == "all")].sort_values(["cohort", "anchor", "mode", "mean_excess"]).round(3).to_string())
for (c, a, mode), o in G[G.donors == "all"].groupby(["cohort", "anchor", "mode"]):
    print(c, a, mode, "top / bottom excess")
    o = o.sort_values("mean_excess")
    print(pd.concat([o.head(20), o.tail(30)]).round(3).to_string())
for (a, mode), o in X2.groupby(["anchor", "mode"]):
    o = o.assign(both=np.minimum(o.mean_excess_flex, o.mean_excess_snrna)).sort_values("both", ascending=False)
    print(a, mode, "concordant up in both cohorts; spearman", round(stats.spearmanr(o.mean_excess_flex, o.mean_excess_snrna, nan_policy="omit")[0], 3))
    print(o.head(30).round(3).to_string())
    o = o.assign(both=np.maximum(o.mean_excess_flex, o.mean_excess_snrna)).sort_values("both")
    print(o.head(20).round(3).to_string())
