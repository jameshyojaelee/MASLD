"""DIG1 job 52 (exploratory, post hoc to job 51). Job 51 found that the MASLD rise of hepatocyte CCL20 sits
in the most injured quarter of hepatocytes. In mouse MASH, RELB and SOX9 activate CCL20 transcription
(doi 10.1016/j.metabol.2025.156301; there CCL20 came mainly from progenitor and cholangiocyte-like cells).
Here: within single human livers, are the hepatocytes that carry CCL20 RNA the RELB-high and SOX9-high
ones, beyond what the shared injury score explains?
Cohorts: GSE298719 Flex hepatocytes (13 donors) and GSE256398 snRNA hepatocytes (17 donors, job 45L labels).
Statistic per donor: rho(anchor, partner) on ranks residualised on log UMI, log genes, central and
portal zonation scores, the long-gene share (job 47 A2) and the KRT8+SQSTM1 injury score ("beyond
injury"), minus the mean rho of the anchor with 200 random partners matched to the partner on hepatocyte
mean x detection deciles x genomic-span tertile (listed, score and zonation genes excluded; seed 20261024).
A donor counts when anchor and partner are each detected in >= 10 of its hepatocytes.
Primary family (Holm over 4), exact two-sided sign-flip over eligible donors:
  C1 Flex CCL20-RELB, C2 Flex CCL20-SOX9, C3 snRNA CCL20-RELB, C4 snRNA CCL20-SOX9.
Descriptive: other pairs (CCL20 with NFKB2, EPCAM, CD24, KRT7, KRT19; CXCL10 and CXCL9 with RELB, SOX9,
STAT1, IRF1 in Flex); without the injury score; with a cholangiocyte score (KRT7, KRT19) added; MASLD
donors only; co-detection of RELB, SOX9, KRT7, KRT19 in CCL20+ vs CCL20- hepatocytes; share of CCL20
UMIs by lineage (Flex all lineages; snRNA from V3 lineage pseudobulk).
Usage: python -I 52_ccl20_relb_sox9.py <new_out_dir>
"""
import sys
from itertools import product
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
V5 = ROOT / "Analysis/Hypotheses/V5-20261021T122723Z"
V3PB = ROOT / "Analysis/Hypotheses/V3-20261008T015239Z/pseudobulk"
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
N_NULL = 200
SCORE = ["KRT8", "SQSTM1"]
CHOL = ["KRT7", "KRT19"]
ZON = {"central": ["CYP2E1", "GLUL", "CYP1A2", "OAT", "CYP3A4"], "portal": ["CPS1", "HAL", "ASS1", "SDS", "PCK1", "ALDOB"]}
PAIRS = [("CCL20", p) for p in ["RELB", "SOX9", "NFKB2", "EPCAM", "CD24", "KRT7", "KRT19"]] + \
        [(a, p) for a in ["CXCL10", "CXCL9"] for p in ["RELB", "SOX9", "STAT1", "IRF1"]]
PRIMARY = {("flex", "CCL20", "RELB"): "C1", ("flex", "CCL20", "SOX9"): "C2", ("snrna", "CCL20", "RELB"): "C3", ("snrna", "CCL20", "SOX9"): "C4"}
GL = pd.read_csv(RUN / "out_length_confound/gene_lengths_ensembl98.tsv.gz", sep="\t")
by_id = GL.set_index("gene_id").span
by_sym = GL.sort_values("span", ascending=False).drop_duplicates("gene_name").set_index("gene_name").span


def resid_rank(M, Z):
    R = stats.rankdata(M, axis=0)
    return R - Z @ np.linalg.lstsq(Z, R, rcond=None)[0]


def cor_cols(a, M):
    a = (a - a.mean()) / a.std(); sd = M.std(0); sd[sd == 0] = np.nan
    return ((M - M.mean(0)) / sd * a[:, None]).mean(0)


def signflip(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) < 3:
        return np.nan
    S = np.array(list(product([1, -1], repeat=len(v))), dtype=np.int8)
    return np.mean(np.abs((S * v).mean(1)) >= abs(v.mean()) - 1e-12)


def cohort(tag, X, gnames, glen, donor, dis_of, seed):
    rng = np.random.default_rng(seed)
    X = sparse.csr_matrix(X); genes = pd.Index(gnames)
    umi = np.asarray(X.sum(axis=1)).ravel(); ngene = np.asarray((X > 0).sum(axis=1)).ravel()
    L = X.multiply(1e4 / umi[:, None]).tocsr(); L.data = np.log1p(L.data)
    mean_all = np.asarray(L.mean(axis=0)).ravel(); det_all = np.asarray((X > 0).mean(axis=0)).ravel()
    pairs = [(a, p) for a, p in PAIRS if a in genes and p in genes]
    listed = set(x for ap in pairs for x in ap)
    zon = {k: [x for x in v if x in genes] for k, v in ZON.items()}
    excl = listed | set(SCORE) | set(CHOL) | set(x for v in zon.values() for x in v)
    keep = np.where((det_all >= 0.01) & np.isfinite(glen))[0]; kg = np.array(gnames)[keep]
    binm = pd.qcut(pd.Series(mean_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
    bind = pd.qcut(pd.Series(det_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
    bint = pd.qcut(pd.Series(np.log10(glen[keep])).rank(method="first"), 3, labels=False).to_numpy()
    lq = pd.qcut(pd.Series(np.log10(glen[keep])).rank(method="first"), 5, labels=False).to_numpy()
    long_idx = keep[(lq == 4) & ~np.isin(kg, list(excl))]
    long_share = np.asarray(X[:, long_idx].sum(axis=1)).ravel() / umi
    kbin = dict(zip(kg, (binm * 10 + bind) * 3 + bint))
    pool = pd.Series(kbin); pool = pool[~pool.index.isin(list(excl))]
    pool = {b: list(ix) for b, ix in pool.groupby(pool).groups.items()}

    def bin_of(g):  # partner may be below the 1% filter or listed; place it by its own mean, detection and span
        i = genes.get_loc(g)
        rm = min(int((mean_all[keep] < mean_all[i]).mean() * 10), 9); rd = min(int((det_all[keep] < det_all[i]).mean() * 10), 9)
        rt = min(int((np.log10(glen[keep]) < np.log10(glen[i])).mean() * 3), 2) if np.isfinite(glen[i]) else 1
        b = (rm * 10 + rd) * 3 + rt
        if b in pool and len(pool[b]) >= 10:
            return b
        ks = [k for k in pool if len(pool[k]) >= 10]
        return min(ks, key=lambda k: abs(k // 30 - rm) + abs((k // 3) % 10 - rd) + abs(k % 3 - rt))
    nulls = {p: [pool[bin_of(p)][j] for j in rng.integers(0, 10 ** 9, N_NULL) % len(pool[bin_of(p)])] for p in sorted(set(p for _, p in pairs))}
    need = sorted(excl | set(x for v in nulls.values() for x in v))
    Ld = L[:, [genes.get_loc(x) for x in need]].toarray(); gi = {x: i for i, x in enumerate(need)}
    Xd = (X[:, [genes.get_loc(x) for x in need]] > 0).toarray()
    rows, codet = [], []
    for d in sorted(np.unique(donor)):
        m = donor == d; nd = int(m.sum())
        D = Ld[m].astype(np.float64); sd = D.std(0); sd[sd == 0] = np.nan
        Zg = (D - D.mean(0)) / sd
        sc = lambda gs: np.nanmean(Zg[:, [gi[x] for x in gs if x in gi]], 1)
        inj = sc(SCORE); chol = sc([x for x in CHOL if x in gi])
        base = np.column_stack([np.ones(nd), stats.rankdata(np.log(umi[m])), stats.rankdata(np.log(ngene[m])),
                                stats.rankdata(np.nan_to_num(sc(zon["central"]))), stats.rankdata(np.nan_to_num(sc(zon["portal"]))),
                                stats.rankdata(long_share[m])])
        adjs = {"no_injury": base, "beyond_injury": np.column_stack([base, stats.rankdata(inj)]),
                "beyond_injury+chol": np.column_stack([base, stats.rankdata(inj), stats.rankdata(np.nan_to_num(chol))])}
        det = Xd[m]
        for a, p in pairs:
            na, npn = int(det[:, gi[a]].sum()), int(det[:, gi[p]].sum())
            ok = na >= 10 and npn >= 10
            for adj, Zr in adjs.items():
                if not ok or (adj == "beyond_injury+chol" and p in CHOL):
                    rows.append(dict(cohort=tag, donor=d, disease=dis_of[d], anchor=a, partner=p, adjust=adj, n_cells=nd,
                                     n_anchor_pos=na, n_partner_pos=npn, eligible=False)); continue
                ra = resid_rank(D[:, [gi[a]]], Zr)[:, 0]; rp = resid_rank(D[:, [gi[p]]], Zr)[:, 0]
                Rn = resid_rank(D[:, [gi[x] for x in nulls[p]]], Zr)
                r = np.corrcoef(ra, rp)[0, 1]; nr = cor_cols(ra, Rn)
                rows.append(dict(cohort=tag, donor=d, disease=dis_of[d], anchor=a, partner=p, adjust=adj, n_cells=nd, n_anchor_pos=na,
                                 n_partner_pos=npn, eligible=True, rho=r, null_mean=np.nanmean(nr), null_sd=np.nanstd(nr)))
        if "CCL20" in gi:
            pos = det[:, gi["CCL20"]]
            for p in ["RELB", "SOX9", "KRT7", "KRT19", "EPCAM", "NFKB2", "CD24"]:
                if p in gi:
                    codet.append(dict(cohort=tag, donor=d, disease=dis_of[d], gene=p, n_ccl20_pos=int(pos.sum()),
                                      frac_in_ccl20_pos=det[pos, gi[p]].mean() if pos.any() else np.nan, frac_in_ccl20_neg=det[~pos, gi[p]].mean(),
                                      umi_ccl20_pos=np.median(umi[m][pos]) if pos.any() else np.nan, umi_ccl20_neg=np.median(umi[m][~pos])))
        print(tag, d, nd, flush=True)
    return pd.DataFrame(rows), pd.DataFrame(codet)


# Flex
A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad")
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str}); dq = dq[dq.included]
A = A[A.obs.gsm.isin(dq.gsm).to_numpy()].copy()
fg = list(A.var_names)
lin = []
Xa = sparse.csr_matrix(A.X)
for g in ["CCL20", "CXCL10", "CXCL9", "RELB", "SOX9"]:
    j = fg.index(g)
    for (gsm, lab), ix in A.obs.groupby(["gsm", "label"], observed=True).indices.items():
        lin.append(dict(cohort="flex", gene=g, donor=gsm, disease=dict(zip(dq.gsm, dq.disease))[gsm], lineage=lab, n_cells=len(ix),
                        umi=float(Xa[ix, j].sum()), frac_pos=float((Xa[ix, j] > 0).mean())))
H = (A.obs.label == "Hepatocytes").to_numpy()
wf, cf = cohort("flex", Xa[H], fg, by_sym.reindex(fg).to_numpy(float), A.obs.gsm.astype(str).to_numpy()[H], dict(zip(dq.gsm, dq.disease)), 20261024)
del A, Xa
# snRNA
Hs = ad.read_h5ad(RUN / "out_snrna_labels/gse256398_hepatocyte_counts.h5ad")
md = pd.read_csv(V3PB / "donor_metadata.tsv", sep="\t")
sg = list(Hs.var_names); sid = Hs.var.gene_ids.astype(str).to_numpy()
ws, cs = cohort("snrna", Hs.X, sg, by_id.reindex(sid).to_numpy(float), Hs.obs.donor.astype(str).to_numpy(), dict(zip(md.donor.astype(str), md.disease)), 20261024)
ccl20_id = sid[sg.index("CCL20")]
del Hs
for f in sorted(V3PB.glob("pseudobulk_*.tsv.gz")):
    pb = pd.read_csv(f, sep="\t", index_col=0)
    if ccl20_id not in pb.index:
        continue
    row = pb.loc[ccl20_id].drop(labels=[c for c in ["symbol"] if c in pb.columns])
    for d, v in row.items():
        lin.append(dict(cohort="snrna", gene="CCL20", donor=d, disease=dict(zip(md.donor.astype(str), md.disease)).get(str(d), np.nan),
                        lineage=f.name.replace("pseudobulk_", "").replace(".tsv.gz", ""), n_cells=np.nan, umi=float(v), frac_pos=np.nan))

W = pd.concat([wf, ws], ignore_index=True); W["excess"] = W.rho - W.null_mean
W.to_csv(OUT / "donor_pair_coupling.tsv", sep="\t", index=False)
pd.concat([cf, cs], ignore_index=True).to_csv(OUT / "ccl20_codetection.tsv", sep="\t", index=False)
LIN = pd.DataFrame(lin); LIN.to_csv(OUT / "lineage_umi.tsv", sep="\t", index=False)
summ = []
E = W[W.eligible == True]
for (c, a, p, adj), o in E.groupby(["cohort", "anchor", "partner", "adjust"]):
    for grp, oo in [("all", o), ("masld", o[o.disease == 1])]:
        summ.append(dict(cohort=c, anchor=a, partner=p, adjust=adj, donors=grp, n_donors=len(oo), mean_rho=oo.rho.mean(), mean_null=oo.null_mean.mean(),
                         mean_excess=oo.excess.mean(), n_positive=int((oo.excess > 0).sum()), p_signflip=signflip(oo.excess)))
S = pd.DataFrame(summ); S.to_csv(OUT / "summary_signflip.tsv", sep="\t", index=False)
prim = S[(S.adjust == "beyond_injury") & (S.donors == "all")].copy()
prim["id"] = [PRIMARY.get((c, a, p)) for c, a, p in zip(prim.cohort, prim.anchor, prim.partner)]
prim = prim[prim.id.notna()].sort_values("id").reset_index(drop=True)
o = np.argsort(prim.p_signflip.to_numpy()); ph = np.maximum.accumulate(prim.p_signflip.to_numpy()[o] * (len(o) - np.arange(len(o))))
prim.loc[prim.index[o], "p_holm"] = np.minimum(1, ph)
prim.to_csv(OUT / "primary_holm.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 300)
print(prim.round(4).to_string()); print(S.round(3).to_string())
C = pd.concat([cf, cs]); print(C.groupby(["cohort", "disease", "gene"])[["n_ccl20_pos", "frac_in_ccl20_pos", "frac_in_ccl20_neg"]].agg({"n_ccl20_pos": "sum", "frac_in_ccl20_pos": "mean", "frac_in_ccl20_neg": "mean"}).round(3).to_string())
tot = LIN.groupby(["cohort", "gene", "disease", "lineage"]).umi.sum().reset_index()
tot["share"] = tot.umi / tot.groupby(["cohort", "gene", "disease"]).umi.transform("sum")
print(tot[tot.share > 0.01].sort_values(["cohort", "gene", "disease", "share"], ascending=[True, True, True, False]).round(3).to_string())
