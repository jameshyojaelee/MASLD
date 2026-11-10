"""CONF1 C5/C1b step 1 (PRESPEC_C5_C1b.md). GSE289173 hepatocytes.
Writes donor hepatocyte pseudobulk (all genes) and donor table for C5, then runs C1b with the DIG1 job 56
cohort() method unchanged and the frozen NFKB set.
Usage: python -I c5_gse289173_extract_c1b.py <new_out_dir>
"""
import sys
from itertools import product
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.spatial import cKDTree

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
H5 = ROOT / "data/public_r4/GSE289173/derived_H90-20261015T120337Z/gse289173_qc_nuclei.h5ad"
SETS_TSV = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z/out_atlas_state_replication/gene_sets.tsv"
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
K_MATCH, N_RAND, MIN_POS, N_FLIP, MIN_ELIG = 10, 20, 10, 200_000, 6
SEED_ANCHOR, SEED_FLIP = 20261109, 20261008
STAGE = {"ctrl": 0, "MASL": 1, "eMASH": 2, "aMASH": 3}
SCORE = ["KRT8", "SQSTM1"]
ZON = {"central": ["CYP2E1", "GLUL", "CYP1A2", "OAT", "CYP3A4"], "portal": ["CPS1", "HAL", "ASS1", "SDS", "PCK1", "ALDOB"]}
ANCHORS = ["CCL20"]
CORE = ["CXCL1", "CXCL8", "BIRC3", "TNFAIP3", "ICAM1", "ZC3H12A", "RELB", "NFKBIA", "NFKB2", "MAP3K14", "UBD", "FCAMR"]
LISTED = sorted(set(CORE) | {"CCL20", "CXCL10", "CXCL9", "CXCL11", "SOX9", "CX3CL1", "CSF1", "CCL2", "CXCL2", "LCN2", "NFKBIZ", "IL32",
                             "TNFRSF12A", "STAT1", "IRF1", "GBP1", "EPCAM", "KRT7", "KRT19", "SPP1", "BICC1", "PKHD1", "DCDC2", "FGFR2",
                             "CDH6", "ALB", "CYP3A4", "PTPRC", "CD68", "CD3E", "PECAM1", "COL1A1"})
st = pd.read_csv(SETS_TSV, sep="\t")
SETS = {r.set: r.genes.split(",") for r in st.itertuples()}
assert len(SETS["NFKB"]) == 64
pd.DataFrame([dict(set=k, n=len(v), genes=",".join(v)) for k, v in SETS.items()]).to_csv(OUT / "gene_sets_used.tsv", sep="\t", index=False)


def signflip(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) < 3:
        return np.nan
    if len(v) <= 20:
        S = np.array(list(product([1, -1], repeat=len(v))), dtype=np.int8)
    else:
        S = np.where(np.random.default_rng(SEED_FLIP).random((N_FLIP, len(v))) < 0.5, 1, -1).astype(np.int8)
    return np.mean(np.abs((S * v).mean(1)) >= abs(v.mean()) - 1e-12)


def pct(x):
    return (stats.rankdata(x) - 0.5) / len(x)


def cat(o, k):
    g = o[k]; c = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]]); return c[g["codes"][:]]


# ---- read hepatocytes
with h5py.File(H5, "r") as f:
    o = f["obs"]
    gsm = cat(o, "gsm"); grp = cat(o, "group"); lab = cat(o, "label"); sex = cat(o, "sex")
    gnames = [x.decode() for x in f["var"]["_index"][:]]
    rows = np.where(lab == "Hepatocytes")[0]
    Xr = f["X"]; indptr = Xr["indptr"][:]
    blocks = []
    for b0 in range(0, len(indptr) - 1, 50_000):
        b1 = min(b0 + 50_000, len(indptr) - 1)
        sel = rows[(rows >= b0) & (rows < b1)]
        if len(sel) == 0:
            continue
        lo, hi = indptr[b0], indptr[b1]
        B = sparse.csr_matrix((Xr["data"][lo:hi], Xr["indices"][lo:hi], indptr[b0:b1 + 1] - lo), shape=(b1 - b0, len(gnames)))
        blocks.append(B[sel - b0])
        print("block", b0, len(sel), flush=True)
    X = sparse.vstack(blocks).tocsr(); del blocks
assert np.allclose(X.data, np.round(X.data)), "X is not integer counts"
donor = gsm[rows]; dgrp = grp[rows]; dsex = sex[rows]
print("hepatocytes", X.shape, flush=True)

# ---- C5 inputs: donor pseudobulk
dons = sorted(np.unique(donor))
pb = np.column_stack([np.asarray(X[donor == d].sum(axis=0)).ravel() for d in dons]).astype(np.int64)
pd.DataFrame(pb, index=gnames, columns=dons).rename_axis("symbol").reset_index().to_csv(OUT / "pseudobulk_hep.tsv.gz", sep="\t", index=False)
meta = pd.DataFrame(dict(donor=dons,
                         group=[dgrp[donor == d][0] for d in dons],
                         female=[int(dsex[donor == d][0] == "female") for d in dons],
                         n_hep=[int((donor == d).sum()) for d in dons]))
for d in dons:
    assert len(set(dgrp[donor == d])) == 1 and len(set(dsex[donor == d])) == 1
meta["stage_ord"] = meta.group.map(STAGE)
meta.to_csv(OUT / "donors.tsv", sep="\t", index=False)


# ---- C1b: DIG1 job 56 cohort() unchanged except seed and the sets/stage mapping above
def cohort(X, gnames, donor, stage_of, seed):
    rng = np.random.default_rng(seed)
    genes = {}
    for i, g in enumerate(gnames):
        genes.setdefault(g, i)
    umi = np.asarray(X.sum(axis=1)).ravel(); ngene = np.asarray((X > 0).sum(axis=1)).ravel()
    det_all = np.asarray((X > 0).mean(axis=0)).ravel()
    zon = {k: [x for x in v if x in genes] for k, v in ZON.items()}
    setg = set(x for v in SETS.values() for x in v)
    excl = set(LISTED) | set(SCORE) | setg | set(x for v in zon.values() for x in v)
    gsel = np.where((det_all >= 0.01) | np.isin(gnames, sorted(set(LISTED) | setg)))[0]
    anchors = {}
    for a in [x for x in ANCHORS if x in genes]:
        da = det_all[genes[a]]
        cand = sorted(set(g for g, d in zip(gnames, det_all) if 0.67 * da <= d <= 1.5 * da and g not in excl and not g.startswith("ENSG")))
        anchors[a] = list(rng.choice(cand, size=min(N_RAND, len(cand)), replace=False))
    side = sorted(set(SCORE) | set(x for v in zon.values() for x in v))
    Ls = np.log1p(X[:, [genes[g] for g in side]].toarray() * (1e4 / umi[:, None])); si = {g: i for i, g in enumerate(side)}
    gene_rows, frac_rows = [], []
    for d in sorted(np.unique(donor)):
        m = np.where(donor == d)[0]; Xd = X[m]
        frac_rows.append(dict(donor=d, stage=stage_of[d], n_hep=len(m), n_ccl20_pos=int((Xd[:, genes["CCL20"]] > 0).sum()),
                              frac_CCL20=float((Xd[:, genes["CCL20"]] > 0).mean())))
        if len(m) < 6 * MIN_POS:
            continue
        S = Ls[m]; sd = S.std(0); sd[sd == 0] = np.nan; Zs = (S - S.mean(0)) / sd
        sc = lambda gs: np.nan_to_num(np.nanmean(Zs[:, [si[g] for g in gs]], 1))
        inj = sc(SCORE)
        Bm = np.column_stack([np.ones(len(m)), stats.rankdata(np.log(umi[m])), stats.rankdata(np.log(ngene[m]))])
        inj_r = inj - Bm @ np.linalg.lstsq(Bm, inj, rcond=None)[0]
        feats = {"depth": np.column_stack([pct(np.log(umi[m])), pct(np.log(ngene[m]))]),
                 "depth+injury": np.column_stack([pct(np.log(umi[m])), pct(np.log(ngene[m])), pct(inj_r)])}
        Xg = Xd[:, gsel].tocsr()
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
                    gene_rows.append(dict(anchor=a, anc=anc, mode=mode, donor=d, n_pos=len(P), lfc=lfc.astype(np.float32)))
    return gene_rows, pd.DataFrame(frac_rows), anchors, np.array(gnames)[gsel]


stage_of = dict(zip(meta.donor, meta.stage_ord))
gene_rows, FR, AN, gsym = cohort(X, gnames, donor, stage_of, SEED_ANCHOR)
pd.DataFrame([dict(anchor=a, random_anchors=",".join(r)) for a, r in AN.items()]).to_csv(OUT / "random_anchors.tsv", sep="\t", index=False)
GR = pd.DataFrame(gene_rows)
elig = sorted(GR[(GR.anc == "CCL20") & (GR["mode"] == "depth+injury")].donor) if len(GR) else []
FR["eligible"] = FR.donor.isin(elig)
FR = FR.merge(meta[["donor", "group"]], on="donor")
FR.to_csv(OUT / "c1b_detection_gate.tsv", sep="\t", index=False)
print("C1b eligible donors", len(elig), FR[FR.eligible].group.value_counts().to_dict(), flush=True)

gi = {g: i for i, g in enumerate(gsym)}
if len(elig) < MIN_ELIG:
    pd.DataFrame([dict(test="C1b", n_eligible=len(elig), decision="not_testable_detection_gate")]).to_csv(
        OUT / "c1b_primary.tsv", sep="\t", index=False)
else:
    setrows, generows = [], []
    for mode, grp_ in GR.groupby("mode"):
        dd = sorted(grp_[grp_.anc == "CCL20"].donor)
        A_ = np.vstack([grp_[(grp_.anc == "CCL20") & (grp_.donor == d)].lfc.iloc[0] for d in dd])
        Nl = np.vstack([np.mean(np.vstack(grp_[(grp_.anc != "CCL20") & (grp_.donor == d)].lfc.tolist()), 0)
                        if ((grp_.anc != "CCL20") & (grp_.donor == d)).any() else np.full(len(gsym), np.nan) for d in dd])
        Ex = A_ - Nl
        for sname, sg in SETS.items():
            idx = [gi[g] for g in sg if g in gi]
            for d, row in zip(dd, Ex[:, idx]):
                setrows.append(dict(mode=mode, set=sname, n_set_genes=len(idx), donor=d, stage=stage_of[d], set_excess=float(np.nanmean(row))))
        generows.append(pd.DataFrame(dict(mode=mode, gene=gsym, n_donors=len(dd), mean_excess=np.nanmean(Ex, 0),
                                          n_pos_excess=(Ex > 0).sum(0), in_nfkb=np.isin(gsym, SETS["NFKB"]))))
    S = pd.DataFrame(setrows); S.to_csv(OUT / "c1b_donor_set_excess.tsv", sep="\t", index=False)
    pd.concat(generows).to_csv(OUT / "c1b_gene_excess.tsv.gz", sep="\t", index=False)
    rows_ = []
    for (mode, s), o_ in S.groupby(["mode", "set"]):
        v = o_.set_excess.to_numpy()
        rho, prho = stats.spearmanr(o_.stage, v) if o_.stage.nunique() > 1 else (np.nan, np.nan)
        rows_.append(dict(mode=mode, set=s, n_set_genes=o_.n_set_genes.iloc[0], n_donors=len(o_), mean_excess=np.nanmean(v),
                          n_pos=int((v > 0).sum()), p_signflip=signflip(v), spearman_stage=rho, p_stage=prho))
    SS = pd.DataFrame(rows_); SS.to_csv(OUT / "c1b_summary.tsv", sep="\t", index=False)
    prim = SS[(SS["set"] == "NFKB") & (SS["mode"] == "depth+injury")].copy()
    prim["decision"] = np.where((prim.mean_excess > 0) & (prim.p_signflip < 0.05), "replicated", "not_replicated")
    prim.to_csv(OUT / "c1b_primary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(prim.round(4).to_string()); print(SS.round(3).to_string())

with open(OUT / "pip_freeze.txt", "w") as fh:
    import importlib.metadata as im
    fh.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
