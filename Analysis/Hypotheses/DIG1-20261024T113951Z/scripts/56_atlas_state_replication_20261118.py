"""DIG1 job 56 (2026-11-18 deduplicated unique_counts rerun; N03).
Replaces scalesc_human_annotated_celltypist.h5ad (1.23M run-duplicated rows) and Aug-13 stage metadata
with NUCLEUS-FIX-20261107/unique_counts/unique_counts.h5ad (650,083 unique cells) and donor_metadata.tsv,
strictly filtering exclude != True & exclude_stage_analysis != True (drops GSE202379_P98; uses D11=MASL, D14=MASH),
and reads obs['donor'] directly instead of run-level mapping files.
Usage: python -I 56_atlas_state_replication_20261118.py <new_out_dir>
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
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
NF = ROOT / "Analysis/Hypotheses/NUCLEUS-FIX-20261107/unique_counts"
ATLAS = NF / "unique_counts.h5ad"
DONOR_META = NF / "metadata/donor_metadata.tsv"
J53 = RUN / "out_chemokine_hepatocyte_state"
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
K_MATCH, N_RAND, MIN_POS, N_FLIP = 10, 20, 10, 200_000
COHORTS = {"GSE202379": 20261026, "GSE244832": 20261027}
STAGE = {"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2, "Cirrhosis": 3}
SCORE = ["KRT8", "SQSTM1"]
ZON = {"central": ["CYP2E1", "GLUL", "CYP1A2", "OAT", "CYP3A4"], "portal": ["CPS1", "HAL", "ASS1", "SDS", "PCK1", "ALDOB"]}
ANCHORS = ["CCL20", "SOX9", "CXCL10"]
CORE = ["CXCL1", "CXCL8", "BIRC3", "TNFAIP3", "ICAM1", "ZC3H12A", "RELB", "NFKBIA", "NFKB2", "MAP3K14", "UBD", "FCAMR"]
LISTED = sorted(set(CORE) | {"CCL20", "CXCL10", "CXCL9", "CXCL11", "SOX9", "CX3CL1", "CSF1", "CCL2", "CXCL2", "LCN2", "NFKBIZ", "IL32",
                             "TNFRSF12A", "STAT1", "IRF1", "GBP1", "EPCAM", "KRT7", "KRT19", "SPP1", "BICC1", "PKHD1", "DCDC2", "FGFR2",
                             "CDH6", "ALB", "CYP3A4", "PTPRC", "CD68", "CD3E", "PECAM1", "COL1A1"})
TCELL_MARK = ["CCR6", "SLC4A10", "KLRB1", "RORC", "IL23R", "IL17A", "CXCR3", "CXCR6", "IFNG", "TNF"]
SEED_FLIP = 20261008

# gene sets fixed from job 53
X53 = pd.read_csv(J53 / "cross_cohort_excess.tsv.gz", sep="\t")
G53 = pd.read_csv(J53 / "gene_lfc_matched.tsv.gz", sep="\t")


def both_up(a):
    o = X53[(X53.anchor == a) & (X53["mode"] == "depth+injury")]
    return sorted(set(o[(o.mean_excess_flex > 0.5) & (o.mean_excess_snrna > 0.5)].gene) - {a})


g10 = G53[(G53.cohort == "flex") & (G53.anchor == "CXCL10") & (G53["mode"] == "depth+injury") & (G53.donors == "all")]
SETS = {"NFKB": both_up("CCL20"), "BILIARY": both_up("SOX9"),
        "IFN": sorted(set(g10[(g10.q_bh < 0.1) & (g10.mean_excess > 0.5)].gene) - {"CXCL10"}), "CORE": CORE}
pd.DataFrame([dict(set=k, n=len(v), genes=",".join(v)) for k, v in SETS.items()]).to_csv(OUT / "gene_sets.tsv", sep="\t", index=False)
print({k: len(v) for k, v in SETS.items()}, flush=True)


def signflip(v, rng=None):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) < 3:
        return np.nan
    if len(v) <= 20:
        S = np.array(list(product([1, -1], repeat=len(v))), dtype=np.int8)
    else:
        S = np.where(np.random.default_rng(SEED_FLIP).random((N_FLIP, len(v))) < 0.5, 1, -1).astype(np.int8)
    return float(np.mean(np.abs((S * v).mean(1)) >= abs(v.mean()) - 1e-12))


def signflip_matrix(E):
    """Two-sided sign-flip p per column of a donors x genes matrix (NaN donors dropped per column)."""
    n = E.shape[0]
    if n <= 20:
        S = np.array(list(product([1, -1], repeat=n)), dtype=np.float32)
    else:
        S = np.where(np.random.default_rng(SEED_FLIP).random((20_000, n)) < 0.5, 1.0, -1.0).astype(np.float32)
    p = np.full(E.shape[1], np.nan)
    ok = np.isfinite(E).all(0)
    Eo = np.nan_to_num(E[:, ok]).astype(np.float32)
    obs = np.abs(Eo.mean(0))
    cnt = np.zeros(Eo.shape[1])
    for i in range(0, len(S), 2000):
        cnt += (np.abs(S[i:i + 2000] @ Eo / n) >= obs - 1e-7).sum(0)
    p[ok] = cnt / len(S)
    for j in np.where(~ok)[0]:
        p[j] = signflip(E[:, j])
    return p


def bh(p):
    p = np.asarray(p, float); out = np.full(len(p), np.nan); ok = np.isfinite(p)
    if ok.any():
        q = p[ok]; o = np.argsort(q); r = q[o] * len(q) / (np.arange(len(q)) + 1)
        out[np.where(ok)[0][o]] = np.minimum(1, np.minimum.accumulate(r[::-1])[::-1])
    return out


def pct(x):
    return (stats.rankdata(x) - 0.5) / len(x)


def obs_str(o, k):
    g = o[k]
    if isinstance(g, h5py.Group) and "categories" in g:
        c = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]])
        codes = g["codes"][:]
        return np.where(codes >= 0, c[np.maximum(codes, 0)], "")
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g[:]])


def flag_true(s):
    return s.astype(str).str.strip().str.lower().eq("true")


# ---- donors from NUCLEUS-FIX donor_metadata.tsv
meta = pd.read_csv(DONOR_META, sep="\t", dtype=str)
eligible = meta[
    meta.dataset.isin(COHORTS)
    & (~flag_true(meta.exclude))
    & (~flag_true(meta.exclude_stage_analysis))
    & meta.disease_stage_coarse.isin(STAGE)
].copy()
stage_of = {r.donor: STAGE[r.disease_stage_coarse] for r in eligible.itertuples(index=False)}
print("Eligible donors by dataset:", eligible.groupby(["dataset", "disease_stage_coarse"]).size().to_dict(), flush=True)
assert "GSE202379_P98" not in stage_of
assert stage_of["GSE244832_D11"] == STAGE["Steatosis"]
assert stage_of["GSE244832_D14"] == STAGE["Steatohepatitis"]
eligible.to_csv(OUT / "donors_eligible.tsv", sep="\t", index=False)

# ---- read raw counts for eligible donors/cohorts in row blocks
with h5py.File(ATLAS, "r") as f:
    o = f["obs"]
    ds_all = obs_str(o, "dataset")
    donor_all = obs_str(o, "donor")
    ct_all = obs_str(o, "cell_type")
    prep_all = obs_str(o, "preparation_method")
    gnames = [x.decode() if isinstance(x, bytes) else str(x) for x in f["raw"]["var"]["_index"][:]]
    in_cohort = np.isin(ds_all, list(COHORTS))
    keep = in_cohort & np.isin(donor_all, list(stage_of))
    print("Dropped cells in COHORTS (excluded donors e.g. P98):", int((in_cohort & ~keep).sum()), flush=True)
    print("Retained prep methods:", pd.Series(prep_all[keep]).value_counts().to_dict(), flush=True)
    assert not pd.Series(prep_all[keep]).isin(["cd45_positive", "cd45_negative", "sorted"]).any(), "Unexpected sorted prep"
    rows = np.where(keep)[0]
    Xr = f["raw"]["X"]; indptr = Xr["indptr"][:]
    blocks = []
    for b0 in range(0, len(indptr) - 1, 100_000):
        b1 = min(b0 + 100_000, len(indptr) - 1)
        sel = rows[(rows >= b0) & (rows < b1)]
        if len(sel) == 0:
            continue
        lo, hi = indptr[b0], indptr[b1]
        B = sparse.csr_matrix((Xr["data"][lo:hi].astype(np.float32), Xr["indices"][lo:hi], indptr[b0:b1 + 1] - lo),
                              shape=(b1 - b0, len(gnames)))
        blocks.append(B[sel - b0])
        print("block", b0, len(sel), flush=True)
    X = sparse.vstack(blocks).tocsr(); del blocks
cell_ds = ds_all[rows]; cell_ct = ct_all[rows]; cell_donor = donor_all[rows].astype(str)
assert np.allclose(X.data, np.round(X.data)), "raw X is not integer counts"
print("cells", X.shape, "unique_donors", len(np.unique(cell_donor)), flush=True)


def cohort(tag, X, gnames, donor, seed):
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
    gene_rows, cell_rows, frac_rows = [], [], []
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
        tot_d = np.asarray(Xd.sum(axis=1)).ravel()
        tot_g = np.asarray(Xd.sum(axis=0)).ravel()
        frac_rows.append(dict(cohort=tag, donor=d, stage=stage_of[d], n_hep=len(m), median_umi=float(np.median(umi[m])),
                              **{f"n_{a}_pos": int((Xd[:, genes[a]] > 0).sum()) for a in ANCHORS if a in genes},
                              **{f"frac_{a}": float((Xd[:, genes[a]] > 0).mean()) for a in ANCHORS if a in genes},
                              **{f"cpm_{g}": float(tot_g[genes[g]] / tot_g.sum() * 1e6) for g in ["CCL20", "CXCL1", "CXCL8", "SOX9", "CXCL10"] if g in genes}))
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
                    gene_rows.append(dict(anchor=a, anc=anc, mode=mode, donor=d, lfc=lfc.astype(np.float32)))
                    cell_rows.append(dict(cohort=tag, donor=d, stage=stage_of[d], anchor=a, anc=anc, is_anchor=anc == a, mode=mode, n_pos=len(P),
                                          zon_diff=zdiff[P].mean() - (w @ zdiff) / w.sum(), inj_diff=inj_r[P].mean() - (w @ inj_r) / w.sum(),
                                          logumi_diff=np.log(umi[m][P]).mean() - (w @ np.log(umi[m])) / w.sum()))
        print(tag, d, len(m), flush=True)
    gsym = np.array(gnames)[gsel]
    gi = {g: i for i, g in enumerate(gsym)}
    out, setrows = [], []
    GR = pd.DataFrame(gene_rows)
    for (a, mode), grp in GR.groupby(["anchor", "mode"]):
        dons = sorted(grp[grp.anc == a].donor)
        A_ = np.vstack([grp[(grp.anc == a) & (grp.donor == d)].lfc.iloc[0] for d in dons])
        Nl = np.vstack([np.mean(np.vstack(grp[(grp.anc != a) & (grp.donor == d)].lfc.tolist()), 0) if ((grp.anc != a) & (grp.donor == d)).any()
                        else np.full(len(gsel), np.nan) for d in dons])
        Ex = A_ - Nl
        for sname, sg in SETS.items():
            idx = [gi[g] for g in sg if g in gi]
            for d, row in zip(dons, Ex[:, idx]):
                setrows.append(dict(cohort=tag, anchor=a, mode=mode, set=sname, n_set_genes=len(idx), donor=d, stage=stage_of[d],
                                    set_excess=float(np.nanmean(row))))
        p = signflip_matrix(Ex)
        out.append(pd.DataFrame(dict(cohort=tag, anchor=a, mode=mode, gene=gsym, n_donors=len(dons), mean_lfc=A_.mean(0),
                                     mean_null_lfc=np.nanmean(Nl, 0), mean_excess=np.nanmean(Ex, 0), n_pos_excess=(Ex > 0).sum(0),
                                     p_signflip=p, q_bh=bh(p), listed=np.isin(gsym, LISTED))))
    return (pd.DataFrame(cell_rows), pd.concat(out, ignore_index=True), pd.DataFrame(setrows), pd.DataFrame(frac_rows),
            anchors)


Cs, Gs, Ss, Fs, An = [], [], [], [], {}
for tag, seed in COHORTS.items():
    m = (cell_ds == tag) & (cell_ct == "Hepatocytes")
    c, g, s, fr, an = cohort(tag, X[np.where(m)[0]], gnames, cell_donor[m], seed)
    Cs.append(c); Gs.append(g); Ss.append(s); Fs.append(fr); An[tag] = an
C = pd.concat(Cs, ignore_index=True); G = pd.concat(Gs, ignore_index=True); S = pd.concat(Ss, ignore_index=True); F = pd.concat(Fs, ignore_index=True)
C.to_csv(OUT / "donor_anchor_zonation_injury.tsv", sep="\t", index=False)
G.to_csv(OUT / "gene_lfc_matched.tsv.gz", sep="\t", index=False)
S.to_csv(OUT / "donor_set_excess.tsv", sep="\t", index=False)
F.to_csv(OUT / "donor_ccl20_fractions.tsv", sep="\t", index=False)
pd.DataFrame([dict(cohort=c, anchor=a, random_anchors=",".join(r)) for c, an in An.items() for a, r in an.items()]).to_csv(
    OUT / "random_anchors.tsv", sep="\t", index=False)

# set summaries
rows = []
for (c, a, mode, s), o in S.groupby(["cohort", "anchor", "mode", "set"]):
    v = o.set_excess.to_numpy()
    rho, prho = stats.spearmanr(o.stage, v) if o.stage.nunique() > 1 else (np.nan, np.nan)
    rows.append(dict(cohort=c, anchor=a, mode=mode, set=s, n_set_genes=o.n_set_genes.iloc[0], n_donors=len(o), mean_excess=np.nanmean(v),
                     n_pos=int((v > 0).sum()), p_signflip=signflip(v), spearman_stage=rho, p_stage=prho,
                     **{f"mean_stage{k}": np.nanmean(o.set_excess[o.stage == k]) if (o.stage == k).any() else np.nan for k in range(4)}))
SS = pd.DataFrame(rows); SS.to_csv(OUT / "summary_set_excess.tsv", sep="\t", index=False)
prim = SS[(SS.anchor == "CCL20") & (SS.set == "NFKB") & (SS["mode"] == "depth+injury")].copy()
prim.insert(0, "id", prim.cohort.map({"GSE202379": "P1", "GSE244832": "P2"}))
oo = np.argsort(prim.p_signflip.to_numpy()); ph = np.maximum.accumulate(prim.p_signflip.to_numpy()[oo] * (len(oo) - np.arange(len(oo))))
prim.loc[prim.index[oo], "p_holm"] = np.minimum(1, ph)
prim.to_csv(OUT / "primary_holm.tsv", sep="\t", index=False)

# zonation / injury summaries
zr = []
for (c, a, mode), o in C.groupby(["cohort", "anchor", "mode"]):
    for var in ["zon_diff", "inj_diff", "logumi_diff"]:
        piv = o.pivot_table(index="donor", columns="is_anchor", values=var, aggfunc="mean")
        if True not in piv.columns:
            continue
        piv = piv.dropna(subset=[True]); ex = piv[True] - piv.get(False, np.nan)
        zr.append(dict(cohort=c, anchor=a, mode=mode, var=var, n_donors=len(piv), mean_anchor=piv[True].mean(), mean_excess=ex.mean(),
                       n_pos=int((ex > 0).sum()), p_signflip=signflip(ex.to_numpy())))
pd.DataFrame(zr).to_csv(OUT / "summary_zonation_injury.tsv", sep="\t", index=False)

# donor-level immune links
genes = {}
for i, g in enumerate(gnames):
    genes.setdefault(g, i)
irows = []
for d in sorted(set(F.donor)):
    md = cell_donor == d
    t = np.where(md & (cell_ct == "T cells"))[0]
    mac = np.where(md & np.isin(cell_ct, ["Macrophages", "Mono+mono derived cells"]))[0]
    rec = dict(donor=d, n_cells=int(md.sum()), n_t=len(t), neut_share=float((md & (cell_ct == "Neutrophils")).sum() / md.sum()),
               t_share=float(len(t) / md.sum()))
    if len(t) >= 30:
        Tb = X[t][:, [genes[g] for g in TCELL_MARK]] > 0
        rec.update({f"t_frac_{g}": float(Tb[:, j].mean()) for j, g in enumerate(TCELL_MARK)})
    if len(mac) >= 30:
        tot = np.asarray(X[mac].sum(axis=0)).ravel()
        rec.update({f"mac_cpm_{g}": float(tot[genes[g]] / tot.sum() * 1e6) for g in ["TNF", "IL1B", "IL1A", "OSM", "CXCL8"]})
    irows.append(rec)
I = F.merge(pd.DataFrame(irows), on="donor", how="left")
cc = S[(S.anchor == "CCL20") & (S.set == "NFKB") & (S["mode"] == "depth+injury")][["donor", "set_excess"]].rename(columns={"set_excess": "ccl20_nfkb_excess"})
I = I.merge(cc, on="donor", how="left")
I.to_csv(OUT / "donor_immune_links.tsv", sep="\t", index=False)
lr = []
for c, o in I.groupby("cohort"):
    for x in ["frac_CCL20", "cpm_CCL20"]:
        for y in [k for k in I.columns if k.startswith(("t_frac_", "mac_cpm_")) or k in ("neut_share", "t_share")]:
            oo_ = o[[x, y, "stage"]].dropna()
            if len(oo_) < 8:
                continue
            r, p = stats.spearmanr(oo_[x], oo_[y])
            rx = stats.rankdata(oo_[x]); ry = stats.rankdata(oo_[y]); rs = stats.rankdata(oo_.stage)
            Z = np.column_stack([np.ones(len(rs)), rs])
            ex_ = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]; ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
            lr.append(dict(cohort=c, x=x, y=y, n=len(oo_), spearman=r, p=p, partial_stage=np.corrcoef(ex_, ey)[0, 1]))
L = pd.DataFrame(lr)
if len(L):
    L["q_bh_within_cohort"] = L.groupby("cohort").p.transform(lambda v: bh(v.to_numpy()))
L.to_csv(OUT / "immune_link_correlations.tsv", sep="\t", index=False)

with open(OUT / "pip_freeze.txt", "w") as fh:
    import importlib.metadata as im
    fh.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 400)
print(prim.round(4).to_string()); print(SS.round(3).to_string()); print(pd.DataFrame(zr).round(3).to_string())
print(F.groupby(["cohort", "stage"]).agg(n=("donor", "size"), frac_CCL20=("frac_CCL20", "median"), cpm_CCL20=("cpm_CCL20", "median"),
                                          frac_SOX9=("frac_SOX9", "median")).round(4).to_string())
print(L.round(3).to_string() if len(L) else "no immune links")
