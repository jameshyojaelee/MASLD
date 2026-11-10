"""DIG1 job 49 (exploratory; follows jobs 46-48). CosMx in situ: in intact tissue, do the hepatocytes
with more KRT8/KRT18/SQSTM1 RNA also carry more NF-kB response RNA, IL32 and TNFRSF12A (Fn14), as in
Flex cells and snRNA nuclei (jobs 47, 48)? Are those hepatocytes nearer macrophages?
Built from 38_insitu_heatshock.py: same data (Govaere 4 sections, GSE292268 4 biopsies, Fukumoto 2
samples), QC, lineage labels (job 14b), core-hepatocyte rule, normalisation and matched-set null.
Sets (fixed before reading outcomes; panel genes from the job 48 held-out and seen NF-kB lists that
are on all three panels): NFKB = NFKBIA, NFKB1, SOD2, ICAM1, IER3; MDB = KRT8, KRT18, SQSTM1.
Statistic per unit: rho(set, MDB) with ranks residualised on log total and log n_genes, minus the
mean rho of 200 random panel gene sets matched on core-hepatocyte mean x detection decile bins
(all listed genes excluded from the pool; seed 20261024).
Primary family (Holm over 3), exact sign-flip over evaluable units (>= 500 core hepatocytes):
  J1  NFKB set excess.  J2  IL32 excess.  J3  TNFRSF12A excess.
Descriptive: single genes, leave-one-out NFKB sets, spillover controls (stromal, cholangiocyte,
macrophage genes; negative probes); sensitivities as job 38. Neighbourhood: rho between a test score
and the fraction of each lineage among the 30 nearest cells (both residualised on depth), minus the
same for the matched null sets; macrophage genes in hepatocytes are the spillover reference.
Units are sections or biopsies, treated as separate people [Unverified for Govaere sections].
Usage: python -I 49_insitu_nfkb.py <new_out_dir>
"""
import sys
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
V7 = ROOT / "Analysis/Hypotheses/V7-20261021T235717Z"
DL = V7 / "downloads_untrusted"
GOVAERE_DIR = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external/govaere2026_natgenetics/cosmx")
CONTROL_PREFIXES = ("Negative", "NegPrb", "SystemControl", "FalseCode", "NegativeAdd")
SEED, K_NN, EDGE_MULT, K_WIDE = 0, 6, 2.5, 30
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
rng = np.random.default_rng(20261024)
N_NULL = 200
NFKB = ["NFKBIA", "NFKB1", "SOD2", "ICAM1", "IER3"]
MDB = ["KRT8", "KRT18", "SQSTM1"]
GENES = {"nfkb": NFKB + ["RELA", "BIRC3"], "injury": ["IL32", "TNFRSF12A", "ANXA2", "GDF15", "S100A6", "VIM", "CD44", "LGALS3",
         "CDKN1A", "MT2A", "MT1X", "LCN2", "SPP1", "TIMP1", "HSPB1"], "chemokine": ["CCL2", "CXCL8"], "ap1": ["JUN", "JUNB", "FOS", "ATF3"],
         "liver_function": ["FGG", "APOA1", "GLUL", "CRP"], "ligand": ["TNF", "TNFSF12", "IL1B"],
         "stromal_spill": ["ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "COL3A1"], "cholangiocyte_spill": ["KRT7", "KRT19"],
         "macrophage_spill": ["CD68", "CD163", "C1QA", "C1QB", "C1QC", "LYZ", "MARCO", "CD74", "HLA-DRA"]}
SENS_TESTS = ["NFKB_set", "IL32", "TNFRSF12A", "FGG", "HSPB1", "KRT7", "CD68", "__negprobe__"]
NBH_TESTS = ["NFKB_set", "IL32", "TNFRSF12A", "FGG", "HSPB1", "__MDB__", "CD68", "C1QA", "CD74", "COL1A1", "KRT7", "__negprobe__"]
WIDE_L = ["cholangiocyte", "stellate_fibroblast", "macrophage", "endothelial", "hepatocyte"]


def split_controls(df):
    neg_cols = [c for c in df.columns if str(c).startswith(CONTROL_PREFIXES)]
    negp = [c for c in neg_cols if str(c).lower().startswith("neg")]
    targets = [c for c in df.columns if c not in neg_cols]
    neg = df[negp].sum(axis=1).to_numpy(float) if negp else np.full(len(df), np.nan)
    return df[targets], neg, len(negp)


def load_atomx(expr_path, meta_path, unit_prefix):
    e = pd.read_csv(expr_path)
    e = e[e.cell_ID != 0]
    m = pd.read_csv(meta_path, usecols=["fov", "cell_ID", "CenterX_global_px", "CenterY_global_px"])
    e = e.merge(m, on=["fov", "cell_ID"], how="inner")
    obs = pd.DataFrame({"fov": e.fov.astype(str).values, "x": e.CenterX_global_px.values, "y": e.CenterY_global_px.values,
                        "cell": [f"{unit_prefix}_{f}_{c}" for f, c in zip(e.fov, e.cell_ID)]})
    tg, neg, n_negp = split_controls(e.drop(columns=["fov", "cell_ID", "CenterX_global_px", "CenterY_global_px"]))
    return tg, neg, obs, n_negp


def load(dataset):
    parts = []
    if dataset == "govaere":
        for gsm, name in [("GSM9351704", "Leuven_1"), ("GSM9351705", "Leuven_2"), ("GSM9351706", "Leuven_3"), ("GSM9351707", "Leuven_4")]:
            tg, neg, obs, n_negp = load_atomx(GOVAERE_DIR / f"{gsm}_{name}_exprMat_file.csv.gz",
                                              GOVAERE_DIR / f"{gsm}_{name}_metadata_file.csv.gz", name)
            obs["unit"] = name
            parts.append((tg, neg, obs, n_negp))
    elif dataset == "gse292268":
        for gsm, name in [("GSM8855042", "Slide1"), ("GSM8855043", "Slide2")]:
            d = DL / "gse292268"
            tg, neg, obs, n_negp = load_atomx(d / f"{gsm}_{name}_exprMat_file.csv.gz", d / f"{gsm}_{name}_metadata_file.csv.gz", name)
            fovpos = pd.read_csv(d / f"{gsm}_{name}_fov_positions_file.csv.gz")
            km = KMeans(n_clusters=2, random_state=SEED, n_init=10).fit(fovpos[["X_mm", "Y_mm"]].to_numpy())
            fov_unit = dict(zip(fovpos.FOV.astype(str), [f"{name}_biopsy{k + 1}" for k in km.labels_]))
            obs["unit"] = obs.fov.map(fov_unit)
            parts.append((tg, neg, obs, n_negp))
    elif dataset == "fukumoto":
        d = DL / "fukumoto"
        for name, ef, mf in [("sample1", "CosMx_expression_data_MASH_sample1_NT.csv", "CosMx_metadata_MASH_sample1_NT.csv"),
                             ("sample2", "CosMx_expression_data_MASH_sample_2_NT.csv", "CosMx_metadata_MASH_sample_2_NT.csv")]:
            e = pd.read_csv(d / ef, index_col=0)
            m = pd.read_csv(d / mf, index_col=0, usecols=lambda c: c in ("Unnamed: 0", "", "fov", "x_slide_mm", "y_slide_mm"))
            m = m.reindex(e.index)
            obs = pd.DataFrame({"fov": m.fov.astype(str).values, "x": m.x_slide_mm.values, "y": m.y_slide_mm.values,
                                "cell": [f"{name}_{c}" for c in e.index], "unit": name})
            tg, neg, n_negp = split_controls(e)
            parts.append((tg, neg, obs, n_negp))
    genes = sorted(set.intersection(*[set(p[0].columns) for p in parts]))
    X = sparse.vstack([sparse.csr_matrix(p[0][genes].to_numpy(np.float32)) for p in parts]).tocsr()
    obs = pd.concat([p[2] for p in parts], ignore_index=True)
    obs["neg_counts"] = np.concatenate([p[1] for p in parts])
    obs["fov"] = obs.unit + "_" + obs.fov
    return X, obs, pd.Index(genes)



def signflip(v):
    v = np.asarray(v, float)
    S = np.array(list(product([1, -1], repeat=len(v))))
    return np.mean(np.abs((S * v).mean(1)) >= abs(v.mean()) - 1e-12)


def resid_rank(M, Z):
    R = np.apply_along_axis(stats.rankdata, 0, M)
    return R - Z @ np.linalg.lstsq(Z, R, rcond=None)[0]


def cor_cols(a, M):
    a = (a - a.mean()) / a.std(); sd = M.std(0); sd[sd == 0] = np.nan
    return ((M - M.mean(0)) / sd * a[:, None]).mean(0)


rows, nbh, lm, nbhc = [], [], [], []
for dataset in ["govaere", "gse292268", "fukumoto"]:
    X, obs, genes = load(dataset)
    obs["total"] = np.asarray(X.sum(axis=1)).ravel()
    obs["n_genes"] = np.asarray((X > 0).sum(axis=1)).ravel()
    keep = ((obs.total >= 20) & (obs.n_genes >= 10)).to_numpy()
    X, obs = X[keep], obs[keep].reset_index(drop=True)
    labs = pd.read_csv(RUN / "out_insitu_detection" / dataset / "cell_labels_with_id.tsv.gz", sep="\t")
    assert len(labs) == len(obs) and (labs.cell.to_numpy() == obs.cell.to_numpy()).all()
    assert np.allclose(labs.total.to_numpy(), obs.total.to_numpy())
    obs["lineage"] = labs.lineage.to_numpy()
    lin = obs.lineage.to_numpy()
    scale = np.median(obs.total) / obs.total.to_numpy()
    Xn = (sparse.diags(scale) @ X).tocsr(); Xn.data = np.log1p(Xn.data)
    obs["neg_norm"] = np.log1p(obs.neg_counts.to_numpy() * scale)
    gi = {g: i for i, g in enumerate(genes)}
    core = np.zeros(len(obs), bool)
    wide = {l: np.full(len(obs), np.nan) for l in WIDE_L}
    for fov, idx in obs.groupby("fov").indices.items():
        if len(idx) <= K_WIDE:
            continue
        xy = obs.loc[idx, ["x", "y"]].to_numpy(float)
        d, j = cKDTree(xy).query(xy, k=K_WIDE + 1)
        d, j = d[:, 1:], j[:, 1:]
        med = np.median(d[:, 0])
        ok6 = d[:, :K_NN] <= EDGE_MULT * med
        nl = lin[idx][j]
        core[idx] = ((nl[:, :K_NN] == "hepatocyte") | ~ok6).all(axis=1) & ok6.any(axis=1)
        for l in WIDE_L:
            wide[l][idx] = (nl == l).mean(axis=1)
    # lineage means and detection, per unit
    show = [g for gs in GENES.values() for g in gs if g in gi] + [g for g in MDB if g in gi]
    for unit, uidx in obs.groupby("unit").indices.items():
        for l in np.unique(lin[uidx]):
            li = uidx[lin[uidx] == l]
            sub = Xn[li][:, [gi[g] for g in show]].toarray()
            lm.append(dict(dataset=dataset, unit=unit, lineage=l, n=len(li), **{f"mean_{g}": v for g, v in zip(show, sub.mean(0))},
                           **{f"det_{g}": v for g, v in zip(show, (sub > 0).mean(0))}, det_negprobe=(obs.neg_counts.to_numpy()[li] > 0).mean()))
    excl = set(show) | set(MDB)
    for unit, uidx in obs.groupby("unit").indices.items():
        hep_all = uidx[lin[uidx] == "hepatocyte"]
        for subset, hidx in [("core", hep_all[core[hep_all]]), ("all", hep_all)]:
            if len(hidx) < 500:
                rows.append(dict(dataset=dataset, unit=unit, subset=subset, n_cells=len(hidx), evaluable=False)); continue
            M = Xn[hidx].toarray()
            mu, det = M.mean(0), (M > 0).mean(0)
            cand = np.array([g not in excl and det[gi[g]] >= 0.01 for g in genes])
            cg = np.array(genes)[cand]
            bm = pd.qcut(pd.Series(mu[cand]).rank(method="first"), 10, labels=False).to_numpy()
            bd = pd.qcut(pd.Series(det[cand]).rank(method="first"), 10, labels=False).to_numpy()
            binof = dict(zip(cg, zip(bm, bd)))
            pool = pd.Series(list(zip(bm, bd)), index=cg)
            pool = {b: list(ix) for b, ix in pool.groupby(pool).groups.items()}

            def nearest_bin(g):  # test genes are excluded from the pool; place them by their own mean and detection
                rm = (mu[cand] < mu[gi[g]]).mean(); rd = (det[cand] < det[gi[g]]).mean()
                b = (min(int(rm * 10), 9), min(int(rd * 10), 9))
                if b in pool and len(pool[b]) >= 10:
                    return b
                ks = [k for k in pool if len(pool[k]) >= 10]
                return min(ks, key=lambda k: abs(k[0] - b[0]) + abs(k[1] - b[1]))
            sd = M.std(0); sd[sd == 0] = np.nan
            Zm = (M - mu) / sd
            Zc = np.column_stack([np.ones(len(hidx)), stats.rankdata(np.log(obs.total.to_numpy()[hidx])), stats.rankdata(np.log(obs.n_genes.to_numpy()[hidx]))])
            W = np.column_stack([wide[l][hidx] for l in WIDE_L[:-1]])
            Zw = np.column_stack([Zc, W[:, W.std(0) > 0]]) if np.isfinite(W).all() else None
            sc = lambda gs: np.nanmean(Zm[:, [gi[g] for g in gs]], 1)
            mdb_v = {"MDB": sc([g for g in MDB if g in gi]), "SQSTM1_only": Zm[:, gi["SQSTM1"]]}
            tests = {"NFKB_set": [g for g in NFKB if g in gi]}
            tests.update({f"NFKB_minus_{g}": [h for h in NFKB if h != g and h in gi] for g in NFKB})
            tests.update({g: [g] for g in show if g not in MDB})
            tests["__negprobe__"] = None
            neg_v = obs.neg_norm.to_numpy()[hidx]

            def neg_bin():
                rm = (mu[cand] < neg_v.mean()).mean(); rd = (det[cand] < (neg_v > 0).mean()).mean()
                b = (min(int(rm * 10), 9), min(int(rd * 10), 9))
                if b in pool and len(pool[b]) >= 10:
                    return b
                return min([k for k in pool if len(pool[k]) >= 10], key=lambda k: abs(k[0] - b[0]) + abs(k[1] - b[1]))
            Wd = {l: wide[l][hidx] for l in WIDE_L[:-1]}
            nbh_ok = subset == "core" and all(np.isfinite(w).all() and w.std() > 0 for w in Wd.values())
            Wr = {l: resid_rank(w[:, None], Zc)[:, 0] for l, w in Wd.items()} if nbh_ok else {}

            def nbh_rows(tname, ro, NR):
                for l, wl in Wr.items():
                    nbhc.append(dict(dataset=dataset, unit=unit, test=tname, lineage=l, n_cells=len(hidx),
                                     rho=np.corrcoef(wl, ro)[0, 1], null_mean=np.nanmean(cor_cols(wl, NR))))
            for adj, Zr in [("depth", Zc), ("depth+neighbourhood", Zw)]:
                if Zr is None:
                    continue
                for mname, mv in mdb_v.items():
                    if subset == "all" and (adj != "depth" or mname != "MDB"):
                        continue
                    main = adj == "depth" and mname == "MDB"
                    ri = resid_rank(mv[:, None], Zr)[:, 0]
                    for tname, gs in tests.items():
                        if not main and tname not in SENS_TESTS:
                            continue
                        if gs is None:
                            v = neg_v; size = 1; bins_used = [neg_bin()]
                        else:
                            v = sc(gs); size = len(gs); bins_used = [nearest_bin(g) for g in gs]
                        ro = resid_rank(v[:, None], Zr)[:, 0]
                        r = np.corrcoef(ri, ro)[0, 1]
                        NS = np.column_stack([np.nanmean(Zm[:, [gi[rng.choice(pool[b])] for b in bins_used]], 1) for _ in range(N_NULL)])
                        NR = resid_rank(NS, Zr)
                        nr = cor_cols(ri, NR)
                        rows.append(dict(dataset=dataset, unit=unit, subset=subset, adjust=adj, mdb=mname, test=tname, set_size=size,
                                         n_cells=len(hidx), evaluable=True, rho=r, null_mean=np.nanmean(nr), null_sd=np.nanstd(nr),
                                         mean_test=np.nanmean(M[:, [gi[g] for g in gs]]) if gs else np.nan,
                                         det_test=np.mean(M[:, [gi[g] for g in gs]] > 0) if gs else (v > 0).mean()))
                        if main and nbh_ok and tname in NBH_TESTS:
                            nbh_rows(tname, ro, NR)
                    if main and nbh_ok:
                        bins_used = [nearest_bin(g) for g in MDB if g in gi]
                        NS = np.column_stack([np.nanmean(Zm[:, [gi[rng.choice(pool[b])] for b in bins_used]], 1) for _ in range(N_NULL)])
                        nbh_rows("__MDB__", resid_rank(mdb_v["MDB"][:, None], Zc)[:, 0], resid_rank(NS, Zc))
            if subset == "core":
                hs = sc(tests["NFKB_set"])
                top = hs >= np.quantile(hs, 0.90)
                for l in WIDE_L:
                    w = wide[l][hidx]
                    nbh.append(dict(dataset=dataset, unit=unit, lineage=l, n_top=int(top.sum()), frac_top=np.nanmean(w[top]), frac_rest=np.nanmean(w[~top])))
                # whole-hepatocyte neighbourhood view: core and non-core hepatocytes together
                Ma = Xn[hep_all][:, [gi[g] for g in tests["NFKB_set"]]].toarray()
                hsa = np.nanmean((Ma - Ma.mean(0)) / Ma.std(0), 1)
                topa = hsa >= np.quantile(hsa, 0.90)
                for l in WIDE_L:
                    w = wide[l][hep_all]
                    nbh.append(dict(dataset=dataset, unit=unit, lineage=l + "__allhep", n_top=int(topa.sum()), frac_top=np.nanmean(w[topa]), frac_rest=np.nanmean(w[~topa])))
            print(dataset, unit, subset, len(hidx), flush=True)
    del X, Xn, obs

R = pd.DataFrame(rows)
R["excess"] = R.rho - R.null_mean
R.to_csv(OUT / "unit_coupling.tsv", sep="\t", index=False)
pd.DataFrame(nbh).to_csv(OUT / "neighbourhood_of_nfkb_high.tsv", sep="\t", index=False)
pd.DataFrame(lm).to_csv(OUT / "lineage_means_detection.tsv", sep="\t", index=False)
NB = pd.DataFrame(nbhc)
NB["excess"] = NB.rho - NB.null_mean
NB.to_csv(OUT / "unit_neighbourhood_coupling.tsv", sep="\t", index=False)
E = R[R.evaluable == True]
summ = []
for (subset, adj, mname, t), o in E.groupby(["subset", "adjust", "mdb", "test"]):
    summ.append(dict(subset=subset, adjust=adj, mdb=mname, test=t, n_units=len(o), mean_rho=o.rho.mean(), mean_null=o.null_mean.mean(),
                     mean_excess=o.excess.mean(), n_positive=int((o.excess > 0).sum()), p_signflip=signflip(o.excess.to_numpy())))
summ = pd.DataFrame(summ)
summ["bh"] = summ.groupby(["subset", "adjust", "mdb"]).p_signflip.transform(stats.false_discovery_control)
summ.to_csv(OUT / "summary_signflip.tsv", sep="\t", index=False)
nsum = pd.DataFrame([dict(test=t, lineage=l, n_units=len(o), mean_rho=o.rho.mean(), mean_null=o.null_mean.mean(), mean_excess=o.excess.mean(),
                          n_positive=int((o.excess > 0).sum()), p_signflip=signflip(o.excess.to_numpy()))
                     for (t, l), o in NB.groupby(["test", "lineage"])])
nsum["bh"] = stats.false_discovery_control(nsum.p_signflip)
nsum.to_csv(OUT / "summary_neighbourhood_signflip.tsv", sep="\t", index=False)
base = summ[(summ.subset == "core") & (summ.adjust == "depth") & (summ.mdb == "MDB")].set_index("test")
prim = pd.DataFrame([dict(id=i, test=f"core hepatocytes: {t} vs MDB, excess", est=base.loc[t, "mean_excess"], n_units=base.loc[t, "n_units"],
                          n_positive=base.loc[t, "n_positive"], p=base.loc[t, "p_signflip"])
                     for i, t in [("J1", "NFKB_set"), ("J2", "IL32"), ("J3", "TNFRSF12A")]])
o = np.argsort(prim.p.to_numpy()); ph = np.maximum.accumulate(prim.p.to_numpy()[o] * (len(o) - np.arange(len(o))))
prim.loc[prim.index[o], "p_holm"] = np.minimum(1, ph)
prim.to_csv(OUT / "primary_holm.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
print(prim.round(4).to_string())
print(summ.round(3).to_string())
print(E[(E.subset == "core") & (E.adjust == "depth") & (E.mdb == "MDB")].pivot(index="test", columns="unit", values="excess").round(3).to_string())
print(nsum.round(3).to_string())
nb = pd.DataFrame(nbh); nb["diff"] = nb.frac_top - nb.frac_rest
print(nb.pivot(index="lineage", columns="unit", values="diff").round(3).to_string())
