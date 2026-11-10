"""DIG1 job 38 (exploratory; follows jobs 36b, 37b). CosMx in situ: are the hepatocytes that carry
more heat-shock RNA the same hepatocytes that carry more KRT8/KRT18/SQSTM1 RNA, in intact tissue
without dissociation? Does HSPB1 go the other way, as a published study of ballooned hepatocytes
implies (doi 10.1038/srep22528)?
Data: Govaere (4 sections), GSE292268 (4 biopsies), Fukumoto (2 samples); 10 units. Loading and QC
as 03_insitu_position.py; lineage labels reused from out_insitu_detection/<dataset>/
cell_labels_with_id.tsv.gz (job 14b, same rules). Core hepatocytes = all 6 graph neighbours are
hepatocytes (V7 rule), to limit segmentation spillover.
Values: counts scaled to the median total, log1p, z within unit over core hepatocytes.
Sets (fixed): HSF1 = HSPA1A/B, BAG3, HSP90AA1 (on all panels); MDB = KRT8, KRT18, SQSTM1.
Statistic per unit: rho(set, MDB) with ranks residualised on log total and log n_genes, minus the
mean rho of 200 random panel gene sets of the same size matched on core-hepatocyte mean x detection
decile bins (test and control genes excluded; seed 20261024). Single genes the same, with 200
random single genes from the same bin.
Primary family (Holm over 2), exact sign-flip over the 10 units (1,024 patterns), two-sided:
  I1  HSF1 set excess.
  I2  HSPB1 excess.
Descriptive: per gene and unit (HSPA1A/B BAG3 HSP90AA1 HSP90AB1 HSPB1 CRYAB TUBB; spillover
controls ACTA2 TAGLN MYL9 COL1A1 COL1A2 COL3A1 KRT7 KRT19; negative probes); sensitivities: all
hepatocytes, MDB = SQSTM1 only, and wide-neighbourhood composition (30 nearest cells) added to the
residualisation; neighbourhood composition of the top 10% HSF1-score core hepatocytes vs the rest;
lineage means and detection.
Units are sections or biopsies, treated as separate people [Unverified for Govaere sections].
Usage: python -I 38_insitu_heatshock.py <new_out_dir>
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
HSF1 = ["HSPA1A/B", "BAG3", "HSP90AA1"]
MDB = ["KRT8", "KRT18", "SQSTM1"]
GENES = {"heat_shock": ["HSPA1A/B", "BAG3", "HSP90AA1", "HSP90AB1", "HSPB1", "CRYAB"], "cytoskeleton": ["TUBB"],
         "stromal_spill": ["ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "COL3A1"], "cholangiocyte_spill": ["KRT7", "KRT19"]}
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


rows, nbh, lm = [], [], []
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
            tests = {"HSF1_set": [g for g in HSF1 if g in gi]}
            tests.update({g: [g] for g in show if g not in MDB})
            tests["__negprobe__"] = None
            for adj, Zr in [("depth", Zc), ("depth+neighbourhood", Zw)]:
                if Zr is None:
                    continue
                for mname, mv in mdb_v.items():
                    if subset == "all" and (adj != "depth" or mname != "MDB"):
                        continue
                    ri = resid_rank(mv[:, None], Zr)[:, 0]
                    for tname, gs in tests.items():
                        if gs is None:
                            v = obs.neg_norm.to_numpy()[hidx]; size = 1
                            mneg, dneg = v.mean(), (v > 0).mean()
                            rm = (mu[cand] < mneg).mean(); rd = (det[cand] < dneg).mean()
                            bins_used = [(min(int(rm * 10), 9), min(int(rd * 10), 9))]
                            bins_used = [b if b in pool and len(pool[b]) >= 10 else min([k for k in pool if len(pool[k]) >= 10], key=lambda k: abs(k[0] - b[0]) + abs(k[1] - b[1])) for b in bins_used]
                        else:
                            v = sc(gs); size = len(gs); bins_used = [nearest_bin(g) for g in gs]
                        ro = resid_rank(v[:, None], Zr)[:, 0]
                        r = np.corrcoef(ri, ro)[0, 1]
                        NS = np.column_stack([np.nanmean(Zm[:, [gi[rng.choice(pool[b])] for b in bins_used]], 1) for _ in range(N_NULL)])
                        nr = cor_cols(ri, resid_rank(NS, Zr))
                        rows.append(dict(dataset=dataset, unit=unit, subset=subset, adjust=adj, mdb=mname, test=tname, set_size=size,
                                         n_cells=len(hidx), evaluable=True, rho=r, null_mean=np.nanmean(nr), null_sd=np.nanstd(nr),
                                         mean_test=np.nanmean(M[:, [gi[g] for g in gs]]) if gs else np.nan,
                                         det_test=np.mean(M[:, [gi[g] for g in gs]] > 0) if gs else (v > 0).mean()))
            if subset == "core":
                hs = sc(tests["HSF1_set"])
                top = hs >= np.quantile(hs, 0.90)
                for l in WIDE_L:
                    w = wide[l][hidx]
                    nbh.append(dict(dataset=dataset, unit=unit, lineage=l, n_top=int(top.sum()), frac_top=np.nanmean(w[top]), frac_rest=np.nanmean(w[~top])))
                # whole-hepatocyte neighbourhood view: core and non-core hepatocytes together
                Ma = Xn[hep_all][:, [gi[g] for g in tests["HSF1_set"]]].toarray()
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
pd.DataFrame(nbh).to_csv(OUT / "neighbourhood_of_hsf1_high.tsv", sep="\t", index=False)
pd.DataFrame(lm).to_csv(OUT / "lineage_means_detection.tsv", sep="\t", index=False)
E = R[R.evaluable == True]
summ = []
for (subset, adj, mname, t), o in E.groupby(["subset", "adjust", "mdb", "test"]):
    summ.append(dict(subset=subset, adjust=adj, mdb=mname, test=t, n_units=len(o), mean_rho=o.rho.mean(), mean_null=o.null_mean.mean(),
                     mean_excess=o.excess.mean(), n_positive=int((o.excess > 0).sum()), p_signflip=signflip(o.excess.to_numpy())))
summ = pd.DataFrame(summ)
summ["bh"] = summ.groupby(["subset", "adjust", "mdb"]).p_signflip.transform(stats.false_discovery_control)
summ.to_csv(OUT / "summary_signflip.tsv", sep="\t", index=False)
base = summ[(summ.subset == "core") & (summ.adjust == "depth") & (summ.mdb == "MDB")].set_index("test")
prim = pd.DataFrame([dict(id="I1", test="core hepatocytes: HSF1 set vs MDB, excess", est=base.loc["HSF1_set", "mean_excess"],
                          n_units=base.loc["HSF1_set", "n_units"], n_positive=base.loc["HSF1_set", "n_positive"], p=base.loc["HSF1_set", "p_signflip"]),
                     dict(id="I2", test="core hepatocytes: HSPB1 vs MDB, excess", est=base.loc["HSPB1", "mean_excess"],
                          n_units=base.loc["HSPB1", "n_units"], n_positive=base.loc["HSPB1", "n_positive"], p=base.loc["HSPB1", "p_signflip"])])
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
nb = pd.DataFrame(nbh); nb["diff"] = nb.frac_top - nb.frac_rest
print(nb.pivot(index="lineage", columns="unit", values="diff").round(3).to_string())
