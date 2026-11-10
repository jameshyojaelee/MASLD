"""DIG1 job 54 (exploratory, post hoc to job 53). In intact tissue, are the two injured-hepatocyte states of
job 53 (an NF-kB chemokine state marked by CCL20, a biliary-like state marked by SOX9, and the CXCL10
interferon state) coherent within hepatocytes, and where do they sit?
Data, QC, labels, core-hepatocyte rule and matched-set null as job 49 (Govaere 4 sections, GSE292268 4
biopsies, Fukumoto 2 samples; units with >= 500 core hepatocytes).
State sets (on all three panels; chosen from job 53 cross-cohort genes before any in situ data were read):
  NFKB_CHEM = CXCL8, ICAM1, NFKBIA, IER3, CX3CL1, CSF1   (up in CCL20+ hepatocytes in Flex and snRNA)
  BILIARY   = EPCAM, FGFR2, CD24, SPP1, KRT23, ANXA4    (up in SOX9+ hepatocytes in Flex and snRNA)
  IFN       = CXCL9, STAT1                              (up in CXCL10+ Flex hepatocytes)
Coupling per unit (core hepatocytes): rho(anchor, set) on ranks residualised on log total, log n_genes and
the KRT8/KRT18/SQSTM1 score ("beyond injury"), minus the mean rho with 200 random panel sets matched on
core-hepatocyte mean x detection deciles (listed genes excluded; seed 20261024).
Location per unit (core hepatocytes): rho(score, fraction of a lineage among the 30 nearest cells), both
residualised on depth, minus the same for the matched null; "t_nk|portal" also residualises on the
cholangiocyte and stellate/fibroblast fractions.
Detection gate: a unit enters a primary test only if, in its core hepatocytes, the anchor and at least half
of the set genes are detected more often than the summed negative probes. Fewer than 3 gated units =
untestable.
Primary families, exact two-sided sign-flip over gated units:
  S (Holm over 3): S1 CCL20 x NFKB_CHEM, S2 SOX9 x BILIARY, S3 CXCL10 x IFN (beyond injury).
  L (Holm over 3): L1 BILIARY vs cholangiocyte fraction, L2 NFKB_CHEM vs cholangiocyte fraction,
                   L3 IFN vs t_nk|portal fraction (depth adjusted).
Descriptive: cross couplings (CCL20 x BILIARY, SOX9 x NFKB_CHEM, ...), negative-probe and spillover
anchors (KRT7, KRT19 x BILIARY; CD68 x NFKB_CHEM; CD3E x IFN), all lineages, beyond-injury location,
Fukumoto alone.
Usage: python -I 54_insitu_two_states.py <new_out_dir>
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
MDB = ["KRT8", "KRT18", "SQSTM1"]
STATE = {"NFKB_CHEM": ["CXCL8", "ICAM1", "NFKBIA", "IER3", "CX3CL1", "CSF1"],
         "BILIARY": ["EPCAM", "FGFR2", "CD24", "SPP1", "KRT23", "ANXA4"],
         "IFN": ["CXCL9", "STAT1"]}
ANCH = ["CCL20", "SOX9", "CXCL10"]
SPILL = ["KRT7", "KRT19", "CD68", "C1QA", "COL1A1", "CD3E"]
COUPLE = [(a, s) for a in ANCH for s in STATE] + [("CCL20", "CXCL10"), ("CCL20", "SOX9"), ("SOX9", "CXCL10")] + \
         [("__negprobe__", s) for s in STATE] + [("KRT7", "BILIARY"), ("KRT19", "BILIARY"), ("CD68", "NFKB_CHEM"), ("CD3E", "IFN")]
LOC_TESTS = ANCH + list(STATE) + ["__MDB__", "__negprobe__"] + SPILL
LIN = ["cholangiocyte", "stellate_fibroblast", "macrophage", "t_nk", "endothelial", "b_plasma"]
PRIM_S = {"S1": ("CCL20", "NFKB_CHEM"), "S2": ("SOX9", "BILIARY"), "S3": ("CXCL10", "IFN")}
PRIM_L = {"L1": ("BILIARY", "cholangiocyte"), "L2": ("NFKB_CHEM", "cholangiocyte"), "L3": ("IFN", "t_nk|portal")}


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




def holm(p):
    p = np.asarray(p, float); out = np.full(len(p), np.nan); ok = np.where(np.isfinite(p))[0]
    if len(ok):
        o = ok[np.argsort(p[ok])]; ph = np.maximum.accumulate(p[o] * (len(o) - np.arange(len(o))))
        out[o] = np.minimum(1, ph)
    return out


coup, loc, det_rows = [], [], []
for dataset in ["govaere", "gse292268", "fukumoto"]:
    X, obs, genes = load(dataset)
    obs["total"] = np.asarray(X.sum(axis=1)).ravel()
    obs["n_genes"] = np.asarray((X > 0).sum(axis=1)).ravel()
    keep = ((obs.total >= 20) & (obs.n_genes >= 10)).to_numpy()
    X, obs = X[keep], obs[keep].reset_index(drop=True)
    labs = pd.read_csv(RUN / "out_insitu_detection" / dataset / "cell_labels_with_id.tsv.gz", sep="\t")
    assert len(labs) == len(obs) and (labs.cell.to_numpy() == obs.cell.to_numpy()).all()
    assert np.allclose(labs.total.to_numpy(), obs.total.to_numpy())
    lin = labs.lineage.to_numpy()
    scale = np.median(obs.total) / obs.total.to_numpy()
    Xn = (sparse.diags(scale) @ X).tocsr(); Xn.data = np.log1p(Xn.data)
    neg_norm = np.log1p(obs.neg_counts.to_numpy() * scale)
    gi = {g: i for i, g in enumerate(genes)}
    core = np.zeros(len(obs), bool)
    wide = {l: np.full(len(obs), np.nan) for l in LIN}
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
        for l in LIN:
            wide[l][idx] = (nl == l).mean(axis=1)
    listed = [g for g in ANCH + sum(STATE.values(), []) + SPILL + MDB if g in gi]
    excl = set(listed)
    for unit, uidx in obs.groupby("unit").indices.items():
        neg_det_unit = {}
        for l in np.unique(lin[uidx]):
            li = uidx[lin[uidx] == l]
            sub = Xn[li][:, [gi[g] for g in listed]].toarray()
            for g, dv, mv in zip(listed, (sub > 0).mean(0), sub.mean(0)):
                det_rows.append(dict(dataset=dataset, unit=unit, cells=l, n=len(li), gene=g, det=dv, mean=mv))
            det_rows.append(dict(dataset=dataset, unit=unit, cells=l, n=len(li), gene="__negprobe__", det=(obs.neg_counts.to_numpy()[li] > 0).mean(),
                                 mean=neg_norm[li].mean()))
        hidx = uidx[(lin[uidx] == "hepatocyte") & core[uidx]]
        if len(hidx) < 500:
            continue
        M = Xn[hidx].toarray()
        mu, det = M.mean(0), (M > 0).mean(0)
        negv = neg_norm[hidx]; neg_det = (negv > 0).mean()
        for g in listed:
            det_rows.append(dict(dataset=dataset, unit=unit, cells="core_hepatocyte", n=len(hidx), gene=g, det=det[gi[g]], mean=mu[gi[g]]))
        det_rows.append(dict(dataset=dataset, unit=unit, cells="core_hepatocyte", n=len(hidx), gene="__negprobe__", det=neg_det, mean=negv.mean()))
        cand = np.array([g not in excl and det[gi[g]] >= 0.01 for g in genes])
        cg = np.array(genes)[cand]
        bm = pd.qcut(pd.Series(mu[cand]).rank(method="first"), 10, labels=False).to_numpy()
        bd = pd.qcut(pd.Series(det[cand]).rank(method="first"), 10, labels=False).to_numpy()
        pool = pd.Series(list(zip(bm, bd)), index=cg)
        pool = {b: list(ix) for b, ix in pool.groupby(pool).groups.items()}

        def bin_of(m_, d_):
            b = (min(int((mu[cand] < m_).mean() * 10), 9), min(int((det[cand] < d_).mean() * 10), 9))
            if b in pool and len(pool[b]) >= 10:
                return b
            return min([k for k in pool if len(pool[k]) >= 10], key=lambda k: abs(k[0] - b[0]) + abs(k[1] - b[1]))
        sd = M.std(0); sd[sd == 0] = np.nan
        Zm = np.nan_to_num((M - mu) / sd)
        Zc = np.column_stack([np.ones(len(hidx)), stats.rankdata(np.log(obs.total.to_numpy()[hidx])), stats.rankdata(np.log(obs.n_genes.to_numpy()[hidx]))])
        mdb = np.mean(Zm[:, [gi[g] for g in MDB if g in gi]], 1)
        Zi = np.column_stack([Zc, stats.rankdata(mdb)])
        W = {l: wide[l][hidx] for l in LIN}
        Zn = np.column_stack([Zi] + [W[l] for l in LIN if np.isfinite(W[l]).all() and W[l].std() > 0])

        def vec(name):  # score vector and the bins its null sets are drawn from
            if name == "__negprobe__":
                return negv, [bin_of(negv.mean(), (negv > 0).mean())], True
            if name == "__MDB__":
                gs = [g for g in MDB if g in gi]
            elif name in STATE:
                gs = [g for g in STATE[name] if g in gi]
            else:
                gs = [name] if name in gi else []
            if not gs:
                return None, None, False
            gate = np.mean([det[gi[g]] > neg_det for g in gs]) >= 0.5
            return np.mean(Zm[:, [gi[g] for g in gs]], 1), [bin_of(mu[gi[g]], det[gi[g]]) for g in gs], gate

        def null_scores(bins):
            return np.column_stack([np.mean(Zm[:, [gi[rng.choice(pool[b])] for b in bins]], 1) for _ in range(N_NULL)])
        for a, s in COUPLE:
            va, _, ga = vec(a); vs, bs, gs_ = vec(s)
            if va is None or vs is None:
                continue
            NS = null_scores(bs)
            for adj, Zr in [("depth", Zc), ("beyond_injury", Zi), ("beyond_injury+nbh", Zn)]:
                ra = resid_rank(va[:, None], Zr)[:, 0]; rs = resid_rank(vs[:, None], Zr)[:, 0]
                nr = cor_cols(ra, resid_rank(NS, Zr))
                coup.append(dict(dataset=dataset, unit=unit, anchor=a, partner=s, adjust=adj, n_cells=len(hidx), gate=bool(ga and gs_),
                                 rho=np.corrcoef(ra, rs)[0, 1], null_mean=np.nanmean(nr), null_sd=np.nanstd(nr)))
        Wok = {l: w for l, w in W.items() if np.isfinite(w).all() and w.std() > 0}
        portal = [W[l] for l in ["cholangiocyte", "stellate_fibroblast"] if l in Wok]
        for t in LOC_TESTS:
            vt, bt, gt = vec(t)
            if vt is None:
                continue
            NS = null_scores(bt)
            for adj, Zr in [("depth", Zc), ("beyond_injury", Zi)]:
                rt = resid_rank(vt[:, None], Zr)[:, 0]; RN = resid_rank(NS, Zr)
                targets = {l: (w, Zr) for l, w in Wok.items()}
                if "t_nk" in Wok and portal:
                    targets["t_nk|portal"] = (Wok["t_nk"], np.column_stack([Zr] + portal))
                for l, (w, Zw) in targets.items():
                    if l == "t_nk|portal":
                        rt2 = resid_rank(vt[:, None], Zw)[:, 0]; RN2 = resid_rank(NS, Zw)
                    else:
                        rt2, RN2 = rt, RN
                    wl = resid_rank(w[:, None], Zw)[:, 0]
                    loc.append(dict(dataset=dataset, unit=unit, test=t, lineage=l, adjust=adj, n_cells=len(hidx), gate=bool(gt),
                                    rho=np.corrcoef(wl, rt2)[0, 1], null_mean=np.nanmean(cor_cols(wl, RN2))))
        print(dataset, unit, len(hidx), flush=True)
    del X, Xn, obs

C = pd.DataFrame(coup); C["excess"] = C.rho - C.null_mean
L = pd.DataFrame(loc); L["excess"] = L.rho - L.null_mean
D = pd.DataFrame(det_rows)
C.to_csv(OUT / "unit_coupling.tsv", sep="\t", index=False)
L.to_csv(OUT / "unit_location.tsv", sep="\t", index=False)
D.to_csv(OUT / "detection_by_unit.tsv", sep="\t", index=False)


def sflip(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    return signflip(v) if len(v) >= 3 else np.nan


def summarise(T, keys):
    rows = []
    for k, o in T.groupby(keys):
        for grp, oo in [("all", o), ("gated", o[o.gate]), ("fukumoto", o[o.dataset == "fukumoto"])]:
            rows.append(dict(zip(keys, k), units=grp, n_units=len(oo), mean_excess=oo.excess.mean(), n_positive=int((oo.excess > 0).sum()),
                             p_signflip=sflip(oo.excess) if grp != "fukumoto" else np.nan))
    return pd.DataFrame(rows)


SC = summarise(C, ["anchor", "partner", "adjust"]); SC.to_csv(OUT / "summary_coupling.tsv", sep="\t", index=False)
SL = summarise(L, ["test", "lineage", "adjust"]); SL.to_csv(OUT / "summary_location.tsv", sep="\t", index=False)
prim = []
for i, (a, s) in PRIM_S.items():
    r = SC[(SC.anchor == a) & (SC.partner == s) & (SC.adjust == "beyond_injury") & (SC.units == "gated")]
    prim.append(dict(id=i, test=f"{a} x {s}, beyond injury", **(r.iloc[0][["n_units", "mean_excess", "n_positive", "p_signflip"]].to_dict() if len(r) else {})))
for i, (t, l) in PRIM_L.items():
    r = SL[(SL.test == t) & (SL.lineage == l) & (SL.adjust == "depth") & (SL.units == "gated")]
    prim.append(dict(id=i, test=f"{t} vs {l} fraction", **(r.iloc[0][["n_units", "mean_excess", "n_positive", "p_signflip"]].to_dict() if len(r) else {})))
P = pd.DataFrame(prim)
P["p_holm"] = np.nan
for fam in ["S", "L"]:
    m = P.id.str.startswith(fam).to_numpy()
    P.loc[m, "p_holm"] = holm(P.loc[m, "p_signflip"].to_numpy())
P.to_csv(OUT / "primary_holm.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 600)
print(P.round(4).to_string())
dc = D[D.cells == "core_hepatocyte"].pivot_table(index="gene", columns="unit", values="det")
print(dc.round(3).to_string())
print(SC[SC.units != "fukumoto"].round(3).to_string())
print(SC[SC.units == "fukumoto"].round(3).to_string())
print(C[C.adjust == "beyond_injury"].pivot_table(index=["anchor", "partner"], columns="unit", values="excess").round(3).to_string())
print(SL[(SL.adjust == "depth")].round(3).to_string())
print(L[(L.adjust == "depth")].pivot_table(index=["test", "lineage"], columns="unit", values="excess").round(3).to_string())
