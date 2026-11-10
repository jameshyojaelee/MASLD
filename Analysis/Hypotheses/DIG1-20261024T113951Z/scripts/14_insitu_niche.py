"""DIG1 job 14 (exploratory, descriptive). CosMx niche test: are hepatocytes whose fibroblast
neighbours carry more FGF7 higher in SOX9, EPCAM, KRT8, CDKN1A and FGFR2, beyond neighbour
composition, lobular position and depth? Question from jobs 1 and 12: fibroblast FGF7 rises with
MASLD in both single-cell cohorts and tracks hepatocyte SOX9 across donors (Flex rho 0.79, V3 0.62).
Mouse precedent: Takase et al. 2013 Genes Dev (doi:10.1101/gad.204776.112).

Loading, QC and lineage labels are lines 12-138 of 03_insitu_position.py, unchanged (V7 rules).
Neighbourhood: 30 nearest cells within FOV (cap 8 x median first-neighbour distance). Exposure for
ligand L from sender S: mean log1p normalised L over S cells in the neighbourhood, in hepatocytes
with >= 1 S neighbour. Statistic per unit: partial Spearman of exposure with the hepatocyte
outcome; ranks residualised on log1p counts of S, cholangiocyte, macrophage and endothelial
neighbours, hepatocyte log total counts, mean log total of the S neighbours, log distance to the
nearest cholangiocyte and lobular coordinate (when the unit has >= 50 cholangiocytes), with FOV
fixed effects. Null: the same statistic for up to 200 other S-enriched panel genes (S mean >= 2x
hepatocyte mean, S detection >= 5%); null_pct = fraction of null genes below the observed rho.
Subsets: all eligible hepatocytes and core hepatocytes (all six nearest neighbours hepatocytes),
which removes direct spillover. Spillover outcomes (KRT7, KRT19; COL1A1, COL1A2, LUM, DCN) and
hepatocyte reference genes (APOA1, SERPINA1) are reported alongside.
Read-out fixed before running: FGF7 -> SOX9 and FGFR2 in core hepatocytes; the link is supported
only if rho > 0 and null_pct >= 0.95 in most units with >= 200 eligible core hepatocytes, and the
cholangiocyte spillover genes do not show the same pattern. Units are sections or biopsies, not
donors with disease labels; nothing here is a disease contrast.
Usage: python -I 14_insitu_niche.py {govaere|fukumoto|gse292268}
"""
import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
V7 = ROOT / "Analysis/Hypotheses/V7-20261021T235717Z"
DL = V7 / "downloads_untrusted"
GOVAERE_DIR = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external/govaere2026_natgenetics/cosmx")
HOT = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
V3_SCRIPT = ROOT / "scripts/analysis/hypotheses_2026q4/V3_program_holdout/annotate_pseudobulk.py"
STROMAL_ECM = "hotspot_hepatocytes_f05c535ae5bbc0b9"
DUCTULAR = "hotspot_hepatocytes_48f39dd4d817a10e"
CONTROL_PREFIXES = ("Negative", "NegPrb", "SystemControl", "FalseCode", "NegativeAdd")
SEED, K_NN, EDGE_MULT, K_WIDE = 0, 6, 2.5, 30
GENES = {"member": ["FGF13", "MIR4435-2HG", "ITGAV", "EPHA3"],
         "cholangiocyte_spill": ["KRT7", "KRT19", "EPCAM", "SOX9", "CD24", "SPP1", "ANXA4"],
         "stromal_spill": ["COL1A1", "COL1A2", "COL3A1", "LUM", "DCN"],
         "hepatocyte_ref": ["APOA1", "SERPINA1", "CRP"]}

dataset = sys.argv[1]
OUT = RUN / "out_insitu_niche" / dataset
OUT.mkdir(parents=True, exist_ok=False)


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


X, obs, genes = load(dataset)
obs["total"] = np.asarray(X.sum(axis=1)).ravel()
obs["n_genes"] = np.asarray((X > 0).sum(axis=1)).ravel()
keep = ((obs.total >= 20) & (obs.n_genes >= 10)).to_numpy()
X, obs = X[keep], obs[keep].reset_index(drop=True)
mem = pd.read_csv(HOT / "program_membership_v2.tsv", sep="\t")
mem = mem[mem.mapped_symbol_status == "gencode_v49_unique_symbol_confirmed"]
members_conf = set(mem[mem.program_uid.isin([STROMAL_ECM, DUCTULAR])].mapped_symbol)

spec = importlib.util.spec_from_file_location("v3_annotate", V3_SCRIPT)
v3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v3)
A = ad.AnnData(X=X.copy(), obs=obs.set_index("cell"), var=pd.DataFrame(index=genes))
sc.pp.normalize_total(A, target_sum=float(np.median(obs.total)))
sc.pp.log1p(A)
A.var["use_pca"] = ~A.var_names.isin(members_conf)
sc.tl.pca(A, n_comps=30, mask_var="use_pca", random_state=SEED)
sc.pp.neighbors(A, n_neighbors=15, random_state=SEED)
sc.tl.leiden(A, resolution=1.0, random_state=SEED, flavor="igraph", n_iterations=2, directed=False)
markers = {k: [g for g in v if g in genes and g not in members_conf] for k, v in v3.MARKERS.items()}
markers = {k: v for k, v in markers.items() if v}
for k, v in markers.items():
    sc.tl.score_genes(A, v, score_name=f"mk_{k}", random_state=SEED)
cols = [f"mk_{k}" for k in markers]
cl = A.obs.groupby("leiden", observed=True)[cols].mean()
zc = (cl[cols] - cl[cols].mean()) / cl[cols].std(ddof=1)
cl["lineage"] = zc.idxmax(axis=1).str.replace("mk_", "", regex=False).where(zc.max(axis=1) >= 0.5, "unassigned")
obs["lineage"] = A.obs.leiden.map(cl.lineage).astype(str).values
del A
ct = pd.crosstab(obs.unit, obs.lineage)
ct.to_csv(OUT / "cells_per_unit_lineage.tsv", sep="\t")
v7_ct = pd.read_csv(V7 / "outputs" / dataset / "cells_per_unit_lineage.tsv", sep="\t", index_col=0)
(ct.reindex_like(v7_ct).fillna(0) - v7_ct).to_csv(OUT / "label_count_difference_vs_V7.tsv", sep="\t")

scale = np.median(obs.total) / obs.total.to_numpy()
Xn = (sparse.diags(scale) @ X).tocsr()
gi = {g: i for i, g in enumerate(genes)}
obs["neg_norm"] = obs.neg_counts.to_numpy() * scale
lin = obs.lineage.to_numpy()

# ---- job 14 niche analysis ----
obs[["unit", "fov", "x", "y", "lineage", "total"]].to_csv(OUT / "cell_labels.tsv.gz", sep="\t", index=False)
PAIRS = [("FGF7", "stellate_fibroblast"), ("HGF", "stellate_fibroblast"), ("CXCL12", "stellate_fibroblast"),
         ("TNFSF12", "macrophage")]
OUTCOMES = {"plasticity": ["SOX9", "EPCAM", "KRT8", "CDKN1A", "AKR1B10"],
            "receptor": ["FGFR2", "MET", "ACKR3", "CXCR4", "TNFRSF12A"],
            "cholangiocyte_spill": ["KRT7", "KRT19"],
            "stromal_spill": ["COL1A1", "COL1A2", "LUM", "DCN"],
            "hepatocyte_ref": ["APOA1", "SERPINA1"]}
NB_CAP, MIN_HEP, N_NULL_MAX = 8.0, 200, 200
rng = np.random.default_rng(20261024)
logX = Xn.copy()
logX.data = np.log1p(logX.data)
logtot = np.log(obs.total.to_numpy(float))

nb = np.full((len(obs), K_WIDE), -1, dtype=np.int64)
core = np.zeros(len(obs), bool)
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
    jj = idx[j]
    jj[d > NB_CAP * med] = -1
    nb[idx] = jj


def count_nb(rows, lineage):
    n = nb[rows]
    return ((n >= 0) & (lin[np.where(n >= 0, n, 0)] == lineage)).sum(axis=1)


def resid(M, C, groups):
    """Residualise columns of M on C with group fixed effects (within-group demeaning)."""
    M = M - pd.DataFrame(M).groupby(groups).transform("mean").to_numpy()
    C = C - pd.DataFrame(C).groupby(groups).transform("mean").to_numpy()
    beta = np.linalg.lstsq(C, M, rcond=None)[0]
    return M - C @ beta


def rank_cols(M):
    return np.column_stack([stats.rankdata(M[:, k]) for k in range(M.shape[1])])


hep_mean = np.asarray(logX[lin == "hepatocyte"].mean(axis=0)).ravel()
outs = [(g, role) for role, gs in OUTCOMES.items() for g in gs if g in gi]
pd.DataFrame(outs, columns=["gene", "role"]).to_csv(OUT / "outcomes_used.tsv", sep="\t", index=False)
gl = Xn[:, gi["GLUL"]].toarray().ravel() if "GLUL" in gi else None
rows = []
for unit, uidx in obs.groupby("unit").indices.items():
    hep = uidx[lin[uidx] == "hepatocyte"]
    chol = uidx[lin[uidx] == "cholangiocyte"]
    xy = obs[["x", "y"]].to_numpy(float)
    covs_pos = []
    if len(chol) >= 50:
        dp = cKDTree(xy[chol]).query(xy[hep], k=1)[0]
        central = hep[gl[hep] >= np.quantile(gl[hep], 0.95)]
        dc = cKDTree(xy[central]).query(xy[hep], k=1)[0]
        covs_pos = [np.log(dp + 1e-6), dp / (dp + dc + 1e-12)]
    for lig, S in PAIRS:
        if lig not in gi:
            continue
        nS = count_nb(hep, S)
        el = nS >= 1
        h = hep[el]
        if len(h) < MIN_HEP:
            rows.append(dict(unit=unit, ligand=lig, sender=S, n=len(h), evaluable=False))
            continue
        n = nb[h]
        m = (n >= 0) & (lin[np.where(n >= 0, n, 0)] == S)
        r_i = np.repeat(np.arange(len(h)), m.sum(axis=1))
        W = sparse.csr_matrix((1.0 / m.sum(axis=1)[r_i], (r_i, n[m])), shape=(len(h), len(obs)))
        Sm = lin == S
        s_mean = np.asarray(logX[Sm].mean(axis=0)).ravel()
        s_det = np.asarray((logX[Sm] > 0).mean(axis=0)).ravel()
        excl = {lig} | {g for g, _ in outs} | {p[0] for p in PAIRS}
        cand = [g for g in genes if s_mean[gi[g]] >= 2 * hep_mean[gi[g]] and s_det[gi[g]] >= 0.05 and g not in excl]
        if len(cand) > N_NULL_MAX:
            cand = list(rng.choice(cand, N_NULL_MAX, replace=False))
        E = (W @ logX[:, [gi[lig]] + [gi[g] for g in cand]]).toarray()
        C = np.column_stack([np.log1p(nS[el]), np.log1p(count_nb(h, "cholangiocyte")), np.log1p(count_nb(h, "macrophage")),
                             np.log1p(count_nb(h, "endothelial")), logtot[h], W @ logtot] + [c[el] for c in covs_pos])
        Y = logX[h][:, [gi[g] for g, _ in outs]].toarray()
        for subset, keep_s in [("all", np.ones(len(h), bool)), ("core", core[h])]:
            if keep_s.sum() < MIN_HEP:
                rows.append(dict(unit=unit, ligand=lig, sender=S, subset=subset, n=int(keep_s.sum()), evaluable=False))
                continue
            grp = obs.fov.to_numpy()[h][keep_s]
            Re = resid(rank_cols(E[keep_s]), C[keep_s], grp)
            Ry = resid(rank_cols(Y[keep_s]), C[keep_s], grp)
            Re /= np.linalg.norm(Re, axis=0, keepdims=True) + 1e-12
            Ry /= np.linalg.norm(Ry, axis=0, keepdims=True) + 1e-12
            R = Re.T @ Ry  # (1 + n_null) x outcomes
            for k, (g, role) in enumerate(outs):
                null = R[1:, k]
                rows.append(dict(unit=unit, ligand=lig, sender=S, subset=subset, evaluable=True, n=int(keep_s.sum()),
                                 position_covariates=bool(covs_pos), outcome=g, role=role,
                                 outcome_detect=float((Y[keep_s, k] > 0).mean()), rho=R[0, k],
                                 null_median=float(np.median(null)), null_pct=float((null < R[0, k]).mean()), n_null=len(null)))
res = pd.DataFrame(rows)
res.to_csv(OUT / "niche_partial_spearman.tsv", sep="\t", index=False)
sel = res[(res.evaluable == True) & (res.ligand == "FGF7")]
print(sel[["unit", "subset", "n", "outcome", "role", "rho", "null_median", "null_pct"]].round(3).to_string(index=False), flush=True)
