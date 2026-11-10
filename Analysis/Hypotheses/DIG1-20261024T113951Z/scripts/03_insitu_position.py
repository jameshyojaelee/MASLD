"""DIG1 job 3 (exploratory; EXPLORATORY_PLAN.md question A5). CosMx in situ position of
hepatocytes that express Ductular injury members (FGF13, MIR4435-2HG, ITGAV) and EPHA3.

Loading, QC and lineage labels are copied from V7 01_phase1.py (lines 40-147, unchanged rules).
Lobular coordinate per hepatocyte: L = d_portal / (d_portal + d_central), with cholangiocytes as the
portal anchor and the top 5% GLUL hepatocytes of the unit as the central anchor (GLUL is the only
zonation marker on these panels). Core hepatocytes (all graph neighbours hepatocytes) remove
direct segmentation spillover; spillover-control genes test that. Sections or biopsies are units;
everything is descriptive.
Usage: python -I 03_insitu_position.py {govaere|fukumoto|gse292268}
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
OUT = RUN / "out_insitu" / dataset
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

# core hepatocytes (V7 rule) and wide-neighbourhood composition, within FOV
core = np.zeros(len(obs), bool)
wide = {l: np.full(len(obs), np.nan) for l in ["cholangiocyte", "stellate_fibroblast", "macrophage", "endothelial", "hepatocyte"]}
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
    for l in wide:
        wide[l][idx] = (nl == l).mean(axis=1)

gl = Xn[:, gi["GLUL"]].toarray().ravel()
gene_list = [(g, role) for role, gs in GENES.items() for g in gs if g in gi]
pd.DataFrame(gene_list, columns=["gene", "role"]).to_csv(OUT / "genes_used.tsv", sep="\t", index=False)
rows, prof, nb = [], [], []
for unit, uidx in obs.groupby("unit").indices.items():
    hep = uidx[lin[uidx] == "hepatocyte"]
    chol = uidx[lin[uidx] == "cholangiocyte"]
    if len(hep) < 500 or len(chol) < 50:
        rows.append(dict(unit=unit, n_hep=len(hep), n_chol=len(chol), evaluable=False))
        continue
    central = hep[gl[hep] >= np.quantile(gl[hep], 0.95)]
    xy = obs[["x", "y"]].to_numpy(float)
    dp = cKDTree(xy[chol]).query(xy[hep], k=1)[0]
    dc = cKDTree(xy[central]).query(xy[hep], k=1)[0]
    L = dp / (dp + dc + 1e-12)
    is_core = core[hep]
    for g, role in gene_list + [("__negprobe__", "negprobe")]:
        v = obs.neg_norm.to_numpy()[hep] if g == "__negprobe__" else Xn[hep, gi[g]].toarray().ravel()
        if np.isnan(v).all():
            continue
        for subset, m in [("all", np.ones(len(hep), bool)), ("core", is_core)]:
            r = stats.spearmanr(L[m], v[m])[0]
            q = pd.qcut(L[m], 5, labels=False, duplicates="drop")
            means = pd.Series(v[m]).groupby(q).mean()
            rows.append(dict(unit=unit, n_hep=len(hep), n_chol=len(chol), n_core=int(is_core.sum()), evaluable=True,
                             gene=g, role=role, subset=subset, rho_L=r, mean=v[m].mean(), detect=(v[m] > 0).mean(),
                             log2_portal_over_central=np.log2((means.iloc[0] + 0.01) / (means.iloc[-1] + 0.01))))
            prof += [dict(unit=unit, gene=g, role=role, subset=subset, L_quintile=int(k), mean=val, rel=val / (v[m].mean() + 1e-9))
                     for k, val in means.items()]
        # neighbourhood of member-high hepatocytes: top 10% of the gene within unit, all hepatocytes
        if role in ("member", "cholangiocyte_spill", "hepatocyte_ref"):
            top = v >= np.quantile(v, 0.90) if (v > 0).mean() > 0.10 else v > 0
            for l in wide:
                w = wide[l][hep]
                nb.append(dict(unit=unit, gene=g, role=role, neighbour_lineage=l, n_top=int(top.sum()),
                               frac_top=np.nanmean(w[top]), frac_rest=np.nanmean(w[~top]),
                               L_top=np.median(L[top]), L_rest=np.median(L[~top])))
pd.DataFrame(rows).to_csv(OUT / "lobular_gradient_by_unit.tsv", sep="\t", index=False)
pd.DataFrame(prof).to_csv(OUT / "lobular_profile_quintiles.tsv", sep="\t", index=False)
pd.DataFrame(nb).to_csv(OUT / "neighbourhood_of_high_hepatocytes.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).query("evaluable == True and subset == 'core'")[["unit", "gene", "role", "rho_L", "log2_portal_over_central"]].to_string(index=False), flush=True)
