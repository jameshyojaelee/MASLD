"""DIG1 job 14b (exploratory). Power check for job 14: detection and mean of the ligands and
outcome genes by lineage and unit on each CosMx panel, with the same loading, QC and labels
(lines 12-138 of 03_insitu_position.py). Also saves per-cell IDs with labels for reuse.
Usage: python -I 14b_insitu_detection.py {govaere|fukumoto|gse292268}
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
OUT = RUN / "out_insitu_detection" / dataset
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

obs[["cell", "unit", "fov", "lineage", "total"]].to_csv(OUT / "cell_labels_with_id.tsv.gz", sep="\t", index=False)
G = ["FGF7", "HGF", "CXCL12", "TNFSF12", "SOX9", "EPCAM", "KRT8", "CDKN1A", "AKR1B10", "FGFR2", "MET", "ACKR3", "CXCR4",
     "TNFRSF12A", "KRT7", "KRT19", "COL1A1", "LUM", "DCN", "APOA1", "SERPINA1", "GLUL"]
rows = []
for (unit, l), idx in obs.groupby(["unit", "lineage"]).indices.items():
    for g in G:
        if g not in gi:
            continue
        v = X[idx, gi[g]].toarray().ravel()
        rows.append(dict(unit=unit, lineage=l, n_cells=len(idx), gene=g, frac_detected=float((v > 0).mean()), mean_counts=float(v.mean())))
    rows.append(dict(unit=unit, lineage=l, n_cells=len(idx), gene="__negprobe__", frac_detected=float((obs.neg_counts.to_numpy()[idx] > 0).mean()),
                     mean_counts=float(np.nanmean(obs.neg_counts.to_numpy()[idx]))))
d = pd.DataFrame(rows)
d.to_csv(OUT / "detection_by_unit_lineage.tsv", sep="\t", index=False)
print(d[d.lineage.isin(["stellate_fibroblast", "hepatocyte", "cholangiocyte", "macrophage"])].pivot_table(
    index="gene", columns=["lineage"], values="frac_detected", aggfunc="median").round(3).to_string(), flush=True)
