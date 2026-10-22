"""V7 phase 1: CosMx observability rows and in situ lineage of program member genes.
Follows PRESPEC_V7_phase1.md. One dataset per run.

Usage: python -I 01_phase1.py {govaere|fukumoto|gse292268}
Writes outputs/<dataset>/ (refuses to overwrite).
"""
import gzip
import hashlib
import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/V7-20261021T235717Z"
DL = RUN / "downloads_untrusted"
GOVAERE_DIR = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external/govaere2026_natgenetics/cosmx")
HOT = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
MEMBERSHIP, REGISTRY = HOT / "program_membership_v2.tsv", HOT / "program_registry_v2.tsv"
V3_SCRIPT = ROOT / "scripts/analysis/hypotheses_2026q4/V3_program_holdout/annotate_pseudobulk.py"
STROMAL_ECM = "hotspot_hepatocytes_f05c535ae5bbc0b9"
DUCTULAR = "hotspot_hepatocytes_48f39dd4d817a10e"
CONTROLS = {"stellate_fibroblast": ["COL1A1", "COL1A2", "COL3A1", "LUM", "PDGFRB", "RGS5"],
            "cholangiocyte": ["KRT7", "KRT19", "EPCAM", "SOX9"]}
PROGRAM_REF = {STROMAL_ECM: "stellate_fibroblast", DUCTULAR: "cholangiocyte"}
LINEAGE_OF = {"hepatocytes": "hepatocyte", "fibroblasts": "stellate_fibroblast", "cholangiocytes": "cholangiocyte",
              "macrophages": "macrophage", "tcells": "t_nk"}
CONTROL_PREFIXES = ("Negative", "NegPrb", "SystemControl", "FalseCode", "NegativeAdd")
SEED, EPS, DET_MIN, MIN_GENES, MIN_WEIGHT, MIN_LINEAGE = 0, 0.01, 0.01, 8, 0.80, 100
MIN_HEP, MIN_REF, K_NN, EDGE_MULT = 200, 50, 6, 2.5


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


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
            fovpos.assign(unit=fovpos.FOV.astype(str).map(fov_unit)).to_csv(OUT / f"fov_split_{name}.tsv", sep="\t", index=False)
            obs["unit"] = obs.fov.map(fov_unit)
            parts.append((tg, neg, obs, n_negp))
    elif dataset == "fukumoto":
        d = DL / "fukumoto"
        for name, ef, mf in [("sample1", "CosMx_expression_data_MASH_sample1_NT.csv", "CosMx_metadata_MASH_sample1_NT.csv"),
                             ("sample2", "CosMx_expression_data_MASH_sample_2_NT.csv", "CosMx_metadata_MASH_sample_2_NT.csv")]:
            e = pd.read_csv(d / ef, index_col=0)
            m = pd.read_csv(d / mf, index_col=0, usecols=lambda c: c in ("Unnamed: 0", "", "fov", "x_slide_mm", "y_slide_mm"))
            m = m.reindex(e.index)
            if m.x_slide_mm.isna().any():
                sys.exit(f"STOP: {name} expression cells missing from metadata")
            obs = pd.DataFrame({"fov": m.fov.astype(str).values, "x": m.x_slide_mm.values, "y": m.y_slide_mm.values,
                                "cell": [f"{name}_{c}" for c in e.index], "unit": name})
            tg, neg, n_negp = split_controls(e)
            parts.append((tg, neg, obs, n_negp))
    else:
        sys.exit("unknown dataset")
    genes = sorted(set.intersection(*[set(p[0].columns) for p in parts]))
    X = sparse.vstack([sparse.csr_matrix(p[0][genes].to_numpy(np.float32)) for p in parts]).tocsr()
    obs = pd.concat([p[2] for p in parts], ignore_index=True)
    obs["neg_counts"] = np.concatenate([p[1] for p in parts])
    obs["fov"] = obs.unit + "_" + obs.fov
    return X, obs, pd.Index(genes), [p[3] for p in parts]


dataset = sys.argv[1]
OUT = RUN / "outputs" / dataset
assert sha256(RUN / "PRESPEC_V7_phase1.md") == (RUN / "PRESPEC_V7_phase1.md.sha256").read_text().split()[0], "prespec hash"
assert sha256(MEMBERSHIP).startswith("feec3fc9"), "membership hash"
OUT.mkdir(parents=True, exist_ok=False)

X, obs, genes, n_negp = load(dataset)
obs["total"] = np.asarray(X.sum(axis=1)).ravel()
obs["n_genes"] = np.asarray((X > 0).sum(axis=1)).ravel()
keep = ((obs.total >= 20) & (obs.n_genes >= 10)).to_numpy()
X, obs = X[keep], obs[keep].reset_index(drop=True)
print(f"{dataset}: {X.shape[0]} cells after QC, {len(genes)} targets, negative probes per part {n_negp}", flush=True)

mem = pd.read_csv(MEMBERSHIP, sep="\t")
mem = mem[mem.mapped_symbol_status == "gencode_v49_unique_symbol_confirmed"]
members_conf = set(mem[mem.program_uid.isin([STROMAL_ECM, DUCTULAR])].mapped_symbol)

# ---- labels ------------------------------------------------------------------------------
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
# argmax of marker scores standardized across clusters: raw means favour abundant hepatocyte
# transcripts, which spill into every CosMx cell (see DEVIATIONS.md)
zc = (cl[cols] - cl[cols].mean()) / cl[cols].std(ddof=1)
cl["lineage"] = zc.idxmax(axis=1).str.replace("mk_", "", regex=False).where(zc.max(axis=1) >= 0.5, "unassigned")
cl["z_max"] = zc.max(axis=1)
cl["n_cells"] = A.obs.leiden.value_counts().reindex(cl.index).values
cl.to_csv(OUT / "cluster_lineage_assignment.tsv", sep="\t")
obs["lineage"] = A.obs.leiden.map(cl.lineage).astype(str).values
del A
pd.Series({k: ";".join(v) for k, v in markers.items()}).to_csv(OUT / "markers_used.tsv", sep="\t", header=False)
pd.crosstab(obs.unit, obs.lineage).to_csv(OUT / "cells_per_unit_lineage.tsv", sep="\t")

# ---- estimand 1: observability rows -------------------------------------------------------
reg = pd.read_csv(REGISTRY, sep="\t", usecols=["cell_type", "module", "program_uid", "module_name"])
det_by_lin = {}
for lin in set(LINEAGE_OF.values()):
    m = (obs.lineage == lin).to_numpy()
    det_by_lin[lin] = (pd.Series(np.asarray((X[m] > 0).mean(axis=0)).ravel(), index=genes) if m.sum() else pd.Series(0.0, index=genes), int(m.sum()))
rows = []
for _, r in reg.iterrows():
    lin = LINEAGE_OF[r.cell_type]
    det, n_lin = det_by_lin[lin]
    w = mem[mem.program_uid == r.program_uid][["mapped_symbol", "original_l1_weight"]]
    on = w[w.mapped_symbol.isin(genes)]
    dd = on[det.reindex(on.mapped_symbol).to_numpy() >= DET_MIN]
    state = "lineage not captured" if n_lin < MIN_LINEAGE else (
        "testable" if len(dd) >= MIN_GENES and dd.original_l1_weight.sum() >= MIN_WEIGHT else "not measured")
    rows.append(dict(dataset=dataset, program_uid=r.program_uid, cell_type=r.cell_type, module=r.module,
                     module_name=r.module_name, lineage_cells=n_lin, n_members=len(w), panel_genes=len(on),
                     panel_weight=on.original_l1_weight.sum(), detected_genes=len(dd),
                     detected_weight=dd.original_l1_weight.sum(), state=state))
pd.DataFrame(rows).to_csv(OUT / "observability_117_programs.tsv", sep="\t", index=False)

# ---- estimand 2: in situ lineage ----------------------------------------------------------
scale = np.median(obs.total) / obs.total.to_numpy()
Xn = (sparse.diags(scale) @ X).tocsr()
gi = {g: i for i, g in enumerate(genes)}
obs["neg_norm"] = obs.neg_counts.to_numpy() * scale

# neighbour graph within FOV: k nearest, edges at most EDGE_MULT x FOV median nearest-neighbour distance
nbr_ref = {lin: np.zeros(len(obs), bool) for lin in CONTROLS}
nbr_all_hep = np.zeros(len(obs), bool)
has_edge = np.zeros(len(obs), bool)
lin_arr = obs.lineage.to_numpy()
for fov, idx in obs.groupby("fov").indices.items():
    if len(idx) <= K_NN:
        continue
    xy = obs.loc[idx, ["x", "y"]].to_numpy(float)
    d, j = cKDTree(xy).query(xy, k=K_NN + 1)
    d, j = d[:, 1:], j[:, 1:]
    ok = d <= EDGE_MULT * np.median(d[:, 0])
    nl = lin_arr[idx][j]
    has_edge[idx] = ok.any(axis=1)
    for lin in CONTROLS:
        nbr_ref[lin][idx] = ((nl == lin) & ok).any(axis=1)
    nbr_all_hep[idx] = ((nl == "hepatocyte") | ~ok).all(axis=1) & ok.any(axis=1)

res = []
for uid, ref in PROGRAM_REF.items():
    mg = sorted(g for g in mem[mem.program_uid == uid].mapped_symbol if g in gi)
    cg = [g for g in CONTROLS[ref] if g in gi]
    for unit, uidx in obs.groupby("unit").indices.items():
        lin_u = lin_arr[uidx]
        hep = uidx[lin_u == "hepatocyte"]
        refc = uidx[lin_u == ref]
        core = hep[nbr_all_hep[hep]]
        inter = hep[nbr_ref[ref][hep] & has_edge[hep]]
        evaluable = len(hep) >= MIN_HEP and len(refc) >= MIN_REF
        for g in mg + cg + ["__negprobe__"]:
            if g == "__negprobe__":
                v = obs.neg_norm.to_numpy()
                if np.isnan(v).all():
                    continue
                get = lambda ii: float(np.nanmean(v[ii])) if len(ii) else np.nan
            else:
                col = Xn[:, gi[g]].toarray().ravel()
                get = lambda ii, col=col: float(col[ii].mean()) if len(ii) else np.nan
            mh, mr, mc, mi = get(hep), get(refc), get(core), get(inter)
            res.append(dict(dataset=dataset, unit=unit, program_uid=uid, reference=ref,
                            gene=g, role="member" if g in mg else ("negprobe" if g == "__negprobe__" else "control"),
                            n_hep=len(hep), n_ref=len(refc), n_core=len(core), n_interface=len(inter), evaluable=evaluable,
                            mean_hep=mh, mean_ref=mr, mean_core=mc, mean_interface=mi,
                            R=np.log2((mh + EPS) / (mr + EPS)), C=np.log2((mc + EPS) / (mi + EPS))))
res = pd.DataFrame(res)
ctrl = res[res.role == "control"].groupby(["unit", "program_uid"]).agg(R_ctrl=("R", "median"), C_ctrl=("C", "median")).reset_index()
res = res.merge(ctrl, on=["unit", "program_uid"], how="left")
res["E"] = res.R - res.R_ctrl
res["core_excess"] = res.C - res.C_ctrl
res.to_csv(OUT / "insitu_lineage_by_unit.tsv", sep="\t", index=False)

ev = res[(res.role == "member") & res.evaluable]
summ = ev.groupby(["program_uid", "gene"]).agg(
    n_units=("unit", "size"), median_E=("E", "median"), median_core_excess=("core_excess", "median"),
    n_beyond=("E", lambda s: int(((s > 1) & (ev.loc[s.index, "core_excess"] > 0)).sum())),
    n_spill=("E", lambda s: int(((s <= 0) & (ev.loc[s.index, "core_excess"] <= 0)).sum()))).reset_index()
summ["reading"] = np.where(summ.n_beyond > summ.n_units / 2, "hepatocyte expression beyond spillover",
                           np.where(summ.n_spill > summ.n_units / 2, "consistent with spillover", "mixed"))
summ["source_owned"] = (dataset == "gse292268") & (summ.gene == "IGFBP7")
summ.to_csv(OUT / "insitu_lineage_summary.tsv", sep="\t", index=False)
print(summ.to_string(index=False), flush=True)
print(pd.DataFrame(rows).state.value_counts().to_string(), flush=True)
