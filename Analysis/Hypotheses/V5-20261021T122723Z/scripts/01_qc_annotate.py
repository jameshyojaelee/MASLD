"""V5 job 1: G0, cell QC, CellTypist and marker labels, G1 (PC1), G2, G3, PC2,
D1-D3, D5-D9. Computes no program score. Follows PRESPEC_V5.md.

Usage: python -I 01_qc_annotate.py
"""
import hashlib
import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
from scipy import sparse, stats

sys.path.insert(0, str(Path(__file__).parent))
import v5_common as C

out = C.OUT
if any(out.iterdir()):
    sys.exit(f"refusing to write into non-empty {out}")
C.DERIVED.mkdir(exist_ok=True)
QC_H5AD = C.DERIVED / "gse298719_qc_cells.h5ad"
if QC_H5AD.exists():
    sys.exit(f"refusing to overwrite {QC_H5AD}")
print(f"prespec sha256 {C.check_prespec()}", flush=True)
if C.sha256(C.CT_MODEL) != C.CT_MODEL_SHA:
    sys.exit("STOP: CellTypist model hash differs")

# ---- G0 identity -----------------------------------------------------------------
design = C.load_design()
g0 = []
feat_ref = None
for r in design.itertuples():
    pre = C.MAT / f"{r.gsm}_{r.title}_"
    files = [Path(f"{pre}{s}") for s in ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")]
    present = all(f.exists() for f in files)
    content = hashlib.sha256("\n".join(C.read_gz_lines(files[2])).encode()).hexdigest() if present else ""
    feat_ref = feat_ref or content
    g0.append(dict(gsm=r.gsm, title=r.title, files_present=present, features_content_sha256=content,
                   features_identical=content == feat_ref))
g0 = pd.DataFrame(g0)
feat = pd.read_csv(f"{C.MAT}/{design.gsm[0]}_{design.title[0]}_features.tsv.gz", sep="\t", header=None,
                   names=["id", "symbol", "type"])
g0_pass = bool(g0.files_present.all() and g0.features_identical.all() and (feat.type == "Gene Expression").all()
               and len(feat) == 18082)
g0.to_csv(out / "gate_G0_identity.tsv", sep="\t", index=False)
if not g0_pass:
    (out / "STOP_G0.txt").write_text("untestable: deposited matrices do not match the metadata\n")
    sys.exit("STOP G0")
print(f"G0 pass: 13 samples, {len(feat)} features, identical feature files", flush=True)

ens = feat.id.str.split(".").str[0].to_numpy()
uniq = pd.Index(pd.unique(feat.symbol.astype(str)))
col = uniq.get_indexer(feat.symbol.astype(str))
S = sparse.csr_matrix((np.ones(len(feat), np.float32), (np.arange(len(feat)), col)), shape=(len(feat), len(uniq)))
print(f"unique symbols {len(uniq)} (duplicates summed)", flush=True)

mats, mats_ens, obs = [], [], []
for r in design.itertuples():
    pre = C.MAT / f"{r.gsm}_{r.title}_"
    m = scipy.io.mmread(f"{pre}matrix.mtx.gz").tocsr().T.tocsr().astype(np.float32)
    if m.shape[1] != len(feat):
        sys.exit(f"STOP G0: {r.gsm} matrix has {m.shape[1]} features")
    bc = C.read_gz_lines(Path(f"{pre}barcodes.tsv.gz"))
    if len(bc) != m.shape[0]:
        sys.exit(f"STOP G0: {r.gsm} barcode count differs from matrix rows")
    mats_ens.append(m)
    mats.append((m @ S).tocsr())
    obs.append(pd.DataFrame({"gsm": r.gsm, "barcode": bc}))
X = sparse.vstack(mats).tocsr()
XE = sparse.vstack(mats_ens).tocsr()
del mats, mats_ens
obs = pd.concat(obs, ignore_index=True)
raw_cells = obs.groupby("gsm").size()
print(f"raw cells {X.shape[0]:,} x {X.shape[1]:,} symbols", flush=True)

# ---- QC: gene filter, then cell filter ---------------------------------------------
gkeep = np.asarray(X.sum(axis=0)).ravel() >= 3
X = X[:, gkeep].tocsr()
genes = uniq[gkeep]
gi = {g: i for i, g in enumerate(genes)}
ekeep = np.asarray(XE.sum(axis=0)).ravel() >= 3
XE = XE[:, ekeep].tocsr()
ens_kept = ens[ekeep]
umi = np.asarray(X.sum(axis=1)).ravel()
obs["umi"] = umi
obs["n_genes"] = np.diff(X.indptr)
mt_idx = [gi[g] for g in genes if g.startswith("MT-")]
obs["pct_mt"] = 100 * np.asarray(X[:, mt_idx].sum(axis=1)).ravel() / np.maximum(umi, 1)
obs["log10_umi"] = np.log10(np.maximum(umi, 1))
med = obs.groupby("gsm").log10_umi.transform("median")
mad = obs.groupby("gsm").log10_umi.transform(lambda v: stats.median_abs_deviation(v, scale="normal"))
obs["umi_upper"] = 10 ** (med + 5 * mad)
qc = (obs.umi >= 200) & (obs.log10_umi <= med + 5 * mad) & (obs.pct_mt < 25)
obs["qc_h90_rule"] = (obs.umi >= 200) & (obs.umi <= 10000)
h90_counts = obs.groupby("gsm").qc_h90_rule.sum()
X = X[qc.to_numpy()].tocsr()
XE = XE[qc.to_numpy()].tocsr()
obs = obs[qc].reset_index(drop=True).merge(design, on="gsm", how="left")
print(f"genes kept {gkeep.sum():,}; cells passing QC {len(obs):,}", flush=True)

# ---- G1 / PC1 sex from Y probes --------------------------------------------------
codes = pd.Categorical(obs.gsm, categories=design.gsm).codes
ind = sparse.csr_matrix((np.ones(len(obs)), (codes, np.arange(len(obs)))), shape=(len(design), len(obs)))
tot = np.asarray(ind @ X.sum(axis=1)).ravel()
yp = [g for g in C.Y_PROBES if g in gi]
ycpm = np.asarray((ind @ X[:, [gi[g] for g in yp]]).sum(axis=1)).ravel() / tot * 1e6
pc1 = design[["gsm", "id", "sex"]].copy()
pc1["y_cpm_sum"] = ycpm
pc1["y_probes"] = ",".join(yp)
lv = np.log10(ycpm + 1)
order = np.sort(lv)
gap = int(np.argmax(np.diff(order)))
cut = (order[gap] + order[gap + 1]) / 2
pc1["log10_y"] = lv
pc1["sex_call"] = np.where(lv > cut, "male", "female")
pc1["concordant"] = pc1.sex_call == pc1.sex
pc1.to_csv(out / "D5_PC1_sex.tsv", sep="\t", index=False)
n_disc = int((~pc1.concordant).sum())
g1 = n_disc < 2
print(f"PC1: {13 - n_disc} of 13 concordant (cut log10 {cut:.3f})", flush=True)

# ---- CellTypist ----------------------------------------------------------------
import scanpy as sc
import celltypist
A = ad.AnnData(X=X.copy(), obs=obs.set_index(obs.gsm + "_" + obs.barcode), var=pd.DataFrame(index=genes))
sc.pp.normalize_total(A, target_sum=1e4)
sc.pp.log1p(A)
model = celltypist.models.Model.load(model=str(C.CT_MODEL))
np.random.seed(C.SEED)
pred = celltypist.annotate(A, model=model, majority_voting=True)
obs["label"] = pred.predicted_labels["majority_voting"].astype(str).values
obs["label_raw"] = pred.predicted_labels["predicted_labels"].astype(str).values
obs["conf"] = pred.probability_matrix.max(axis=1).values
print("CellTypist:", obs.label.value_counts().to_dict(), flush=True)

# ---- marker labels (sensitivity; program genes removed from markers) ---------------
spec = importlib.util.spec_from_file_location("v3_annotate", C.V3_SCRIPT)
v3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v3)
mem = C.membership()
prog_genes = set()
for uid in (C.PRIMARY, C.SECONDARY, C.NC1):
    prog_genes |= set(C.weights(mem, uid).mapped_symbol)
markers = {k: [g for g in v if g in gi and g not in prog_genes] for k, v in v3.MARKERS.items()}
pd.Series({k: ";".join(v) for k, v in markers.items()}).to_csv(out / "markers_used.tsv", sep="\t", header=False)
A.layers["counts"] = X.copy()
A.obs["donor"] = obs.gsm.values
sc.pp.highly_variable_genes(A, n_top_genes=2000, flavor="seurat_v3", batch_key="donor", layer="counts")
sc.tl.pca(A, n_comps=30, mask_var="highly_variable", random_state=v3.SEED)
sc.pp.neighbors(A, n_neighbors=15, random_state=v3.SEED)
sc.tl.leiden(A, resolution=1.0, random_state=v3.SEED, flavor="igraph", n_iterations=2, directed=False)
for k, gl in markers.items():
    sc.tl.score_genes(A, gl, score_name=f"mk_{k}", random_state=v3.SEED)
cols = [f"mk_{k}" for k in markers]
cl = A.obs.groupby("leiden", observed=True)[cols].mean()
cl["lineage"] = cl[cols].idxmax(axis=1).str.replace("mk_", "", regex=False)
cl.to_csv(out / "marker_cluster_assignment.tsv", sep="\t")
obs["marker_lineage"] = A.obs["leiden"].map(cl["lineage"]).astype(str).values
del A
agree = pd.crosstab(obs.label, obs.marker_lineage)
agree.to_csv(out / "D2_celltypist_vs_marker.tsv", sep="\t")

# ---- lineage counts, G2 ---------------------------------------------------------
lin = pd.crosstab(obs.gsm, obs.label).reindex(design.gsm).fillna(0).astype(int)
lin.to_csv(out / "lineage_counts_per_donor.tsv", sep="\t")
dq = design.copy()
dq["raw_cells"] = raw_cells.reindex(dq.gsm).values
dq["cells_qc"] = obs.groupby("gsm").size().reindex(dq.gsm).fillna(0).astype(int).values
dq["cells_h90_umi_rule"] = h90_counts.reindex(dq.gsm).values
dq["umi_upper_bound"] = obs.groupby("gsm").umi_upper.first().reindex(dq.gsm).values
dq["n_hep"] = lin.get("Hepatocytes", pd.Series(0, index=lin.index)).reindex(dq.gsm).values
dq["n_hep_marker"] = obs[obs.marker_lineage == "hepatocyte"].groupby("gsm").size().reindex(dq.gsm).fillna(0).astype(int).values
hep = obs[obs.label == "Hepatocytes"]
for c, f in (("umi", "median_umi_hep"), ("n_genes", "median_genes_hep"), ("pct_mt", "median_pct_mt_hep")):
    dq[f] = hep.groupby("gsm")[c].median().reindex(dq.gsm).values
dq["pc1_concordant"] = pc1.concordant.values
dq["excluded_reason"] = ""
if g1:
    dq.loc[~dq.pc1_concordant, "excluded_reason"] = "PC1 sex discordant"
dq.loc[(dq.n_hep < C.MIN_CELLS) & (dq.excluded_reason == ""), "excluded_reason"] = "fewer than 30 hepatocytes"
dq["included"] = dq.excluded_reason == ""
inc = dq[dq.included]
n_h, n_m = int((inc.disease == 0).sum()), int((inc.disease == 1).sum())
g2 = n_h >= 5 and n_m >= 5
dq.to_csv(out / "donor_qc.tsv", sep="\t", index=False)
dq.groupby("group").agg(donors=("gsm", "size"), raw_cells=("raw_cells", "median"), cells_qc=("cells_qc", "median"),
                        n_hep_median=("n_hep", "median"), n_hep_min=("n_hep", "min"),
                        median_umi_hep=("median_umi_hep", "median"), included=("included", "sum")) \
    .to_csv(out / "D1_qc_by_group.tsv", sep="\t")
frac = pd.crosstab(obs.gsm, obs.label, normalize="index").reindex(design.gsm)
frac.to_csv(out / "D2_label_fractions_per_donor.tsv", sep="\t")
frac.join(design.set_index("gsm").group).groupby("group").mean().to_csv(out / "D2_label_fractions_by_group.tsv", sep="\t")
hep.groupby("group").conf.median().to_csv(out / "D2_hep_conf_by_group.tsv", sep="\t")

# ---- G3 / D3 coverage on included donors ----------------------------------------
inc_mask = obs.gsm.isin(inc.gsm).to_numpy()
det_by_lin = {}
for ct, lab in C.LINEAGE_LABEL.items():
    m = inc_mask & (obs.label == lab).to_numpy()
    det_by_lin[ct] = np.asarray((X[m] > 0).mean(axis=0)).ravel() if m.sum() else np.zeros(len(genes))
cov_rows, gene_rows = [], []
for uid, grp in mem.groupby("program_uid"):
    ct = grp.cell_type.iloc[0]
    w = C.weights(mem, uid)
    pres = w.mapped_symbol.isin(gi)
    det = w.mapped_symbol.map(lambda g: det_by_lin[ct][gi[g]] if g in gi else np.nan)
    scored = pres & (det >= C.DET_MIN)
    cov_rows.append(dict(program_uid=uid, cell_type=ct, module=int(grp.module.iloc[0]), n_confirmed=len(w),
                         n_present=int(pres.sum()), present_weight=float(w.original_l1_weight[pres].sum()),
                         n_scored=int(scored.sum()), retained_weight=float(w.original_l1_weight[scored].sum()),
                         lineage_cells=int((inc_mask & (obs.label == C.LINEAGE_LABEL[ct]).to_numpy()).sum())))
    if uid in (C.PRIMARY, C.SECONDARY, C.NC1):
        for g, wt, p, dd, s in zip(w.mapped_symbol, w.original_l1_weight, pres, det, scored):
            gene_rows.append(dict(program_uid=uid, gene=g, l1_weight=wt, on_axis=p, detection=dd, scored=s))
cov = pd.DataFrame(cov_rows)
cov["passes_g3_rule"] = (cov.n_scored >= C.MIN_GENES) & (cov.retained_weight >= C.MIN_WEIGHT)
cov.to_csv(out / "D3_coverage.tsv", sep="\t", index=False)
pd.DataFrame(gene_rows).to_csv(out / "D3_confirmatory_gene_detection.tsv", sep="\t", index=False)
row = cov.set_index("program_uid")
g3p, g3s, g3n = (bool(row.loc[u, "passes_g3_rule"]) for u in (C.PRIMARY, C.SECONDARY, C.NC1))

# ---- PC2 collagen in fibroblasts ---------------------------------------------------
fm = (obs.label == "Fibroblasts").to_numpy()
fcodes = pd.Categorical(obs.gsm[fm], categories=design.gsm).codes
find = sparse.csr_matrix((np.ones(fm.sum()), (fcodes, np.arange(fm.sum()))), shape=(len(design), fm.sum()))
Xf = X[fm]
col_c = np.asarray((find @ Xf[:, [gi["COL1A1"], gi["COL1A2"]]]).sum(axis=1)).ravel()
tot_f = np.asarray((find @ Xf).sum(axis=1)).ravel()
pc2 = design[["gsm", "id", "group", "disease", "sex_f"]].copy()
pc2["n_fib"] = np.asarray(find.sum(axis=1)).ravel().astype(int)
pc2["log2_cpm_col1"] = np.log2(np.divide(col_c, tot_f, out=np.zeros_like(col_c), where=tot_f > 0) * 1e6 + 1)
pc2.to_csv(out / "D6_PC2_collagen.tsv", sep="\t", index=False)
sub = pc2[(pc2.n_fib >= 20) & pc2.gsm.isin(inc.gsm)]
res = dict(n_healthy=int((sub.disease == 0).sum()), n_masld=int((sub.disease == 1).sum()))
if res["n_healthy"] >= 3 and res["n_masld"] >= 3:
    Xd = np.column_stack([np.ones(len(sub)), sub.disease, sub.sex_f]).astype(float)
    b, se, sig, df, _ = C.hc3(sub.log2_cpm_col1.to_numpy(float), Xd)
    res.update(beta=b[1], hc3_se=se[1], hc3_p=2 * stats.t.sf(abs(b[1] / se[1]), df), std_diff=b[1] / sig, df=df,
               readable=True, prediction_met=bool(b[1] > 0 and b[1] / sig > 0.5))
else:
    res.update(readable=False, prediction_met=False)
pd.DataFrame([res]).to_csv(out / "PC2_collagen_result.tsv", sep="\t", index=False)
print("PC2:", res, flush=True)

# ---- D7 ambient proxies in hepatocytes -------------------------------------------
d7_genes = [g for g in ["KRT7", "KRT19", "EPCAM", "PTPRC", "PECAM1"] if g in gi and g not in prog_genes]
hm = (obs.label == "Hepatocytes").to_numpy()
hob = obs[hm].reset_index(drop=True)
norm = C.log_normalize(X[hm][:, [gi[g] for g in d7_genes]].tocsr(), hob.umi.to_numpy(float))
hcodes = pd.Categorical(hob.gsm, categories=design.gsm).codes
hind = sparse.csr_matrix((np.ones(len(hob)), (hcodes, np.arange(len(hob)))), shape=(len(design), len(hob)))
nh = np.maximum(np.asarray(hind.sum(axis=1)).ravel(), 1)
d7 = pd.DataFrame(np.asarray((hind @ norm).todense()) / nh[:, None], columns=d7_genes)
d7.insert(0, "gsm", design.gsm.values)
d7.insert(1, "group", design.group.values)
d7.to_csv(out / "D7_ambient_proxies_per_donor.tsv", sep="\t", index=False)

# ---- D8 leverage, D9 balance and VIF (design only) ---------------------------------
inc = inc.copy()
inc["log10_n_hep"] = np.log10(inc.n_hep)
Xd = np.column_stack([np.ones(len(inc)), inc.disease, inc.sex_f, inc.log10_n_hep]).astype(float)
inc["leverage"] = np.einsum("ij,jk,ik->i", Xd, np.linalg.inv(Xd.T @ Xd), Xd)
inc[["gsm", "id", "group", "sex", "plate", "n_hep", "log10_n_hep", "leverage"]].to_csv(out / "D8_leverage.tsv", sep="\t", index=False)
Z = Xd[:, [0, 2, 3]]
yv = inc.disease.to_numpy(float)
r = yv - Z @ np.linalg.lstsq(Z, yv, rcond=None)[0]
vif = float(1 / (r @ r / ((yv - yv.mean()) ** 2).sum()))
inc.groupby("group").agg(n=("gsm", "size"), female=("sex_f", "sum"), age=("age", "mean"), bmi=("bmi", "mean"),
                         log10_n_hep=("log10_n_hep", "mean"), median_umi_hep=("median_umi_hep", "median")) \
    .to_csv(out / "D9_balance.tsv", sep="\t")
pd.crosstab(inc.plate, inc.group).to_csv(out / "D9_plate_by_group.tsv", sep="\t")
primary_model = "score ~ disease + sex + log10(n_hep)" if vif <= 5 else "score ~ disease"

summ = dict(G0_pass=g0_pass, G1_pass=g1, n_discordant_sex=n_disc, G2_pass=g2, n_healthy=n_h, n_masld=n_m,
            G3_primary_pass=g3p, primary_scored=int(row.loc[C.PRIMARY, "n_scored"]),
            primary_weight=float(row.loc[C.PRIMARY, "retained_weight"]),
            G3_secondary_pass=g3s, secondary_scored=int(row.loc[C.SECONDARY, "n_scored"]),
            secondary_weight=float(row.loc[C.SECONDARY, "retained_weight"]),
            NC1_coverage_ok=g3n, NC1_scored=int(row.loc[C.NC1, "n_scored"]), NC1_weight=float(row.loc[C.NC1, "retained_weight"]),
            VIF_disease=vif, primary_model=primary_model, n_included=len(inc),
            hep_cells_included=int(inc.n_hep.sum()), programs_passing_g3_rule=int(cov.passes_g3_rule.sum()))
pd.Series(summ).to_csv(out / "gate_G0_G3.tsv", sep="\t", header=["value"])
print(pd.Series(summ).to_string(), flush=True)

# ---- save QC cells (symbol counts) and Ensembl counts for jobs 2-3 ------------------
A = ad.AnnData(X=X, obs=obs.astype({"gsm": str}).set_index(obs.gsm + "_" + obs.barcode), var=pd.DataFrame(index=genes))
A.write_h5ad(QC_H5AD)
AE = ad.AnnData(X=XE, obs=pd.DataFrame(index=A.obs_names), var=pd.DataFrame(index=pd.Index(ens_kept, name="gene_id")))
AE.write_h5ad(C.DERIVED / "gse298719_qc_cells_ensembl.h5ad")
print(f"wrote {QC_H5AD}", flush=True)
for k, v, msg in (("G1", g1, "untestable: sample labels cannot be linked to matrices"),
                  ("G2", g2, "untestable: too few donors"), ("G3", g3p, "untestable: program not measured")):
    if not v:
        (out / f"STOP_{k}.txt").write_text(msg + "\n")
        sys.exit(f"STOP {k}: {msg}")
print("G0-G3 pass")
