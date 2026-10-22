"""H90 job 6301: nucleus QC, strip, CellTypist, lineage counts, coverage, PC1, PC2,
D1-D3, D5-D9, gates G1-G3. Computes no program score."""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
from scipy import sparse, stats

sys.path.insert(0, str(Path(__file__).parent))
import h90_common as C

C.check_inputs()
out = C.OUT
C.DERIVED.mkdir(parents=True, exist_ok=True)
QC_H5AD = C.DERIVED / "gse289173_qc_nuclei.h5ad"

design = pd.read_csv(C.WRITE / "donor_design.tsv", sep="\t").sort_values("gsm").reset_index(drop=True)
design["sex_f"] = (design.sex == "female").astype(int)
feat = pd.read_csv(C.FEATURES_FL514, sep="\t", header=None, names=["id", "symbol", "type"])
if not (feat.type == "Gene Expression").all():
    sys.exit("unexpected non-gene features")
uniq = pd.Index(pd.unique(feat.symbol.astype(str)))
col = uniq.get_indexer(feat.symbol.astype(str))
S = sparse.csr_matrix((np.ones(len(feat), np.float32), (np.arange(len(feat)), col)), shape=(len(feat), len(uniq)))
print(f"features {len(feat)}, unique symbols {len(uniq)} (duplicates summed by symbol)", flush=True)

mats, obs = [], []
for r in design.itertuples():
    pre = f"{r.gsm}_{r.title}_"
    m = scipy.io.mmread(C.RAW / f"{pre}matrix.mtx.gz").tocsr().T.tocsr().astype(np.float32)
    bc = C.read_gz_lines(C.RAW / f"{pre}barcodes.tsv.gz")
    mats.append((m @ S).tocsr())
    obs.append(pd.DataFrame({"gsm": r.gsm, "barcode": bc}))
X = sparse.vstack(mats).tocsr()
del mats
obs = pd.concat(obs, ignore_index=True)
print(f"raw nuclei {X.shape[0]:,} x {X.shape[1]:,} symbols, nnz {X.nnz:,}", flush=True)
filt_bc = obs.groupby("gsm").size()

# A8 order: gene filter (>= 3 counts over raw nuclei) then nucleus filter (200..10,000 UMI)
gsum = np.asarray(X.sum(axis=0)).ravel()
gkeep = gsum >= 3
X = X[:, gkeep].tocsr()
genes = uniq[gkeep]
umi = np.asarray(X.sum(axis=1)).ravel()
ngen = np.diff(X.indptr)
obs["umi"] = umi
obs["n_genes"] = ngen
qc = (umi >= 200) & (umi <= 10000)
obs["qc_umi"] = qc
obs["qc_gene_rule"] = (ngen >= 200) & (umi <= 10000)
gene_rule_counts = obs.groupby("gsm").qc_gene_rule.sum()
print(f"genes kept (>=3 counts): {gkeep.sum():,}; nuclei passing UMI rule: {qc.sum():,}", flush=True)
X = X[qc].tocsr()
obs = obs[qc].reset_index(drop=True)
obs = obs.merge(design[["gsm", "title", "group", "sex", "sex_f", "series", "stage"]], on="gsm", how="left")

# ---- PC1 sex genes (unstripped counts, all QC nuclei) ----------------------------
gi = {g: i for i, g in enumerate(genes)}
ygenes = ["RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "USP9Y"]
codes = pd.Categorical(obs.gsm, categories=design.gsm).codes
ind = sparse.csr_matrix((np.ones(len(obs)), (codes, np.arange(len(obs)))), shape=(len(design), len(obs)))
tot = np.asarray(ind @ X.sum(axis=1)).ravel()
def dsum(gl):
    idx = [gi[g] for g in gl if g in gi]
    if not idx:
        return np.zeros(len(design))
    return np.asarray((ind @ X[:, idx]).sum(axis=1)).ravel()
xist = dsum(["XIST"]) / tot * 1e6
ysum = dsum(ygenes) / tot * 1e6
pc1 = design[["gsm", "title", "group", "sex"]].copy()
pc1["xist_cpm"] = xist
pc1["y_cpm_sum"] = ysum
pc1["y_genes_present"] = ",".join(g for g in ygenes if g in gi)
pc1["log2_ratio"] = np.log2((xist + 1) / (ysum + 1))
pc1["sex_call"] = np.where(pc1.log2_ratio > 0, "female", "male")
pc1["concordant"] = pc1.sex_call == pc1.sex
pc1.to_csv(out / "D5_PC1_sex.tsv", sep="\t", index=False)
n_disc = int((~pc1.concordant).sum())
print(f"PC1: {48 - n_disc} of 48 concordant", flush=True)

# ---- CellTypist ------------------------------------------------------------------
import scanpy as sc
import celltypist
A = ad.AnnData(X=X.copy(), obs=obs.set_index(obs.gsm + "_" + obs.barcode), var=pd.DataFrame(index=genes))
sc.pp.normalize_total(A, target_sum=1e4)
sc.pp.log1p(A)
model = celltypist.models.Model.load(model=str(C.CT_MODEL))
np.random.seed(C.SEED)
pred = celltypist.annotate(A, model=model, majority_voting=True)
pl = pred.predicted_labels
obs["label"] = pl["majority_voting"].astype(str).values
obs["label_raw"] = pl["predicted_labels"].astype(str).values
obs["conf"] = pred.probability_matrix.max(axis=1).values
del A
print("CellTypist labels:", obs.label.value_counts().to_dict(), flush=True)

# ---- library size on the shared stripped gene set ----------------------------
lib_genes = (out / "library_genes.txt").read_text().split()
lib_idx = np.array([gi[g] for g in lib_genes if g in gi])
obs["lib"] = np.asarray(X[:, lib_idx].sum(axis=1)).ravel()
print(f"library genes present after gene filter: {len(lib_idx)} of {len(lib_genes)}", flush=True)

# ---- lineage counts and inclusion ------------------------------------------------
lin = pd.crosstab(obs.gsm, obs.label).reindex(design.gsm).fillna(0).astype(int)
lin.to_csv(out / "lineage_counts_per_donor.tsv", sep="\t")
dq = design[["gsm", "title", "group", "stage", "sex", "series"]].copy()
dq["filtered_barcodes"] = filt_bc.reindex(dq.gsm).values
qc_all = obs.groupby("gsm").size()
dq["nuclei_qc_umi"] = qc_all.reindex(dq.gsm).fillna(0).astype(int).values
dq["n_hep"] = lin.get("Hepatocytes", pd.Series(0, index=lin.index)).reindex(dq.gsm).values
hep = obs[obs.label == "Hepatocytes"]
dq["median_umi_hep"] = hep.groupby("gsm").umi.median().reindex(dq.gsm).values
dq["median_genes_hep"] = hep.groupby("gsm").n_genes.median().reindex(dq.gsm).values
dq["median_lib_hep"] = hep.groupby("gsm").lib.median().reindex(dq.gsm).values
dq["pc1_concordant"] = pc1.concordant.values
dq["excluded_reason"] = ""
if n_disc >= 3:
    g1 = False
else:
    g1 = True
    dq.loc[~dq.pc1_concordant, "excluded_reason"] = "PC1 sex discordant"
dq.loc[(dq.n_hep < C.MIN_NUCLEI) & (dq.excluded_reason == ""), "excluded_reason"] = "fewer than 30 hepatocyte nuclei"
dq["included"] = dq.excluded_reason == ""
inc = dq[dq.included]
n_ctrl = int((inc.group == "ctrl").sum()); n_masl = int((inc.group == "MASL").sum()); n_mash = int(inc.group.isin(["eMASH", "aMASH"]).sum())
g2 = n_ctrl >= 18 and n_masl >= 5 and n_mash >= 15

# ---- coverage (D3, G3) -------------------------------------------------------------
mem = C.program_table()
cov_rows = []
det_by_lin = {}
for ct, lab in C.LINEAGE_LABEL.items():
    m = (obs.label == lab).to_numpy()
    if m.sum() == 0:
        det_by_lin[ct] = np.zeros(len(genes))
    else:
        det_by_lin[ct] = np.asarray((X[m] > 0).mean(axis=0)).ravel()
for (ct, mod), _ in mem.groupby(["cell_type", "module"]):
    w = C.program_weights(mem, ct, mod)
    pres = w.mapped_symbol.isin(gi)
    detd = w.mapped_symbol.map(lambda g: det_by_lin[ct][gi[g]] >= 0.01 if g in gi else False)
    meas = w[pres & detd]
    cov_rows.append(dict(cell_type=ct, module=mod, n_mapped=len(w), n_present=int(pres.sum()),
                         n_scored=len(meas), retained_weight=float(meas.original_l1_weight.sum()),
                         lineage_nuclei=int((obs.label == C.LINEAGE_LABEL[ct]).sum()),
                         testable=bool(len(meas) >= 8 and meas.original_l1_weight.sum() >= 0.20)))
cov = pd.DataFrame(cov_rows)
cov.to_csv(out / "D3_coverage.tsv", sep="\t", index=False)
def g3(ct, mod):
    r = cov[(cov.cell_type == ct) & (cov.module == mod)].iloc[0]
    return bool(r.n_scored >= 8 and r.retained_weight >= 0.80), r
g3p, r20 = g3(*C.PRIMARY)
g3s, r8 = g3(*C.SECONDARY)
g3n, r2 = g3(*C.NC1)

# ---- D1 / D2 ------------------------------------------------------------------------
# the 200-gene rule (>= 200 genes, <= 10,000 UMI) beside the UMI rule, counted on all raw nuclei
dq["nuclei_qc_gene_rule"] = gene_rule_counts.reindex(dq.gsm).fillna(0).astype(int).values
dq.to_csv(out / "donor_qc.tsv", sep="\t", index=False)
d1 = dq.groupby("group").agg(donors=("gsm", "size"), filtered_barcodes_median=("filtered_barcodes", "median"),
                             nuclei_qc_median=("nuclei_qc_umi", "median"), n_hep_median=("n_hep", "median"),
                             n_hep_min=("n_hep", "min"), median_umi_hep=("median_umi_hep", "median"),
                             median_genes_hep=("median_genes_hep", "median"), included=("included", "sum"))
d1.to_csv(out / "D1_qc_by_group.tsv", sep="\t")
frac = pd.crosstab(obs.gsm, obs.label, normalize="index").reindex(design.gsm)
frac.to_csv(out / "D2_label_fractions_per_donor.tsv", sep="\t")
fg = frac.join(design.set_index("gsm").group).groupby("group").mean()
fg.to_csv(out / "D2_label_fractions_by_group.tsv", sep="\t")
confg = hep.merge(design[["gsm"]], on="gsm").groupby("group").conf.median()
confg.to_csv(out / "D2_hep_conf_by_group.tsv", sep="\t")

# ---- PC2 collagen in fibroblast nuclei ----------------------------------------------
def hc3(y, Xd):
    XtX = np.linalg.inv(Xd.T @ Xd)
    b = XtX @ Xd.T @ y
    e = y - Xd @ b
    h = np.einsum("ij,jk,ik->i", Xd, XtX, Xd)
    meat = Xd.T @ (Xd * (e ** 2 / (1 - h) ** 2)[:, None])
    V = XtX @ meat @ XtX
    df = len(y) - Xd.shape[1]
    sigma = np.sqrt((e ** 2).sum() / df)
    return b, np.sqrt(np.diag(V)), sigma, df, h

fib = obs[obs.label == "Fibroblasts"]
fcodes = pd.Categorical(fib.gsm, categories=design.gsm).codes
find = sparse.csr_matrix((np.ones(len(fib)), (fcodes, np.arange(len(fib)))), shape=(len(design), len(fib)))
Xf = X[(obs.label == "Fibroblasts").to_numpy()]
col_c = np.asarray((find @ Xf[:, [gi["COL1A1"], gi["COL1A2"]]]).sum(axis=1)).ravel()
tot_f = np.asarray((find @ Xf).sum(axis=1)).ravel()
pc2 = design[["gsm", "title", "group", "sex", "sex_f"]].copy()
pc2["n_fib"] = np.asarray(find.sum(axis=1)).ravel().astype(int)
pc2["col1_counts"] = col_c
pc2["fib_total_umi"] = tot_f
pc2["log2_cpm_col1"] = np.log2(np.divide(col_c, tot_f, out=np.zeros_like(col_c), where=tot_f > 0) * 1e6 + 1)
pc2.to_csv(out / "D6_PC2_collagen.tsv", sep="\t", index=False)
sub = pc2[(pc2.n_fib >= 20) & pc2.group.isin(["MASL", "aMASH"])]
n_m, n_a = int((sub.group == "MASL").sum()), int((sub.group == "aMASH").sum())
pc2_res = dict(n_MASL=n_m, n_aMASH=n_a)
if n_m >= 4 and n_a >= 6:
    Xd = np.column_stack([np.ones(len(sub)), (sub.group == "aMASH").astype(float), sub.sex_f.astype(float)])
    b, se, sig, df, _ = hc3(sub.log2_cpm_col1.to_numpy(float), Xd)
    pc2_res.update(beta_aMASH_vs_MASL=b[1], hc3_se=se[1], hc3_t=b[1] / se[1], hc3_p=2 * stats.t.sf(abs(b[1] / se[1]), df),
                   std_diff=b[1] / sig, df=df, readable=True, prediction_met=bool(b[1] > 0 and b[1] / sig > 0.5))
else:
    pc2_res.update(readable=False, prediction_met=False)
pd.DataFrame([pc2_res]).to_csv(out / "PC2_collagen_result.tsv", sep="\t", index=False)
print("PC2:", pc2_res, flush=True)

# ---- D7 ambient / doublet proxies in hepatocyte nuclei ----------------------------
prog_genes = set()
for k in (C.PRIMARY, C.SECONDARY, C.NC1):
    prog_genes |= set(C.program_weights(mem, *k).mapped_symbol)
d7_genes = [g for g in ["KRT7", "KRT19", "EPCAM", "PTPRC", "PECAM1"] if g not in prog_genes and g in gi]
hm = (obs.label == "Hepatocytes").to_numpy()
Xh = X[hm]
hob = obs[hm].reset_index(drop=True)
norm_d7 = C.log_normalize(Xh[:, [gi[g] for g in d7_genes]].tocsr(), hob.lib.to_numpy(float))
hcodes = pd.Categorical(hob.gsm, categories=design.gsm).codes
hind = sparse.csr_matrix((np.ones(len(hob)), (hcodes, np.arange(len(hob)))), shape=(len(design), len(hob)))
nh = np.asarray(hind.sum(axis=1)).ravel()
d7 = pd.DataFrame(np.asarray((hind @ norm_d7).todense()) / np.maximum(nh, 1)[:, None], columns=d7_genes)
d7.insert(0, "gsm", design.gsm.values)
d7.insert(1, "group", design.group.values)
krt_any = np.asarray((Xh[:, [gi["KRT7"], gi["KRT19"]]] > 0).sum(axis=1)).ravel() > 0
d7["frac_hep_any_KRT7_KRT19"] = np.asarray(hind @ krt_any.astype(float)).ravel() / np.maximum(nh, 1)
d7.to_csv(out / "D7_ambient_proxies_per_donor.tsv", sep="\t", index=False)
d7.drop(columns="gsm").groupby("group").mean().to_csv(out / "D7_ambient_proxies_by_group.tsv", sep="\t")

# ---- D8 leverage, D9 balance, realized MDE ---------------------------------------
inc = dq[dq.included].copy()
inc["log10_n_hep"] = np.log10(inc.n_hep)
inc["sex_f"] = (inc.sex == "female").astype(int)
Xd = np.column_stack([np.ones(len(inc)), inc.stage, inc.sex_f, inc.log10_n_hep]).astype(float)
h = np.einsum("ij,jk,ik->i", Xd, np.linalg.inv(Xd.T @ Xd), Xd)
inc["leverage"] = h
inc[["gsm", "title", "group", "stage", "sex", "n_hep", "log10_n_hep", "leverage"]].to_csv(out / "D8_leverage.tsv", sep="\t", index=False)
r_stage_n = float(np.corrcoef(inc.stage, inc.log10_n_hep)[0, 1])
Z = Xd[:, [0, 2, 3]]
res = inc.stage.to_numpy(float) - Z @ np.linalg.lstsq(Z, inc.stage.to_numpy(float), rcond=None)[0]
r2_stage = 1 - (res ** 2).sum() / ((inc.stage - inc.stage.mean()) ** 2).sum()
vif = 1 / (1 - r2_stage)
n = len(inc)
dfm = n - 4
sx = np.sqrt(((inc.stage - inc.stage.mean()) ** 2).sum())
tc = stats.t.ppf(0.975, dfm)
from scipy.optimize import brentq
def _pow(bb):  # two-sided power; the lower tail underflows to NaN in scipy at large ncp, so it is set to 0
    lo = stats.nct.cdf(-tc, dfm, bb * sx)
    return stats.nct.sf(tc, dfm, bb * sx) + (0.0 if not np.isfinite(lo) else lo)
mde_nocov = brentq(lambda bb: _pow(bb) - 0.8, 0.01, 3)
d9 = inc.groupby("stage").agg(n=("gsm", "size"), log10_n_hep_mean=("log10_n_hep", "mean"),
                              median_umi_hep=("median_umi_hep", "median"), female=("sex_f", "sum"))
d9.to_csv(out / "D9_balance.tsv", sep="\t")
summ = dict(G1_pass=g1, n_discordant_sex=n_disc, discordant=";".join(pc1.title[~pc1.concordant]),
            G2_pass=g2, n_ctrl=n_ctrl, n_MASL=n_masl, n_MASH=n_mash, n_eMASH=int((inc.group == "eMASH").sum()), n_aMASH=int((inc.group == "aMASH").sum()),
            G3_primary_pass=g3p, primary_scored=int(r20.n_scored), primary_weight=float(r20.retained_weight),
            G3_secondary_pass=g3s, secondary_scored=int(r8.n_scored), secondary_weight=float(r8.retained_weight),
            NC1_scored=int(r2.n_scored), NC1_weight=float(r2.retained_weight), NC1_coverage_ok=g3n,
            r_stage_log10nhep=r_stage_n, VIF_stage=vif, MDE_realized_nocov=mde_nocov, MDE_realized_x_sqrtVIF=mde_nocov * np.sqrt(vif),
            gate_MDE_x_sqrtVIF=0.443347672664669 * np.sqrt(vif), n_included=n,
            hep_nuclei_total_included=int(inc.n_hep.sum()))
pd.Series(summ).to_csv(out / "gate_G1_G3.tsv", sep="\t", header=["value"])
print(pd.Series(summ).to_string(), flush=True)

# ---- save QC nuclei for job 6302 ---------------------------------------------------
if QC_H5AD.exists():
    sys.exit(f"refusing to overwrite {QC_H5AD}")
A = ad.AnnData(X=X.astype(np.float32), obs=obs.astype({"gsm": str}).set_index(obs.gsm + "_" + obs.barcode), var=pd.DataFrame(index=genes))
A.write_h5ad(QC_H5AD)
print(f"wrote {QC_H5AD}", flush=True)
for k, v in (("G1", g1), ("G2", g2), ("G3", g3p)):
    if not v:
        msg = {"G1": "untestable: sample labels cannot be linked to matrices",
               "G2": "untestable: too few donors", "G3": "untestable: program not measured"}[k]
        (out / f"STOP_{k}.txt").write_text(msg + "\n")
        sys.exit(f"STOP {k}: {msg}")
print("G1-G3 pass")
