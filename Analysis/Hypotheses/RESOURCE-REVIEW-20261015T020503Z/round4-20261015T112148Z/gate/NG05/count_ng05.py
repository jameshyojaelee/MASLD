"""NG05 (L04-03) feasibility counts. Gene-level observability per Catalog class.
Counting and a power check only; not the analysis."""
import gzip, glob, json, os
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import norm

np.random.seed(20260923)
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = os.path.dirname(os.path.abspath(__file__))
cat = pd.read_parquet(f"{ROOT}/RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z/passport_gene_index.parquet")
cat = cat[["ensembl_id", "symbol", "primary_evidence_class"]].copy()
cat["primary_evidence_class"] = cat["primary_evidence_class"].astype(str)
gen = pd.read_csv(f"{ROOT}/data/gencode_v49_gene_metadata.tsv.gz", sep="\t")
cat = cat.merge(gen[["ensembl_base", "gene_biotype"]], left_on="ensembl_id", right_on="ensembl_base", how="left")
deg = pd.read_csv(f"{ROOT}/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
deg["ens"] = deg["gene"].str.split(".").str[0]
cat = cat.merge(deg[["ens", "AveExpr"]], left_on="ensembl_id", right_on="ens", how="left")

def symset(path, col="Genes"):
    t = pd.read_csv(path, sep="\t", usecols=[col])
    s = set()
    for v in t[col].dropna().astype(str):
        s.update(x.strip() for x in v.split(";"))
    return s

def symset_frac(path, thr=0.5):
    t = pd.read_csv(path, sep="\t")
    q = t.iloc[:, 3:]
    frac = q.notna().mean(axis=1)
    s = set()
    for v, f in zip(t["Genes"].astype(str), frac):
        if f >= thr:
            s.update(x.strip() for x in v.split(";"))
    return s, q.shape[1]

liver = symset(f"{ROOT}/data/PXD051911/liver_protein_quant.txt")
plasma = symset(f"{ROOT}/data/PXD051911/plasma_protein_quant.txt")
liver50, n_liver = symset_frac(f"{ROOT}/data/PXD051911/liver_protein_quant.txt")
plasma50, n_plasma = symset_frac(f"{ROOT}/data/PXD051911/plasma_protein_quant.txt")
with gzip.open(f"{ROOT}/data/external/govaere2026_natgenetics/cosmx/GSM9351704_Leuven_1_exprMat_file.csv.gz", "rt") as fh:
    cols = fh.readline().strip().split(",")
cosmx = {c for c in cols if not any(k in c for k in ("NegPrb", "FalseCode", "SystemControl", "fov", "cell_ID", "cell"))}

# single-cell lineage arm: saved per-lineage pseudobulk (columns are sequencing samples/runs)
sc_det = {}
sc_info = {}
for f in sorted(glob.glob(f"{ROOT}/Analysis/SingleCell/results_gpu_v2/pseudobulk/*_pseudobulk.csv")):
    lin = os.path.basename(f).replace("_pseudobulk.csv", "")
    m = pd.read_csv(f, index_col=0)
    sc_info[lin] = m.shape[1]
    frac = (m >= 10).mean(axis=1)
    for g, v in frac.items():
        sc_det[g] = max(sc_det.get(g, 0.0), v)
sc = {g for g, v in sc_det.items() if v >= 0.5}

cat["cosmx"] = cat["symbol"].isin(cosmx).astype(int)
cat["liver_prot"] = cat["symbol"].isin(liver).astype(int)
cat["liver_prot50"] = cat["symbol"].isin(liver50).astype(int)
cat["plasma_prot"] = cat["symbol"].isin(plasma).astype(int)
cat["plasma_prot50"] = cat["symbol"].isin(plasma50).astype(int)
cat["sc_lineage"] = cat["symbol"].isin(sc).astype(int)
cat["sc_symbol_present"] = cat["symbol"].isin(sc_det.keys()).astype(int)

assays = ["cosmx", "liver_prot", "liver_prot50", "plasma_prot", "plasma_prot50", "sc_lineage"]
res = {"n_liver_samples": n_liver, "n_plasma_samples": n_plasma, "cosmx_panel_genes": len(cosmx),
       "liver_symbols": len(liver), "plasma_symbols": len(plasma), "sc_lineage_columns": sc_info,
       "sc_genes_detected": len(sc)}
cat.groupby("primary_evidence_class")[assays + ["sc_symbol_present"]].agg(["sum", "mean"]).to_csv(f"{OUT}/coverage_by_class.tsv", sep="\t")
cat.groupby("primary_evidence_class")["AveExpr"].describe().to_csv(f"{OUT}/aveexpr_by_class.tsv", sep="\t")
cat.groupby(["primary_evidence_class", "gene_biotype"]).size().unstack(fill_value=0).loc[:, lambda d: d.sum() > 20].to_csv(f"{OUT}/biotype_by_class.tsv", sep="\t")

# analysis set: anchored vs established, with AveExpr
a = cat[cat["primary_evidence_class"].isin(["genetically_anchored", "established_state_associated"])].copy()
a["anch"] = (a["primary_evidence_class"] == "genetically_anchored").astype(int)
a["pc"] = (a["gene_biotype"] == "protein_coding").astype(int)
res["n_by_class"] = a["primary_evidence_class"].value_counts().to_dict()
res["n_missing_aveexpr"] = a.groupby("primary_evidence_class")["AveExpr"].apply(lambda s: int(s.isna().sum())).to_dict()
res["pc_by_class"] = a.groupby("primary_evidence_class")["pc"].sum().to_dict()
a = a.dropna(subset=["AveExpr"])
# common support: share of each class within the other's 5-95% AveExpr range
q = a.groupby("anch")["AveExpr"].quantile([0.05, 0.95]).unstack()
res["support"] = {
    "anch_in_est_5_95": float(((a.anch == 1) & a.AveExpr.between(q.loc[0, 0.05], q.loc[0, 0.95])).sum() / (a.anch == 1).sum()),
    "est_in_anch_5_95": float(((a.anch == 0) & a.AveExpr.between(q.loc[1, 0.05], q.loc[1, 0.95])).sum() / (a.anch == 0).sum()),
}
a["dec"] = pd.qcut(a["AveExpr"], 10, labels=False)
dec_tab = a.groupby(["dec", "anch"]).size().unstack(fill_value=0)
dec_tab.to_csv(f"{OUT}/decile_by_class.tsv", sep="\t")
z = norm.ppf(1 - 0.0125 / 2) + norm.ppf(0.8)
res["models"] = {}
for y in assays:
    r = {}
    try:
        m0 = smf.logit(f"{y} ~ anch", data=a).fit(disp=0)
        r["unadj_OR"] = float(np.exp(m0.params["anch"]))
        m = smf.logit(f"{y} ~ anch + bs(AveExpr, df=4) + pc", data=a).fit(disp=0, maxiter=200)
        b, se = m.params["anch"], m.bse["anch"]
        r.update(adj_OR=float(np.exp(b)), ci=[float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))], p=float(m.pvalues["anch"]),
                 se=float(se), MDE_OR_holm0125=float(np.exp(z * se)), MDE_OR_lower=float(np.exp(-z * se)),
                 aveexpr_slope_sign_unadj=float(smf.logit(f"{y} ~ AveExpr", data=a).fit(disp=0).params["AveExpr"]),
                 n_events={int(k): int(v) for k, v in a.groupby("anch")[y].sum().items()})
    except Exception as e:
        r["error"] = repr(e)
    res["models"][y] = r
# positive control on unresolved protein-coding genes: detection by AveExpr decile
u = cat[(cat.primary_evidence_class == "unresolved") & (cat.gene_biotype == "protein_coding")].dropna(subset=["AveExpr"]).copy()
u["dec"] = pd.qcut(u["AveExpr"], 10, labels=False)
pcdec = u.groupby("dec")[["cosmx", "liver_prot", "plasma_prot", "sc_lineage"]].mean()
pcdec.to_csv(f"{OUT}/posctrl_detection_by_decile_unresolved_pc.tsv", sep="\t")
res["posctrl_n_unresolved_pc"] = len(u)
json.dump(res, open(f"{OUT}/gate_counts.json", "w"), indent=1, default=str)
cat.to_csv(f"{OUT}/per_gene_observability.tsv.gz", sep="\t", index=False)
print(json.dumps(res, indent=1, default=str))
