"""Feasibility checks for candidate L03-04 (login node, small files only)."""
import numpy as np
import pandas as pd
from scipy.linalg import helmert
from scipy.stats import norm

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
P = ROOT + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/"
NEW = ROOT + "Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz"
OLD = ROOT + "Analysis/Deconvolution/results/{g}/{g}_bayesprism_proportions.tsv"
MEMB = ROOT + "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/program_membership_v2.tsv"
DS = ["GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"]

# 1. calls
deg = pd.read_csv(P + "deg_results.csv")
calls = deg[(deg.padj < 0.05) & (deg.logFC.abs() > 0.5)]
print("genes in fit", len(deg))
print("calls", len(calls), "up", (calls.logFC > 0).sum(), "down", (calls.logFC < 0).sum())
print("calls SE median", round(calls.SE.median(), 4), "IQR", calls.SE.quantile([.25, .75]).round(4).tolist())
pstar = deg.loc[deg.padj < 0.05, "P.Value"].max()
print("n padj<0.05 (any FC)", (deg.padj < 0.05).sum(), "largest raw p with padj<0.05", pstar)

# 2. design table
dz = pd.read_csv(P + "model_design.tsv", sep="\t").set_index("sample_id")
w = pd.read_csv(P + "lvqw_sample_weights.tsv", sep="\t").set_index("sample_id").sample_weight
new = pd.read_csv(NEW, sep="\t").set_index("sample_id")
ct = [c for c in new.columns if c != "dataset"]
print("new fractions", new.shape, "cell types", len(ct))
old_parts = []
for g in DS:
    o = pd.read_csv(OLD.format(g=g), sep="\t", index_col=0)
    old_parts.append(o)
old = pd.concat(old_parts)
print("earlier fraction rows (all 5 files)", old.shape, "columns identical to new", list(old.columns) == ct)
ids = dz.index
print("774 ids in new", ids.isin(new.index).sum(), "in earlier", ids.isin(old.index).sum())
old = old.loc[ids, ct]
new = new.loc[ids, ct]

def ilr(f):
    return np.log(np.maximum(f.to_numpy(float), 1e-6)) @ (-helmert(16).T)

for name, f in [("rebuilt", new), ("earlier", old)]:
    print(name, "rowsum range", f.sum(1).min().round(6), f.sum(1).max().round(6))

# 3. per cell type agreement between inputs across 774
agree = pd.DataFrame({
    "pearson_logfloor": [np.corrcoef(np.log(np.maximum(new[c], 1e-6)), np.log(np.maximum(old[c], 1e-6)))[0, 1] for c in ct],
    "mean_rebuilt": new.mean().values, "mean_earlier": old.mean().values}, index=ct)
print(agree.round(3).to_string())

# 4. R2 of disease on ILR within dataset + sex, weighted and unweighted -> SE inflation
X0 = dz.drop(columns=["group_binaryDisease"]).to_numpy(float)
y = dz["group_binaryDisease"].to_numpy(float)

def resid_var(X, y, wt):
    sw = np.sqrt(wt)
    b, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
    r = (y - X @ b) * sw
    return (r ** 2).sum()

for name, f in [("rebuilt", new), ("earlier", old)]:
    Z = ilr(f)
    X1 = np.hstack([X0, Z])
    for wname, wt in [("unweighted", np.ones(len(y))), ("lvqw", w.loc[ids].to_numpy(float))]:
        rv0 = resid_var(X0, y, wt); rv1 = resid_var(X1, y, wt)
        infl = np.sqrt(rv0 / rv1)
        print(f"{name} {wname}: partial R2 of group on 15 ILR | dataset+sex = {1 - rv1 / rv0:.3f}; SE inflation sqrt(VIF ratio) = {infl:.3f}")
        if wname == "lvqw":
            for tag, p in [("BH cut from current fit", pstar), ("Bonferroni 24115", 0.05 / 24115)]:
                k = norm.ppf(1 - p / 2) + norm.ppf(0.8)
                print(f"   per-gene MDE at {tag} (two-sided p {p:.2e}): unadjusted {k * calls.SE.median():.3f}, adjusted {k * calls.SE.median() * infl:.3f} log2 units")

# 5. cohort-level structure: groups per dataset
print(pd.crosstab(dz.index.map(new.assign(ds=pd.read_csv(NEW, sep='\t').set_index('sample_id').dataset).ds), y))

# 6. positive control genes
m = pd.read_csv(MEMB, sep="\t")
ecm = m[m.program_uid == "hotspot_hepatocytes_f05c535ae5bbc0b9"]
print("Stromal ECM genes", ecm.mapped_symbol.nunique())
sym_calls = set(calls.symbol.astype(str))
in_fit = ecm.mapped_symbol.isin(set(deg.symbol.astype(str))).sum()
hit = ecm[ecm.mapped_symbol.isin(sym_calls)]
print("Stromal ECM genes in fit", in_fit, "in 1854 calls", len(hit))
cc = calls.set_index("symbol").loc[hit.mapped_symbol, ["logFC", "SE", "padj"]]
print(cc.round(3).to_string())
