"""n by cohort and group for every arm, and design-only SE inflation per arm (no expression data, no adjusted outcome)."""
import numpy as np, pandas as pd
from scipy.linalg import helmert
from scipy.stats import norm
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
P = ROOT + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/"
NEW = ROOT + "Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz"
OLD = ROOT + "Analysis/Deconvolution/results/{g}/{g}_bayesprism_proportions.tsv"
MAN = ROOT + "Analysis/Deconvolution/corrected_774_20261009T141331Z/bulk/sample_manifest.tsv"
SM = ROOT + "figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/stage_extension_sample_manifest.tsv"
MEMB = ROOT + "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/program_membership_v2.tsv"
DS = ["GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"]
dz = pd.read_csv(P + "model_design.tsv", sep="\t").set_index("sample_id")
w = pd.read_csv(P + "lvqw_sample_weights.tsv", sep="\t").set_index("sample_id").sample_weight.loc[dz.index]
man = pd.read_csv(MAN, sep="\t").set_index("sample_id").loc[dz.index]
new = pd.read_csv(NEW, sep="\t").set_index("sample_id")
ct = [c for c in new.columns if c != "dataset"]
new = new.loc[dz.index, ct]
old = pd.concat([pd.read_csv(OLD.format(g=g), sep="\t", index_col=0) for g in DS]).loc[dz.index, ct]
ds = man.dataset
grp = dz.group_binaryDisease
print("units", len(dz), "unique analysis_unit_id", man.analysis_unit_id.nunique())
print(pd.crosstab(ds, grp.map({0: "control", 1: "disease"}), margins=True).to_string())
main = ds.isin(["GSE130970", "GSE162694", "GSE213621"]); sep = ds.isin(["GSE126848", "GSE135251"])
print("DR-main units", int(main.sum()), "control", int((main & (grp == 0)).sum()), "disease", int((main & (grp == 1)).sum()))
print("DR-sep units", int(sep.sum()), "control", int((sep & (grp == 0)).sum()), "disease", int((sep & (grp == 1)).sum()))
deg = pd.read_csv(P + "deg_results.csv")
calls = deg[(deg.padj < 0.05) & (deg.logFC.abs() > 0.5)]
print("genes", len(deg), "calls", len(calls), "up", int((calls.logFC > 0).sum()), "down", int((calls.logFC < 0).sum()))
print("calls |logFC| in (0.5,0.6]", int((calls.logFC.abs() <= 0.6).sum()))
m = pd.read_csv(MEMB, sep="\t")
ecm = m[m.program_uid == "hotspot_hepatocytes_f05c535ae5bbc0b9"].source_gene.unique()
print("Stromal ECM source genes", len(ecm), "in calls", int(calls.symbol.isin(ecm).sum()), sorted(calls.symbol[calls.symbol.isin(ecm)]))
sm = pd.read_csv(SM, sep="\t"); s = sm[sm.contrast == "F4_vs_F0"]
print("F4_vs_F0:", pd.crosstab(s.dataset, s.group).to_string().replace("\n", " | "))
# design-only SE inflation per arm, sample weights from the stored fit
LIN7 = {"hep_epi": ["Hepatocytes", "Cholangiocytes"], "endo": ["Endothelial cells"], "fib": ["Fibroblasts"],
        "mac_mono": ["Macrophages", "Mono+mono derived cells"], "oth_my": ["cDC1s", "cDC2s", "pDCs", "Neutrophils", "Basophils"],
        "T_NK": ["T cells", "Circulating NK/NKT", "Resident NK"], "B_pl": ["B cells", "Plasma cells"]}
DR5 = {"hep_epi": LIN7["hep_epi"], "stromal": LIN7["fib"] + LIN7["endo"], "mac_mono": LIN7["mac_mono"],
       "oth_my": LIN7["oth_my"], "lymph": LIN7["T_NK"] + LIN7["B_pl"]}
DR4 = {"hep_epi": LIN7["hep_epi"], "mac_mono": LIN7["mac_mono"], "oth_my": LIN7["oth_my"],
       "npc_nonmy": LIN7["fib"] + LIN7["endo"] + LIN7["T_NK"] + LIN7["B_pl"]}
def amalg(f, parts):
    a = np.column_stack([f[v].sum(1) for v in parts.values()]); return a / a.sum(1, keepdims=True)
def ilr(x, floor):
    x = np.maximum(np.asarray(x, float), floor); D = x.shape[1]
    return np.log(x) @ (-helmert(D).T)
X0 = dz.drop(columns=["group_binaryDisease"]).to_numpy(float); y = grp.to_numpy(float); wt = w.to_numpy(float)
def infl(X0, Z, y, wt):
    def rv(X):
        sw = np.sqrt(wt); b, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None); r = (y - X @ b) * sw; return r @ r
    return np.sqrt(rv(X0) / rv(np.hstack([X0, Z])))
def unscaled(X, wt):  # unscaled SE of the group coefficient (last column) under sample weights
    Xw = X * np.sqrt(wt)[:, None]; return np.sqrt(np.linalg.inv(Xw.T @ Xw)[-1, -1])
cse = calls.SE.median()
print(f"\nmedian call SE {cse:.4f}")
rows = []
for nm, fr in [("rebuilt", new), ("earlier", old)]:
    arms = {"A16": (ilr(fr.to_numpy(), 1e-6), np.ones(len(y), bool)),
            "C7": (ilr(amalg(fr, LIN7), 1e-6), np.ones(len(y), bool)),
            "D1": (ilr(fr.to_numpy(), 1e-4), np.ones(len(y), bool)),
            "DRmain": (ilr(amalg(fr, DR5), 1e-6), main.to_numpy()),
            "DRsep": (ilr(amalg(fr, DR4), 1e-6), sep.to_numpy())}
    for arm, (Z, keep) in arms.items():
        Xk = X0[keep][:, np.abs(X0[keep]).sum(0) > 0]
        k = infl(Xk, Z[keep], y[keep], wt[keep])
        su = unscaled(np.column_stack([Xk, y[keep]]), wt[keep]) / unscaled(np.column_stack([X0, y]), wt)
        se_u = cse * su
        rows.append([nm, arm, int(keep.sum()), round(k, 3), round(su, 3), round(se_u * k, 4),
                     round(2.80 * se_u * np.sqrt(k**2 - 1), 3), round((norm.ppf(1 - 0.025 / 1854) + norm.ppf(0.8)) * se_u * np.sqrt(k**2 - 1), 3),
                     round(se_u * np.sqrt(k**2 - 0.75), 4)])
t = pd.DataFrame(rows, columns=["input", "arm", "n_units", "SE_inflation_adj", "SE_unadj_vs_774", "median_SE_adj",
                                "MDE_attenuation_nominal", "MDE_attenuation_bonf1854", "SE_keep_contrast"])
print(t.to_string(index=False))
