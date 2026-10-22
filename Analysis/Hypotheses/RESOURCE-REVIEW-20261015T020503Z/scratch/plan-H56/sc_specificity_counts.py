"""Count the 1,854 bulk calls by single-cell lineage specificity (donor-collapsed lineage profiles).
Reads only single-cell pseudobulk counts and the stored unadjusted call list; no adjusted fit."""
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
SC = ROOT / "Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk"
DEG = ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
PAIR = {ds: ROOT / f"data/{ds}/metadata/donor_pairing.csv" for ds in ("GSE244832", "GSE185477", "GSE202379", "GSE136103")}
LIN = ["Hepatocytes", "Fibroblasts", "Cholangiocytes", "Macrophages", "T_cells", "Endothelial_cells"]
s2d = {}
for ds, fp in PAIR.items():
    dp = pd.read_csv(fp, dtype=str)
    for _, r in dp.iterrows():
        for s in str(r["rna_srrs"]).split(";"):
            if s and s != "nan":
                s2d[s.strip()] = f"{ds}_{r['donor_id']}"
prof, nd = {}, {}
for ln in LIN:
    m = pd.read_csv(SC / f"pseudobulk_{ln}.tsv.gz", sep="\t", index_col=0)
    m = m.groupby(level=0).sum()
    cols = pd.Series(m.columns).map(s2d).fillna(pd.Series(m.columns)).values
    m = m.T.groupby(cols).sum().T
    cpm = m.div(m.sum(axis=0).replace(0, np.nan), axis=1) * 1e6
    prof[ln] = np.log2(cpm + 1).mean(axis=1)
    nd[ln] = m.shape[1]
print("donors per lineage profile:", nd)
P = pd.DataFrame(prof).dropna(how="all")
Z = P.sub(P.mean(1), axis=0).div(P.std(1).replace(0, np.nan), axis=0)
srt = np.sort(P.values, axis=1)
top_gap = pd.Series(srt[:, -1] - srt[:, -2], index=P.index)
rng = pd.Series(srt[:, -1] - srt[:, 0], index=P.index)
top = P.idxmax(1)
deg = pd.read_csv(DEG)
calls = deg[(deg.padj < 0.05) & (deg.logFC.abs() > 0.5)].copy()
sym_n = deg.symbol.value_counts()
calls["sym_unique"] = calls.symbol.map(sym_n).eq(1)
calls["in_sc"] = calls.symbol.isin(P.index) & calls.sym_unique
print("calls", len(calls), "unique symbol", int(calls.sym_unique.sum()), "found in sc profiles", int(calls.in_sc.sum()))
c = calls[calls.in_sc].set_index("symbol")
zmax = Z.loc[c.index].max(1)
cls = pd.Series("other", index=c.index)
spec = (zmax >= 1.5) & (top_gap.loc[c.index] >= 1)
cls[spec] = "specific_" + top.loc[c.index][spec]
uni = (rng.loc[c.index] < 1) & (P.loc[c.index].min(1) >= 1)
cls[uni & ~spec] = "uniform"
c["cls"] = cls
c["dir"] = np.where(c.logFC > 0, "up", "down")
print(pd.crosstab(c.cls, c.dir, margins=True).to_string())
print("\nlow-specificity alternatives (counts only, among calls found in sc profiles):")
for r_cut in (1.0, 1.5, 2.0):
    u = (rng.loc[c.index] < r_cut) & (P.loc[c.index].min(1) >= 1)
    print(f"  range < {r_cut} log2 and min >= 1: {int(u.sum())} (up {int((u & (c.logFC > 0)).sum())}, down {int((u & (c.logFC < 0)).sum())})")
for zc in (0.8, 1.0, 1.2):
    u = (zmax < zc)
    print(f"  max z < {zc}: {int(u.sum())}")
u15 = (rng.loc[c.index] < 1.5) & (P.loc[c.index].min(1) >= 1) & ~spec
print(f"\nFINAL uniform class (range < 1.5 log2, min >= 1, not specific): {int(u15.sum())} (up {int((u15 & (c.logFC > 0)).sum())}, down {int((u15 & (c.logFC < 0)).sum())})")
print("FINAL stromal-specific (fibroblasts or endothelial):", int(c.cls.isin(["specific_Fibroblasts", "specific_Endothelial_cells"]).sum()),
      "; hepatocyte-specific:", int((c.cls == "specific_Hepatocytes").sum()),
      "; immune-specific (macrophages or T cells):", int(c.cls.isin(["specific_Macrophages", "specific_T_cells"]).sum()),
      "; cholangiocyte-specific:", int((c.cls == "specific_Cholangiocytes").sum()))
ecm12 = ["C7","CACNA1C","CLSTN2","COL4A1","COL4A2","CRISPLD2","EPHA3","LAMA2","NAV3","PLCXD3","ROR2","SVEP1"]
print("Stromal ECM calls by class:", c.reindex(ecm12).cls.value_counts(dropna=False).to_dict())
