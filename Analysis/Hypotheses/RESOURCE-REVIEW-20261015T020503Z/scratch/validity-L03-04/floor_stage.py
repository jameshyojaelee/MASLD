import pandas as pd, numpy as np
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
new = pd.read_csv(ROOT + "Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz", sep="\t").set_index("sample_id")
man = pd.read_csv(ROOT + "Analysis/Deconvolution/corrected_774_20261009T141331Z/bulk/sample_manifest.tsv", sep="\t").set_index("sample_id").loc[new.index]
old = pd.concat([pd.read_csv(ROOT + f"Analysis/Deconvolution/results/{g}/{g}_bayesprism_proportions.tsv", sep="\t", index_col=0) for g in sorted(new.dataset.unique())]).loc[new.index]
g = man.dataset == "GSE135251"
print("GSE135251 fibroblast fraction < 1e-6 by fibrosis stage (rebuilt / earlier / n)")
for s in sorted(man.fibrosis_stage[g].dropna().unique()):
    i = g & (man.fibrosis_stage == s)
    print(f"  F{int(s)}: {(new.Fibroblasts[i] < 1e-6).sum()} / {(old.Fibroblasts[i] < 1e-6).sum()} / {i.sum()}")
for nm, fr in [("rebuilt", new), ("earlier", old)]:
    strom = fr["Fibroblasts"] + fr["Endothelial cells"]
    print(f"{nm}: fibroblast+endothelial < 1e-4 in GSE135251: {(strom[g] < 1e-4).sum()} of {g.sum()}; "
          f"median fibroblast fraction GSE135251 {fr.Fibroblasts[g].median():.2e}, other cohorts {fr.Fibroblasts[~g].median():.2e}")
# correlation of rebuilt vs earlier log10 fibroblast within GSE135251 among units above floor in both
both = g & (new.Fibroblasts > 1e-6) & (old.Fibroblasts > 1e-6)
print("GSE135251 units above floor in both inputs:", both.sum(),
      "Pearson log10 fibroblast rebuilt vs earlier:", round(np.corrcoef(np.log10(new.Fibroblasts[both]), np.log10(old.Fibroblasts[both]))[0,1], 3))
