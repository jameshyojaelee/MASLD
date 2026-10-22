"""Per-cohort 5th percentile of each fraction part and of candidate compartments (inputs only, no outcome)."""
import pandas as pd, numpy as np
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
NEW = ROOT + "Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz"
OLD = ROOT + "Analysis/Deconvolution/results/{g}/{g}_bayesprism_proportions.tsv"
DS = ["GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"]
new = pd.read_csv(NEW, sep="\t").set_index("sample_id")
ds = new.dataset
ct = [c for c in new.columns if c != "dataset"]
old = pd.concat([pd.read_csv(OLD.format(g=g), sep="\t", index_col=0) for g in DS]).loc[new.index, ct]
new = new[ct]
comp = {
 "hepatocytes": ["Hepatocytes"],
 "cholangiocytes": ["Cholangiocytes"],
 "hepatic_epithelium": ["Hepatocytes", "Cholangiocytes"],
 "fibroblasts": ["Fibroblasts"],
 "endothelial": ["Endothelial cells"],
 "stromal_vascular": ["Fibroblasts", "Endothelial cells"],
 "macrophage_monocyte": ["Macrophages", "Mono+mono derived cells"],
 "other_myeloid": ["cDC1s", "cDC2s", "pDCs", "Neutrophils", "Basophils"],
 "all_myeloid": ["Macrophages", "Mono+mono derived cells", "cDC1s", "cDC2s", "pDCs", "Neutrophils", "Basophils"],
 "T_NK": ["T cells", "Circulating NK/NKT", "Resident NK"],
 "B_plasma": ["B cells", "Plasma cells"],
 "lymphoid": ["T cells", "Circulating NK/NKT", "Resident NK", "B cells", "Plasma cells"],
 "nonparenchymal_nonmyeloid": ["Fibroblasts", "Endothelial cells", "T cells", "Circulating NK/NKT", "Resident NK", "B cells", "Plasma cells"],
}
rows = []
for nm, fr in [("rebuilt", new), ("earlier", old)]:
    for c in ct:
        rows.append([nm, c] + [fr.loc[ds == g, c].quantile(0.05) for g in DS])
    for k, v in comp.items():
        s = fr[v].sum(1)
        rows.append([nm, "C:" + k] + [s[ds == g].quantile(0.05) for g in DS])
t = pd.DataFrame(rows, columns=["input", "part"] + DS)
t["min_p05"] = t[DS].min(1)
t["min_p05_excl_GSE135251"] = t[[g for g in DS if g != "GSE135251"]].min(1)
pd.set_option("display.width", 250)
fmt = lambda x: f"{x:.1e}"
print(t.to_string(index=False, float_format=fmt))
print("\nunits per cohort:", ds.value_counts().to_dict())
