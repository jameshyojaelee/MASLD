import pandas as pd
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
new = pd.read_csv(ROOT + "Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz", sep="\t").set_index("sample_id")
man = pd.read_csv(ROOT + "Analysis/Deconvolution/corrected_774_20261009T141331Z/bulk/sample_manifest.tsv", sep="\t").set_index("sample_id").loc[new.index]
old = pd.concat([pd.read_csv(ROOT + f"Analysis/Deconvolution/results/{g}/{g}_bayesprism_proportions.tsv", sep="\t", index_col=0) for g in sorted(new.dataset.unique())]).loc[new.index]
for c in ["Fibroblasts", "Cholangiocytes", "B cells", "Neutrophils"]:
    t = pd.DataFrame({"rebuilt": (new[c] < 1e-6).groupby([man.dataset, man.group_binary]).sum(),
                      "earlier": (old[c] < 1e-6).groupby([man.dataset, man.group_binary]).sum(),
                      "n": man.groupby(["dataset", "group_binary"]).size()})
    print(c, "units below 1e-6 by cohort x group"); print(t.to_string()); print()
