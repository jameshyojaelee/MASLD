import pandas as pd
r = pd.read_csv("resource_identity_roster.tsv", sep="\t", keep_default_na=False, na_values=[""])
b = pd.read_csv("bulk_roster.tsv", sep="\t", keep_default_na=False, na_values=[""])
p = pd.read_csv("full_filtered_related_pairs.tsv", sep="\t")
p.columns = [c.lstrip("#") for c in p.columns]
ds = dict(zip(r.sample_id, r.dataset))
sel = set(b.loc[b.primary_selected.astype(str)=="True","sample_id"])
print("bulk primary selected:", len(sel))
p["d1"] = p.IID1.map(ds); p["d2"] = p.IID2.map(ds)
p["in1"] = p.IID1.isin(sel); p["in2"] = p.IID2.isin(sel)
p["cross"] = p.d1 != p.d2
print(p.groupby(["cross"]).size())
print(p[p.cross][["IID1","IID2","d1","d2","KINSHIP","NSNP","in1","in2"]].to_string())
g = p[(p.d1=="GSE213621")|(p.d2=="GSE213621")]
print("GSE213621 pairs:", len(g)); print(g[["IID1","IID2","d1","d2","KINSHIP","in1","in2"]].to_string())
print(r[r.dataset=="GSE213621"].groupby(["passes_min_sites","identity_qualified"]).size())
print(r.groupby("dataset").autosomal_called.describe())
