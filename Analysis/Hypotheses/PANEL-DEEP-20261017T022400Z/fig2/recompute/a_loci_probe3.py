import pandas as pd, numpy as np
a = pd.read_csv("out_D/D_gated_instances_with_clusters.tsv.gz", sep="\t", dtype={"chr": str})
a["anc5"]=a.ancestry.replace({"CSA":"SAS"})
def ncl(df,w):
    n=0
    for c,g in df.groupby("chr"):
        p=np.sort(g.pos.unique()); n+=1+int(np.sum(np.diff(p)>w))
    return n
for w in [0, 5e5, 1e6]:
    t = {(tr,an): ncl(g,w) for (tr,an),g in a.groupby(["trait","anc5"])}
    tt = {tr: ncl(g,w) for tr,g in a.groupby("trait")}
    print(int(w), t); print(" per trait", tt, sum(tt.values()), "global", ncl(a,w))
# distinct study-locus per trait-ancestry
print(a.groupby(["trait","anc5"]).locus.nunique().to_dict())
print(a.groupby(["trait"]).locus.nunique().to_dict(), a.locus.nunique())
