import pandas as pd
S="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/candidates/panel-deep-fixes-20261017T035633Z/fig2/source_tables/"
for p in ["2H","2G"]:
    h=pd.read_csv(S+p+".tsv",sep="\t")
    print(p, h.track.value_counts().to_dict())
    for t in ["gwas","eqtl"]:
        d=h[h.track==t].sort_values("p").head(3)
        print(t, d[["position","p","value","mark","in_cs"]].to_dict("records"))
    pip=h[h.track.str.contains("pip",case=False)]
    print("pip rows",len(pip), pip.sort_values("value",ascending=False).head(4)[["position","value","in_cs","mark"]].to_dict("records"))
    print(h[h.track=="gene"][["label","lfc","padj","state"]].to_dict("records") if "gene" in set(h.track) else h[h.track.str.contains("gene")][["label","lfc","padj","state"]].drop_duplicates().to_dict("records"))
