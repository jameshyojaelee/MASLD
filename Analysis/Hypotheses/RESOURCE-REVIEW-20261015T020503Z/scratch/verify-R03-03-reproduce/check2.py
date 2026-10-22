import pandas as pd
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
M=R+"RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/"
S=R+"figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/"
w=pd.read_csv(M+"lvqw_sample_weights.tsv",sep="\t")
f5=pd.read_csv(S+"five_cohort_sample_manifest.tsv",sep="\t")
x=w.merge(f5[["sample_id","dataset"]],on="sample_id",how="left")
print("weights rows",len(w),"unmatched",x.dataset.isna().sum())
s=x.groupby("dataset").sample_weight.agg(["count","sum"]); s["share"]=s["sum"]/s["sum"].sum(); print(s)
cols=["contrast","FDR","n_cohorts","n_reference","n_comparison","bh_family_size"]
g=pd.read_csv(S+"stage_extension_all_gene_results.tsv",sep="\t",usecols=cols)
t=g.groupby("contrast").agg(n=("FDR","size"),sig=("FDR",lambda v:(v<0.05).sum()),ncoh=("n_cohorts","first"),nref=("n_reference","first"),ncmp=("n_comparison","first"))
print(t)
