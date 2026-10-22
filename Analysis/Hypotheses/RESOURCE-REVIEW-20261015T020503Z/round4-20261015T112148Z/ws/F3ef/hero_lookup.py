import pandas as pd, glob, os
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"; RUN=ROOT+"/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z"
OUT=ROOT+"/figures/candidates/review-r3-fixes-20261015T112148Z/F3ef"
ext=pd.read_csv(ROOT+"/figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv",sep="\t",usecols=["gene_id_versioned","gene_name","logFC","FDR","contrast","axis"])
H=ext[ext.gene_name.isin(["GNMT","MAT1A","CYP2C19"])&ext.axis.isin(["fibrosis","fibrosis_adjacent"])]
ids=dict(zip(H.gene_id_versioned,H.gene_name))
rows=[dict(gene=r.gene_name,fit="released",contrast=r.contrast,logFC=r.logFC,q=r.FDR) for r in H.itertuples()]
def add(fit,con,p,sep="\t"):
    g=pd.read_csv(p,sep=sep); g=g[g.gene_id_versioned.isin(ids)]
    for r in g.itertuples(): rows.append(dict(gene=ids[r.gene_id_versioned],fit=fit,contrast=con,logFC=r.logFC,q=r.FDR))
for f in glob.glob(RUN+"/jobs/4211/fits/*.tsv.gz"):
    lab,con=os.path.basename(f).replace(".tsv.gz","").split("__"); add("v4211_"+lab,con,f)
for f in glob.glob(RUN+"/jobs/1721/fits/*.genes.tsv.gz"):
    lab,con=os.path.basename(f).replace(".genes.tsv.gz","").split("__"); add("r1721_"+lab,con,f)
a=pd.read_csv(ROOT+"/figures/candidates/f0-reference-arms-r8-20261008T193309Z/analysis/f0_arms/f0_arm_results.tsv.gz",sep="\t")
a=a[a.gene_id_versioned.isin(ids)]
for r in a.itertuples(): rows.append(dict(gene=ids[r.gene_id_versioned],fit="arm_"+r.arm,contrast=r.contrast,logFC=r.logFC,q=r.FDR))
d=pd.DataFrame(rows); d.to_csv(OUT+"/hero_genes_stage_sensitivity.tsv",sep="\t",index=False)
pd.set_option("display.width",250)
for c in ["F0_to_F1","F1_vs_F0","F1_to_F2","F2_vs_F0","F3_vs_F0","F4_vs_F0"]:
    x=d[d.contrast==c].pivot_table(index="fit",columns="gene",values="q"); print(c); print(x.applymap(lambda v:f"{v:.3g}").to_string())
