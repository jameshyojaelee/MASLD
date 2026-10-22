# Audit H80: which tier-1/2 gene-state x summary-fit groups hold the 145 rds refit on 2026-08-15/16
# (run 2026-10-15 on the login node; output in refit_mtime_check.out)
import glob, os, datetime, pandas as pd
FM="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/eqtl_susie_polyfun"
W="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H80/gene_state_by_polyfun_fit.tsv"
rows=[(os.path.basename(p)[:-10], datetime.date.fromtimestamp(os.path.getmtime(p)).isoformat()) for p in glob.glob(FM+"/chr*/ENSG*_susie.rds")]
r=pd.DataFrame(rows,columns=["ensembl","mtime"]); late=r[r.mtime>="2026-08-15"]
g=pd.read_csv(W,sep="\t"); x=late.merge(g[["ensembl","susie_state_t12","fit"]],on="ensembl",how="left")
print("rds with mtime >= 2026-08-15:",len(late)); print(pd.crosstab(x.susie_state_t12.fillna("not_in_gene_table"),x.fit.fillna("NA")).to_string())
