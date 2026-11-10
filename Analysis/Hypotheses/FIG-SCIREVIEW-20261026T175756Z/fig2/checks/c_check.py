# Fig 2C: recount classes from gene-level table; re-check lead rule against COLOC r2 gene table
import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
c=pd.read_csv(R+"figures/candidates/panel-deep-fixes-20261017T035633Z/fig2/source_tables/2C_gene_level.tsv",sep="\t")
print("genes",c.ensembl.nunique(),"distinct leads",c.top_snp.nunique())
print(c.class_panel.value_counts().to_dict())
print(c[c.symbol.isin(['FABP1','RORA','SAMM50','MAU2','PNPLA3','TM6SF2','HSD17B13','GCKR'])][['symbol','gwas_name','top_snp','s4','a4','class_panel','coding_gene_symbols']].to_string())
g=pd.read_csv(R+"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
u=g[(g.coloc_best_susie_pp4>0.5)|(g.coloc_best_abf_pp4>0.5)]
print("union",len(u),"in 2C",u.ensembl.isin(c.ensembl).sum())
mx=np.fmax(c.s4.fillna(-1),c.a4.fillna(-1)); um=u.set_index('ensembl')
best=np.fmax(um.coloc_best_susie_pp4.fillna(-1),um.coloc_best_abf_pp4.fillna(-1))
cc=c.set_index('ensembl'); print("lead PP.H4 equals gene max:", np.isclose(np.fmax(cc.s4.fillna(-1),cc.a4.fillna(-1)), best.reindex(cc.index)).mean())
