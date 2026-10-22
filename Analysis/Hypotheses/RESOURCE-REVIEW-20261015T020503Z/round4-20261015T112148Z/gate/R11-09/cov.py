import pandas as pd
C='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet'
d=pd.read_parquet(C)
print(d.columns.tolist()); print(len(d), d.symbol.nunique())
p=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv',sep='\t',keep_default_na=False,na_values=[''])
print('rows',len(p),'uniq',p.gene.nunique(),'multi-gene ids',p.gene.str.contains('[;,|]').sum())
ps=set(p.gene)
print('protein genes not in catalog symbols:',len(ps-set(d.symbol)), sorted(ps-set(d.symbol))[:30])
d['m']=d.symbol.isin(ps)
print(d.groupby('primary_evidence_class').m.agg(['size','sum','mean']))
known=['ALB','APOA1','APOB','CYP2E1','CYP3A4','FASN','PNPLA3','HSD17B13','TM6SF2','MTARC1','GCKR','PLIN2','COL1A1','ACTA2','LUM','APOE']
print({g:(g in ps) for g in known})
neg=['INS','GCG','MYH7','TNNT2','NEUROD1','SYP','OLIG2','TTN','KRT1','CD8A']
print({g:(g in ps) for g in neg})
