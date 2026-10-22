import pandas as pd
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/'
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
s=pd.read_csv(R+'susie_coloc_signal_pairs.tsv.gz',sep='\t')
print(s.shape, s.eligible.value_counts().to_dict(), s.gene.nunique())
el=s[s.eligible==True]
gm=el.groupby('ensembl').agg(h4=('PP.H4','max'),h3=('PP.H3','max'))
print('genes with eligible pairs',len(gm),'maxH4>0.5',(gm.h4>0.5).sum(),'maxH3>0.8 & maxH4<=0.5',((gm.h3>0.8)&(gm.h4<=0.5)).sum(),'maxH3>0.9 & maxH4<0.1',((gm.h3>0.9)&(gm.h4<0.1)).sum())
e=pd.read_parquet(D+'passport_evidence_long.parquet'); ge=e[e.evidence_domain=='genetics'].set_index('ensembl_id')
gm.index=gm.index.str.split('.').str[0]
j=gm.join(ge[['call_state','symbol']],how='left')
x=j[(j.h3>0.8)&(j.h4<=0.5)]
print('H3>0.8 & H4<=0.5 genes by Catalog genetics call', x.call_state.value_counts(dropna=False).to_dict())
print(x.sort_values('h3',ascending=False).head(10))
print('gwas in pairs', s.gwas.nunique())
