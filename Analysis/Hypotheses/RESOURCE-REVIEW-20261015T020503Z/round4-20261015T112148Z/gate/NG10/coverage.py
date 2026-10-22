import pandas as pd, numpy as np
from scipy.stats import binomtest
C='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
K='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/public_r4/biorxiv_2026.08.26.747338/media-1.xlsx'
E='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/elison2025_sc_chromatin/elison2025_sc_chromatin_calls.csv'
g=pd.read_parquet(C+'passport_gene_index.parquet')
e=pd.read_parquet(C+'passport_evidence_long.parquet')
bulk=e[e.evidence_domain=='transcriptomics'][['ensembl_id','estimate','q_value','call_state','testability_state']]
print('catalog classes', g.primary_evidence_class.value_counts().to_dict())
sym=g.groupby('symbol').ensembl_id.nunique(); print('symbols mapping >1 ensembl', (sym>1).sum())
g1=g.drop_duplicates('symbol',keep=False)
for s in ['NAS','Fibrosis']:
    k=pd.read_csv(f'kumak_{s}.tsv',sep='\t')
    k['dir']=np.where(k.total_upregulated>k.total_downregulated,'up',np.where(k.total_upregulated<k.total_downregulated,'down','tie'))
    m=k.merge(g1,left_on='Gene',right_on='symbol',how='left')
    print(f'\n{s}: n={len(k)} unique={k.Gene.nunique()} matched={m.ensembl_id.notna().sum()} dir={k.dir.value_counts().to_dict()}')
    print(m.primary_evidence_class.value_counts(dropna=False).to_dict())
    mb=m.merge(bulk,on='ensembl_id',how='inner')
    t=mb[(mb.dir!='tie')&mb.estimate.notna()]
    agree=((t.dir=='up')&(t.estimate>0))|((t.dir=='down')&(t.estimate<0))
    print('bulk tested with direction',len(t),'sign agree',agree.sum(), binomtest(int(agree.sum()),len(t),0.5,alternative='greater').pvalue)
    sig=t[t.q_value<0.05]; ag=((sig.dir=='up')&(sig.estimate>0))|((sig.dir=='down')&(sig.estimate<0))
    print('q<0.05 subset',len(sig),'agree',ag.sum())
    for gene in ['COL1A1','COL3A1','THBS2','FBLN5','PDGFA','LAMA2','LAMA3','GNMT','MAT1A','CYP2C19']:
        r=m[m.Gene==gene]
        if len(r): print(' ',gene,r.dir.iloc[0],r.primary_evidence_class.iloc[0])
el=pd.read_csv(E); m=el.merge(g1,left_on='gene',right_on='symbol',how='left')
print('\nElison n',len(el),'matched',m.ensembl_id.notna().sum()); print(m[['gene','elison_finding_type','primary_evidence_class']].to_string())
