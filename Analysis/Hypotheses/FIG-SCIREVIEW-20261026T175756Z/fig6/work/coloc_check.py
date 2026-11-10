import pandas as pd, numpy as np
C='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/'
t=pd.read_csv(C+'tier12/gene_level_coloc_tier12.csv')
print('tier12 gene rows',len(t)); print(t.susie_state_t12.value_counts())
print('SuSiE>0.5',(t.coloc_best_susie_pp4>0.5).sum())
tiers=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t')
t12=set(tiers[tiers.tier.isin([1,2])].study_name); print('n tier1/2 studies',len(t12), (tiers.tier==1).sum(), (tiers.tier==2).sum())
genes={'ENSG00000112164':'GLP1R','ENSG00000156510':'HKDC1','ENSG00000186205':'MTARC1','ENSG00000151090':'THRB'}
rows=[]
for ch in pd.read_csv(C+'susie_coloc_all_gwas.csv',usecols=['gwas_name','gene','ensembl','PP.H4.abf','PP.H4.susie','PP.H3.susie','method','n_cs_pairs_eligible'],chunksize=200000):
    rows.append(ch[ch.ensembl.isin(genes)])
a=pd.concat(rows); a=a[a.gwas_name.isin(t12)]
for e,s in genes.items():
    x=a[a.ensembl==e]
    fin=x['PP.H4.susie'].notna()
    print(f'\n{s}: tier1/2 records {len(x)}, finite SuSiE {fin.sum()}, SuSiE>0.5 {(x["PP.H4.susie"]>0.5).sum()}')
    if fin.any():
        b=x.loc[x['PP.H4.susie'].idxmax()]; print('  best SuSiE',b['PP.H4.susie'],b.gwas_name)
    b=x.loc[x['PP.H4.abf'].idxmax()]; print('  best ABF',b['PP.H4.abf'],b.gwas_name)
    print(x[fin][['gwas_name','PP.H4.susie','PP.H3.susie','PP.H4.abf']].to_string())
    print('  methods',x.method.value_counts().to_dict())
print(t[t.ensembl.isin(genes)].to_string())
