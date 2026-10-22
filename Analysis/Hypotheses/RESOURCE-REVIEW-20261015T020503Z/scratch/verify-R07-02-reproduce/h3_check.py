import pandas as pd
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/'
C='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
s=pd.read_csv(R+'susie_coloc_signal_pairs.tsv.gz',sep='\t')
cmp_=pd.read_csv(R+'qc/eligibility_rule_comparison.tsv',sep='\t')
tier=cmp_[['gwas_name','tier']].drop_duplicates()
print('gwas tiers in comparison file:', tier.tier.value_counts().to_dict())
s=s.merge(tier,left_on='gwas',right_on='gwas_name',how='left')
print('signal pairs tier:', s.tier.value_counts(dropna=False).to_dict())
el=s[s.eligible==True]
def summarize(df,label):
    g=df.groupby(['gene','ensembl']).agg(maxH3=('PP.H3','max'),maxH4=('PP.H4','max')).reset_index()
    a=g[(g.maxH3>0.8)&(g.maxH4<=0.5)]
    b=g[(g.maxH3>0.9)&(g.maxH4<0.1)]
    print(label,'genes:',g.shape[0],' H3>0.8&maxH4<=0.5:',a.shape[0],' H3>0.9&maxH4<0.1:',b.shape[0])
    return g,a,b
# variant 1: per gene over all eligible pairs across all studies
g1,a1,b1=summarize(el,'V1 all studies eligible pairs')
# variant 2: tier 1/2 only
g2,a2,b2=summarize(el[el.tier.isin([1,2])],'V2 tier12 eligible pairs')
# variant 3: per gene-study (any study with H3>0.8 and maxH4 in that study <=0.5) -> count genes
gs=el.groupby(['gene','ensembl','gwas']).agg(maxH3=('PP.H3','max'),maxH4=('PP.H4','max')).reset_index()
v3=gs[(gs.maxH3>0.8)&(gs.maxH4<=0.5)]
print('V3 gene-study level genes:',v3.gene.nunique())
e=pd.read_parquet(C+'passport_evidence_long.parquet')
gen=e[e.evidence_domain=='genetics'][['ensembl_id','symbol','call_state','testability_state','estimate']]
for lab,a in [('V1a',a1),('V1b',b1),('V2a',a2),('V2b',b2)]:
    m=a.merge(gen,left_on='ensembl',right_on='ensembl_id',how='left')
    print(lab, m.call_state.value_counts(dropna=False).to_dict())
m=a1.merge(gen,left_on='ensembl',right_on='ensembl_id',how='left')
print(m[['gene','maxH3','maxH4','call_state','testability_state','estimate']].sort_values('maxH3',ascending=False).to_string())
mb=b1.merge(gen,left_on='ensembl',right_on='ensembl_id',how='left')
print('V1b genes:',sorted(mb.gene.tolist()))
# indeterminate genetics overall and how many have no SuSiE record
print('genetics indeterminate total:', (gen.call_state=='indeterminate').sum())
