import pandas as pd
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/'
tier=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t')
t=dict(zip(tier.study_name,tier.tier))
sp=pd.read_csv(R+'susie_coloc_signal_pairs.tsv.gz',sep='\t')
print(sp.shape, sp.eligible.value_counts().to_dict())
sp['tier']=sp.gwas.map(t)
sp=sp[sp.tier.isin([1,2])]
sp=sp[sp.eligible.astype(str).str.lower().isin(['true','1'])] if sp.eligible.dtype!=bool else sp[sp.eligible]
sp['hit']=sp['PP.H4']>0.5
g=sp.groupby('gwas').agg(tier=('tier','first'),pairs=('hit','size'),hitpairs=('hit','sum'),genes_hit=('ensembl',lambda s: s[sp.loc[s.index,'hit']].nunique()))
loc=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/corrected_locus_coordinates.tsv',sep='\t')
cs=loc.groupby('study').size().rename('gated_cs')
g=g.join(cs,how='outer'); g['tier']=g.index.map(t)
g=g[g.tier.isin([1,2])].fillna(0)
print(g.sort_values(['tier','gated_cs'],ascending=[True,False]).to_string())
for k,d in g.groupby('tier'):
    print('tier',k,'studies',len(d),'gated_cs',int(d.gated_cs.sum()),'studies_with_pairs',int((d.pairs>0).sum()),'eligible_pairs',int(d.pairs.sum()),'hit_pairs',int(d.hitpairs.sum()))
gl=pd.read_csv(R+'gene_level_coloc_tier12.csv')
h=gl[gl.coloc_best_susie_pp4>0.5]
print('genes susie>0.5',len(h)); print(h.driving_tier.value_counts().to_dict())
# genes with any tier1 hit pair
hit=sp[sp.hit]
g1=set(hit[hit.tier==1].ensembl); g2=set(hit[hit.tier==2].ensembl)
print('genes with tier1 hit',len(g1),'tier2 hit',len(g2),'both',len(g1&g2),'union',len(g1|g2))
from scipy.stats import fisher_exact
for k,d in sp.groupby('tier'):
    gg=d.groupby('ensembl').hit.any()
    print('tier',k,'genes tested',len(gg),'genes hit',int(gg.sum()),'frac',round(gg.mean(),4),'pair frac',round(d.hit.mean(),4))
print('fisher pairs (nominal, pairs not independent)',fisher_exact([[30,171-30],[810,3844-810]]))
