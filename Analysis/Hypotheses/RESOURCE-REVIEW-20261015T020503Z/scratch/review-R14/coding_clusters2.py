import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np
d=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/corrected_locus_coordinates.tsv',sep='\t')
d['coding']=d.protein_altering_pip_mass+d.canonical_splice_pip_mass
for sc,sub in d.groupby('trait_scope'):
    print(sc,'credible-set-weighted coding frac',round(sub.coding.mean()*100,3),'n',len(sub))
x=d[(d.trait_scope=='tier1_direct_masld_pdff')&(d.chrom.isin([11,1,14,2,19,4,17]))]
cols=['study','trait','chrom','pos','n_members','protein_altering_pip_mass','canonical_splice_pip_mass','coding']
for (c,lo,hi) in [(11,800000,850000),(1,220900000,221000000),(14,94800000,94900000),(2,27700000,27800000),(19,45400000,45420000),(4,100200000,100500000),(17,64200000,64250000)]:
    print(x[(x.chrom==c)&(x.pos>=lo)&(x.pos<=hi)][cols].to_string(index=False))
reg=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/snapshot/GWAS/finemapping/config/gwas_registry.tsv',sep='\t')
print(reg.columns.tolist())
