import pandas as pd
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
D=R+'/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L08'
t12=pd.read_csv(R+'/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv')
gen=set(t12.index[t12.susie_state_t12=='positive'].astype(str))
print('primary SuSiE genes',len(gen), sorted(gen)[:5])
deg=pd.read_csv(R+'/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv')
can=deg[(deg.padj<0.05)&(deg.logFC.abs()>0.5)]
up=set(can[can.logFC>0].symbol); dn=set(can[can.logFC<0].symbol)
tested=set(deg.symbol.dropna())
p=pd.read_csv(D+'/scdvp/scdvp_table3_healthy.tsv',sep='\t')
p['g']=p.protein.astype(str).str.split(';')
z=p.explode('g').drop_duplicates('g').set_index('g')
prot=set(z.index)
print('scDVP proteins',len(p),'symbols',len(prot),'in bulk tested universe',len(prot&tested))
for name,s in [('canonical_up',up),('canonical_down',dn),('genetic_436',gen),('bulk_tested_background',tested)]:
    zz=z.loc[sorted(s&prot)]
    print(name,'n_with_protein',len(zz),'zonated',(zz.zonated=='Yes').sum(),'portal',((zz.coefficient>0)&(zz.zonated=='Yes')).sum(),'central',((zz.coefficient<0)&(zz.zonated=='Yes')).sum(),'median_coef %.3f'%zz.coefficient.median())
