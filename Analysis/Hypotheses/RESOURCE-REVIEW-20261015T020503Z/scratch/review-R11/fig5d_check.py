import pandas as pd
P='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv'
d=pd.read_csv(P,sep='\t',keep_default_na=False,na_values=[''])
pairs=d.drop_duplicates(['gwas_name','gene'])
print('pairs',len(pairs),'genes',pairs.gene.nunique(), pairs.trait_class.value_counts().to_dict())
print('min mapped posterior', d.mapped_posterior_mass.min())
# per-pair, max over lineages of any_accessible_mass ; how many pairs have <1% in every lineage
m=d.groupby(['gwas_name','gene','trait_class']).any_accessible_mass.max().reset_index()
for thr in [0.001,0.01,0.05,0.5]:
    print(f'pairs with max-over-lineages accessible mass < {thr}:', (m.any_accessible_mass<thr).sum(),'/',len(m), m[m.any_accessible_mass<thr].trait_class.value_counts().to_dict())
