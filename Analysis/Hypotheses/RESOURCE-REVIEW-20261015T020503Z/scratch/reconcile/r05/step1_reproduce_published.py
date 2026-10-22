# Re-aggregate published node-level stage/continuum betas into 43 system values; compare with Fig 4H source.
import pandas as pd, numpy as np
from scipy.stats import spearmanr
ROOT='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
nodes=pd.read_csv(ROOT+'figures/candidates/continuum-774-panels-20261011T152106Z/molecular_systems_atlas/overlays/all_node_stage_continuum_overlays.tsv.gz',sep='\t')
print(nodes.shape, nodes.component_role.value_counts(dropna=False).to_dict(), 'NA community', nodes.community_id.isna().sum())
def agg(df,col):
    d=df[df.community_id.notna()]
    fam=d.groupby(['community_id','collection'])[col].median().dropna()
    return fam.groupby('community_id').mean()
s=agg(nodes,'stage_beta'); c=agg(nodes,'continuum_beta')
h=pd.read_csv(ROOT+'figures/main/fig4_singlecell_programs/source_tables/corrected_774/fig4h_single_outline/fig4h_source.tsv',sep='\t').set_index('community_id')
print('n systems',len(s),'max|diff| stage',(s-h.stage_beta).abs().max(),'max|diff| cont',(c-h.continuum_beta).abs().max())
print('published rho',spearmanr(h.stage_beta,h.continuum_beta)[0],'same dir',int((np.sign(h.stage_beta)==np.sign(h.continuum_beta)).sum()),'median stage',h.stage_beta.median())
