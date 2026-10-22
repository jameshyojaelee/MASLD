# Re-aggregate job 2171 feature-level stage fits to 43 systems; compare with Fig 4H source continuum axis.
import pandas as pd, numpy as np
from scipy.stats import spearmanr
ROOT='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
J='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/2171/out/'
res=pd.read_csv(J+'feature_stage_meta_three_rosters.tsv.gz',sep='\t')
nodes=pd.read_csv(ROOT+'figures/candidates/continuum-774-panels-20261011T152106Z/molecular_systems_atlas/overlays/all_node_stage_continuum_overlays.tsv.gz',sep='\t')
h=pd.read_csv(ROOT+'figures/main/fig4_singlecell_programs/source_tables/corrected_774/fig4h_single_outline/fig4h_source.tsv',sep='\t').set_index('community_id')
pub=nodes[['feature_id','collection','community_id']]
def agg(df,col):
    d=df[df.community_id.notna()]
    return d.groupby(['community_id','collection'])[col].median().dropna().groupby('community_id').mean()
for roster in ['published','disease_only','disease_only_restd']:
    r=res[res.roster==roster].merge(pub,on=['feature_id','collection'])
    s=agg(r,'beta').loc[h.index]; c=h.continuum_beta
    print(roster,'n_units',int(r.n_GSE162694.iloc[0]),int(r.n_GSE213621.iloc[0]),'rho %.4f'%spearmanr(s,c)[0],'same_dir',int((np.sign(s)==np.sign(c)).sum()),'of',len(s),
          'max|diff| vs panel stage_beta' if roster=='published' else '', (s-h.stage_beta).abs().max() if roster=='published' else '')
