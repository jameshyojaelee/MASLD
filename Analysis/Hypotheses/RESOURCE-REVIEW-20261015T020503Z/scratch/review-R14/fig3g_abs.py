import pandas as pd, numpy as np
from scipy.stats import spearmanr, mannwhitneyu
d=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/main/fig3_bulk_transcriptomics/source_tables/fig3g_raw_treat_interface.tsv',sep='\t')
d=d[np.isfinite(d.logFC)&np.isfinite(d.susie_pp4)]
print('n',len(d), d.state.value_counts().to_dict())
print('signed rho',spearmanr(d.susie_pp4,d.logFC))
print('abs rho',spearmanr(d.susie_pp4,d.logFC.abs()))
hi=d.susie_pp4>0.5
print('median |logFC| pp4>0.5',d.logFC[hi].abs().median(),'n',hi.sum(),' pp4<=0.5',d.logFC[~hi].abs().median(),'n',(~hi).sum())
print('MWU |logFC|',mannwhitneyu(d.logFC[hi].abs(),d.logFC[~hi].abs()))
print('up/down among pp4>0.5',(d.logFC[hi]>0).sum(),(d.logFC[hi]<0).sum())
