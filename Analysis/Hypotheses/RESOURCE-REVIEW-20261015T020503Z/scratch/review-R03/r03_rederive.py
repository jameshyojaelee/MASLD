# Review R03 re-derivations (read-only; small tables). Run from repo root.
import pandas as pd, numpy as np
from scipy.stats import fisher_exact, spearmanr, mannwhitneyu
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
print('== 3C per-cohort calls, union, tiers, direction')
d=pd.read_csv(R+'figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv',sep='\t')
c=d[(d.FDR<0.05)&(d.logFC.abs()>0.5)]
print(c.groupby('dataset').size().to_dict(), 'memberships',len(c))
n=c.groupby('gene_id_versioned').dataset.nunique(); conc=c.groupby('gene_id_versioned').logFC.apply(lambda x:(x>0).all() or (x<0).all())
print('union',len(n), n.value_counts().sort_index().to_dict())
print('discordant by tier', (~conc).groupby(n).sum().to_dict(), 'multi',int((n>=2).sum()),'discordant multi',int((~conc[n>=2]).sum()))
disc=conc.index[~conc]; inv=c[c.gene_id_versioned.isin(disc)].groupby('gene_id_versioned').dataset.apply(lambda x:'GSE126848' in set(x))
print('discordant involving GSE126848',int(inv.sum()),'| GSE126848 unique',int(((n==1)&c.groupby('gene_id_versioned').dataset.first().eq('GSE126848')).sum()))
print('== Identity screen: depth exclusion by cohort and group')
r=pd.read_csv(R+'figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/all_screened_roster.tsv',sep='\t',keep_default_na=False,na_values=[''])
r=r[r.dataset.isin(['GSE126848','GSE130970','GSE135251','GSE162694','GSE213621'])]
r['depth_excl']=r.exclusion_reason.eq('insufficient_filtered_identity_sites')
print(pd.crosstab([r.dataset,r.group_binary],r.depth_excl).to_string())
g=r[r.dataset=='GSE126848']; print('GSE126848 Fisher', fisher_exact(pd.crosstab(g.group_binary,g.depth_excl).values))
print('GSE126848 full_called median ctrl/dis', g.groupby('group_binary').full_called.median().to_dict(), 'MWU p', mannwhitneyu(g[g.group_binary=='Control'].full_called,g[g.group_binary=='Disease'].full_called).pvalue)
k=g[g.group_binary=='Disease']; print('GSE126848 disease retained by condition', pd.crosstab(k.condition,~k.depth_excl).to_dict())
print('== 3F dependence on mixed F0')
f=pd.read_csv(R+'figures/candidates/f0-reference-arms-r8-20261008T193309Z/analysis/f0_arms/f0_arm_results.tsv.gz',sep='\t')
m=pd.read_csv(R+'figures/main/fig3_bulk_transcriptomics/source_tables/fig3f_stage_genetics_matrix.tsv',sep='\t')
w=m.pivot_table(index='gene_id_versioned',columns='contrast',values='FDR')
b=f[(f.arm=='B_disease_only_F0')&(f.contrast=='F0_to_F1')].set_index('gene_id_versioned').FDR
w['B']=w.index.map(b)
print('genes',len(w),'F0->F1 marks',int((w.F0_to_F1<0.05).sum()),'kept disease-only F0',int(((w.F0_to_F1<0.05)&(w.B<0.05)).sum()))
only=(w.F0_to_F1<0.05)&(w[['F1_to_F2','F2_to_F3','F3_to_F4']].min(axis=1)>=0.05)
print('rows entering only via F0->F1',int(only.sum()),'and lost with disease-only F0',int((only&(w.B>=0.05)).sum()), 'q range', w.B[only&(w.B>=0.05)].min().round(4), w.B[only&(w.B>=0.05)].max().round(4))
print('== 3G')
t=pd.read_csv(R+'figures/main/fig3_bulk_transcriptomics/source_tables/fig3g_raw_treat_interface.tsv',sep='\t')
print(len(t), t.state.value_counts().to_dict(), 'signed', spearmanr(t.logFC,t.susie_pp4), 'abs', spearmanr(t.logFC.abs(),t.susie_pp4))
print('== Stage roster composition')
s=pd.read_csv(R+'figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/stage_extension_sample_manifest.tsv',sep='\t')
fm=pd.read_csv(R+'figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/five_cohort_sample_manifest.tsv',sep='\t',keep_default_na=False)
x=s.merge(fm[['sample_id','group_binary']],on='sample_id')
for ax in ['fibrosis','NAS']:
    y=x[x.axis==ax].drop_duplicates('sample_id'); print(ax, pd.crosstab(y.group,y.group_binary).to_dict())
print('stage cohorts', sorted(s.dataset.unique()))
