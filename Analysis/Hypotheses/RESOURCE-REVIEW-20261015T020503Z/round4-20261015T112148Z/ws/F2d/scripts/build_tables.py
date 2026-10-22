import pandas as pd, sys
H='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H58-20261015T112423Z/jobs/5800/H58-20261015T113519Z_22252009/'
M='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/main/fig2_genetics/panels/Fig2D_corrected_coding_interval.tsv'
W,F=sys.argv[1],sys.argv[2]
a=pd.read_csv(M,sep='\t').iloc[0]
C=pd.read_csv(H+'family_C.tsv',sep='\t').set_index('member')
R=pd.read_csv(H+'descriptive_refits.tsv',sep='\t').set_index('label')
rows=[dict(row='Pooled, all ancestries (adopted)',strata='all',n_clusters_direct=19,n_clusters_enzyme=228,fraction_direct=a.fraction_direct,fraction_enzyme=a.fraction_enzyme,
  difference=a.difference,lo95=a.diff_boot_lo,hi95=a.diff_boot_hi,ratio=a.ratio,source='adopted fig2_family_multiplicity.tsv, seed 20260817',h58_call='not part of H58 family')]
lab={'C1':'EUR+AMR, fixed weights 16/26, 10/26','C2':'EUR only','C3':'AMR only'}
for m in ['C1','C2','C3']:
    r=C.loc[m]; rows.append(dict(row=lab[m],strata=r.strata,n_clusters_direct=None,n_clusters_enzyme=None,fraction_direct=r.fraction_direct,fraction_enzyme=r.fraction_enzyme,
     difference=r.difference,lo95=r.lo95,hi95=r.hi95,lo98_333=r.lo98_333,hi98_333=r.hi98_333,ratio=r.ratio,source='H58 family_C.tsv, seed 20260923, 10,000 draws',h58_call=r.call))
r=R.loc['nonEUR_pooled']; rows.append(dict(row='All non-EUR pooled',strata='AFR+AMR+EAS+SAS',n_clusters_direct=int(r.n_clusters_direct),n_clusters_enzyme=int(r.n_clusters_enzyme),
  fraction_direct=r.fraction_direct,fraction_enzyme=r.fraction_enzyme,difference=r.difference,lo95=r.diff_boot_lo,hi95=r.diff_boot_hi,ratio=r.ratio,source='H58 descriptive_refits.tsv',h58_call='descriptive, no call'))
cn={'C2':(16,14),'C3':(10,61)}
t=pd.DataFrame(rows); t.to_csv(F+'/Fig2D_coding_by_ancestry_source.tsv',sep='\t',index=False); print(t.to_string())
s=pd.read_csv(W+'/out/fig2d_trait_subclass_split.tsv',sep='\t'); s=s[s.enzyme_subset=='ENZ']
s['row']=s.direct_subset.map({'DIRECT':'All direct (NAFLD, NASH, PDFF)','PDFF':'PDFF only (3 EUR studies)','NAFLD_NASH':'NAFLD and NASH only'})
s.to_csv(F+'/Fig2D_direct_subclass_source.tsv',sep='\t',index=False); print(s.to_string())
x=pd.read_csv(H+'context_descriptive.tsv',sep='\t')
x=x[x.label.isin(['credible-set-weighted repaired','credible-set-weighted adopted'])]
x.to_csv(F+'/Fig2_companion_noncoding_context_repaired_source.tsv',sep='\t',index=False); print(x.to_string())
c=pd.read_csv(H+'FIG2A_cells.tsv',sep='\t'); c[['trait','ancestry']]=c['trait|ancestry'].str.split('|',expand=True)
c.to_csv(F+'/Fig2A_finemapped_support_by_cell_source.tsv',sep='\t',index=False)
