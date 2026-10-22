# Light read: per-gene genetic observability states from COLOC r2 ABF posteriors (tier 1-2 main GWAS)
import pandas as pd, numpy as np
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
tier=pd.read_csv(R+'GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t')
main=set(tier.loc[tier.placement=='main','study_name'])
cols=['gwas_name','gene','ensembl','PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf','PP.H4.susie','method']
d=pd.read_csv(R+'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv',usecols=cols)
print('rows',len(d),'gwas',d.gwas_name.nunique(),'genes',d.ensembl.nunique())
d=d[d.gwas_name.isin(main)]
print('main rows',len(d),'main gwas',d.gwas_name.nunique(),'genes',d.ensembl.nunique())
d['gwas_sig']=d['PP.H1.abf']+d['PP.H3.abf']+d['PP.H4.abf']
d['eqtl_sig']=d['PP.H2.abf']+d['PP.H3.abf']+d['PP.H4.abf']
g=d.groupby('ensembl').agg(gene=('gene','first'),max_h4=('PP.H4.abf','max'),max_h4s=('PP.H4.susie','max'),
  max_gwas=('gwas_sig','max'),max_eqtl=('eqtl_sig','max'),max_h1=('PP.H1.abf','max'),max_h3=('PP.H3.abf','max'))
# state rule (prespecified here, descriptive): 
def st(r):
    if (r.max_h4s>0.5) : return 'A_shared_susie'
    if r.max_h4>0.5: return 'B_shared_abf_only'
    if r.max_gwas>=0.8 and r.max_eqtl>=0.8: return 'C_gwas_and_eqtl_distinct_or_unresolved'
    if r.max_gwas>=0.8: return 'D_gwas_signal_no_eqtl'
    if r.max_eqtl>=0.8: return 'E_eqtl_only_no_gwas'
    return 'F_neither_signal'
g['state']=g.apply(st,axis=1)
print(g.state.value_counts().sort_index())
for s in ['PNPLA3','TM6SF2','MTARC1','HSD17B13','MBOAT7','APOB','GCKR','TRIB1','MTTP','ERLIN1','SERPINA1','PCSK7','ADH1B','GPAM','THRB','HKDC1','TMC4','GLP1R']:
    x=g[g.gene==s]
    print(s, x[['max_h4','max_h4s','max_gwas','max_eqtl','max_h1','max_h3','state']].round(3).to_dict('records'))
g.to_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L03/c1_gene_states_main_tier.tsv',sep='\t')
