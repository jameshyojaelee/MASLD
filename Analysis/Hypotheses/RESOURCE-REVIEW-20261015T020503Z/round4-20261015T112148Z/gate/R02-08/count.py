import pandas as pd, numpy as np
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv'
cols=['gwas_name','gene','ensembl','chr','PP.H4.abf','PP.H4.susie','n_snps','n_palindromic_kept','n_palindromic_af_tested','n_palindromic_af_dropped','method','top_snp','ancestry']
m=pd.read_csv(R,usecols=cols)
print('all rows',len(m),'strata',m.gwas_name.nunique())
p=m[m['PP.H4.susie']>0.5]
print('susie-positive rows',len(p),'genes',p.gene.nunique())
u=p[p.gwas_name.isin(['UKBB_ALT','UKBB_AST','UKBB_GGT'])]
print('UKBB positive rows',len(u),'genes',u.gene.nunique())
ug=set(u.gene); og=set(p[~p.gwas_name.isin(['UKBB_ALT','UKBB_AST','UKBB_GGT'])].gene)
print('UKBB-positive genes also positive in another stratum',len(ug&og),'UKBB-only genes',len(ug-og))
g=p.groupby('gwas_name')[['n_snps','n_palindromic_kept','n_palindromic_af_tested','n_palindromic_af_dropped']].sum()
g['rows']=p.groupby('gwas_name').size()
g['drop_frac_of_tested']=g.n_palindromic_af_dropped/g.n_palindromic_af_tested
print(g.sort_values('rows',ascending=False).to_string())
eur=g[g.index.str.contains('EUR')|g.index.str.startswith('FinnGen')]
print('EUR AF-tested pooled dropped/tested', eur.n_palindromic_af_dropped.sum(), eur.n_palindromic_af_tested.sum())
a=m.groupby('gwas_name')[['n_palindromic_af_tested','n_palindromic_af_dropped']].sum()
a['f']=a.n_palindromic_af_dropped/a.n_palindromic_af_tested
print('all rows, AF drop fraction by stratum'); print(a[a.n_palindromic_af_tested>0].sort_values('f').to_string())
for gene in ['PNPLA3','TM6SF2','GCKR','HSD17B13','MBOAT7','APOE']:
    print(m[(m.gene==gene)&m.gwas_name.str.startswith('UKBB_')][['gwas_name','gene','PP.H4.abf','PP.H4.susie','top_snp','n_palindromic_kept']].to_string())
u.to_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/R02-08/ukbb_positive_rows.tsv',sep='\t',index=False)
