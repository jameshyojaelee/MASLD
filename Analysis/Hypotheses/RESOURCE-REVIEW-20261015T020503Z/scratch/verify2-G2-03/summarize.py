import pandas as pd, numpy as np, glob, os
J='Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/4211'
rel=pd.read_csv('figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/stage_all_gene_results.tsv',sep='\t',usecols=['gene_id_versioned','gene_id_base','gene_name','logFC','P.Value','FDR','contrast'])
print('released genome-wide BH<0.05 per contrast:'); print(rel.groupby('contrast').FDR.apply(lambda x:(x<0.05).sum()))
fits={os.path.basename(f).replace('.tsv.gz',''):pd.read_csv(f,sep='\t') for f in glob.glob(J+'/fits/*.tsv.gz')}
b=fits['base__F0_to_F1'].merge(rel[rel.contrast=='F0_to_F1'],on='gene_id_versioned',suffixes=('','_rel'))
print('base reproduction F0_to_F1: n',len(b),'max|dlogFC|',np.abs(b.logFC-b.logFC_rel).max(),'max|dlog10P|',np.abs(np.log10(b['P.Value'])-np.log10(b['P.Value_rel'])).max(),'sig agree',((b.FDR<0.05)==(b.FDR_rel<0.05)).all())
genes=pd.read_csv('GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv')
pos=set(genes.loc[genes.coloc_best_susie_pp4>0.5,'ensembl']); print('PP.H4>0.5 genes',len(pos))
rel['base']=rel.gene_id_versioned.str.split('.').str[0]
name=rel.drop_duplicates('base').set_index('base').gene_name
print('in stage table',len(pos & set(rel.base)))
# released panel
m=pd.read_csv('figures/candidates/fig3f-no-pip-774-20261012T011128Z/source_tables/fig3f_stage_genetics_matrix.tsv',sep='\t')
relmarks=set(zip(m.loc[m.FDR<0.05,'gene_id_base'],m.loc[m.FDR<0.05,'contrast']))
relgenes=set(m.gene_id_base)
f01=sorted(g for g,c in relmarks if c=='F0_to_F1'); print('released F0_to_F1 marks',len(f01))
for lab in ['no130970','lib','gini','nf','joint']:
    f=fits[lab+'__F0_to_F1']; f['base']=f.gene_id_versioned.str.split('.').str[0]
    sig=set(f.loc[f.FDR<0.05,'base'])
    surv=[g for g in f01 if g in sig]
    print(f'{lab}: F0_to_F1 marks surviving BH<0.05: {len(surv)} of {len(f01)}', [name[g] for g in surv] if len(surv)<=10 else '')
# joint reselection
J4=pd.concat([fits['joint__'+c].assign(base=lambda x:x.gene_id_versioned.str.split('.').str[0]) for c in ['F0_to_F1','F1_to_F2','F2_to_F3','F3_to_F4']])
J4=J4[J4.base.isin(pos)]
sel=set(J4.loc[J4.FDR<0.05,'base'])
jm=J4[J4.base.isin(sel)]
print('joint reselection: genes',len(sel),'effects',len(jm),'marks',(jm.FDR<0.05).sum(), jm[jm.FDR<0.05].groupby('contrast').size().to_dict())
stay=sel&relgenes; leave=relgenes-sel; enter=sel-relgenes
print('stay',len(stay),'leave',len(leave),'enter',len(enter))
print('leave:',sorted(name[g] for g in leave)); print('enter:',sorted(name[g] for g in enter))
jmarks=set(zip(jm.loc[jm.FDR<0.05,'base'],jm.loc[jm.FDR<0.05,'contrast']))
print('released marks staying under joint',len(relmarks&jmarks),'of',len(relmarks))
print('genome-wide joint BH<0.05:',{k:(v.FDR<0.05).sum() for k,v in fits.items()})
# per-cohort composition of F0_to_F1 contrast
r=pd.read_csv(J+'/stage_roster_covariates.tsv',sep='\t')
print(r[r.fib_group.isin(['F0','F1'])].groupby(['dataset','fib_group']).size().unstack())
print(r.groupby(['dataset','fib_group'])[['c_lib','c_gini','c_nf']].mean().round(2).loc[(slice(None),['F0','F1']),:])
