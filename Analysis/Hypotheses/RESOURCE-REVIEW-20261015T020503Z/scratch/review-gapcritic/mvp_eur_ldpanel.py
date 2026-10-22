import pandas as pd
RUN='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/'
m=pd.read_pickle(RUN+'scratch/review-R02/r2_main.pkl')
p=m[m['PP.H4.susie']>0.5].copy()
p['ldp']=p.ld_panel.fillna('<NA>').replace('', '<empty>')
print(p.groupby(['ancestry','ldp']).agg(rows=('ensembl','size'),genes=('ensembl','nunique')).to_string())
mvp=p.gwas_name.str.startswith('MVP_')&p.gwas_name.str.endswith('_EUR')
g=p.groupby('ensembl').apply(lambda d:(d.gwas_name.str.startswith('MVP_')&d.gwas_name.str.endswith('_EUR')).all())
print('genes positive ONLY in MVP EUR strata:',int(g.sum()),'of',p.ensembl.nunique())
fin=m[m['PP.H4.susie'].notna()]
fin['ldp']=fin.ld_panel.fillna('<NA>').replace('','<empty>')
print('finite SuSiE rows by ancestry x ld_panel'); print(fin.groupby(['ancestry','ldp']).size().to_string())
print('MVP EUR finite SuSiE rows by study:',fin[fin.gwas_name.str.startswith('MVP_')&fin.gwas_name.str.endswith('_EUR')].groupby(['gwas_name','ldp']).size().to_dict())
