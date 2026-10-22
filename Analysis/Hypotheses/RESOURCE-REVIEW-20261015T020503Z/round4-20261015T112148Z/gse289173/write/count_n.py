# Count GSE289173 donors by group, sex and control series from the GEO series matrix only (no expression read).
import gzip, pandas as pd, collections
f='Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gse289173/gate/geo_meta/GSE289173_series_matrix.txt.gz'
rows={}
for l in gzip.open(f,'rt'):
    if l.startswith('!Sample_title') or l.startswith('!Sample_geo_accession'):
        rows[l.split('\t')[0]]=[x.strip().strip('"') for x in l.rstrip('\n').split('\t')[1:]]
    if l.startswith('!Sample_characteristics_ch1'):
        v=[x.strip().strip('"') for x in l.rstrip('\n').split('\t')[1:]]
        k=v[0].split(':')[0]; rows[k]=[x.split(': ',1)[1] for x in v]
d=pd.DataFrame({'gsm':rows['!Sample_geo_accession'],'title':rows['!Sample_title'],'group':rows['disease group'],'sex':rows['Sex']})
d['series']=d.title.str.startswith('ctrl_').map({True:'ctrl_xx',False:'FL'})
d['stage']=d.group.map({'ctrl':0,'MASL':1,'eMASH':2,'aMASH':2})
print('donors',len(d),'unique titles',d.title.nunique(),'unique gsm',d.gsm.nunique())
print(pd.crosstab(d.group,[d.series]).to_string())
print(pd.crosstab(d.group,d.sex).to_string())
print(pd.crosstab(d.stage,d.series).to_string())
d.to_csv('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gse289173/write/donor_design.tsv',sep='\t',index=False)
