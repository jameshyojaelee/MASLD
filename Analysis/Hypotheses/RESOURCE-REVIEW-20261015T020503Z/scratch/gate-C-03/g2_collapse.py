import pandas as pd, numpy as np
D='Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-obs-20260926T155651Z/coloc_direction/full-20260930T120450Z/direction_by_trait.tsv'
t=pd.read_csv(D,sep='\t',low_memory=False,keep_default_na=False,na_values=[''])
r=t[t.release=='candidate_r2_20260924T053244Z'].copy()
r['dirrow']=r.direction_sign.isin([1,-1])
def isindel(a,b): return (a.astype(str).str.len()>1)|(b.astype(str).str.len()>1)
r['indel1']=isindel(r.gwas_allele1,r.gwas_allele2); r['indel2']=isindel(r.eqtl_ea_hit2,r.eqtl_nea_hit2)
r['indel_any']=r.indel1|r.indel2
d=r[r.dirrow]
print('directional rows',len(d),'genes',d.ensembl.nunique())
print('directional rows with indel GWAS lead',int(d.indel1.sum()),'indel eQTL lead',int(d.indel2.sum()),'either',int(d.indel_any.sum()),'genes either',d[d.indel_any].ensembl.nunique())
def collapse(g, gtex):
    dd=g[g.dirrow]
    if len(dd)==0: return 'untestable'
    if gtex and (dd.check_gtex_liver_same_sign_hit2==0).any(): return 'ambiguous'
    s=dd.direction_sign.unique()
    return ('call_up' if s[0]>0 else 'call_down') if len(s)==1 else 'ambiguous'
a=r.groupby('ensembl').apply(lambda g: collapse(g,False))
b=r.groupby('ensembl').apply(lambda g: collapse(g,True))
print('no GTEx rule',a.value_counts().to_dict())
print('with GTEx rule',b.value_counts().to_dict())
# genes whose call rests on indel-lead rows
calls=b[b.isin(['call_up','call_down'])].index
gi=d[d.ensembl.isin(calls)].groupby('ensembl').indel_any.agg(['any','all'])
print('called genes with >=1 indel directional row',int(gi['any'].sum()),'called genes with all directional rows indel',int(gi['all'].sum()))
ggt=d[d.ensembl.isin(calls)].groupby('ensembl').trait.apply(lambda s: set(s)=={'GGT'})
print('called genes resting on GGT only',int(ggt.sum()))
enz=d[d.ensembl.isin(calls)].groupby('ensembl').tier_label.apply(lambda s: set(s)=={'liver_enzyme'})
print('called genes resting only on enzyme traits',int(enz.sum()),'of',len(calls))
b.rename('dir').to_csv('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-C-03/gene_dir_h26rule.tsv',sep='\t')
gi.to_csv('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-C-03/called_gene_indel.tsv',sep='\t')
