import pandas as pd, numpy as np
e=pd.read_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L10-04/eur_tier12_pairs.tsv",sep="\t"); p=pd.read_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L10-04/platelet_eur_pairs.tsv",sep="\t")
print("liver arm pairs",len(e),"studies",e.gwas.nunique()); print(e.groupby('gwas').hit1.nunique().sort_values(ascending=False).to_string())
def clus(df,gap=5e5):
    d=df[['gwas','hit1']].drop_duplicates().copy()
    d['chr']=d.hit1.str.split(':').str[0]; d['pos']=d.hit1.str.split(':').str[1].astype(int)
    d=d.sort_values(['chr','pos']); cid=[];c=-1;prev=(None,None)
    for ch,pos in zip(d.chr,d.pos):
        if ch!=prev[0] or pos-prev[1]>gap: c+=1
        cid.append(c); prev=(ch,pos)
    d['locus']=cid; return d
d=clus(e)
ns=d.groupby('locus').gwas.nunique(); nl=d.groupby('locus').hit1.nunique()
print("liver loci",d.locus.nunique()); print("studies per locus: mean %.2f median %d"%(ns.mean(),ns.median())); print(ns.value_counts().sort_index().to_string())
print("distinct lead positions per locus: mean %.2f max %d"%(nl.mean(),nl.max()))
dp=clus(p); print("platelet loci (before liver exclusion)",dp.locus.nunique(),"leads per locus mean %.2f"%dp.groupby('locus').hit1.nunique().mean())
pl=pd.read_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3841/pilot_locus_max_pp4.tsv",sep="\t")
print(pl.groupby('study').agg(n=('hit1','size'),med_peaks=('n_peaks','median'),q25=('n_peaks',lambda x:x.quantile(.25)),q75=('n_peaks',lambda x:x.quantile(.75)),gt05=('best_pp4',lambda x:(x>0.5).sum()),gt08=('best_pp4',lambda x:(x>0.8).sum())).to_string())
