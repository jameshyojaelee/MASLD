import pandas as pd, numpy as np, time
e=pd.read_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L10-04/eur_tier12_pairs.tsv",sep="\t"); p=pd.read_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L10-04/platelet_eur_pairs.tsv",sep="\t")
def pos(df):
    d=df[['gwas','hit1']].drop_duplicates().copy(); d['chr']=d.hit1.str.split(':').str[0]; d['pos']=d.hit1.str.split(':').str[1].astype(int); return d
le=pos(e); lp=pos(p)
# platelet-only: >500kb from every liver lead
def far(r):
    l=le[le.chr==r.chr]; return (np.abs(l.pos-r.pos)>5e5).all()
lp['platelet_only']=lp.apply(far,axis=1)
print("platelet leads",len(lp),"platelet-only",lp.platelet_only.sum())
sets={'UKBB_GGT':le[le.gwas=='UKBB_GGT'],'MVP_ALT_EUR':le[le.gwas=='MVP_ALT_EUR'],'MVP_Platelet_EUR':lp[lp.platelet_only]}
rows=[]
for st,L in sets.items():
    t=time.time()
    s=pd.read_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/data/sumstats/%s_reformatted_hg19.tsv"%st,sep="\t",usecols=[0,1,6],dtype={0:str},low_memory=False)
    s.columns=['chr','pos','p']
    for ch,Lc in L.groupby('chr'):
        sc=s[s.chr==ch].sort_values('pos'); P=sc.pos.values; PV=sc.p.values
        for _,r in Lc.iterrows():
            i0,i1=np.searchsorted(P,r.pos-5e5),np.searchsorted(P,r.pos+5e5,side='right')
            j=np.searchsorted(P,r.pos); lead=PV[j] if j<len(P) and P[j]==r.pos else np.nan
            w=PV[i0:i1]; rows.append(dict(study=st,hit1=r.hit1,lead_p=lead,win_minp=np.nanmin(w),n_gws_5e8=int((w<5e-8).sum())))
    print(st,"read+scan %.0fs"%(time.time()-t)); del s
o=pd.DataFrame(rows)
o['gap']=-np.log10(o.win_minp.clip(lower=1e-300))+np.log10(o.lead_p.clip(lower=1e-300))
o.to_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/validity-L10-04/window_signal_strength.tsv",sep="\t",index=False)
for st,g in o.groupby('study'):
    print(st,"leads",len(g),"lead missing",g.lead_p.isna().sum(),"median -log10 lead p %.1f"%(-np.log10(g.lead_p.clip(lower=1e-300))).median(),
      "| window has a signal >=3 log10 units stronger than lead: %d (%.0f%%)"%((g.gap>=3).sum(),100*(g.gap>=3).mean()),
      "| lead within 0.5 of window top: %d"%((g.gap<0.5).sum()))
