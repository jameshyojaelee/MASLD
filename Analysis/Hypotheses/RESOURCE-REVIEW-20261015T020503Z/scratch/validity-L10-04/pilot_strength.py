import pandas as pd, numpy as np
from scipy.stats import spearmanr, mannwhitneyu
pl=pd.read_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3841/pilot_locus_max_pp4.tsv",sep="\t")
pl['chr']=pl.hit1.str.split(':').str[0]; pl['pos']=pl.hit1.str.split(':').str[1].astype(int)
out=[]
for (st,ch),g in pl.groupby(['study','chr']):
    f="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3841/%s_chr%s.tsv"%(st,ch)
    try: s=pd.read_csv(f,sep="\t",header=None,usecols=[1,6],names=['pos','p'])
    except FileNotFoundError: print("missing",f); continue
    s=s.groupby('pos').p.min()
    for _,r in g.iterrows():
        w=s[(s.index>=r.pos-5e5)&(s.index<=r.pos+5e5)]
        out.append(dict(study=st,hit1=r.hit1,lead_p=s.get(r.pos,np.nan),win_minp=w.min(),best_pp4=r.best_pp4,n_peaks=r.n_peaks))
o=pd.DataFrame(out); o['mlog_win']=-np.log10(o.win_minp.clip(lower=1e-300)); o['mlog_lead']=-np.log10(o.lead_p.clip(lower=1e-300))
print(o.groupby('study')[['mlog_lead','mlog_win','best_pp4','n_peaks']].median().to_string()); print(o.groupby('study').size())
for st,g in o.groupby('study'):
    print(st,"spearman(-log10 window min p, best_pp4)",spearmanr(g.mlog_win,g.best_pp4,nan_policy='omit'))
print("all: spearman",spearmanr(o.mlog_win,o.best_pp4,nan_policy='omit'))
print("all: spearman npeaks vs pp4",spearmanr(o.n_peaks,o.best_pp4,nan_policy='omit'))
o.to_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/validity-L10-04/pilot_strength.tsv",sep="\t",index=False)
