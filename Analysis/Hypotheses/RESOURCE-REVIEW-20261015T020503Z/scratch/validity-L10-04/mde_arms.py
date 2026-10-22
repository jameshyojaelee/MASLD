import pandas as pd, numpy as np
o=pd.read_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/validity-L10-04/window_signal_strength.tsv",sep="\t")
o['chr']=o.hit1.str.split(':').str[0]; o['pos']=o.hit1.str.split(':').str[1].astype(int)
def nloci(d,gap=5e5):
    d=d.sort_values(['chr','pos']); c=0; prev=(None,None)
    for ch,pos in zip(d.chr,d.pos):
        if ch!=prev[0] or pos-prev[1]>gap: c+=1
        prev=(ch,pos)
    return c
from math import asin, sqrt, sin
def mde(n1,n2,p0):
    h=(1.959964+0.841621)*sqrt(1/n1+1/n2)
    p1=sin(asin(sqrt(p0))+h/2)**2
    return p1-p0
k={st:nloci(g) for st,g in o.groupby('study')}
print("loci at 500 kb:",k)
for st in ['UKBB_GGT','MVP_ALT_EUR']:
    for p0 in [0.16,0.30,0.46]:
        print(st,"vs platelet-only (%d vs %d loci) control rate %.2f -> MDE %.3f"%(k[st],k['MVP_Platelet_EUR'],p0,mde(k[st],k['MVP_Platelet_EUR'],p0)))
print("pooled 279 vs 127, p0 0.16 -> %.3f"%mde(279,127,0.16))
for st,g in o.groupby('study'):
    print(st,"median GWS variants in window",g.n_gws_5e8.median(),"q75",g.n_gws_5e8.quantile(.75))
