import csv
from collections import defaultdict
root='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
MAIN={r['study_name'] for r in csv.DictReader(open(root+'GWAS/finemapping/config/gwas_trait_tier.tsv'),delimiter='\t') if r['placement']=='main'}
def trait(s):
    for t in ('NAFLD','NASH','PDFF','ALT','AST','GGT'):
        if t in s.upper(): return t
win={}
for r in csv.DictReader(open(root+'GWAS/finemapping/results/credible_sets.csv')):
    s=r['study']
    if s not in MAIN: continue
    p=int(float(r['position'])); k=(s,r['locus'])
    w=win.get(k)
    if w is None: win[k]=[r['chromosome'],p,p,trait(s)]
    else: w[1]=min(w[1],p); w[2]=max(w[2],p)
def merge(iv,pad=250000):
    iv=sorted((a-pad,b+pad) for a,b in iv); n=0; end=None
    for a,b in iv:
        if end is not None and a<=end: end=max(end,b)
        else: n+=1; end=b
    return n
bt=defaultdict(list); bc=defaultdict(list)
for c,a,b,t in win.values(): bt[(t,c)].append((a,b)); bc[c].append((a,b))
per=defaultdict(int)
for (t,c),iv in bt.items(): per[t]+=merge(iv)
print('per-trait',dict(per),'sum',sum(per.values()),'global',sum(merge(iv) for iv in bc.values()))
