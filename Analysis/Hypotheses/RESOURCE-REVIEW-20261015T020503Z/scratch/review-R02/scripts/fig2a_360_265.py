import pandas as pd
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
tier=pd.read_csv(B+'GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t'); main=tier[tier.placement=='main']
tr=dict(zip(main.study_name,main.trait))
cs=pd.read_csv(B+'GWAS/finemapping/results/credible_sets.csv',usecols=['chromosome','position','study','locus'],dtype={'locus':str})
cs=cs[cs.study.isin(tr)]
w=cs.groupby(['study','locus']).agg(ch=('chromosome','first'),lo=('position','min'),hi=('position','max')).reset_index()
w['trait']=w.study.map(tr)
def nmerge(d,pad=250000):
    iv=sorted(zip(d.lo-pad,d.hi+pad)); n=0; end=None
    for a,b in iv:
        if end is not None and a<=end: end=max(end,b)
        else: n+=1; end=b
    return n
within=sum(nmerge(d) for _,d in w.groupby(['trait','ch']))
glob=sum(nmerge(d) for _,d in w.groupby('ch'))
print('within-trait loci',within,'global unique',glob)
print(w.groupby('trait').apply(lambda d: sum(nmerge(x) for _,x in d.groupby('ch'))).to_dict())
