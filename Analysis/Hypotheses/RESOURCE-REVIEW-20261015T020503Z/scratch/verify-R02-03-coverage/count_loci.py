import pandas as pd, sys
from collections import defaultdict
root='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
tier=pd.read_csv(root+'GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t')
MAIN=set(tier.loc[tier.placement=='main','study_name'])
reg=pd.read_csv(root+'GWAS/finemapping/config/gwas_registry.tsv',sep='\t')
reg=reg[reg.study_name.isin(MAIN)]
print('n MAIN',len(MAIN),'n reg rows',len(reg))
cs=pd.read_csv(root+'GWAS/finemapping/results/credible_sets.csv',usecols=['chromosome','position','study','locus','trait','recommended_pip','max_pip'],dtype={'locus':str,'chromosome':str})
cs=cs[cs.study.isin(MAIN)]
per=cs.groupby('study').locus.nunique()
out=[]
for s in sorted(MAIN):
    out.append((s, int(per.get(s,0)), reg.loc[reg.study_name==s,'ancestry'].iloc[0] if (reg.study_name==s).any() else 'NA'))
for s,n,a in out: print(s,a,n)
print('studies with 0 credible-set loci:',[s for s,n,a in out if n==0])
print('total (study,locus) pairs', sum(n for s,n,a in out))
