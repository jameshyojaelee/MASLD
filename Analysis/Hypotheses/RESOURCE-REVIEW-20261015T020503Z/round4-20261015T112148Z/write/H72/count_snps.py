import pandas as pd
W='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H72'
b=pd.read_csv(W+'/best_enzyme_row.tsv',sep='\t')
m=pd.read_csv(W+'/mde_rows.tsv',sep='\t'); g=m[m.kind=='gene']
for lab,x in [('PP>=0.5',b[b.pp>=0.5]),('PP>=0.8',b[b.pp>=0.8])]:
    u=x.drop_duplicates(['trait','top_snp'])
    print(lab,'genes',len(x),'unique trait x top_snp',len(u),dict(u.trait.value_counts()))
    def clus(df,w=500000):
        n=0
        for (t,c),d in df.assign(c=df.top_snp.str.split(':').str[0],p=df.top_snp.str.split(':').str[1].astype(int)).groupby(['trait','c']):
            ps=sorted(d.p); n+=1+sum(1 for a,bb in zip(ps,ps[1:]) if bb-a>w)
        return n
    print(lab,'500kb single-linkage clusters (per trait)',clus(u))
x=b[b.pp>=0.5]; u=x.drop_duplicates(['trait','top_snp'])
print('PP>=0.5 unique SNPs present in mde_rows (both sexes usable):',len(g),'; genes whose SNP is usable:',x.merge(g[['trait','id']].rename(columns={'id':'ensembl'}),how='left').shape)
