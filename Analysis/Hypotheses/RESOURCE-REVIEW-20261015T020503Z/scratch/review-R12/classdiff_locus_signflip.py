# Locus-clustered sign-flip for the post hoc class-difference comparison (281 pairs).
import numpy as np, pandas as pd
root='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
P=pd.read_csv(root+'figures/candidates/bulk-context-774-20261014T140904Z_22236311/class_difference/results/pairs_singlecell_hep_specificity.tsv',sep='\t')
L=pd.read_csv(root+'Analysis/Hypotheses/H10-20261007T210859Z/loci.tsv',sep='\t')
m={}
for _,r in L.iterrows():
    for e in str(r.ensembl).split(','):
        m[e.strip().split('.')[0]]=r.locus_id
P['d']=P.y_genetic-P.y_state
P['locus']=P.g_genetic.str.split('.').str[0].map(m)
print('pairs',len(P),'mean d %.6f'%P.d.mean(),'mapped to locus',P.locus.notna().sum(),'distinct loci',P.locus.nunique())
P['locus']=P.locus.fillna(P.g_genetic)  # unmapped genes are their own cluster
rng=np.random.default_rng(20261015); B=100000
obs=P.d.mean()
# pair-level sign flip (reproduce)
d=P.d.values; s=rng.choice([-1,1],size=(B,len(d)))
null=(s*d).mean(1); print('pair-level sign-flip p %.6f'%((1+(np.abs(null)>=abs(obs)).sum())/(B+1)))
# locus-level sign flip: flip all pairs of a locus together; statistic = mean over pairs
g=P.groupby('locus').d.sum(); n=len(P)
s=rng.choice([-1,1],size=(B,len(g)))
null=(s*g.values).sum(1)/n
print('n clusters',len(g),'locus-level sign-flip p %.6f'%((1+(np.abs(null)>=abs(obs)).sum())/(B+1)))
# cluster sizes
print('cluster size distribution', P.groupby('locus').size().value_counts().sort_index().to_dict())
# locus-mean statistic (each locus weight 1)
gm=P.groupby('locus').d.mean()
s=rng.choice([-1,1],size=(B,len(gm))); null=(s*gm.values).mean(1)
print('locus-mean statistic %.4f'%gm.mean(),'p %.6f'%((1+(np.abs(null)>=abs(gm.mean())).sum())/(B+1)))
# state genes reused across pairs?
print('distinct state genes',P.g_state.nunique(),'distinct genetic genes',P.g_genetic.nunique())
