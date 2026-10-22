import pandas as pd
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/'
c=pd.read_csv(D+'corrected_locus_coordinates.tsv',sep='\t')
ex=pd.read_csv(D+'excluded_credible_sets.tsv',sep='\t')
c=c[~c.credible_set_uid.isin(ex.credible_set_uid)]
c=c[c.pip_sum_gate_passed.astype(str).str.lower()=='true']
c['cohort']=c.study.str.split('_').str[0]
g=c.groupby(['trait_scope','cohort','ancestry'])
t=g.agg(n_cs=('credible_set_uid','size'),n_loci=('locus','nunique'),
        med_members=('n_members','median'),q90_members=('n_members',lambda x:x.quantile(.9)),
        coding_frac=('protein_altering_pip_mass','mean'))
t['cs_per_locus']=t.n_cs/t.n_loci
pd.set_option('display.width',200)
print(t.round(3))
# loci with >=5 CS
lc=c.groupby(['study','locus']).size()
print('loci with >=5 CS by study:'); print(lc[lc>=5].groupby(level=0).size())
print('max CS per locus:',lc.max())
