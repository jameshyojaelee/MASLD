import pandas as pd
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/'
c=pd.read_csv(D+'corrected_locus_coordinates.tsv',sep='\t')
ex=pd.read_csv(D+'excluded_credible_sets.tsv',sep='\t')
c=c[~c.credible_set_uid.isin(ex.credible_set_uid)]
c=c[c.pip_sum_gate_passed.astype(str).str.lower()=='true']
print('gated instances by scope', c.groupby('trait_scope').size().to_dict())
c['k']=c.groupby(['study','locus']).credible_set_uid.transform('size')
c['cohort']=c.study.str.split('_').str[0]
pd.set_option('display.width',250); pd.set_option('display.max_columns',20)
t=c.groupby(['trait_scope','cohort','ancestry']).agg(n_cs=('credible_set_uid','size'),
   n_loci=('locus','nunique'), loci_at_cap=('k',lambda x:(x==10).sum()/10),
   cs_in_cap_loci=('k',lambda x:(x==10).sum()), med_members=('n_members','median'),
   coding_mean=('protein_altering_pip_mass','mean'))
print(t.round(4).to_string())
e=c[c.trait_scope=='tier2_liver_enzyme']
for lab,sub in [('enzyme all',e),('enzyme loci<10 CS',e[e.k<10]),('enzyme loci<=3 CS',e[e.k<=3]),('enzyme loci==1 CS',e[e.k==1])]:
    m=sub.protein_altering_pip_mass+sub.canonical_splice_pip_mass
    print(lab,'n_cs',len(sub),'coding+splice mean',round(m.mean(),4))
d=c[c.trait_scope!='tier2_liver_enzyme']
m=d.protein_altering_pip_mass+d.canonical_splice_pip_mass
print('direct all n_cs',len(d),'coding+splice mean',round(m.mean(),4))
