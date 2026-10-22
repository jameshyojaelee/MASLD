import pandas as pd
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/'
c=pd.read_csv(D+'credible_sets.csv',low_memory=False)
print('shape',c.shape,'studies',c.study.nunique(),'loci(study-level str)',c[['study','locus']].drop_duplicates().shape[0])
print('susie_cs values',c.susie_cs.value_counts(dropna=False).head(8).to_dict())
cs=c[c.susie_cs>0]
print('rows in SuSiE CS',len(cs),'unique variants',cs.variant_id.nunique(),'study-locus-cs',cs.groupby(['study','locus','susie_cs']).ngroups)
print('unique chr:pos',cs[['chromosome','position']].drop_duplicates().shape[0])
print('susie_converged',cs.susie_converged.value_counts(dropna=False).to_dict(),'susie_reliable',cs.susie_reliable.value_counts(dropna=False).to_dict())
# pip within CS
g=cs.groupby(['study','locus','susie_cs']).agg(n=('variant_id','size'),pipsum=('susie_pip','sum'),maxpip=('susie_pip','max'),chrom=('chromosome','first'))
print('CS size quantiles',g.n.quantile([0,.25,.5,.75,.9,.99,1]).to_dict())
print('CS pip sum quantiles',g.pipsum.quantile([0,.05,.5,.95,1]).round(3).to_dict())
print('CS per study'); print(g.reset_index().groupby('study').size().describe().round(1).to_dict())
print('ancestry',c.drop_duplicates('study').ancestry.value_counts().to_dict())
print('traits',c.drop_duplicates('study').trait.value_counts().to_dict())
