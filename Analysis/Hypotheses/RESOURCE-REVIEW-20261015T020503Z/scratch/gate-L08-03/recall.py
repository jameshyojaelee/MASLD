import pandas as pd, sys
O=sys.argv[1]
GW='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/'
G='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L08/gtex_v10/'
t12=pd.read_csv(GW+'susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv')
pos=set(t12.loc[t12.susie_state_t12=='positive','ensembl'])
e=pd.read_csv(O+'/observed_pairs_eQTLs.tsv',sep='\t')
hit=set(e.loc[e.clpp>=0.1,'ens']); hitr=set(e.loc[(e.clpp>=0.1)&e.resource,'ens'])
eq=pd.read_parquet(G+'Liver.v10.eQTLs.SuSiE_summary.parquet',columns=['phenotype_id']); eqg=set(eq.phenotype_id.str.extract(r'(ENSG\d+)')[0])
print('COLOC r2 SuSiE-positive genes',len(pos),'| with any GTEx v10 liver eQTL CS',len(pos&eqg),'| with eQTL CLPP>=0.1 (any GWAS CS)',len(pos&hit),'| (resource CS)',len(pos&hitr))
for cls in ['sQTLs','apaQTLs']:
    a=pd.read_parquet(G+'Liver.v10.%s.SuSiE_summary.parquet'%cls,columns=['phenotype_id']); ag=set(a.phenotype_id.str.extract(r'(ENSG\d+)')[0])
    d=pd.read_csv(O+'/observed_pairs_%s.tsv'%cls,sep='\t'); h=set(d.loc[d.clpp>=0.1,'ens'])
    print(cls,'COLOC-positive genes with any GTEx v10 liver CS',len(pos&ag),'| with CLPP>=0.1',len(pos&h))
d=pd.read_csv(O+'/observed_pairs_sQTLs.tsv',sep='\t'); h=d[(d.clpp>=0.1)]
eh=set(e.loc[e.clpp>=0.1,'cs_uid'])
b=h[~h.cs_uid.isin(eh)&~h.ens.isin(pos)]
print('sQTL beyond-eQTL rows (full universe):'); print(b[['cs_uid','qcs','clpp','resource']].to_string())
