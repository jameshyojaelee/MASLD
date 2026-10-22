import pandas as pd, sys
O=sys.argv[1]
GW='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/'
t12=pd.read_csv(GW+'susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv')
pos=set(t12.loc[t12.susie_state_t12=='positive','ensembl'])
st=t12.set_index('ensembl').susie_state_t12
e=pd.read_csv(O+'/observed_pairs_eQTLs.tsv',sep='\t')
eh=e[e.clpp>=0.1]
eqcs=set(eh.cs_uid); eqgene_cs=set(zip(eh.cs_uid,eh.ens))
print('eQTL CLPP>=0.1 genes (all):',sorted(set(e.loc[e.clpp>=0.1].ens)))
gname={}
for cls in ['eQTLs','sQTLs','apaQTLs']:
    a=pd.read_parquet('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L08/gtex_v10/Liver.v10.%s.SuSiE_summary.parquet'%cls,columns=['phenotype_id','gene_name'])
    a['ens']=a.phenotype_id.str.extract(r'(ENSG\d+)')[0]; gname.update(dict(zip(a.ens,a.gene_name)))
for cls in ['eQTLs','sQTLs','apaQTLs']:
    d=pd.read_csv(O+'/observed_pairs_%s.tsv'%cls,sep='\t')
    h=d[d.clpp>=0.1].copy()
    h['gene']=h.ens.map(gname); h['study']=h.cs_uid.str.split('|').str[0]
    h['coloc_r2_state']=h.ens.map(st).fillna('not_in_tier12')
    h['gwas_cs_has_eqtl_hit']=h.cs_uid.isin(eqcs)
    h['same_gene_eqtl_hit_same_cs']=[(c,g) in eqgene_cs for c,g in zip(h.cs_uid,h.ens)]
    g=h.groupby('gene').agg(n_gwas_cs=('cs_uid','nunique'),n_resource_cs=('resource','sum'),max_clpp=('clpp','max'),
        coloc_r2=('coloc_r2_state','first'),any_cs_without_eqtl_hit=('gwas_cs_has_eqtl_hit',lambda x:(~x).sum()),
        same_gene_eqtl=('same_gene_eqtl_hit_same_cs','sum'),studies=('study',lambda x:';'.join(sorted(set(x)))[:90]))
    print('\n==',cls,'GWAS CS>=0.1:',h.cs_uid.nunique(),'resource:',h[h.resource].cs_uid.nunique()); print(g.sort_values('n_gwas_cs',ascending=False).to_string())
s=pd.read_csv(O+'/observed_pairs_smallRNA.tsv',sep='\t'); print('\n== smallRNA lead-in-CS'); print(s.to_string())
