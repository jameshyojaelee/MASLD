import pandas as pd
O='test5'
GW='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/'
t12=pd.read_csv(GW+'susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv').set_index('ensembl')
genes={'PNPLA3':'ENSG00000100344','EXOC3L4':'ENSG00000205436','TM4SF4':'ENSG00000169903','C12orf43':'ENSG00000157895','CYP2A6':'ENSG00000255974','RORA':'ENSG00000069667'}
for g,e in genes.items():
    st=t12.susie_state_t12.get(e,'not_in_tier12'); pp=t12.coloc_best_susie_pp4.get(e)
    row=[g,e,'coloc_r2_state='+str(st),'best_susie_pp4='+str(pp)]
    for cls in ['eQTLs','sQTLs','apaQTLs']:
        d=pd.read_csv(O+'/observed_pairs_%s.tsv'%cls,sep='\t'); d=d[d.ens==e]
        row.append('%s max=%s nCS>=0.1=%d'%(cls, round(d.clpp.max(),4) if len(d) else 'NA', d[d.clpp>=0.1].cs_uid.nunique()))
    print(' | '.join(row))
# what is max_pip in the substrate vs per-CS susie_pip at PNPLA3 sQTL lead
s=pd.read_csv(GW+'seqfunc/variant_substrate_hg38.tsv',sep='\t',low_memory=False)
print('substrate rows',len(s),'unique hg19 ids',s.variant_id_hg19.nunique(),'genes',s.gene.nunique())
c=pd.read_csv(GW+'credible_sets.csv',low_memory=False)
x=c[c.variant_id.isin(s[s.gene=='PNPLA3'].variant_id_hg19)]
m=x.groupby('variant_id').agg(max_maxpip=('max_pip','max'),max_susie=('susie_pip','max'),max_carma=('carma_pip','max'),max_recommended=('recommended_pip','max')).sort_values('max_maxpip',ascending=False).head(5)
print(m.to_string())
print(s[s.gene=='PNPLA3'].sort_values('max_pip',ascending=False).head(5)[['variant_id_hg19','max_pip','pos_hg38']].to_string())
