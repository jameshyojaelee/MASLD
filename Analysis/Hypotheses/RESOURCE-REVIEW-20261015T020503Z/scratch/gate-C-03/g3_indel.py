import pandas as pd
D='Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-obs-20260926T155651Z/coloc_direction/full-20260930T120450Z/direction_by_trait.tsv'
t=pd.read_csv(D,sep='\t',low_memory=False,keep_default_na=False,na_values=[''])
r=t[(t.release=='candidate_r2_20260924T053244Z')&t.direction_sign.isin([1,-1])].copy()
O=pd.read_csv('GWAS/finemapping/results/alphagenome_atlas/p0-dbsnp-orientation-20260914T210732Z/tables/dbsnp_indel_orientation.tsv',sep='\t',low_memory=False,keep_default_na=False,na_values=[''])
O['key']=O.hg19_chrom.str.replace('chr','')+':'+O.hg19_pos.astype(str)
O['ak']=O.apply(lambda x: x.key+':'+'/'.join(sorted([str(x.source_allele1),str(x.source_allele2)])),axis=1)
print('verdict_vs_current',O.verdict_vs_current.value_counts().to_dict())
rows=[]
for side,pos,a1,a2 in [('gwas','hit1','gwas_allele1','gwas_allele2'),('eqtl','hit2','eqtl_ea_hit2','eqtl_nea_hit2')]:
    x=r[(r[a1].astype(str).str.len()>1)|(r[a2].astype(str).str.len()>1)].copy()
    x['ak']=x[pos].astype(str)+':'+x.apply(lambda z: '/'.join(sorted([str(z[a1]),str(z[a2])])),axis=1)
    m=x.merge(O[['ak','verdict_vs_current','current_orientation_ambiguous','dbsnp_hg19_state','assemblies_agree']],on='ak',how='left')
    print(side,'indel directional rows',len(x),'matched in dbSNP table',int(m.verdict_vs_current.notna().sum()))
    print(' verdicts',m.verdict_vs_current.value_counts(dropna=False).to_dict())
    print(' current_orientation_ambiguous',m.current_orientation_ambiguous.value_counts(dropna=False).to_dict())
