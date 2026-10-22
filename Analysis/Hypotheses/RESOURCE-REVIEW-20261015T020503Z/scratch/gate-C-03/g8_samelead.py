import pandas as pd
D='Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-obs-20260926T155651Z/coloc_direction/full-20260930T120450Z/direction_by_trait.tsv'
t=pd.read_csv(D,sep='\t',low_memory=False,keep_default_na=False,na_values=[''])
r=t[(t.release=='candidate_r2_20260924T053244Z')&t.direction_sign.isin([1,-1])]
print('directional r2 rows',len(r),'same_lead',r.same_lead.value_counts(dropna=False).to_dict())
d=pd.read_csv("Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-C-03/gene_dir_h26rule.tsv",sep="\t")
calls=d[d.dir.isin(['call_up','call_down'])].ensembl
x=r[r.ensembl.isin(calls)].groupby('ensembl').same_lead.apply(lambda s: s.astype(str).str.upper().eq('TRUE').any())
print('called genes with >=1 same-lead row',int(x.sum()),'called genes resting only on LD-joined rows',int((~x).sum()),'of',len(x))
