import pyarrow.parquet as pq, pyarrow.compute as pc, pandas as pd
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
f=pq.ParquetFile(B+'Analysis/SingleCell/results_gpu_v2/hotspot_modules/cell_scores_all.parquet')
vals={'SRR22677587':[],'SRR22677595':[]}
for i in range(f.metadata.num_row_groups):
    t=f.read_row_group(i,columns=['cell_id','module','score'])
    t=t.filter(pc.equal(t['module'],'hepatocytes__20'))
    d=t.to_pandas(); d['s']=d.cell_id.str.split('_',n=1).str[0]
    for k in vals: vals[k]+=d.loc[d.s==k,'score'].tolist()
p=pd.read_csv(B+'figures/candidates/source-diagnosis-programs-20261014T161500Z_22237556/hotspot/donor_program_scores_primary.tsv',sep='\t')
for k,v in vals.items():
    r=p[(p.donor==k)&(p.cell_type=='hepatocytes')&(p.module==20)]
    print(k,'n_cells',len(v),'parquet mean',sum(v)/len(v),'table score',r.score.values, r.disease_stage_coarse.values, r.stage_ordinal.values)
