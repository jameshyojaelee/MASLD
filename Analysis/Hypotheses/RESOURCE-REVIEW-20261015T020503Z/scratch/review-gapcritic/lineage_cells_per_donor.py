# Cells behind each donor x lineage score in the Fig 4C primary fit (module 1 of each lineage as a proxy for all modules of that lineage).
import pyarrow.parquet as pq, pyarrow.compute as pc, pandas as pd, collections
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
f=pq.ParquetFile(B+'Analysis/SingleCell/results_gpu_v2/hotspot_modules/cell_scores_all.parquet')
cnt=collections.Counter()
for i in range(f.metadata.num_row_groups):
    t=f.read_row_group(i,columns=['cell_id','module','cell_type'])
    t=t.filter(pc.match_substring_regex(t['module'],r'__1$'))
    d=t.to_pandas(); d['s']=d.cell_id.str.split('_',n=1).str[0]
    for (ct,s),n in d.groupby(['cell_type','s']).size().items(): cnt[(ct,s)]+=n
c=pd.Series(cnt).rename_axis(['cell_type','sample']).reset_index(name='n_cells')
p=pd.read_csv(B+'figures/candidates/source-diagnosis-programs-20261014T161500Z_22237556/hotspot/donor_program_scores_primary.tsv',sep='\t',usecols=['donor','cell_type','module','dataset','stage_ordinal'])
p=p[p.module==1].drop(columns='module')
p=p[p.stage_ordinal.notna()]
m=p.merge(c,left_on=['donor','cell_type'],right_on=['cell_type','sample'][::-1][::-1] if False else ['sample','cell_type'],how='left')
m['run_level_key']=m.n_cells.notna()
pd.set_option('display.width',200)
print('donor x lineage rows with stage (fit roster proxy):',len(m))
print('rows whose donor id is a single run (count exact):',int(m.run_level_key.sum()))
s=m[m.run_level_key]
print(s.assign(lt50=s.n_cells<50,lt20=s.n_cells<20).groupby('cell_type').agg(donors=('donor','nunique'),lt50=('lt50','sum'),lt20=('lt20','sum'),min_cells=('n_cells','min')).to_string())
print(s[s.n_cells<50].sort_values(['cell_type','n_cells'])[['cell_type','donor','dataset','stage_ordinal','n_cells']].to_string())
m.to_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-gapcritic/donor_lineage_cell_counts.tsv',sep='\t',index=False)
