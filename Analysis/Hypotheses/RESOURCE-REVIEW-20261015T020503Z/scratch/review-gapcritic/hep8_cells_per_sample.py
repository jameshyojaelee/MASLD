# Count cells contributing to the hepatocytes__8 (Ductular injury) donor score, per atlas sample (run).
import pyarrow.parquet as pq, pyarrow.compute as pc, collections
f=pq.ParquetFile('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/SingleCell/results_gpu_v2/hotspot_modules/cell_scores_all.parquet')
cnt=collections.Counter()
for i in range(f.metadata.num_row_groups):
    t=f.read_row_group(i,columns=['cell_id','module'])
    t=t.filter(pc.equal(t['module'],'hepatocytes__8'))
    for cid in t['cell_id'].to_pylist():
        cnt[cid.split('_',1)[0]]+=1
items=sorted(cnt.items(),key=lambda x:x[1])
print('samples with hepatocytes__8 scores:',len(items),' total cells:',sum(cnt.values()))
for s in ['SRR22677587','SRR22677595']: print(s,cnt.get(s,0))
print('samples with <50 cells:',sum(1 for k,v in items if v<50))
print('smallest 15:',items[:15])
