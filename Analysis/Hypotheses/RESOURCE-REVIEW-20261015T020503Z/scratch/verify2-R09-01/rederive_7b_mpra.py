import csv, json
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/'
rows=[r for r in csv.DictReader(open(B+'fig7_panels_rebuild-20260901T012245Z/tables/panel_7b_gains.tsv'),delimiter='\t') if 'MPRA' in r['task']]
print('MPRA rows in panel 7B table:',len(rows))
print('baselines named:',sorted({r['baseline'] for r in rows}))
print('source runs:',sorted({r['source'].split('/')[3] for r in rows}))
cls={'above':0,'below':0,'crosses':0}
for r in rows:
    lo,hi=float(r['lo']),float(r['hi'])
    c='above' if lo>0 else ('below' if hi<0 else 'crosses'); cls[c]+=1
    print(f"  {r['model']:<60} gain {float(r['gain']):+.4f} [{lo:+.4f},{hi:+.4f}] {c}")
print('interval classes:',cls)
h={(r['model_id'],r['head_id']):r for r in csv.DictReader(open(B+'model-cpu-train-605-21099008/evaluation/model_head_summary.tsv'),delimiter='\t')}
hy=float(h[('hyenadna','delta_ridge')]['primary_score']); ai=float(h[('available_simple_controls','allele_identity_ridge')]['primary_score'])
print('605 hyenadna delta_ridge',hy,'allele_identity_ridge',ai,'gain recomputed',hy-ai)
print('605 control rows family:',sorted({r['family'] for k,r in h.items() if r['model_id']=='available_simple_controls'}))
for run in ['model-cpu-train-604-21097189','model-cpu-train-605-21099008']:
    d=json.load(open(B+run+'/evaluation/receipt.json'))
    print(run,'mandatory_task_native_baselines_complete=',d['mandatory_task_native_baselines_complete'],'missing=',d['missing_baselines'],'control=',d['available_control_reference']['head_id'])
