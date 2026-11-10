#!/usr/bin/env python3
"""Audit reference variant identities and retained COLOC variant exports."""
import concurrent.futures
import csv
import gzip
import hashlib
import json
from pathlib import Path
import pandas as pd

RUN = Path(__file__).resolve().parents[1]
ROOT = RUN.parents[2]
FM = ROOT / 'GWAS/finemapping'
OUT = RUN / 'audit'
OUT.mkdir(exist_ok=False)
PANELS = ['polyfun_eur', '1kg_afr_gram', '1kg_amr_gram', '1kg_eas_gram', '1kg_sas_gram', '1kg_eur']


def inspect(task):
    panel, path = task
    names = ['chr', 'rsid', 'cm', 'pos', 'a1', 'a2']
    d = pd.read_csv(path, sep=r'\s+', header=None, names=names, dtype={'a1':str,'a2':str})
    d['row_1based'] = range(1, len(d)+1)
    duplicated = d[d.duplicated(['chr','pos'], keep=False)].copy()
    rows=[]
    for (chrom,pos), g in duplicated.groupby(['chr','pos'], sort=False):
        pairs = {tuple(sorted((a,b))) for a,b in zip(g.a1,g.a2)}
        ordered = set(zip(g.a1,g.a2))
        kind = ('multiple_allele_pairs' if len(pairs)>1 else
                'opposite_order_same_allele_pair' if len(ordered)>1 else 'exact_allele_duplicate')
        for r in g.to_dict('records'):
            rows.append(dict(panel=panel, block=str(path), collision_kind=kind,
                             n_at_position=len(g), **r))
    info=dict(panel=panel,path=str(path),variants=len(d),duplicate_position_rows=len(duplicated),
              duplicate_positions=duplicated[['chr','pos']].drop_duplicates().shape[0],
              sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    return info,rows


tasks=[]; missing=[]
for panel in PANELS:
    folder=FM/'data/ld_ref'/panel
    blocks=pd.read_csv(folder/'approx_LD_blocks.txt',sep=r'\s+')
    for r in blocks.itertuples(index=False):
        prefix=f'{r.start}.{r.stop}'
        path=folder/f'chr{r.chr}'/prefix/(prefix+'.bim')
        if path.exists():tasks.append((panel,path))
        else:missing.append(dict(panel=panel,path=str(path)))
infos=[];collisions=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    for i,(info,rows) in enumerate(pool.map(inspect,tasks),1):
        infos.append(info);collisions.extend(rows)
        if i%200==0:print(f'BIM blocks checked {i}/{len(tasks)}',flush=True)
pd.DataFrame(infos).to_csv(OUT/'LD_blocks.tsv',sep='\t',index=False)
pd.DataFrame(collisions).to_csv(OUT/'LD_position_collisions.tsv.gz',sep='\t',index=False)
pd.DataFrame(missing,columns=['panel','path']).to_csv(OUT/'missing_LD_blocks.tsv',sep='\t',index=False)
pos=pd.DataFrame(collisions)[['panel','chr','pos']].drop_duplicates()
pos.to_csv(OUT/'LD_ambiguous_positions.tsv',sep='\t',index=False)

r2=FM/'results/susie_coloc_r2_20260924T053244Z'
cols=['gwas_name','gene','ensembl','chr','PP.H4.susie','method','ancestry','lambda_s_locus']
m=pd.read_csv(r2/'susie_coloc_all_gwas.csv',usecols=cols)
assert not m.duplicated(['gwas_name','ensembl','chr']).any()
assert m.gwas_name.nunique()==50
m.to_csv(OUT/'all_50_study_records.tsv.gz',sep='\t',index=False)
m[m['PP.H4.susie'].gt(.5)].to_csv(OUT/'all_positive_gene_studies.tsv',sep='\t',index=False)

replay=ROOT/'Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/genetics'
plan=pd.read_csv(replay/'prepared/replay_plan.tsv',sep='\t',dtype=str)
exports=[]
for i,r in enumerate(plan.itertuples(index=False),1):
    p=replay/'replay_execution'/f'batch_{r.batch_id}'/'exports'/r.gwas_name/f'chr{r.chr}'/r.ensembl/'variant_posteriors.tsv.gz'
    rec=dict(study=r.gwas_name,ensembl=r.ensembl,chromosome=r.chr,path=str(p),available=p.exists())
    if p.exists():
        d=pd.read_csv(p,sep='\t')
        keys=[c for c in ['signal_pair_index','variant_id'] if c in d]
        rec['rows']=len(d);rec['columns']=';'.join(d.columns)
        rec['duplicate_export_keys']=int(d.duplicated(keys).sum()) if len(keys)==2 else None
    exports.append(rec)
    if i%200==0:print(f'Posterior exports checked {i}/{len(plan)}',flush=True)
pd.DataFrame(exports).to_csv(OUT/'posterior_export_inventory.tsv',sep='\t',index=False)
summary=dict(study_records=len(m),studies=m.gwas_name.nunique(),native_genes=m.ensembl.nunique(),
             positive_gene_studies=int(m['PP.H4.susie'].gt(.5).sum()),
             LD_blocks_checked=len(tasks),LD_blocks_missing=len(missing),
             collision_positions_by_panel=pos.groupby('panel').size().to_dict(),
             posterior_exports_planned=len(plan),posterior_exports_present=sum(x['available'] for x in exports),
             interpretation='Reference collision exposure is not confirmed GWAS-fit inclusion; posterior exports may omit unshared variants.')
(OUT/'metadata_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
