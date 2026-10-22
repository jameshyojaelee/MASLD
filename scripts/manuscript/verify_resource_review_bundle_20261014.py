#!/usr/bin/env python3
"""Independent standard-library verification after snapshot becomes read-only."""
import csv,hashlib,json,os,re
from pathlib import Path
R=Path(__file__).resolve().parents[2]
name=os.environ.get('RESOURCE_REVIEW_ID')
if name:
 if not re.fullmatch(r'resource-review-[A-Za-z0-9][A-Za-z0-9_.-]*',name) or '..' in name:raise ValueError(name)
 B=R/'figures/candidates'/name
else:
 matches=list((R/'figures/candidates').glob('resource-review-*_'+os.environ['SLURM_JOB_ID']));assert len(matches)==1;B=matches[0]
def h(p):
 z=hashlib.sha256()
 with p.open('rb') as f:
  while True:
   b=f.read(1048576)
   if not b:break
   z.update(b)
 return z.hexdigest()
checks=json.loads((B/'checks.json').read_text());assert checks['bundle_id']==B.name
with (B/'file_inventory.tsv').open() as f:inventory=list(csv.DictReader(f,delimiter='\t'))
assert len(inventory)==checks['files'];assert len({r['copied_path'] for r in inventory})==len(inventory)
for row in inventory:
 p=B/row['copied_path'];assert p.is_relative_to(B/'payload')
 assert h(p)==row['sha256'] and h(Path(row['source_path']))==row['sha256'],p
 assert p.stat().st_size==int(row['bytes']);assert not p.stat().st_mode&0o222
sums=[]
for line in (B/'SHA256SUMS').read_text().splitlines():
 value,rel=line.split('  ',1);p=B/rel;assert p.is_relative_to(B)
 assert h(p)==value,p;sums.append(rel)
allfiles={str(p.relative_to(B)) for p in B.rglob('*') if p.is_file()}
assert set(sums)==allfiles-{'SHA256SUMS'}
assert all(not p.stat().st_mode&0o222 for p in B.rglob('*') if p.is_file())
assert sum(r['role']=='main_panel' for r in inventory)==43
assert checks['main_panel_counts']=={'1':2,'2':9,'3':7,'4':8,'5':6,'6':5,'7':6}
assert checks['catalog_rows']==28233 and sorted(checks['catalog_classes'].values())==[26,400,1129,13127,13551]
external=json.loads((B/'external_sources.json').read_text());scientific=[r for r in external if r['component']=='catalog_source_identity'];assert len(scientific)==55
for row in scientific:assert h(Path(row['source_path']))==row['sha256']
if os.environ.get('RESOURCE_REVIEW_ID'):
 assert all(checks['authority_self_references'].values())
 for n in ['README','PAPER','RESULTS','STATUS','ROADMAP']:assert B.name in (B/'payload/docs'/f'{n}.md').read_text()
context=Path(checks['context_candidate_path'])
copied_context=B/'payload'/context.relative_to(R)/'figure_source_context.tsv'
with copied_context.open() as f:context_rows=list(csv.DictReader(f,delimiter='\t'))
assert len(context_rows)==checks['current_context_source_hashes_checked']
for row in context_rows:assert h(Path(row['source_path']))==row['sha256'],row['source_path']
h17_rows=[r for r in inventory if r['component']=='h17_terminal_verification' and r['role']=='final_terminal_verification'];assert len(h17_rows)==1
h17=json.loads((B/h17_rows[0]['copied_path']).read_text());assert h17['status']=='PASS'
assert h17['resampling_counts']=={'S1_signperm_null.tsv.gz':10000,'S1_matchedset_null.tsv.gz':10000,'negctrl_invivo_label_permutation.tsv':1000,'negctrl_ipsc_label_swap.tsv':100,'bootstrap_S2.tsv':1000}
n3_rows=[r for r in inventory if r['component']=='n3_terminal_verification' and r['role']=='final_terminal_verification'];assert len(n3_rows)==1
n3=json.loads((B/n3_rows[0]['copied_path']).read_text())
for key in ['terminal_raw_fit','terminal_corrected_fit','terminal_ambient_summary','source_only_metadata_changes','native_input_hashes_unchanged','fixed_target_roster_51','tested_families_raw_45_corrected_40','permutations_raw_1000_corrected_1000','seed_20260923_preserved','independent_BH_and_empirical_fallback_checks','fixed_historical_diagnostic_roster_9','controls_preserved_and_pass','correction_coverage_verified']:assert n3.get(key) is True,key
assert checks['review_only'] and not checks['official_release_state_changed'] and not checks['external_published']
assert any('Figure1A' in x for x in checks['incompleteness']);assert checks['compact_figure6_layout_status']=='candidate_not_adopted'
result=dict(bundle_id=B.name,independent_file_sha256_checks=len(inventory),checksum_inventory_checks=len(sums),all_new_files_read_only=True,scientific_input_sha256_checks=55,current_figure_context_source_checks=len(context_rows),main_panel_count=43,catalog_classes_unchanged=True,official_release_state_changed=False,result='PASS')
P=R/'logs'/('resource_review_'+os.environ['SLURM_JOB_ID']+'.independent_check.json');P.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
