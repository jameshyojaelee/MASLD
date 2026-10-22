#!/usr/bin/env python3
"""Bounded review-only snapshot; never changes scientific or release state."""
import csv, hashlib, json, os, shutil, re
from pathlib import Path
from datetime import datetime, timezone
import pandas as pd
R=Path(__file__).resolve().parents[2]
review_id=os.environ.get('RESOURCE_REVIEW_ID') or ('resource-review-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+os.environ['SLURM_JOB_ID'])
if not re.fullmatch(r'resource-review-[A-Za-z0-9][A-Za-z0-9_.-]*',review_id) or '..' in review_id:raise ValueError('RESOURCE_REVIEW_ID must be a safe single resource-review directory name')
B=R/'figures/candidates'/review_id
B.mkdir(exist_ok=False); selected={}; expected={}; external=[]; panels={}
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''): h.update(b)
 return h.hexdigest()
def add(p,c,role='source',sha=None):
 p=Path(p); p=p if p.is_absolute() else R/p
 if not p.is_file(): raise FileNotFoundError(p)
 if p.stat().st_size>256*1024*1024: raise ValueError('Unexpected large selected file: '+str(p))
 selected[p]=(c,role)
 if sha: expected[p]=sha

def rows(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def smalltree(p,c):
 for f in p.rglob('*'):
  if f.is_file() and f.suffix in {'.tsv','.csv','.json','.txt','.md','.parquet','.sha256'} and not any(x in str(f.relative_to(p)) for x in ['before','legacy','archive','pre-','portal_candidate','portal_export','preflight','ambient_failed_']): add(f,c)
for n in ['README','PAPER','RESULTS','STATUS','ROADMAP']:add(R/'docs'/f'{n}.md','manuscript','authority')
for n in ['README.md','METHODS.md','05_figure_legends.md']:add(R/'docs/manuscript'/n,'manuscript')
for f in (R/'docs/manuscript/draft').glob('*.md'):add(f,'manuscript')
add(R/'figures/main/INDEX.md','figures','index')
add(R/'docs/technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23_OUTCOMES.md','technical_outcomes')
D={i:next((R/'figures/main').glob(f'fig{i}_*')) for i in range(1,8)}
for i,d in D.items():
 for f in d.glob('*'):
  if f.is_file() and f.suffix in {'.md','.tsv','.json'}:add(f,f'figure{i}','index')
 if i==1:
  ps=[d/'fig1b_complementary_maps.pdf',d/'fig1c_genetic_trait_scope.pdf']
  for f in d.glob('fig1[bc]*'):
   if f.suffix in {'.tsv','.txt'}:add(f,'figure1')
 elif i==2:
  ps=[]
  for row in rows(d/'panels/Figure2_reindex_manifest.tsv'):
   if row['panel'] in list('ABCDEFGHI'):
    p=d/'panels'/row['canonical_pdf'];ps.append(p);expected[p]=row['sha256']
  for f in (d/'panels').glob('*'):
   if f.is_file() and f.suffix in {'.tsv','.csv','.txt'} and (f.name.startswith(('Fig2','fig2')) or f.name=='Figure2_reindex_manifest.tsv'):add(f,'figure2')
  smalltree(d/'source_tables','figure2')
  for n in ['Fig2_companion_noncoding_annotation_gap.pdf','Fig2_companion_eqtl_target_gene_biotype.pdf']:add(d/'panels'/n,'figure2','companion')
 elif i==5:
  ps=[]
  for row in rows(d/'CANONICAL_MAIN_PANELS.tsv'):
   p=d/row['filename'];ps.append(p);expected[p]=row['sha256']
  smalltree(d/'panels/data','figure5')
  for n in ['composite_mrna_protein_corrected.csv','composite_protein_histology_partial_corrected.csv']:add(d/'data'/n,'figure5')
  add(d/'fig5_molecular_context.pdf','figure5','review_assembly')
 else:
  ps=[]
  for row in rows(d/'CANDIDATE_PANEL_INDEX.tsv'):
   p=Path(row['candidate_relative_path']);p=R/p if str(p).startswith('figures/') else d/p
   if row['callout'].startswith(str(i)):ps.append(p)
   elif i==3 and row['callout']=='unassigned_inset':add(p,'figure3','unassigned_companion')
  idx=d/'CURRENT_CORRECTED_SOURCE_INDEX.tsv'
  if idx.exists():
   for row in rows(idx):
    p=R/row['path']
    if p.suffix!='.pdf':add(p,f'figure{i}',sha=row.get('sha256'))
  if (d/'review_assembly').exists():
   for f in (d/'review_assembly').glob('*'):
    if f.is_file() and f.suffix in {'.pdf','.tsv','.json','.txt','.md'}:add(f,f'figure{i}','review_assembly')
 for p in ps:add(p,f'figure{i}','main_panel')
 panels[str(i)]=len(ps)
assert panels=={'1':2,'2':9,'3':7,'4':8,'5':6,'6':5,'7':6},panels
model_tables=R/'Analysis/MASLD_Model_Benchmark/executions/fig7_panels_rebuild-20260901T012245Z/tables'
for f in model_tables.glob('*'):
 if f.is_file() and '.pre-' not in f.name:
  if f.name in {'panel_7c_exposure.tsv','panel_7f_exposure_grid.tsv'}:add(f,'model_exposure_historical','historical_superseded_by_20261014_source_correction')
  elif f.name=='panel_7a_framework.tsv':add(f,'model_census_labels_historical','historical_labels_superseded_20261014_counts_unchanged')
  elif f.suffix in {'.tsv','.json','.txt'}:add(f,'figure7','unchanged_numeric_or_design_source')
model_corrected=R/'Analysis/MASLD_Model_Benchmark/executions/fig7-source-exposure-correction-20261014'
smalltree(model_corrected/'tables','model_source_correction')
for n in ['README.md','checks.json','source_receipt_checks.json']:add(model_corrected/n,'model_source_correction')
for f in (model_corrected/'scripts').glob('*'):
 if f.is_file() and f.suffix in {'.py','.R','.sbatch'}:add(f,'model_source_correction','executed_source')
smalltree(R/'Analysis/MASLD_Model_Benchmark/executions/encoder-exposure-step4-20260922T175715Z/out','model_exact_source_receipts')
C=R/'figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8'
for f in C.glob('*'):
 if f.is_file():add(f,'catalog_r8')
K=R/'figures/candidates/catalog-r8-completion-20261014T141645Z_22236527'
for n in ['review_completion.json','input_checksums_before.json','input_checksums_after.json']:add(K/n,'catalog_r8')
add(K/'molecular_objects/molecular_object_sidecar.jsonl','catalog_r8','object_sidecar')
for n in ['sidecar_completion.json','worked_entries.tsv']:add(K/'molecular_objects'/n,'catalog_r8')
smalltree(K/'identity_reconciliation_v3','catalog_identity')
before=json.loads((K/'input_checksums_before.json').read_text());after=json.loads((K/'input_checksums_after.json').read_text());assert before==after and len(before)==55
for path,h in before.items():
 p=Path(path);actual=digest(p);assert actual==h,(path,actual,h)
 external.append(dict(component='catalog_source_identity',source_path=path,sha256=h,copied=p in selected,state='55_unchanged_scientific_inputs'))
G=R/'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z'
for n in ['ASSEMBLY_COMPLETE.json','gene_level_coloc.csv','gene_level_coloc_tier12.csv','input_manifest.tsv','output_sha256.tsv','r2_summary.json','sessionInfo.txt']:add(G/n,'genetics_current')
smalltree(G/'summary','genetics_current')
add(R/'GWAS/finemapping/results/coloc_r2_adoption_review_20261014_22236252/review.json','genetics_current')
method_audit=R/'GWAS/finemapping/results/coloc_r2_adoption_review_20261014_22236252/method_availability_22239649'
smalltree(method_audit,'genetics_method_availability')
for f in method_audit.glob('*'):
 if f.is_file():add(f,'genetics_method_availability_precursor','earlier_precursor_superseded_identifier_reconciliation')
current_method_audit=R/'GWAS/finemapping/results/coloc_r2_adoption_review_20261014_22236252/method_availability_22239724'
for f in current_method_audit.glob('*'):
 if f.is_file() and f.suffix in {'.tsv','.json','.txt','.py','.R','.sbatch'}:add(f,'genetics_method_availability_current','current_identified_keys_with_anonymous_row_reconciliation')

U=R/'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z'
for f in U.glob('*'):
 if f.is_file() and f.suffix in {'.json','.tsv','.txt'}:add(f,'bulk_current')
for f in (U/'primary_source_controls').glob('*'):
 if f.is_file() and (f.suffix in {'.csv','.tsv','.json','.txt'} or f.name.endswith('.tsv.gz')):add(f,'bulk_current')
V=R/'figures/candidates/bulk-context-774-20261014T140904Z_22236311'
for n in ['installed_5bc.json','panel4c_bulk_gene_ids.tsv','scientific_checks.txt','native_GSE126848_checks.tsv','inputs_before.sha256']:add(V/n,'bulk_context')
smalltree(V/'proteomics','bulk_context');smalltree(V/'class_difference/results','class_comparison');smalltree(V/'class_difference/prespec','class_comparison')
add(V/'class_difference/VERDICT.md','class_comparison')
smalltree(V/'class_difference/negative-control-20261014T144652Z_22236840','class_comparison_control')
add(R/'Analysis/Multimodal_Program_Projection/config/panel4c_fixed_rows.tsv','bulk_context')
L=R/'GWAS/ldsc/rg_missing_compatible_20261014T150911Z_22237000'
terminal=next(L.glob('*66pair*terminal*.tsv'));add(terminal,'ldsc_candidate');add(L/'CHECKED.txt','ldsc_candidate');add(L/'COMPLETE.txt','ldsc_candidate')
assert len(pd.read_csv(terminal,sep='\t'))==66
# Exposure alternatives remain separate from the six indexed Figure7 panels.
for sub in ['framework7','encoder7']:
 p=R/'figures/candidates/resource-corrections-20261005-v1/figure7_review/resumed_review_v2'/sub
 if p.exists():
  for f in p.rglob('*'):
   if f.is_file() and f.suffix in {'.tsv','.json','.md'}:external.append(dict(component='model_exposure_alternative',source_path=str(f),sha256=digest(f),copied=False,state='candidate_not_selected_main_panels'))
# Explicit current correction companions; prior named sources remain unchanged.
correction=R/'figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014'
for f in correction.glob('*'):
 if f.is_file() and '.before_' not in f.name:add(f,'program_stage_correction')
smalltree(R/'figures/main/fig4_singlecell_programs/source_tables/current_candidate/source_qualified_ccc_20261014','source_qualified_ccc')
for n in ['figs4c_hotspot_stage_heatmap_full.pdf','figs4d_communication_heatmap.pdf','figs4e_program_linked_communication.pdf']:
 p=D[4]/'panels/supplementary'/n
 if p.exists():add(p,'figure4','current_supplementary_companion')
for root in ['Analysis/SingleCell/candidates/source-qualified-CCC-20261014_22237840','Analysis/Hypotheses/H34-stage-source-corrected-20261014_22237596','Analysis/Hypotheses/H17-source-stage-corrected-20261014_22237956','Analysis/Hypotheses/N3-source-diagnosis-20261014']:
 p=R/root
 if p.exists():smalltree(p,'current_source_qualification')
# Compact analysis output products only: exclude raw count matrices/embeddings/RDS.
for root in ['Analysis/Hypotheses/H17-source-stage-corrected-20261014_22237956','Analysis/Hypotheses/N3-source-diagnosis-20261014']:
 p=R/root
 for f in p.rglob('*'):
  if not f.is_file():continue
  if 'ambient_failed_' in str(f.relative_to(p)):continue
  if f.name.endswith(('.tsv.gz','.csv.gz')) and not any(x in f.name.lower() for x in ['donor_gene_celltype_counts','hep_donor_counts','counts_raw','count_matrix','embedding']):add(f,'terminal_resampling_and_effect_products','compact_completed_output_not_raw_matrix')
  if f.suffix in {'.py','.R','.sbatch'}:add(f,'terminal_verification_source','executed_source_or_launcher')
h17_terminal=Path(os.environ.get('H17_TERMINAL_VERIFICATION_ROOT',str(R/'Analysis/Hypotheses/H17-terminal-check-20261014_22239530')))
assert h17_terminal.resolve().is_relative_to(R)
h17_verification=json.loads((h17_terminal/'verification.json').read_text())
assert h17_verification['status']=='PASS'
assert h17_verification['resampling_counts']=={'S1_signperm_null.tsv.gz':10000,'S1_matchedset_null.tsv.gz':10000,'negctrl_invivo_label_permutation.tsv':1000,'negctrl_ipsc_label_swap.tsv':100,'bootstrap_S2.tsv':1000}
add(h17_terminal/'verification.json','h17_terminal_verification','final_terminal_verification')
smalltree(h17_terminal,'h17_terminal_verification')
for f in h17_terminal.glob('*'):
 if f.is_file() and f.suffix in {'.py','.R','.sbatch'}:add(f,'h17_terminal_verification','independent_checker_or_launcher')
n3_terminal=R/'Analysis/Hypotheses/N3-source-diagnosis-20261014/verification.json'
n3_verification=json.loads(n3_terminal.read_text())
n3_required=['terminal_raw_fit','terminal_corrected_fit','terminal_ambient_summary','source_only_metadata_changes','native_input_hashes_unchanged','fixed_target_roster_51','tested_families_raw_45_corrected_40','permutations_raw_1000_corrected_1000','seed_20260923_preserved','independent_BH_and_empirical_fallback_checks','fixed_historical_diagnostic_roster_9','controls_preserved_and_pass','correction_coverage_verified']
assert all(n3_verification.get(k) is True for k in n3_required),n3_verification
add(n3_terminal,'n3_terminal_verification','final_terminal_verification')
context=Path(os.environ['FIGURE_CONTEXT_REVIEW_ROOT'])
assert (context/'checks.json').exists()
context_rows=rows(context/'figure_source_context.tsv')
assert len({r['source_path'] for r in context_rows})==len(context_rows)
for row in context_rows:assert digest(Path(row['source_path']))==row['sha256'],row['source_path']
smalltree(context,'figure_claim_context')
for f in context.glob('*.jsonl'):add(f,'figure_claim_context')
layout=R/'figures/candidates/fig6-compact-publication-20261014T160441Z_22237441'
for n in ['Figure6.pdf','panel_manifest.tsv','panel_identity_checks.tsv','inputs_before.sha256','inputs_unchanged.txt','COMPLETE.txt','sessionInfo.txt']:add(layout/n,'figure6_layout_proposal','candidate_not_adopted')
for root, component, filenames in [
 ('figures/candidates/fig1-embedded-review-proposal-22238522','figure1_layout_proposal',['Figure1.pdf','panel_manifest.tsv','inputs_before.sha256','inputs_unchanged.txt','typography_check.json','assemble_fig1_embedded_review.py']),
 ('figures/candidates/fig1-review-proposal-22238377','figure1d_proposal',['fig1d_analysis_workflow.pdf','fig1d_workflow_source.tsv','render_fig1d_resource.R','assemble_fig1_review.py']),
 ('figures/candidates/fig2-native-review-proposal-20261014_22238526','figure2_layout_proposal',['Figure2_A-I_native_review_PROPOSAL.pdf','CHECKED.txt','identity_font_hash_checks.json','panel_identity_geometry.json','individual_panel_printed_text.txt','unassembled_companion_references.txt','assemble_fig2_native_proposal.py']),
 ('figures/candidates/figure7-embedded-native-review-proposal-22239269','figure7_layout_proposal',['Figure7_embedded_native_review_PROPOSAL.pdf','panel_manifest.tsv','embedded_identity_font_checks.json','source_identity.json','README.txt'])]:
 for name in filenames:add(R/root/name,component,'candidate_layout_not_adopted')
p=Path(os.environ.get('FIGURE1A_PROPOSAL_ROOT',str(R/'figures/candidates/fig1a-resource-proposal-22238205')));smalltree(p,'figure1a_proposal')
for f in p.glob('*'):
 if f.suffix in {'.pdf','.R'}:add(f,'figure1a_proposal','candidate_not_current_indexed_authority')
add(Path(__file__),'collector','executed_source')
add(R/'scripts/manuscript/verify_resource_review_bundle_20261014.py','collector','independent_checker')
add(R/'scripts/manuscript/run_resource_review_bundle_20261014.sbatch','collector','submitted_launcher')
gene=pd.read_parquet(C/'passport_gene_index.parquet');classes=gene['primary_evidence_class'].value_counts().to_dict();assert sorted(classes.values())==[26,400,1129,13127,13551],classes
inventory=[]
for p,(component,role) in sorted(selected.items(),key=lambda x:str(x[0])):
 h=digest(p)
 if p in expected:assert h==expected[p],(p,h,expected[p])
 rel=Path('payload')/p.relative_to(R);dest=B/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
 assert digest(dest)==h and digest(p)==h,p
 inventory.append(dict(component=component,role=role,source_path=str(p),copied_path=str(rel),sha256=h,bytes=p.stat().st_size,state='review_only'))
with (B/'file_inventory.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(inventory[0]),delimiter='\t');w.writeheader();w.writerows(inventory)
(B/'external_sources.json').write_text(json.dumps(external,indent=2)+'\n')
components={c:sum(x['component']==c for x in inventory) for c in sorted({x['component'] for x in inventory})}
authority_refs={n:B.name in (R/'docs'/f'{n}.md').read_text() for n in ['README','PAPER','RESULTS','STATUS','ROADMAP']}
if os.environ.get('RESOURCE_REVIEW_ID'):assert all(authority_refs.values()),authority_refs
checks=dict(bundle_id=B.name,authority_self_references=authority_refs,review_only=True,external_published=False,official_release_state_changed=False,files=len(inventory),components=components,bytes=sum(x['bytes'] for x in inventory),main_panel_counts=panels,main_panel_total=sum(panels.values()),catalog_classes=classes,catalog_rows=len(gene),catalog_scientific_inputs_verified=55,ldsc_pairs=66,ldsc_status='candidate_not_adopted',exposure_alternatives_status='candidate_not_selected_main_panels',context_candidate_path=str(context),current_context_source_hashes_checked=len(context_rows),compact_figure6_layout_status='candidate_not_adopted',figure1_2_7_layout_status='candidate_not_adopted',model_current_exposure_source='source-corrected-20261014 TF5clean2seen Geneformer0clean; earlier Sept7C/7F exposure explicitly historical',source_correction_companions='current_separate_from_named_r8',figure1a_proposal_status='candidate_not_current_indexed_authority_if_supplied',incompleteness=['Figure1A has no indexed PDF; unindexed figure1_JL.png not copied.','Figure6 native-size review assembly is not final journal layout.'],copied_sha256_checks='all_passed')
(B/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
(B/'SHA256SUMS').write_text(''.join(digest(p)+'  '+str(p.relative_to(B))+'\n' for p in sorted(B.rglob('*')) if p.is_file()))
for p in B.rglob('*'):
 if p.is_file():p.chmod(0o444)
print(json.dumps(dict(path=str(B),**checks),indent=2))
