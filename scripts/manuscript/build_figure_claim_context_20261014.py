#!/usr/bin/env python3
"""Read-only source context and exact r8 joins; no new scientific calls."""
import csv, hashlib,json,os,gzip
from pathlib import Path
from datetime import datetime,timezone
import pandas as pd
R=Path(__file__).resolve().parents[2]; O=R/'figures/candidates'/('figure-claim-context-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+os.environ['SLURM_JOB_ID']);O.mkdir(exist_ok=False)
C=R/'figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8';K=R/'figures/candidates/catalog-r8-completion-20261014T141645Z_22236527'
def h(p):
 z=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):z.update(b)
 return z.hexdigest()
gp=C/'passport_gene_index.parquet';ep=C/'passport_evidence_long.parquet';gh=h(gp);eh=h(ep)
g=pd.read_parquet(gp).astype(object).where(pd.notna(pd.read_parquet(gp)), '');e=pd.read_parquet(ep).astype(object).where(pd.notna(pd.read_parquet(ep)), '');assert g.ensembl_id.is_unique
assert sorted(g.primary_evidence_class.value_counts().values)==[26,400,1129,13127,13551]
lookup={str(r['ensembl_id']).split('.')[0]:r for r in g.to_dict('records')};symbols={}
for r in g.to_dict('records'):symbols.setdefault(str(r['symbol']),[]).append(r)
ev={}
for r in e.to_dict('records'):ev.setdefault(str(r['ensembl_id']).split('.')[0],[]).append(r)
side={}
with (K/'molecular_objects/molecular_object_sidecar.jsonl').open() as f:
 for line in f:
  r=json.loads(line);side[r['ensembl_id_base']]=r
vocab=pd.read_csv(C/'passport_controlled_vocabularies.tsv',sep='\t');allowed=set(vocab.loc[vocab.vocabulary_name=='call_state','value']);assert set(e.call_state)<=allowed
B=R/'figures/candidates/resource-review-20261014T155914Z_22237418'; inv=pd.read_csv(B/'file_inventory.tsv',sep='\t').fillna('')
# Refresh exact working source paths from current indices plus diagnosis companion.
records={r['source_path']:r for r in inv.to_dict('records')}
for i in range(1,7):
 d=next((R/'figures/main').glob(f'fig{i}_*')); idx=d/'CURRENT_CORRECTED_SOURCE_INDEX.tsv'
 if idx.exists():
  for rr in pd.read_csv(idx,sep='\t',dtype=str,keep_default_na=False).to_dict('records'):
   p=R/rr['path']; assert h(p)==rr['sha256'], p
   records[str(p)]={'component':f'figure{i}','source_path':str(p),'sha256':rr['sha256']}
companion=R/'figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014'
for p in companion.glob('*'):
 if p.is_file() and '.before_' not in p.name:records[str(p)]={'component':'figure4','source_path':str(p),'sha256':h(p)}
for rr in records.values():
 if rr['component'] in [f'figure{i}' for i in range(1,7)]:rr['sha256']=h(Path(rr['source_path']))
ccc=R/'figures/main/fig4_singlecell_programs/source_tables/current_candidate/source_qualified_ccc_20261014'
for p in ccc.glob('*'):
 if p.is_file():records[str(p)]={'component':'figure4','source_path':str(p),'sha256':h(p)}
inv=pd.DataFrame(records.values())
def panel_context(component,path):
 name=path.name.lower(); text=str(path).lower()
 unit='Source-native unit as explicitly recorded; metadata rows are not observations.'; family='Not applicable to source metadata; no statistical test introduced.'; alt='This source records inputs or descriptive quantities and does not establish a molecular mechanism.'
 if component=='figure1':
  unit='Stable-Ensembl genes for 1B and eligible primary SuSiE genes/trait-support classes for 1C; descriptive counts, not donors.'
  family='1B RNA BH 24,115 and |log2FC|>0.5; genetics retains native SuSiE eligibility. 1C any eligible tier 1/2 SuSiE gene-trait pair. Count display adds no test.'
  alt='Joint eligibility and source phenotype definitions constrain overlap; low overlap does not establish antagonism or complementary mechanisms.'
 elif component=='figure2':
  unit='Gene-study signal pairs or gene-level descriptive counts, according to native source key; GWAS participants are source-specific.'
  family='Eligible tier 1/2 COLOC posterior evidence is retained; census/scatter/ancestry summaries add no new hypothesis family.'
  alt='Colocalization and annotation are inherited-trait evidence; neither specifies a causal gene nor connects disease-state direction.'
  if any(k in name for k in ['architecture','region_group','contrasts']):
   unit='233 physical sentinel clusters; 19 direct / 228 enzyme / 14 shared; shared sampled blocks carry both classes. This is not an asserted LD-independent locus census.'
   family='10,000 physical-cluster bootstrap draws, seed 20260817; 57 correlated sensitivity comparisons; no calibrated trait-class P value.'
   alt='Repeated trait/study ascertainment and physical grouping affect PIP allocation; coding mass is not heritability or gene-causality probability.'
 elif component=='figure3':
  unit='774 conservative genetic-identity clusters, 138 controls / 636 disease across five cohorts; these are not verified independent donors.'
  family='Canonical disease/control BH 24,115 and |log2FC|>0.5; stage/composition families are separately specified by native source.'
  alt='Cross-sectional source eligibility, cohort composition and cell mixture can change apparent disease-state abundance; not longitudinal change.'
  if 'fig3c' in name or 'cohort_deg' in name:
   unit='Gene-by-cohort membership census across five 24,115-gene families; descriptive ribbon widths/counts.';family='Each cohort BH 24,115 and |log2FC|>0.5; union membership census adds no replication P value.';alt='Cohort-specific observability and unequal disease/control eligibility influence shared calls; membership is not independent cohort confirmation.'
  if 'endpoint' in name:
   unit='392 overlapping identity clusters: 113 strict control / 279 advanced, four eligible cohorts.';family='Group-blind descriptive PCA; no endpoint significance or new hypothesis family.';alt='Cohort/sex-adjusted embedding can summarize histology contrast but is not independent validation or donor verification.'
  if 'fig3e' in name:
   family='Four F1–F4 versus F0 BH 24,115-gene families without fold-change filter; F0 contains 103 units including 89 disease biopsies.';alt='Reference disease biopsies and cell composition affect cross-sectional stage differences; counts are not canonical disease/control calls.'
  if 'fig3f' in name:
   unit='205/190/157/94 identity clusters in adjacent stage contrasts; 54 stage-selected primary genetic genes displayed.';family='Four complete 24,115-gene BH families; 54 display genes selected by stage significance, not a new enrichment test.';alt='Outcome-conditioned display and unsigned genetic support do not establish stage-specific inherited mechanisms.'
  if 'fig3g' in name:
   unit='3,776 finite RNA-gene/COLOC pairs; descriptive gene-level association.';family='Spearman nominal P=0.18; canonical RNA BH 24,115/|log2FC|>0.5; no directionally harmonized genetic effect.';alt='Unsigned PP.H4 cannot test risk-allele direction, antagonism or biological complementarity.'
 elif component=='figure4':
  unit='Participants/donors for source inferential models; 117 fixed programs and cells only where descriptive.';family='Source-specific complete prespecified program/system family; stage eligibility, coverage and robustness gates remain distinct.';alt='Cross-sectional stage, source ascertainment and cell mixture may drive score differences; program effects do not establish member-gene effects.'
  if 'source_diagnosis_correction' in text or 'hotspot_stage' in name:
   unit='Primary max 62 donors; ambient max 60; original Healthy/Steatosis/Steatohepatitis rule excludes two GSE174748 cirrhotic participants.';family='Complete 117 program BH family; classical and HC3 families separate. Only BICC1 primary-selected; zero HC3-supported. 116 tested-nonsignificant; zero strict tested-negative.';alt='Leverage and source ascertainment limit stability; nonsignificance is not adequate-test negative evidence.'
  if 'transport' in name or 'tissue_state' in name:
   unit='457 source records: 364 selected units plus 93 GSE174478 records, distinct from 774 unit disease fit.';family='224 tests = 112 eligible programs × fibrosis/NAS axes; composition adjusted arm separate.';alt='Source-record identity, cell composition and linear ordinal coding limit interpretation; associations are not temporal transport or nonlinear disease change.'
  if 'source_qualified_ccc' in text:
   unit='Biological donors matched by native source donor keys; sender/receiver cell counts are covariates, not replicates.';family='26 planned fixed-program/communication tests; complete-family classical and HC3 correction separately; full communication family retains its source-native gates.';alt='Contextual RNA covariation does not establish ligand secretion, receptor activation, or physical cell-cell signalling; source reuse and shared genes may affect association.'
  if 'window' in name or 'trajectory' in name:
   unit='Cohort-window mean fixed program scores, with underlying identity units; windows are not independent patients.';family='Descriptive nine fixed-projection windows; stage/sex-adjusted inferential families are separate companion analyses.';alt='Three construction cohorts overlap score derivation; trajectories do not provide independent support or temporal ordering.'
  if 'system' in name or 'geometry' in name:
   unit='43 fixed systems; 435 stage-complete identity units in two evaluation cohorts where used.';family='Descriptive system-level concordance separate from uniform size-matched, connected-region and strict maxT nulls.';alt='Shared memberships and spatial/system dependence can yield concordance; do not substitute one competitive null for another.'
 elif component=='figure5':
  unit='Assay-native participants; measured genes/programs are descriptive output rows, not biological replicates.';family='Source-native assay family; testability and source dependence retained. No common score or pooled cross-assay P value.';alt='Assay coverage, selection and reuse differ; molecular context does not establish causal mechanism or transcript origin.'
  if 'protein_triage' in name:
   unit='3,054 measured prioritized gene rows; protein participants remain native assay units.';family='Protein originating-assay BH family is retained; 3,054 display denominator is not asserted as its testing family. RNA uses 24,115-gene canonical gate.';alt='July-locked prioritized universe and same protein cohort condition the protein reestimates; bulk RNA comes from separate source cohorts and protein/RNA direction is not independent validation.'
  if 'mrna_protein' in name or 'histology_partial' in name:
   unit='25 fixed selected protein/gene rows; native protein participants underlie estimates.';family='RNA canonical BH 24,115 and |log2FC|>0.5; protein native BH; partial-correlation source family retained, no 25-row discovery test.';alt='Protein selection and protein re-estimation used the same protein cohort; bulk RNA comes from separate source cohorts. These descriptive associations do not independently validate the selected processes.'
  if 'variant_consistent' in name:
   unit='767 gene-study primary pairs, 28 direct / 739 enzyme; static peaks in two snATAC sources / five lineages.';family='Mapped posterior>=0.95 testability; posterior-overlap descriptive states, not new association P values.';alt='Static open-chromatin overlap localizes regulatory DNA but does not assign a target or demonstrate functional regulation.'
  if 'spatial' in name:
   unit='Native participant/section units; Vu arrays lack donor keys and are source-dependent; H&E image is illustrative.';family='Source-native BH-adjusted empirical Moran-excess tests using 9,999 draws; overlay adds no test.';alt='Spot counts are not donor replication, no pathologist region masks establish a named histologic structure, and program support is not IGFBP7/BICC1 gene support.'
 elif component=='figure6':
  unit='28,233 stable gene identities and two qualified gene evidence domains; program/assay contexts remain separate.';family='Descriptive r8 class/route census; no confidence ranking, new P values or gene-call expansion.';alt='Colocalization does not identify a credible-set target or mechanism; unsupported measured evidence differs from untestability.'
 return unit,family,alt
source_rows=[]; joins=[]; exact_rows=[]
units={'figure1':'participants for cohort summaries; genes and studies for descriptive censuses','figure2':'GWAS participants; credible sets/LD regions or gene-study pairs according to source','figure3':'774 genetic-identity units across five cohorts, not independently verified donors','figure4':'cohort participants/donors for inferential contrasts; programs/cells only descriptive','figure5':'assay-native participant units; genes/programs are descriptive or pairwise test units as qualified','figure6':'genes and assay contexts, descriptive Catalog classes'}
for record in inv.to_dict('records'):
 component=record['component']
 if component not in units:continue
 p=Path(record['source_path'])
 if p.suffix=='.pdf':continue
 assert h(p)==record['sha256'],p
 table=None
 if p.name.endswith(('.tsv','.csv','.tsv.gz')):
  try:table=pd.read_csv(p,sep='\t' if '.tsv' in p.name else ',',dtype=str,keep_default_na=False)
  except pd.errors.EmptyDataError:table=pd.DataFrame()
 key=''; keymode='not_gene_keyed'
 if table is not None and 'program_uid' in table:keymode='program_uid_membership_hash'
 elif table is not None and 'cell_type' in table and ('module' in table or 'program' in table):keymode='lineage_module_requires_registry_join'
 if table is not None:
  for candidate in ['ensembl_id','ensembl_id_base','bulk_gene_id','gene_id','ensembl']:
   if candidate in table and table[candidate].astype(str).str.startswith('ENSG').any():key=candidate;keymode='stable_ensembl';break
  if not key:
   for candidate in ['symbol','gene','gene_symbol']:
    if candidate in table and table[candidate].astype(str).isin(symbols).any():key=candidate;keymode='unique_symbol_only';break
 role='gene_rows_with_context' if key else ('aggregate_or_distribution' if table is not None else 'source_metadata')
 name=p.name.lower();claimrule='Descriptive aggregate census/distribution; source counts do not establish gene mechanism or tested-negative calls.'
 if key:claimrule='Source assay values retained exactly; Catalog class is joined context, not a reclassification of this assay or its effect threshold.'
 if 'contrast' in name or 'partial' in name or 'calibration' in name:claimrule='Source-native contrast and uncertainty; inferential unit and hypothesis family remain as specified by the source/legend. No mechanism is established.'
 reason='Applicable at panel/source level; no gene-specific call is created.' if not key else 'Exact native-ID or unique-symbol match to unchanged r8 calls; unmatched/ambiguous rows remain explicit.'
 alternative='Aggregate membership or an association does not connect inherited regulation to disease-state change in the same individual.'
 discrimination='Use only existing Catalog experiment rule for gene rows; aggregate-only summaries do not nominate a gene experiment.'
 qual='Current working source; not an independent validation of source findings.'
 if component=='figure4':qual+=' GSE174748 two cirrhotic participants excluded by unchanged three-stage eligibility; primary max 62 / ambient max 60, only BICC1 primary selected and zero HC3-supported; Liver Atlas excluded from the primary stage-association fits. Named r8 program context is historical where affected; use diagnosis correction companion.'
 if component=='figure5':qual+=' Protein significance and RNA interval-null calls are assay-specific; nonsignificance is not tested_negative. Protein reestimates are selection-conditioned on selection in the same protein cohort where source states so; bulk RNA comes from separate source cohorts.'
 native_unit,hypothesis_family,specific_alternative=panel_context(component,p)
 if table is not None:
  if 'biological_unit' in table:native_unit='; '.join(sorted(set(table.biological_unit.astype(str))-{''})) or native_unit
  if 'multiple_testing_family' in table:hypothesis_family='; '.join(sorted(set(table.multiple_testing_family.astype(str))-{''})) or hypothesis_family
  if 'unresolved_alternative' in table:specific_alternative='; '.join(sorted(set(table.unresolved_alternative.astype(str))-{''})) or specific_alternative
 native_context_columns=';'.join(c for c in ['evidence_state','state_reason','testability_state','testability_reason','biological_unit','multiple_testing_family','unresolved_alternative','claim_boundary','selection_conditioned','interpretation'] if table is not None and c in table)
 source_rows.append(dict(component=component,source_path=str(p),sha256=record['sha256'],source_role=role,native_context_columns=native_context_columns,native_rows=len(table) if table is not None else '',join_key=key,join_basis=keymode,biological_unit=native_unit,hypothesis_family=hypothesis_family,evidence_state_applicability='gene_context_only' if key else 'panel_descriptive_or_source_native',evidence_state_reason=reason,unresolved_alternative=specific_alternative,next_discrimination=discrimination,claim_rule=claimrule,threshold_and_family_basis='Exact source fields and docs/manuscript/05_figure_legends.md; no thresholds recomputed',qualification=qual))
 if table is None:continue
 for ix,r in enumerate(table.to_dict('records')):
  # Full original rows are preserved separately from appended interpretation.
  exact_rows.append(dict(source_path=str(p),source_sha256=record['sha256'],source_row_ordinal=ix+1,native_fields=r))
  if not key:continue
  value=str(r[key]);stable=value.split('.')[0]; matches=[lookup[stable]] if keymode=='stable_ensembl' and stable in lookup else (symbols.get(value,[]) if keymode=='unique_symbol_only' else [])
  matched=len(matches)==1;gr=matches[0] if matched else {};sid=str(gr.get('ensembl_id','')).split('.')[0];so=side.get(sid,{})
  joins.append(dict(source_path=str(p),source_sha256=record['sha256'],source_row_ordinal=ix+1,native_key=value,join_basis=keymode,join_status='exact_unique_match' if matched else 'unmatched_or_ambiguous',ensembl_id=sid,symbol=gr.get('symbol',''),primary_evidence_class=gr.get('primary_evidence_class',''),source_native_values_unchanged=True,call_context=ev.get(sid,[]),reason=gr.get('limitation','') if matched else 'No unique r8 mapping; no call inferred.',unresolved_alternative=so.get('unresolved_alternative','') if matched else '',next_discrimination_rule_id=so.get('recommended_experiment_rule_id',gr.get('next_experiment_rule_id','')) if matched else '',object_context=so))
def writejsonl(p,items):
 with p.open('w') as f:
  for r in items:f.write(json.dumps(r,default=str)+'\n')
with (O/'figure_source_context.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(source_rows[0]),delimiter='\t');w.writeheader();w.writerows(source_rows)
writejsonl(O/'source_rows_exact.jsonl',exact_rows);writejsonl(O/'gene_context_joins.jsonl',joins)
# Independent csv-module readback checks lexical values and row ordinals.
by_source={}
for rr in exact_rows:by_source.setdefault(rr['source_path'],[]).append(rr)
for path,rrs in by_source.items():
 p=Path(path); opener=gzip.open if p.name.endswith('.gz') else open
 with opener(p,'rt') as f:
  original=list(csv.DictReader(f,delimiter='\t' if '.tsv' in p.name else ','))
 assert len(original)==len(rrs),(p,len(original),len(rrs))
 for ordinal,(raw,rr) in enumerate(zip(original,rrs),1):
  raw={k:('' if v is None else v) for k,v in raw.items()}
  assert raw==rr['native_fields'] and rr['source_row_ordinal']==ordinal,(p,ordinal)
# Four worked entries retain exact per-layer calls and class; protein selection fields are untouched.
worked=[rr for rr in exact_rows if rr['source_path'].endswith('/fig6c_worked_entries.tsv')]
assert len(worked)==4
for rr in worked:
 r=rr['native_fields']; gr=symbols[r['symbol']];assert len(gr)==1
 so=side[gr[0]['ensembl_id']];assert r['primary_evidence_class']==gr[0]['primary_evidence_class']
 assert r['genetics_call']==so['genetics_call'] and r['rna_call']==so['rna_call']
selection_rows=[rr for rr in exact_rows if rr['source_path'].endswith('/composite_mrna_protein_corrected.csv')]
assert len(selection_rows)==25 and all(rr['native_fields']['selection_conditioned'].upper()=='TRUE' for rr in selection_rows)
assert h(gp)==gh and h(ep)==eh

registry=pd.read_csv(companion/'program_registry_v2.tsv',sep='\t',dtype=str,keep_default_na=False)
assert len(registry)==117 and registry.program_uid.is_unique and registry.membership_sha256.is_unique
assert registry.primary_selected.str.upper().eq('TRUE').sum()==1
assert registry.hc3_supported.str.upper().eq('TRUE').sum()==0
assert registry.tested_negative.str.upper().eq('TRUE').sum()==0
assert registry.tested_nonsignificant.str.upper().eq('TRUE').sum()==116
assert registry.loc[registry.primary_selected.str.upper().eq('TRUE'),'module_name'].str.contains('Ductular',case=False).all()
old=pd.read_parquet(C/'passport_program_index.parquet').astype(object)
oldctx=pd.read_parquet(C/'passport_program_context.parquet').astype(object)
assert len(old)==117 and old.program_uid.is_unique
joined=old[['program_uid','cell_type','membership_sha256']].merge(registry,on=['program_uid','cell_type','membership_sha256'],how='outer',indicator=True,validate='one_to_one')
assert joined['_merge'].eq('both').all(),joined['_merge'].value_counts().to_dict()
joined.drop(columns='_merge').to_csv(O/'r8_program_identity_current_context_join.tsv',sep='\t',index=False)
reg_uid={r['program_uid']:r for r in registry.to_dict('records')}
reg_pair={(r['cell_type'],r['module']):r for r in registry.to_dict('records')}
program_joins=[]
for rr in exact_rows:
 r=rr['native_fields'];uid=r.get('program_uid','');lineage=r.get('cell_type',r.get('lineage',''));module=r.get('module','')
 if not uid and not (lineage and module):continue
 pr=reg_uid.get(uid) if uid else reg_pair.get((lineage,module))
 if pr is not None and r.get('membership_sha256'):assert r['membership_sha256']==pr['membership_sha256']
 program_joins.append({'source_path':rr['source_path'],'source_sha256':rr['source_sha256'],'source_row_ordinal':rr['source_row_ordinal'],'join_basis':'program_uid' if uid else 'exact_lineage_module_to_registry_uid','join_status':'exact_unique_match' if pr else 'unmatched_no_context_inferred','native_program_uid':uid,'native_lineage':lineage,'native_module':module,'program_uid':pr['program_uid'] if pr else '', 'membership_sha256':pr['membership_sha256'] if pr else '', 'current_stage_context':pr if pr else {},'context_scope':'program_only_no_member_gene_call_expansion','stage_availability_and_coverage':'Independent source-native eligibility/coverage fields retained, not substituted by stage association.'})
writejsonl(O/'program_context_source_joins.jsonl',program_joins)
assert not e.call_state.eq('tested_negative').any()
context_audit={'native_index_has_primary_selected':'primary_selected' in old.columns,'native_index_has_source_stability':'source_stability' in old.columns,'native_context_rows':len(oldctx),'native_context_assays':oldctx.assay.value_counts().to_dict(),'native_context_source_releases':oldctx.source_release_id.value_counts().to_dict(),'program_uids_all_unchanged':True,'membership_hashes_all_unchanged':True,'strict_tested_negative_count':0,'tested_nonsignificant_count':116,'current_registry_sha256':h(companion/'program_registry_v2.tsv'),'current_primary_selected_count':int(registry.primary_selected.str.upper().eq('TRUE').sum()),'current_hc3_supported_count':int(registry.hc3_supported.str.upper().eq('TRUE').sum()),'native_context_policy':'Historical named r8 program context retained; current corrected stage-selection/effects and source-stability available only in explicit companion. No member-gene call expansion.'}
(O/'program_context_audit.json').write_text(json.dumps(context_audit,indent=2)+'\n')
(O/'native_r8_program_index_schema.json').write_text(json.dumps(list(old.columns))+'\n')
(O/'native_r8_program_context_schema.json').write_text(json.dumps(list(oldctx.columns))+'\n')
# Preserve complete historical r8 rows and current source rows for explicit external review.
writejsonl(O/'r8_program_index_historical.jsonl',old.where(pd.notna(old),None).to_dict('records'))
writejsonl(O/'r8_program_context_historical.jsonl',oldctx.where(pd.notna(oldctx),None).to_dict('records'))
registry.to_csv(O/'current_program_context_correction_companion.tsv',sep='\t',index=False)
previous=R/'figures/candidates/figure-claim-context-20261014T164442Z_22237913/current_program_context_correction_companion.tsv'
prior=pd.read_csv(previous,sep='\t',dtype=str,keep_default_na=False)
annotation_fields={'tested_negative','source_native_tested_negative','tested_nonsignificant','registry_state','scoring_direction_agree','selected_hc3_fragile'}
unchanged_fields=[c for c in prior.columns if c in registry and c not in annotation_fields]
assert prior[unchanged_fields].equals(registry[unchanged_fields]),'Numerical/membership/source stability fields changed'
derived_direction=(pd.to_numeric(registry.primary_beta,errors='coerce').gt(0)&pd.to_numeric(registry.equal_run_beta,errors='coerce').gt(0))|(pd.to_numeric(registry.primary_beta,errors='coerce').lt(0)&pd.to_numeric(registry.equal_run_beta,errors='coerce').lt(0))
assert registry.scoring_direction_agree.str.upper().eq('TRUE').equals(derived_direction)
derived_fragile=registry.primary_selected.str.upper().eq('TRUE') & ~registry.hc3_supported.str.upper().eq('TRUE')
assert registry.selected_hc3_fragile.str.upper().eq('TRUE').equals(derived_fragile)
context_audit['changed_derived_annotations']={c:int(prior[c].ne(registry[c]).sum()) for c in prior.columns if c in registry and prior[c].ne(registry[c]).any()}
context_audit['derived_direction_and_fragility_independently_verified']=True
context_audit['numerical_membership_fields_unchanged']=unchanged_fields
context_audit['prior_companion_sha256']=h(previous)
context_audit['annotation_only_change_verified']=True
(O/'program_context_audit.json').write_text(json.dumps(context_audit,indent=2)+'\n')
summary=dict(candidate=True,official_release_state_changed=False,source_files=len(source_rows),source_rows_preserved=len(exact_rows),gene_keyed_rows=len(joins),unique_matched=sum(j['join_status']=='exact_unique_match' for j in joins),unmatched_or_ambiguous=sum(j['join_status']!='exact_unique_match' for j in joins),catalog_gene_rows=len(g),catalog_class_counts=g.primary_evidence_class.value_counts().to_dict(),catalog_gene_sha256=gh,catalog_evidence_sha256=eh,catalog_calls_unchanged=True,program_source_rows=len(program_joins),program_source_unique_matched=sum(j['join_status']=='exact_unique_match' for j in program_joins),program_source_unmatched=sum(j['join_status']!='exact_unique_match' for j in program_joins),tested_negative_calls=0,independent_csv_readback='all source values and ordinals pass',worked_6c_per_layer_calls_verified=4,selection_conditioned_5c_rows_verified=25,allowed_call_states=sorted(allowed),GSE174748_qualification='Source category corrected: GSE174748 two cirrhotic participants excluded by unchanged three-stage eligibility; primary max 62 / ambient max 60; only BICC1 primary selected and zero HC3-supported',limitations=['Source native fields and Catalog classes remain distinct; no inference from nonsignificance.','Aggregate sources carry panel-level claim rules and no false gene-specific mechanism calls.','This candidate adds context only; no new assay confidence ranks, targets, or tested negatives.'])
(O/'checks.json').write_text(json.dumps(summary,indent=2)+'\n');(O/'executed_collector.py').write_bytes(Path(__file__).read_bytes());print(json.dumps(dict(path=str(O),**summary),indent=2))
