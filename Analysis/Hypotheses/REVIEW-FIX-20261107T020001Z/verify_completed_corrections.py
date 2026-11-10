from pathlib import Path
import json,csv,hashlib,subprocess,re,urllib.parse
root=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
run=root/'Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z'
checks={}
def read(rel):return json.loads((run/rel).read_text())
assert read('singlecell/verification.json')['passed'];checks['singlecell_independent_families']=45
changes=list(csv.DictReader((run/'singlecell/result_changes.tsv').open(),delimiter='\t'))
h=next(r for r in changes if r['arm']=='binary' and r['analysis']=='Hepatocytes')
assert [int(h[k]) for k in ['new_donors','new_tested','new_calls','retained_calls','lost_calls','gained_calls']]==[71,29948,90,70,9,20]
checks['hepatocyte_counts_match_manuscript']=True
assert read('genetics/xenium_independent_check.json')['status']=='PASS';checks['xenium_independent_scalar_reconstructions']=96
b=read('bulk/eisa_adjusted_v3/summary.json');assert (b['n_genes'],b['ratio_calls'],b['ratio_calls_abs_log2_ge_0_5'])==(12010,2711,114)
assert abs(b['p_global_call_count']-1/201)<1e-14 and b['null_max']==149 and b['first_null_max_effect_difference']>1
checks['eisa_adjusted_counts_and_valid_null']=True
c=read('catalog_models/outputs/corrected_companion_v2/checks.json');assert c['status']=='PASS' and c['rna_applicability_corrected']==4218 and c['gene_records']==28233
checks['catalog_preserved_named_classes_and_corrected_applicability']=True
c=read('catalog_models/outputs/eligibility_candidate/checks.json');assert c['class_changes']==180 and c['genetics_state_changes']==188 and c['candidate_not_adopted']
assert sum(c['candidate_class_counts'].values())==28233
checks['eligibility_candidate_separate']=True
checks['genetics_independent']=read('genetics/checked_companion_v2/independent_checks.json')
for r in read('working_artwork_install.json'):
 assert hashlib.sha256((root/r['target']).read_bytes()).hexdigest()==r['candidate_sha256']
checks['installed_artifacts_match_checked_candidates']=5
for rel,words in [
 ('figures/main/fig3_bulk_transcriptomics/review_assembly/Figure3.pdf',['SteatoSITE','1,267','3,149']),
 ('figures/main/fig3_bulk_transcriptomics/panels/supplementary/figs3_steatosite_exon_vs_intron.pdf',['114 genes','12,010']),
 ('figures/main/fig7_model_framework/fig7E_verdict_vocabulary_and_refusal_thresholds.pdf',['untestable','6 vs 16'])]:
 text=subprocess.check_output(['pdftotext',str(root/rel),'-'],text=True)
 assert all(w in text for w in words),(rel,words)
 if 'fig7E' in rel:assert 'tested negative' not in text
checks['installed_pdf_text_matches_claims']=True
# Inspect new links only; do not repair unrelated historical link debt.
links=re.compile(r'\[[^\]]*\]\(([^)]+)\)',re.S)
verified=[];fail=[]
for p in [*root.glob('docs/*.md'),*root.glob('docs/manuscript/draft/*.md'),root/'docs/manuscript/METHODS.md',root/'docs/manuscript/05_figure_legends.md',root/'docs/manuscript/README.md',root/'figures/main/INDEX.md',root/'figures/main/fig3_bulk_transcriptomics/README.md',root/'figures/main/fig7_model_framework/README.md',run/'REPORT.md']:
 rel=p.relative_to(root);bp=run/'baseline'/rel
 before=set(links.findall(bp.read_text())) if bp.exists() else set()
 for link in set(links.findall(p.read_text()))-before:
  if re.match(r'^[a-zA-Z][a-zA-Z0-9+.-]*:',link) or link.startswith('#'):continue
  dest=urllib.parse.unquote(link.split('#',1)[0]).strip('<>')
  target=(p.parent/dest).resolve()
  (verified if target.exists() else fail).append({'file':str(rel),'link':link})
assert not fail,fail
checks['new_local_links_resolve']=len(verified)
scope=subprocess.run(['python3','scripts/manuscript/validate_resource_scope.py'],cwd=root,text=True,capture_output=True)
(run/'scope_check.txt').write_text(scope.stdout+scope.stderr);assert scope.returncode==0,scope.stdout+scope.stderr
checks['resource_scope']='PASS'
checkpaths=['docs/README.md','docs/PAPER.md','docs/STATUS.md','docs/RESULTS.md','docs/ROADMAP.md','docs/manuscript/METHODS.md','docs/manuscript/05_figure_legends.md','docs/manuscript/draft/fig2.md','docs/manuscript/draft/fig3.md','docs/manuscript/draft/fig5.md','docs/manuscript/draft/fig5_discussion.md','docs/manuscript/draft/fig7.md','docs/manuscript/draft/supplementary_notes.md','scripts/portal/build_resource_catalog_sidecar.py','scripts/portal/build_evidence_passport_bundle.py','Analysis/SingleCell/scripts/pseudobulk_de.R','Analysis/SingleCell/scripts/run_pseudobulk_de_coarse_stages.R','Analysis/SingleCell/scripts/hotspot_modules/501_run_hotspot.py']
diff=subprocess.run(['git','diff','--check','--',*checkpaths],cwd=root,text=True,capture_output=True)
assert diff.returncode==0,diff.stdout+diff.stderr;checks['scoped_whitespace']='PASS'
checks['status']='PASS';(run/'final_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
(run/'new_link_checks.json').write_text(json.dumps(verified,indent=2)+'\n')
print(json.dumps({k:v for k,v in checks.items() if k!='genetics_independent'},indent=2))
