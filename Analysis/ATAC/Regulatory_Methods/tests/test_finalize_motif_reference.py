import csv,gzip,hashlib,importlib.util,json,tempfile
from pathlib import Path
from unittest.mock import patch
source=Path(__file__).resolve().parents[1]/'scripts/finalize_candidate.py'
spec=importlib.util.spec_from_file_location('finalize',source); f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)
def write_table(p,rows):
 with p.open('w') as h:
  w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
with tempfile.TemporaryDirectory() as tmp:
 root=Path(tmp); execution=root/'execution'; motif=root/'motifs';motif.mkdir(); database=root/'database';database.write_bytes(b'fixture database identity')
 record=dict(database_path=str(database),database_sha256=hashlib.sha256(database.read_bytes()).hexdigest(),jaspar_package='fixture');write_table(motif/'jaspar_database.tsv',[record])
 networks=[];sources=[];coverage=[]
 for cohort in ['GSE244832','GSE281367','GSE296875']:
  for lineage in ['hepatocyte','stellate']:
   network=execution/'networks_full'/f'{cohort}_{lineage}'/'analysis';directory=network/f'{cohort}_{lineage}';directory.mkdir(parents=True);networks.append(network)
   complete=dict(status='completed',pilot=False,interpretation='fixture');(network/'completion.json').write_text(json.dumps(complete));(network/'run_settings.json').write_text('{}')
   write_table(network/'comparisons.tsv',[dict(cohort=cohort,lineage=lineage,lineage_state='source_label',source_lineage=lineage)])
   write_table(directory/'links.tsv',[dict(outer_fold=0,gene_id='gene',peak_id='chr1:1-2',method='distance',rank=1,score=1,score_unit='distance',evidence_state='link',endpoint='rna')])
   common=dict(network_root=str(network),source_scope='development_candidate',source_interpretation='fixture');sources.append(dict(**common,completion=complete))
   available=cohort!='GSE296875'
   if available:
    d=motif/'motif_matches'/cohort/lineage;d.mkdir(parents=True)
    for name in ['matches.mtx','peaks.tsv','motifs.tsv']:(d/name).write_text('fixture')
   coverage.append(dict(**common,cohort=cohort,lineage=lineage,lineage_state='source_label',source_lineage=lineage,endpoint_recorded='True',status='annotated' if available else 'untestable',reason='' if available else 'whole_native_peak_universe_motifscan_not_run_125_peak_pilot_not_global',n_link_records=1,n_family_gene_rows=1 if available else 0,peak_matching='exact_ID_only_no_overlap_borrowing'))
 annotation=root/'annotation';annotation.mkdir();completion=dict(status='completed',exact_ID_only=True,motif_root=str(motif),jaspar_source=record,sources=sources);(annotation/'completion.json').write_text(json.dumps(completion));write_table(annotation/'coverage.tsv',coverage)
 with gzip.open(annotation/'motif_family_gene.tsv.gz','wt') as h:h.write('\t'.join(['network_root','cohort','lineage','gene_id','motif_family','score_unit','evidence_state','tf_activity_inferred','tf_binding_inferred','causal_regulation_inferred','jaspar_database_sha256'])+'\n')
 assert f.validate_motif_annotation(annotation,networks,motif)['annotation']==str(annotation)
 original=f.validate_motif_annotation
 with patch.object(f,'validate_motif_annotation',side_effect=lambda a,n,m:original(a,n,motif)):
  f.assemble(execution,root/'assembled',source.parent,Path('/unused'),False,True,annotation)
 assert (root/'assembled/motif_family_compatibility').is_symlink()
 assert (root/'assembled/motif_family_compatibility').resolve()==annotation
 assert json.loads((root/'assembled/motif_annotation_reference.json').read_text())['completion']==completion
 def rejected():
  try:f.validate_motif_annotation(annotation,networks,motif)
  except (ValueError,FileNotFoundError):return
  raise AssertionError('Invalid annotation accepted')
 completion['sources'][0]['network_root']=str(root/'wrongsource');(annotation/'completion.json').write_text(json.dumps(completion));rejected();completion['sources'][0]['network_root']=str(networks[0])
 completion['status']='running';(annotation/'completion.json').write_text(json.dumps(completion));rejected();completion['status']='completed';(annotation/'completion.json').write_text(json.dumps(completion))
 (annotation/'INVALIDATED_fixture.json').write_text('{}');rejected();(annotation/'INVALIDATED_fixture.json').unlink()
 coverage[-1]['status']='annotated';write_table(annotation/'coverage.tsv',coverage);rejected();coverage[-1]['status']='untestable';write_table(annotation/'coverage.tsv',coverage)
 database.write_bytes(b'changed database');rejected();database.write_bytes(b'fixture database identity')
 (annotation/'motif_family_gene.tsv.gz').unlink();rejected()
 print('PASS: exact six-source reuse/symlink and preserved completion receipt; wrong source, incomplete, invalidated, wrong native evidence state, changed frozen database and missing table rejected')
