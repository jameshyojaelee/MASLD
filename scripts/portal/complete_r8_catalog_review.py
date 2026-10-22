#!/usr/bin/env python3
"""Prepare a versioned local r8 review and molecular-object companion.

Never modifies the selected October 9 bundle or changes a gene call. Historical
class-level assays remain historical; fixed 117-program memberships stay fixed.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import pandas as pd


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();root=args.project_root.resolve();out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    source=root/'RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z'
    bulk=root/'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv'
    coloc=root/'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv'
    raw=root/'Analysis/Multimodal_Program_Projection/candidates/gen-corrected-r8-20261009T142027Z/genetics_context/frozen_evidence_classes.tsv'
    annotation=root/'figures/candidates/resource-corrections-20261005-v1/annotation_resolution_r8'
    certificate=root/'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/qualified_input_certificate.json'
    counts=root/'figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_20261007/assembled_774_primary/counts_raw.rds'
    identity=root/'RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11/annotation/GENCODE_v49/run_19749483'
    scripts=out/'producers';scripts.mkdir()
    for name in ['complete_r8_catalog_review.py','build_resource_catalog_sidecar.py','gene_catalog_v2.py','build_evidence_passport_ui.py','render_evidence_passport_ui.py','run_r8_catalog_completion.sbatch']:
        shutil.copyfile(root/'scripts/portal'/name,scripts/name)
    tracked=[p for p in source.rglob('*') if p.is_file()]+[bulk,coloc,raw,certificate,counts,annotation/'annotation_correspondence.tsv',annotation/'resolution_summary.json',identity/'gene_identity.tsv',identity/'lncrna_genomic_class.tsv']
    before={str(p):sha(p) for p in tracked}
    (out/'input_checksums_before.json').write_text(json.dumps(before,indent=2)+'\n')
    overlay=out/'catalog_local';shutil.copytree(source,overlay)
    genes=pd.read_parquet(source/'passport_gene_index.parquet');evidence=pd.read_parquet(source/'passport_evidence_long.parquet')
    expected={'genetically_anchored':400,'concordant':26,'established_state_associated':1129,'unresolved':13127,'untested':13551}
    assert genes.primary_evidence_class.value_counts().to_dict()==expected
    assert len(genes)==28233 and genes.ensembl_id.is_unique and len(evidence)==56466
    assert not (evidence.call_state=='tested_negative').any()
    # Keep categorical classes and source evidence untouched. Clarify current
    # counting/identity restrictions only in the local presentation companion.
    qualification=' Bulk uses 774 conservative genetic-identity clusters, not verified independent donors. Counting holds 855 annotations untestable; no held annotation is tested. Gene-level COLOC does not establish causality or a signal-to-target link.'
    altered=genes.copy();altered['limitation']=altered.limitation.fillna('').astype(str)+qualification
    altered.to_parquet(overlay/'passport_gene_index.parquet',index=False)
    assert altered.drop(columns='limitation').equals(genes.drop(columns='limitation'))
    assay=pd.read_csv(source/'passport_assay_status.tsv',sep='\t',keep_default_na=False)
    stale=assay.source_input_row_id.str.startswith('MYOJIN_CLASS:')
    # The 28 historical class tests are retained without recomputation and
    # without allowing their statistics to alter current gene classes.
    assert int(stale.sum())==28
    assay.loc[stale,'allowed_wording']='Historical evidence-class assay result; not a test of the selected r8 classes and not a gene-level negative.'
    assay.to_csv(overlay/'passport_assay_status.tsv',sep='\t',index=False)
    dictionary=pd.read_csv(source/'passport_data_dictionary.tsv',sep='\t',keep_default_na=False)
    for column,description in {
        'primary_evidence_class':'Selected r8 descriptive interface class. Untested means not jointly eligible and can coexist with support in one measured assay. Classes do not establish a causal mechanism.',
        'limitation':'Current scientific restriction. Bulk units are conservative genetic-identity clusters, not verified independent donors; 855 held annotations remain untestable.',
    }.items():
        mask=(dictionary.table_name=='passport_gene_index.parquet')&(dictionary.column_name==column)
        assert int(mask.sum())==1;dictionary.loc[mask,'description']=description
    dictionary.to_csv(overlay/'passport_data_dictionary.tsv',sep='\t',index=False)
    d=pd.read_csv(bulk);d['stable_id']=d.gene.str.split('.').str[0]
    c=pd.read_csv(coloc);c['stable_id']=c.ensembl.str.split('.').str[0]
    r=pd.read_csv(raw,sep='\t',keep_default_na=False,dtype=str)
    joint=r[r.joint_testable.str.lower()=='true']
    interface=set(d.stable_id)&set(c.stable_id)
    interface_bulk=d[d.stable_id.isin(interface)]
    categorical=set(joint.gene_symbol);interface_symbols=set(interface_bulk.symbol.astype(str))
    # A symbol-grouped categorical universe and a stable-ID interface cannot
    # share a denominator silently; retain every differing identity explicitly.
    rows=[]
    quarantine=pd.read_csv(source/'passport_identity_quarantine.tsv',sep='\t',dtype=str,keep_default_na=False)
    for universe,symbols in [('categorical_only',categorical-interface_symbols),('interface_only',interface_symbols-categorical)]:
        for symbol in sorted(symbols):
            rr=r[r.gene_symbol==symbol];dd=d[d.symbol.astype(str)==symbol]
            cc=c[c.stable_id.isin(dd.stable_id)];qq=quarantine[quarantine.source_symbol==symbol]
            ids=';'.join(dd.stable_id.astype(str).unique())
            rows.append({'universe':universe,'symbol':symbol,'bulk_stable_ids':ids,'categorical_bulk_ids':';'.join(rr.ensembl_bulk.astype(str).unique()),'categorical_genetic_ids':';'.join(rr.ensembl_genetic.astype(str).unique()),'bulk_rows':len(dd),'coloc_same_stable_id_rows':len(cc),'categorical_rows':len(rr),'joint_categorical_rows':int((rr.joint_testable.str.lower()=='true').sum()),'catalog_quarantine_reason':';'.join(qq.reason_code.astype(str).unique()),'interpretation':'categorical source links differing stable IDs by symbol' if universe=='categorical_only' else 'stable-ID join survives a source-symbol mismatch'})
    pd.DataFrame(rows).to_csv(out/'overlap_universe_difference.tsv',sep='\t',index=False)
    counts_summary={'categorical_source_joint_rows':len(joint),'categorical_unique_symbols':len(categorical),'stable_ID_interface':len(interface),'stable_ID_interface_unique_symbols':len(interface_symbols),'symbol_only_categorical':len(categorical-interface_symbols),'symbol_only_interface':len(interface_symbols-categorical),'catalog_adjudicated_joint_classes':sum(expected.values())-expected['untested'],'source_raw_rows':len(r),'symbol_collisions_raw':int(r.gene_symbol.duplicated().sum()),'interpretation':'The categorical source joins display symbols before GENCODE identity adjudication; the interface joins version-stripped Ensembl IDs. Use 14,704 for stable-ID interface overlap; 14,720 is the symbol-grouped historical categorical denominator. r8 catalog classes retain only adjudicated identities.'}
    assert len(joint)==14720 and len(interface)==14704
    (out/'overlap_universe_summary.json').write_text(json.dumps(counts_summary,indent=2)+'\n')
    # Independently recompute all stored 20 spot checks from source-native
    # thresholds. No RNA evidence is upgraded solely by program membership.
    spot=pd.read_csv(root/'Analysis/Hypotheses/r8-review-20261011T150915Z/spot_check_20.tsv',sep='\t')
    checks=[]
    for row in spot.itertuples(index=False):
        gr=genes[genes.ensembl_id==row.ensembl_id].iloc[0];dd=d[d.stable_id==row.ensembl_id];cc=c[c.stable_id==row.ensembl_id]
        gen=bool(len(cc) and (pd.to_numeric(cc.coloc_best_susie_pp4,errors='coerce')>=.5).any())
        state=bool(len(dd) and ((dd.padj<.05)&(dd.logFC.abs()>.5)).any())
        value='untested' if not(len(dd) and len(cc)) else 'concordant' if gen and state else 'genetically_anchored' if gen else 'established_state_associated' if state else 'unresolved'
        checks.append({'symbol':row.symbol,'ensembl_id':row.ensembl_id,'stored_class':gr.primary_evidence_class,'rederived_class':value,'pass':value==gr.primary_evidence_class})
    assert all(r['pass'] for r in checks)
    pd.DataFrame(checks).to_csv(out/'independent_spot_checks.tsv',sep='\t',index=False)
    for name in ['passport_program_index.parquet','passport_program_membership.parquet']:
        assert sha(source/name)==sha(overlay/name)
    subprocess.run([sys.executable,str(scripts/'build_resource_catalog_sidecar.py'),'--bundle',str(overlay),'--gene-identity',str(identity/'gene_identity.tsv'),'--lncrna-class',str(identity/'lncrna_genomic_class.tsv'),'--output',str(out/'molecular_objects'),'--count-certificate',str(certificate),'--counts',str(counts),'--annotation-correspondence',str(annotation/'annotation_correspondence.tsv')],check=True)
    sys.path.insert(0,str(scripts))
    import build_evidence_passport_ui as ui
    ui.build_ui(overlay,fixture_mode=False)
    subprocess.run([sys.executable,str(scripts/'render_evidence_passport_ui.py'),'--bundle',str(overlay)],check=True)
    meta=json.loads((out/'molecular_objects/sidecar_completion.json').read_text())
    assert meta['primary_class_counts']==expected and meta['object_rows']==28233 and meta['bulk_analysis_unit_count']==774 and meta['n_held_untestable_annotations']==855
    records=[json.loads(line) for line in (out/'molecular_objects/molecular_object_sidecar.jsonl').open()]
    assert all(x['target_link_basis']=='none' and x['credible_set_to_gene_link_status']=='no_authoritative_link' for x in records)
    assert all(x['credible_set_coding_pip_mass'] is None for x in records)
    after={str(p):sha(p) for p in tracked};assert before==after
    (out/'input_checksums_after.json').write_text(json.dumps(after,indent=2)+'\n')
    summary={'status':'PASS','source_bundle':str(source),'local_overlay':str(overlay),'source_modified':False,'external_publication':False,'class_counts':expected,'gene_rows':28233,'evidence_rows':56466,'program_rows':117,'independent_spot_checks_passed':len(checks),'historical_class_assay_rows':int(stale.sum()),'source_files_unchanged':len(before),'source_qualification':'774 conservative genetic-identity clusters; 855 held-untestable annotations; unresolved GSE193066 not upgraded','manual_attestation_created':False,'limits':'No gene-level evidence inferred from program/spatial/static overlap; no causal claim; historical class-assay tests not rerun; local presentation overlay does not rewrite the historical selected bundle.'}
    (out/'review_completion.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))


if __name__=='__main__':main()
