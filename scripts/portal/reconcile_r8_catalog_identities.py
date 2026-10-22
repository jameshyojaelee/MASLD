#!/usr/bin/env python3
"""Explain categorical, native stable-ID and adjudicated Catalog denominators."""
import argparse
import json
from pathlib import Path
import pandas as pd


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();root=a.project_root;a.output.mkdir(parents=True,exist_ok=False)
    source=root/'RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z'
    raw=pd.read_csv(root/'Analysis/Multimodal_Program_Projection/candidates/gen-corrected-r8-20261009T142027Z/genetics_context/frozen_evidence_classes.tsv',sep='\t',dtype=str,keep_default_na=False)
    adj=pd.read_csv(source/'preflight/gen_identity_adjudication_v1/gen_frozen_classes_ensembl_adjudicated.tsv',sep='\t',dtype=str,keep_default_na=False)
    q=pd.read_csv(source/'preflight/gen_identity_adjudication_v1/gen_identity_quarantine.tsv',sep='\t',dtype=str,keep_default_na=False)
    g=pd.read_parquet(source/'passport_gene_index.parquet')
    d=pd.read_csv(root/'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv')
    c=pd.read_csv(root/'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv')
    d['stable_id']=d.gene.str.split('.').str[0];c['stable_id']=c.ensembl.str.split('.').str[0]
    interface=set(d.stable_id)&set(c.stable_id)
    native_symbols=set(d[d.stable_id.isin(interface)].symbol.astype(str))
    raw_joint=raw[raw.joint_testable=='true'];adj_joint=adj[adj.joint_testable=='true']
    aliases=raw_joint[~raw_joint.gene_symbol.isin(native_symbols)].copy()
    annotation=pd.read_csv(root/'RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11/annotation/GENCODE_v49/run_19749483/gene_identity.tsv',sep='\t',dtype=str,keep_default_na=False)
    canonical=annotation[annotation.is_canonical_chromosome=='true']
    known=set(canonical.gene_id_base)
    records=[]
    for row in aliases.itertuples(index=False):
        aa=adj[adj.gene_symbol==row.gene_symbol];qq=q[q.source_symbol==row.gene_symbol]
        records.append({'symbol':row.gene_symbol,'ensembl_bulk':row.ensembl_bulk,'ensembl_genetic':row.ensembl_genetic,'bulk_identity_in_v49':row.ensembl_bulk in known,'genetic_identity_in_v49':row.ensembl_genetic in known,'catalog_id':';'.join(aa.adjudicated_ensembl_id),'identity_resolution_rule':';'.join(aa.identity_resolution_rule),'quarantine_reason':';'.join(qq.reason_code),'raw_primary_genetic':row.primary_genetic,'raw_state_supported':row.established_state_associated,'catalog_class':';'.join(g[g.ensembl_id.isin(aa.adjudicated_ensembl_id)].primary_evidence_class)})
    pd.DataFrame(records).to_csv(a.output/'source_symbol_false_joint_links.tsv',sep='\t',index=False)
    # The excluded native rows are retained as quarantined source records.
    quarantined=raw[raw.gene_symbol.isin(q.source_symbol)].copy()
    quarantined=quarantined.merge(q[['source_symbol','reason_code']],left_on='gene_symbol',right_on='source_symbol',how='left')
    quarantined.to_csv(a.output/'quarantined_source_classes.tsv',sep='\t',index=False)
    catalog_joint=set(g[g.primary_evidence_class!='untested'].ensembl_id)
    diff=[]
    for label,ids in [('catalog_only',catalog_joint-interface),('native_interface_only',interface-catalog_joint)]:
        for ens in sorted(ids):
            gr=g[g.ensembl_id==ens];dd=d[d.stable_id==ens];aa=adj[adj.adjudicated_ensembl_id==ens]
            diff.append({'universe':label,'ensembl_id':ens,'catalog_symbol':';'.join(gr.symbol.astype(str)),'bulk_symbol':';'.join(dd.symbol.astype(str)),'adjudication_rule':';'.join(aa.identity_resolution_rule.astype(str)),'catalog_class':';'.join(gr.primary_evidence_class.astype(str))})
    pd.DataFrame(diff).to_csv(a.output/'adjudicated_vs_native_stable_IDs.tsv',sep='\t',index=False)
    qjoint=quarantined[quarantined.joint_testable=='true']
    summary={'raw_categorical_joint_symbols':len(raw_joint),'native_stable_ID_interface':len(interface),'native_interface_distinct_source_symbols':len(native_symbols),'adjudicated_catalog_joint_IDs':len(catalog_joint),'source_symbol_only_categorical':len(aliases),'source_symbol_only_native':len(native_symbols-set(raw_joint.gene_symbol)),'GEN_quarantined_source_rows':len(q),'GEN_quarantined_joint_rows':len(qjoint),'GEN_quarantined_positive_genetic':int((quarantined.primary_genetic=='true').sum()),'GEN_quarantined_positive_state':int((quarantined.established_state_associated=='true').sum()),'GEN_quarantined_joint_positive_genetic':int((qjoint.primary_genetic=='true').sum()),'GEN_quarantined_joint_positive_state':int((qjoint.established_state_associated=='true').sum()),'catalog_only_native_IDs':len(catalog_joint-interface),'native_only_catalog_IDs':len(interface-catalog_joint),'interpretation':'14720 is the pre-adjudication symbol-grouped universe.14704 is a same-stable-ID bulk/COLOC join.14682 is the Catalog after GENCODE identity adjudication. Native source-symbol comparison is16categorical-only/0interface-only; the previously quoted24/8 cannot be reproduced from these exact source-symbol keys. Multi-ID text must never be treated as one stable gene ID.'}
    assert len(raw_joint)==14720 and len(interface)==14704 and len(catalog_joint)==14682 and len(q)==59
    assert summary['GEN_quarantined_joint_positive_genetic']==0 and summary['GEN_quarantined_joint_positive_state']==0
    summary['adjudicated_source_joint_rows']=len(adj_joint)
    summary['joint_row_change_from_identity_adjudication']=len(raw_joint)-len(adj_joint)
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))


if __name__=='__main__':main()
