#!/usr/bin/env python3
"""Provide local access to unchanged selected calls and molecular-object records."""
import argparse
import html
import json
import os
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reviewed-bundle',type=Path,required=True)
    p.add_argument('--local-interface',type=Path,required=True)
    p.add_argument('--molecular-objects',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    seal=a.reviewed_bundle/'PASS06_VALIDATED';assert seal.is_file()
    metadata=json.loads((a.molecular_objects/'sidecar_completion.json').read_text())
    assert metadata['object_rows']==28233 and metadata['bulk_analysis_unit_count']==774 and metadata['n_held_untestable_annotations']==855
    links=[('Browse the complete Gene Catalog',a.local_interface/'index.html'),('View worked gene entries and molecular objects',a.molecular_objects/'index.html'),('Download gene calls (Parquet)',a.reviewed_bundle/'passport_gene_index.parquet'),('Download separate assay evidence (Parquet)',a.reviewed_bundle/'passport_evidence_long.parquet'),('Download molecular-object records (JSONL)',a.molecular_objects/'molecular_object_sidecar.jsonl'),('Read the gene data dictionary (TSV)',a.reviewed_bundle/'passport_data_dictionary.tsv'),('Read source references and relationships (JSON)',a.reviewed_bundle/'portal_export/passport_source_graph.json'),('Read molecular-object annotation and counting sources (JSON)',a.molecular_objects/'sidecar_completion.json')]
    assert all(path.is_file() for _,path in links)
    page='<meta charset="utf-8"><title>MASLD Gene Catalog</title><style>body{font:16px sans-serif;max-width:900px;margin:3rem auto;line-height:1.5}li{margin:1rem 0}h1{font-size:26px}</style><h1>MASLD Gene Catalog</h1><p>28,233 gene identities. Separate inherited-trait and disease-state evidence, explicit assay eligibility, and proposed discriminating experiments.</p><p>Bulk estimates use 774 conservative genetic-identity clusters, not verified independent donors. Counting retains 855 annotations as untestable. Program and assay context does not create member-gene support, and gene-level colocalization does not establish causality or a signal-to-target link.</p><ul>'
    page+=''.join('<li><a href="'+html.escape(os.path.relpath(path.resolve(),out))+'">'+html.escape(label)+'</a></li>' for label,path in links)
    page+='</ul><p>The selected r8 gene calls and memberships are unchanged. Historical class-assay results remain historical. This interface is local; no external portal has been published.</p>'
    (out/'index.html').write_text(page)
    (out/'source_links.json').write_text(json.dumps({'status':'PASS','all_local_links_exist':len(links),'gene_rows':28233,'selected_source':'gene-catalog-corrected-r8-20261009T142027Z','external_publication':False},indent=2)+'\n')
    print('Local access page:',out/'index.html','8/8 links PASS')


if __name__=='__main__':main()
