#!/usr/bin/env python3
"""Add molecular-object records beside a v1 Catalog without altering its calls.

Provisional review is allowed only with an explicit flag. Final candidates
require the count-and-annotation certificate. Neither mode signs a manual
attestation, publishes a portal, or changes the source passport_* contract.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import html
import json
from pathlib import Path
import pandas as pd
import gene_catalog_v2 as catalog


def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as handle:
        for b in iter(lambda:handle.read(1024*1024),b""):h.update(b)
    return h.hexdigest()


def table(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter="\t"))


def object_record(gene,identity,lncrna,evidence,release,provisional):
    known=identity is not None
    biotype="unknown" if not known else (identity["gene_type"] if identity["gene_type"] in {"protein_coding","lncRNA"} else "other")
    genomic=lncrna.get("lncrna_genomic_class","unknown") if lncrna else "unknown"
    mapping="unique_stable_ensembl" if known and identity["base_id_mapping_status"]=="unique" else "mapping_ambiguous"
    evidence_class=gene["primary_evidence_class"]
    genetic=evidence.get("genetics",{});rna=evidence.get("transcriptomics",{})
    if mapping!="unique_stable_ensembl" or evidence_class=="untested" or biotype not in {"protein_coding","lncRNA"}:
        obj="untestable_nomination";question="missing_decisive_assay"
    elif biotype=="protein_coding":
        obj="protein_coding_gene";question="protein_coding_state_function"
    elif genomic in {"same_strand_exonic_overlap","same_strand_gene_body_overlap","complex"}:
        obj="ambiguous_lncrna_locus";question="lncrna_dna_transcription_or_rna"
    elif genomic in {"unknown","not_applicable"}:
        obj="untestable_nomination";question="missing_decisive_assay"
    else:
        obj="lncrna_transcript"
        question="lncrna_dna_transcription_or_rna" if evidence_class in {"genetically_anchored","concordant"} else "lncrna_rna_product"
    # Unknown genomic classification cannot be presented as a resolved lncRNA object.
    if biotype=="lncRNA" and genomic in {"unknown","not_applicable"}:
        biotype="unknown"
    record={"analysis_release_id":release,"catalog_object_id":release+":"+gene["ensembl_id"]+":"+obj,
        "ensembl_id":identity["gene_id_versioned"] if known else gene["ensembl_id"],
        "ensembl_id_base":gene["ensembl_id"],"symbol":gene["symbol"],"primary_evidence_class":evidence_class,
        "legacy_next_experiment_rule_id":gene["next_experiment_rule_id"],
        "gene_biotype":biotype,"candidate_object":obj,
        "lncrna_genomic_class":genomic if biotype=="lncRNA" else "not_applicable",
        "annotation_release":identity["annotation_release"] if known else "GENCODE v49; identity unresolved",
        "strand_audit_status":"not_audited","mapping_status":mapping,
        "noncoding_dna_context":"unresolved","credible_set_id":"",
        "credible_set_coding_pip_mass":None,"credible_set_noncoding_pip_mass":None,"credible_set_unresolved_pip_mass":None,
        "credible_set_to_gene_link_status":"no_authoritative_link",
        # A gene-level support call alone supplies no signal-pair target link.
        "target_link_basis":"none",
        "assay_applicability":{"bulk_rna":"untestable" if rna.get("testability_state")=="untestable" else "indeterminate",
            "proteomics":"not_applicable" if biotype=="lncRNA" else "indeterminate","chromatin":"indeterminate","spatial":"indeterminate"},
        "open_mechanistic_question":question,"genetics_call":genetic.get("call_state","indeterminate"),
        "rna_call":rna.get("call_state","indeterminate"),"candidate_not_adopted":True,
        "input_state":"provisional_pre_identity_and_count_correction" if provisional else "qualified_candidate",
        "unresolved_alternative":"Gene-level colocalization does not establish a credible-set target, conditional direction, or molecular mechanism.",
        "limitation":gene.get("limitation","")+(" Bulk identity and ALT correspondence remain unresolved." if provisional else "")}
    return catalog.adapt_record(record)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle",type=Path,required=True);p.add_argument("--gene-identity",type=Path,required=True)
    p.add_argument("--lncrna-class",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
    p.add_argument("--provisional-review",action="store_true");p.add_argument("--count-certificate",type=Path)
    args=p.parse_args()
    if not args.provisional_review:
        if not args.count_certificate:raise ValueError("Qualified candidate requires count certificate")
        c=json.loads(args.count_certificate.read_text())
        if not(c.get("annotation_resolution_complete") is True and c.get("n_unresolved_annotations")==0
               and c.get("counts_complete") is True and c.get("counts_sha256")):
            raise ValueError("Identity/counting completion is unresolved")
    args.output.mkdir(parents=True,exist_ok=False)
    identities=defaultdict(list)
    for r in table(args.gene_identity):
        if r["is_canonical_chromosome"]=="true":identities[r["gene_id_base"]].append(r)
    lncrna={r["gene_id_base"]:r for r in table(args.lncrna_class) if r["is_canonical_chromosome"]=="true"}
    index_path=args.bundle/"passport_gene_index.parquet";evidence_path=args.bundle/"passport_evidence_long.parquet"
    genes=pd.read_parquet(index_path,columns=["analysis_release_id","ensembl_id","symbol","primary_evidence_class",
        "next_experiment_rule_id","limitation"]).fillna("").to_dict("records")
    evidence=defaultdict(dict)
    for r in pd.read_parquet(evidence_path,columns=["ensembl_id","evidence_domain","call_state","testability_state"]).fillna("").to_dict("records"):
        evidence[r["ensembl_id"]][r["evidence_domain"]]=r
    release=genes[0]["analysis_release_id"]+":molecular-object-review-v2"
    records=[object_record(g,identities[g["ensembl_id"]][0] if len(identities[g["ensembl_id"]])==1 else None,
                           lncrna.get(g["ensembl_id"]),evidence[g["ensembl_id"]],release,args.provisional_review) for g in genes]
    output=args.output/"molecular_object_sidecar.jsonl"
    with output.open("x") as f:
        for r in records:f.write(json.dumps(r,sort_keys=True,allow_nan=False)+"\n")
    examples=[r for r in records if r["symbol"] in {"THRB","HKDC1","GLP1R","MTARC1"}]
    fields=["symbol","primary_evidence_class","candidate_object","genetics_call","rna_call","open_mechanistic_question",
            "recommended_experiment_rule_id","legacy_next_experiment_rule_id","input_state"]
    with (args.output/"worked_entries.tsv").open("x") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",extrasaction="ignore",lineterminator="\n");w.writeheader();w.writerows(examples)
    metadata={"candidate_not_adopted":True,"public_v1_modified":False,"manual_attestation_created":False,
        "provisional_inputs":args.provisional_review,"gene_rows":len(genes),"object_rows":len(records),
        "primary_class_counts":dict(Counter(g["primary_evidence_class"] for g in genes)),
        "object_route_counts":dict(Counter(r["recommended_experiment_rule_id"] for r in records)),
        "source_bundle":str(args.bundle.resolve()),"source_gene_index_sha256":digest(index_path),
        "source_evidence_sha256":digest(evidence_path),"gene_identity_sha256":digest(args.gene_identity),
        "lncrna_class_sha256":digest(args.lncrna_class),"sidecar_sha256":digest(output),
        "limitation":"No signal-pair target links, PIP masses or cross-assay gene testability were inferred from static overlaps."}
    (args.output/"sidecar_completion.json").write_text(json.dumps(metadata,indent=2)+"\n")
    rows="".join("<tr>"+"".join("<td>"+html.escape(str(r[k]))+"</td>" for k in fields)+"</tr>" for r in examples)
    page='<meta charset="utf-8"><title>MASLD Catalog object review</title><style>body{font:14px sans-serif;margin:2rem}td,th{padding:.5rem;border:1px solid #ccc}table{border-collapse:collapse}</style>'
    page+='<h1>MASLD Gene Catalog molecular-object review</h1><p>Candidate, not adopted. '+("Bulk identity and count correction pending." if args.provisional_review else "Qualified inputs; final scientific review pending.")+'</p>'
    page+='<p><a href="molecular_object_sidecar.jsonl">Download object records</a> | <a href="worked_entries.tsv">Worked entries</a> | <a href="sidecar_completion.json">Input/source state</a></p>'
    page+='<table><tr>'+"".join("<th>"+html.escape(k)+"</th>" for k in fields)+"</tr>"+rows+"</table>"
    (args.output/"index.html").write_text(page)
    print(json.dumps(metadata,sort_keys=True))


if __name__=="__main__":main()
