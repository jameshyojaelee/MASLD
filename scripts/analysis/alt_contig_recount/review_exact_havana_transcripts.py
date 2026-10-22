#!/usr/bin/env python3
"""Check exact manual transcript identities against independent exon coordinates."""
import csv
import gzip
import json
from pathlib import Path
import sys
from resolve_annotation_correspondence import ATTR, digest, alignment_blocks, project, union, length


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


root = Path(sys.argv[1])
annotation = root/"annotation_resolution_r6"
review = root/"bulk_execution_sources/remaining_identity_evidence_r6_review_v2"
include_refseq = len(sys.argv) > 2 and sys.argv[2] == "include-refseq"
out = root/"bulk_execution_sources"/("exact_shared_transcript_r6_review" if include_refseq else "exact_havana_transcript_r6_review")
out.mkdir(exist_ok=False)
summary = json.loads((annotation/"resolution_summary.json").read_text())
gtf = Path(summary["gtf"])
assert digest(gtf) == summary["gtf_sha256"]
links = [r for r in table(review/"remaining_mapped_transcript_links.tsv")
         if r["kind"] == "Havana_transcript" or (include_refseq and r["kind"] == "RefSeq")]
refseq = {}
refseq_path = root/"bulk_execution_sources/followup/gencode49_RefSeq.gz"
if include_refseq:
    with gzip.open(refseq_path, "rt") as handle:
        for line in handle:
            c = line.rstrip().split("\t")
            if len(c) >= 2: refseq.setdefault(c[0], set()).add(c[1])
wanted = {r["gene_id"] for r in links} | {r["primary_gene_id"] for r in links}
models, attrs = {}, {}
with gzip.open(gtf, "rt") as handle:
    for line in handle:
        if line.startswith("#"): continue
        c = line.rstrip().split("\t")
        if len(c) != 9 or c[2] not in {"transcript", "exon"}: continue
        a = dict(ATTR.findall(c[8]))
        if a.get("gene_id") not in wanted: continue
        key = (a["gene_id"], a["transcript_id"])
        if c[2] == "transcript":
            attrs[key] = dict(havana=a.get("havana_transcript", ""), chrom=c[0], strand=c[6])
        else:
            models.setdefault(key, []).append((int(c[3])-1, int(c[4])))
correspondence = {r["gene_id"]: r for r in table(annotation/"annotation_correspondence.tsv")}
rows = []
for link in links:
    source = correspondence[link["gene_id"]]
    target = correspondence[link["primary_gene_id"]]
    alignment_path = Path(source["evidence_source"])
    parent, blocks = alignment_blocks(alignment_path)
    for source_key, source_attr in attrs.items():
        if source_key[0] != source["gene_id"]: continue
        source_accessions = {source_attr["havana"]} if link["kind"] == "Havana_transcript" else refseq.get(source_key[1], set())
        if link["exact_accession"] not in source_accessions: continue
        for target_key, target_attr in attrs.items():
            if target_key[0] != target["gene_id"]: continue
            target_accessions = {target_attr["havana"]} if link["kind"] == "Havana_transcript" else refseq.get(target_key[1], set())
            if link["exact_accession"] not in target_accessions: continue
            exons = union(models[source_key])
            mapped, covered, orientations = project(exons, blocks)
            strand = source_attr["strand"]
            if orientations == {-1}: strand = "-" if strand == "+" else "+"
            target_exons = union(models[target_key])
            rows.append(dict(gene_id=source_key[0], primary_gene_id=target_key[0],
                kind=link["kind"], exact_transcript_accession=link["exact_accession"],
                source_transcript=source_key[1], target_transcript=target_key[1],
                source_exonic_bases=length(exons), mapped_bases=covered,
                exact_complete_exon_correspondence=covered == length(exons) and mapped == target_exons,
                strand_matches=strand == target_attr["strand"] and len(orientations) == 1,
                HGNC_conflict=link["HGNC_conflict"], group_conflict=link["conflicts_with_declared_primary_group"],
                same_gene_biotype=source["gene_biotype"] == target["gene_biotype"],
                n_primary_genes_sharing_accession=link["n_primary_genes_sharing_accession"],
                source_exons=json.dumps(exons), mapped_exons=json.dumps(mapped), primary_exons=json.dumps(target_exons)))
with (out/"exact_transcript_coordinate_checks.tsv").open("x") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
result = dict(status="read_only_no_dispositions_changed", exact_transcript_links=len(links),
    transcript_pairs_checked=len(rows), fully_matching_pairs=[r for r in rows if r["exact_complete_exon_correspondence"] and r["strand_matches"]],
    GTF_sha256=digest(gtf), RefSeq_metadata_sha256=digest(refseq_path) if include_refseq else None)
(out/"transcript_coordinate_review.json").write_text(json.dumps(result, indent=2)+"\n")
print(json.dumps(result, sort_keys=True))
