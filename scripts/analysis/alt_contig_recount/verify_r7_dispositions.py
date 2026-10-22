#!/usr/bin/env python3
"""Independent transition/source checks against the separate isoform review."""
from collections import Counter
import csv
import json
from pathlib import Path
import sys
from resolve_annotation_correspondence import digest


def table(path):
    with path.open() as handle: return list(csv.DictReader(handle, delimiter="\t"))


root = Path(sys.argv[1])
old_path, new_path = root/"annotation_resolution_r6", root/"annotation_resolution_r7"
old = {r["gene_id"]: r for r in table(old_path/"annotation_correspondence.tsv")}
new = {r["gene_id"]: r for r in table(new_path/"annotation_correspondence.tsv")}
assert set(old) == set(new) and len(new) == 86369
assert {gid for gid in old if old[gid]["disposition"] == "excluded_reference"} == {gid for gid in new if new[gid]["disposition"] == "excluded_reference"}
assert new["ENSG00000083067.27"]["disposition"] == "primary_gene"
review = root/"bulk_execution_sources/exact_shared_transcript_r6_review/exact_transcript_coordinate_checks.tsv"
proofs = {(r["gene_id"], r["source_transcript"], r["target_transcript"], r["exact_transcript_accession"]): r
          for r in table(review) if r["exact_complete_exon_correspondence"] == "True" and r["strand_matches"] == "True"
          and r["HGNC_conflict"] == "False" and r["group_conflict"] == "False"
          and r["same_gene_biotype"] == "True" and r["n_primary_genes_sharing_accession"] == "1"}
transitions = []
for gid, row in new.items():
    if row["disposition"] == old[gid]["disposition"]: continue
    assert row["disposition"] in {"equivalent", "duplicated"}
    if old[gid]["disposition"] == "unresolved":
        assert row["transcript_identity_evidence"]
        for serialized in row["transcript_identity_evidence"].split(";"):
            source, target, accessions = serialized.split("|")
            for accession in accessions.split(","):
                proof = proofs[(gid, source, target, accession)]
                if row["disposition"] == "equivalent": assert proof["primary_gene_id"] == row["counting_gene_id"]
    else:
        assert old[gid]["disposition"] == "equivalent" and row["disposition"] == "duplicated"
        assert any(other["disposition"] == "duplicated" and other["chromosome"] == row["chromosome"]
            and any(p["gene_id"] == other["gene_id"] and p["primary_gene_id"] == old[gid]["counting_gene_id"] for p in proofs.values())
            for other in new.values())
    transitions.append(dict(gene_id=gid, previous=old[gid]["disposition"], current=row["disposition"]))
source_records = table(new_path/"source_manifest.tsv")
assert all(digest(r["path"]) == r["sha256"] for r in source_records)
result = dict(status="PASS", independent_complete_RNA_coordinate_evidence_verified=True,
    annotations_preserved=len(new), unresolved=sum(r["disposition"] == "unresolved" for r in new.values()),
    excluded_reference=sum(r["disposition"] == "excluded_reference" for r in new.values()),
    transition_counts=dict(Counter(r["previous"]+"_to_"+r["current"] for r in transitions)),
    transitions=transitions, source_records_rehashed=len(source_records), TRPM3_retained=True,
    primary_counting_annotation_emitted=(new_path/"grouped_counting_annotation.gtf.gz").exists(),
    independent_coordinate_table_sha256=digest(review))
(new_path/"review_checks.json").write_text(json.dumps(result, indent=2)+"\n")
print(json.dumps(result, sort_keys=True))
