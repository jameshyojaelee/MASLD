#!/usr/bin/env python3
"""Independent exact-ID and source-link review; do not invent gene equivalence."""
import csv
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import sys


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    root = Path(sys.argv[1])
    new = root/"annotation_resolution_r6"
    old = root/"annotation_resolution_r5"
    sources = root/"bulk_execution_sources"
    output = Path(sys.argv[2]) if len(sys.argv) > 2 else sources/"remaining_identity_evidence_r6_review"
    output.mkdir(exist_ok=False)
    prior = {r["gene_id"]: r for r in read(old/"annotation_correspondence.tsv")}
    current = {r["gene_id"]: r for r in read(new/"annotation_correspondence.tsv")}
    assert set(prior) == set(current) and len(current) == 86369
    contamination = {gid for gid, r in prior.items() if r["reason"] == "GRC_contamination_reference_policy_pending"}
    false_dup = {r["gene_id"] for r in read(sources/"false_duplication_r5_review/false_duplication_gene_overlaps.tsv")
                 if r["gene_body_wholly_within_false_duplication"] == "True" and r["all_exons_within_false_duplication"] == "True"}
    changed = {gid for gid in current if current[gid]["disposition"] != prior[gid]["disposition"]}
    assert changed == contamination | false_dup and len(contamination) == 7 and len(false_dup) == 56
    assert all(current[gid]["disposition"] == "excluded_reference" and current[gid]["counting_gene_id"] == "" for gid in changed)
    assert current["ENSG00000083067.27"] == prior["ENSG00000083067.27"]  # TRPM3
    remaining = {gid for gid, r in current.items() if r["disposition"] == "unresolved"}
    links = [r for r in read(sources/"unresolved_exact_transcript_source_links.tsv") if r["gene_id"] in remaining]
    by_gene = defaultdict(list)
    for r in links:
        by_gene[r["gene_id"]].append(r)
    primary_members = defaultdict(set)
    membership_path = sources/"ensembl115_all_gtf_alt_allele_membership.tsv"
    members = read(membership_path)
    # The exact group source table contains more than one possible assembly locus.
    for r in members:
        gid = r["gene_id"]
        if gid in current and current[gid]["disposition"] == "primary_gene":
            primary_members[r["alt_allele_group_id"]].add(gid)
    rows = []
    for gid, candidates in sorted(by_gene.items()):
        mapped = [r for r in candidates if r["target_in_official_mapped_exon_candidates"] == "True"]
        group_targets = set().union(*(primary_members[g] for g in current[gid]["ensembl_allele_group_ids"].split(";") if g))
        for r in mapped:
            rows.append(dict(r, remaining_reason=current[gid]["reason"],
                official_group_primary_targets=";".join(sorted(group_targets)),
                conflicts_with_declared_primary_group=bool(group_targets and r["primary_gene_id"] not in group_targets)))
    with (output/"remaining_mapped_transcript_links.tsv").open("x") as handle:
        fields = list(rows[0]) if rows else ["gene_id"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    report = dict(status="PASS_reference_exclusions_source_links_reviewed_no_extra_identity_assigned",
        exact_ID_universe_preserved=len(current), excluded_reference=len(changed),
        contamination_excluded=len(contamination), false_duplication_excluded=len(false_dup), TRPM3_unchanged=True,
        remaining_unresolved=len(remaining), remaining_reason_counts=dict(Counter(current[gid]["reason"] for gid in remaining)),
        remaining_genes_with_mapped_transcript_links=len({r["gene_id"] for r in rows}),
        mapped_link_kinds=dict(Counter(r["kind"] for r in rows)),
        Havana_mapped_links=[r for r in rows if r["kind"] == "Havana_transcript"],
        HGNC_conflicting_links=sum(r["HGNC_conflict"] == "True" for r in rows),
        conflicting_primary_group_links=sum(r["conflicts_with_declared_primary_group"] for r in rows),
        interpretation="Shared external transcript accession alone is not gene identity; official gene-group conflicts are retained.",
        input_sha256={str(p): sha(p) for p in [new/"annotation_correspondence.tsv", old/"annotation_correspondence.tsv", membership_path,
            sources/"unresolved_exact_transcript_source_links.tsv"]})
    (new/"review_checks.json").write_text(json.dumps(report, indent=2)+"\n")
    (output/"remaining_identity_review.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
