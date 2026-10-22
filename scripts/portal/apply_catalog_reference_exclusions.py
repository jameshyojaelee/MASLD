#!/usr/bin/env python3
"""Retain provisional Catalog records while marking reference artifacts non-applicable."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import shutil
import re
from gene_catalog_v2 import canonical_assay_applicability


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--annotation", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.output.exists(): raise ValueError("Refusing to overwrite a Catalog candidate")
    correspondence_path = args.annotation/"annotation_correspondence.tsv"
    summary = json.loads((args.annotation/"resolution_summary.json").read_text())
    assert sha(correspondence_path) == summary["annotation_correspondence_sha256"]
    with correspondence_path.open() as handle:
        correspondence = {r["gene_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
    source_path = args.source/"molecular_object_sidecar.jsonl"
    by_stable_id = {}
    for gene_id in correspondence:
        by_stable_id.setdefault(gene_id.rsplit(".", 1)[0], []).append(gene_id)
    metadata = json.loads((args.source/"sidecar_completion.json").read_text())
    assert metadata["provisional_inputs"] and sha(source_path) == metadata["sidecar_sha256"]
    records = [json.loads(line) for line in source_path.open()]
    baseline_ids = [r["catalog_object_id"] for r in records]
    baseline_classes = Counter(r["primary_evidence_class"] for r in records)
    affected = []
    for r in records:
        identity = correspondence.get(r["ensembl_id"])
        match_basis = "exact_versioned_GENCODE_gene_id"
        # The source Catalog intentionally stores stable IDs without versions
        # where its annotation sidecar is absent. Resolve only that explicit
        # representation within the same release and one exact GTF version.
        if (not identity and re.fullmatch(r"ENSG\d+", r["ensembl_id"])
                and r["annotation_release"].startswith("GENCODE v49")
                and len(by_stable_id.get(r["ensembl_id"], [])) == 1):
            identity = correspondence[by_stable_id[r["ensembl_id"]][0]]
            match_basis = "source_unversioned_stable_id_to_unique_GTF_version_same_GENCODE_release"
        if not identity or identity["disposition"] != "excluded_reference": continue
        r["source_rna_call"] = r["rna_call"]
        r["rna_call"] = "not_applicable"
        applicability = json.loads(r["assay_applicability"])
        applicability["bulk_rna"] = "not_applicable"
        r["assay_applicability"] = canonical_assay_applicability(applicability)
        r["bulk_rna_reference_evidence_state"] = "not_applicable_reference_artifact"
        r["bulk_rna_reference_reason"] = identity["reason"]
        r["bulk_rna_reference_gene_id"] = identity["gene_id"]
        r["bulk_rna_reference_identity_match_basis"] = match_basis
        r["bulk_rna_reference_annotation_sha256"] = summary["annotation_correspondence_sha256"]
        r["limitation"] += " This annotation lies wholly in a GRC invalid reference region; its bulk RNA call is non-applicable. Source evidence class remains historical pending refit."
        affected.append(identity["gene_id"])
    assert baseline_ids == [r["catalog_object_id"] for r in records]
    assert baseline_classes == Counter(r["primary_evidence_class"] for r in records)
    args.output.mkdir(parents=True)
    dest = args.output/"molecular_object_sidecar.jsonl"
    with dest.open("x") as handle:
        for r in records: handle.write(json.dumps(r, sort_keys=True, allow_nan=False)+"\n")
    metadata.update(sidecar_sha256=sha(dest), source_sidecar=str(source_path.resolve()),
        source_sidecar_sha256=sha(source_path), annotation_correspondence_sha256=sha(correspondence_path),
        reference_artifact_gene_records_retained=len(affected), reference_artifact_gene_ids=sorted(affected),
        reference_artifact_annotation_count=summary["excluded_reference_count"],
        primary_classes_historical_unchanged=True,
        identity_join="exact versions or source-unversioned stable IDs with one exact version in the same GENCODE release; never symbols")
    (args.output/"sidecar_completion.json").write_text(json.dumps(metadata, indent=2)+"\n")
    shutil.copyfile(args.source/"worked_entries.tsv", args.output/"worked_entries.tsv")
    page = (args.source/"index.html").read_text()
    page += f"<p>{len(affected)} reference-artifact records are retained with non-applicable bulk RNA evidence. Their source classes remain historical pending the corrected refit.</p>"
    (args.output/"index.html").write_text(page)
    print(json.dumps(dict(status="PASS", object_rows=len(records), unchanged_source_class_counts=dict(baseline_classes),
        reference_artifact_records_retained=len(affected), symbol_join_used=False, public_release_modified=False)))


if __name__ == "__main__": main()
