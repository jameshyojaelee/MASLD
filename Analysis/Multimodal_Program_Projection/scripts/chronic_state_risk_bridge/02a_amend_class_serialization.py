#!/usr/bin/env python3
"""Non-overwriting repair of the sealed corrected-class TSV header."""

from __future__ import annotations

import datetime as dt
import json

from bridge_common import (
    CANDIDATE_ROOT,
    atomic_write_text,
    read_tsv,
    require_validated_seal,
    sha256_file,
    write_tsv,
)


def main() -> None:
    seal = require_validated_seal()
    output_root = CANDIDATE_ROOT / "integrity_correction_v2"
    if output_root.exists():
        raise RuntimeError(f"Refusing to overwrite integrity amendment: {output_root}")
    classes_path = CANDIDATE_ROOT / "corrected_evidence_classes.tsv"
    loci_path = CANDIDATE_ROOT / "corrected_genetic_locus_registry.tsv"
    classes = read_tsv(classes_path)
    loci = read_tsv(loci_path)
    by_gene = {row["gene_symbol"]: row for row in loci}
    correction_fields = [
        "tier12_driving_gwas", "tier12_driving_trait", "tier12_driving_tier",
        "tier12_direct_or_proxy", "tier12_top_snp", "tier12_locus_uid",
        "tier12_susie_pp4", "tier12_abf_pp4", "any_direct_tier1_support",
        "tier1_supporting_gwas", "full_portfolio_driving_gwas",
        "full_portfolio_driving_trait", "driver_corrected",
    ]
    corrected = []
    for row in classes:
        output = dict(row)
        locus = by_gene.get(row["gene_symbol"])
        for field in correction_fields:
            output[field] = locus.get(field, "") if locus else ""
        corrected.append(output)
    fields = list(classes[0]) + correction_fields
    output_path = output_root / "corrected_evidence_classes_v2.tsv"
    write_tsv(output_path, corrected, fields)
    if sum(bool(row["tier12_locus_uid"]) for row in corrected) != 413:
        raise RuntimeError("Corrected class amendment does not contain 413 genetic rows")
    manifest = [
        {
            "artifact": "sealed_header_truncated_class_table",
            "path": str(classes_path.relative_to(CANDIDATE_ROOT)),
            "sha256": sha256_file(classes_path),
            "status": "preserved_superseded_for_downstream_use",
        },
        {
            "artifact": "corrected_evidence_classes_v2",
            "path": str(output_path.relative_to(CANDIDATE_ROOT)),
            "sha256": sha256_file(output_path),
            "status": "active_nonoverwriting_integrity_correction",
        },
        {
            "artifact": "corrected_genetic_locus_registry",
            "path": str(loci_path.relative_to(CANDIDATE_ROOT)),
            "sha256": sha256_file(loci_path),
            "status": "unchanged_authoritative_locus_source",
        },
    ]
    write_tsv(output_root / "integrity_correction_manifest.tsv", manifest, list(manifest[0]))
    amendment = {
        "amendment_id": "plan43-integrity-amendment-01-class-header",
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "seal_sha256": sha256_file(CANDIDATE_ROOT / "SEALED.json"),
        "specification_sha256": seal["specification_sha256"],
        "reason": "the sealed producer inferred the corrected-class header from the first non-genetic row, dropping genetic-only correction columns from serialization",
        "scientific_hypothesis_changed": False,
        "outcome_values_used": False,
        "sealed_artifact_overwritten": False,
        "n_rows": len(corrected),
        "n_genetic_rows_with_corrected_locus": 413,
        "active_table_sha256": sha256_file(output_path),
    }
    atomic_write_text(output_root / "INTEGRITY_AMENDMENT_01.json", json.dumps(amendment, indent=2, sort_keys=True) + "\n")
    print(json.dumps(amendment, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
