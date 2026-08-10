#!/usr/bin/env python3
"""Independently validate a sealed Plan 45 corrected registry v3."""

from __future__ import annotations

import json
from collections import Counter

from relay_common import CANDIDATE_ROOT, read_tsv, sha256_file


def truth(value: str) -> bool:
    return value.strip().lower() == "true"


def main() -> None:
    seal_path = CANDIDATE_ROOT / "REGISTRY_V3_SEALED.json"
    classes_path = CANDIDATE_ROOT / "corrected_evidence_registry_v3.tsv"
    loci_path = CANDIDATE_ROOT / "corrected_primary_genetic_locus_registry_v3.tsv"
    pairs_path = CANDIDATE_ROOT / "corrected_gene_gwas_pair_registry_v3.tsv"
    product_manifest = CANDIDATE_ROOT / "corrected_coloc_product_manifest.tsv"
    required = [seal_path, classes_path, loci_path, pairs_path, product_manifest]
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty registry-v3 artifact: {path}")

    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("status") != "registry_v3_frozen_orientation_pending":
        raise RuntimeError("Registry-v3 seal has an unexpected status")
    if seal.get("orientation_complete") is not False:
        raise RuntimeError("Registry v3 improperly claims completed orientation")
    if seal.get("experimental_target_list_frozen") is not False:
        raise RuntimeError("Registry v3 improperly claims a frozen target list")

    hashes = seal.get("output_sha256", {})
    for name, expected in hashes.items():
        path = CANDIDATE_ROOT / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Registry-v3 output hash mismatch: {name}")

    products = read_tsv(product_manifest)
    product_keys = {(row["study_name"], row["chromosome"]) for row in products}
    if len(products) != 1100 or len(product_keys) != 1100:
        raise RuntimeError("Corrected-coloc manifest must contain 1,100 unique products")
    if any(int(row["size_bytes"]) <= 0 or int(row["n_rows"]) <= 0 for row in products):
        raise RuntimeError("Corrected-coloc manifest contains an empty product")

    classes = read_tsv(classes_path)
    class_counts = Counter(row["static_class"] for row in classes)
    if dict(sorted(class_counts.items())) != seal.get("class_counts"):
        raise RuntimeError("Evidence-class counts do not match the seal")
    allowed = {
        "genetic_only",
        "disease_state_only",
        "convergent",
        "neither",
        "indeterminate_not_jointly_testable",
    }
    if set(class_counts) != allowed:
        raise RuntimeError(f"Unexpected evidence classes: {sorted(set(class_counts) - allowed)}")
    for row in classes:
        joint = truth(row["joint_testable"])
        genetic = truth(row["primary_genetic"])
        state = truth(row["established_state_associated"])
        expected = (
            "indeterminate_not_jointly_testable"
            if not joint
            else "convergent"
            if genetic and state
            else "genetic_only"
            if genetic
            else "disease_state_only"
            if state
            else "neither"
        )
        if row["static_class"] != expected:
            raise RuntimeError(f"Class logic mismatch for {row['gene_symbol']}")
        if genetic and row["orientation_status"] != "requires_targeted_credible_set_pair_export":
            raise RuntimeError(f"Primary genetic row bypassed orientation gate: {row['gene_symbol']}")

    loci = read_tsv(loci_path)
    if len(loci) != int(seal["n_primary_genetic_rows"]):
        raise RuntimeError("Primary-genetic locus-row count does not match the seal")
    by_locus: dict[str, list[dict[str, str]]] = {}
    for row in loci:
        by_locus.setdefault(row["coarse_locus_uid"], []).append(row)
        if row["orientation_status"] != "requires_all_eligible_pair_exports":
            raise RuntimeError(f"Locus bypassed orientation gate: {row['gene_symbol']}")
    if len(by_locus) != int(seal["n_coarse_physical_loci"]):
        raise RuntimeError("Coarse-locus count does not match the seal")
    for locus_uid, rows in by_locus.items():
        representatives = [row for row in rows if truth(row["is_representative"])]
        if len(representatives) != 1:
            raise RuntimeError(f"Locus {locus_uid} has {len(representatives)} representatives")
        representative = representatives[0]["gene_symbol"]
        if any(row["representative_gene"] != representative for row in rows):
            raise RuntimeError(f"Locus {locus_uid} has inconsistent representative labels")

    pairs = read_tsv(pairs_path)
    if len(pairs) != int(seal["n_gene_gwas_pair_rows"]):
        raise RuntimeError("Gene-GWAS pair count does not match the seal")
    pair_keys = {
        (row["gene_symbol"], row["ensembl_id"], row["gwas_name"])
        for row in pairs
    }
    if len(pair_keys) != len(pairs):
        raise RuntimeError("Gene-GWAS pair registry contains duplicate pairs")
    primary_pairs = [row for row in pairs if truth(row["pair_primary_genetic"])]
    if len(primary_pairs) != int(seal["n_primary_gene_gwas_pairs"]):
        raise RuntimeError("Primary gene-GWAS pair count does not match the seal")
    class_by_gene = {row["gene_symbol"]: row["static_class"] for row in classes}
    locus_gene_keys = {
        (row["coarse_locus_uid"], row["ensembl_id"]): row for row in loci
    }
    for row in pairs:
        if row["static_class"] != class_by_gene.get(row["gene_symbol"]):
            raise RuntimeError(f"Pair evidence-class mismatch: {row['gene_symbol']}")
        if truth(row["pair_primary_genetic"]):
            key = (row["coarse_locus_uid"], row["ensembl_id"])
            if key not in locus_gene_keys:
                raise RuntimeError(f"Primary pair lacks a locus-gene summary: {key}")
            expected_representative = locus_gene_keys[key]["is_representative"]
            if row["is_representative_gene_in_locus"] != expected_representative:
                raise RuntimeError(f"Pair representative flag mismatch: {key}")
            if row["orientation_status"] != "requires_targeted_credible_set_pair_export":
                raise RuntimeError(f"Primary pair bypassed orientation gate: {key}")

    print(
        "Plan 45 registry-v3 validation passed: "
        f"products={len(products)}; evidence_rows={len(classes)}; "
        f"primary_rows={len(loci)}; loci={len(by_locus)}; "
        f"gene_gwas_pairs={len(pairs)}; primary_pairs={len(primary_pairs)}; "
        "orientation firewall closed"
    )


if __name__ == "__main__":
    main()
