#!/usr/bin/env python3
"""Assemble immutable release-linked portal products without canonical writes."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from spatial_resource_lib import RESOURCE_RELEASE_ID, SpatialResourceError, sha256_file, write_json, write_tsv


PRODUCTS = (
    ("registry/spatial_dataset_registry.json", "spatial_dataset_registry.json"),
    ("coverage/spatial_gene_context.parquet", "spatial_gene_context.parquet"),
    ("coverage/spatial_program_coverage.parquet", "spatial_program_coverage.parquet"),
    ("effects/spatial_program_effects.parquet", "spatial_program_effects.parquet"),
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-root", type=Path, required=True)
    args = parser.parse_args()
    candidate = args.candidate_root.resolve()
    output = candidate / "portal"
    if (output / "READY").exists():
        raise SpatialResourceError(f"immutable portal candidate already exists: {output}")
    for prerequisite in ("registry/READY", "coverage/READY", "hmsma_label_blind/READY", "effects/READY"):
        if not (candidate / prerequisite).is_file():
            raise SpatialResourceError(f"candidate prerequisite missing: {prerequisite}")
    output.mkdir(parents=True, exist_ok=True)
    manifest = []
    for source_relative, target_name in PRODUCTS:
        source = candidate / source_relative
        target = output / target_name
        shutil.copy2(source, target)
        manifest.append({
            "release_id": RESOURCE_RELEASE_ID,
            "product": target_name,
            "source_relative_path": source_relative,
            "bytes": target.stat().st_size,
            "sha256": sha256_file(target),
        })
    write_json(
        output / "spatial_semantics.json",
        {
            "release_id": RESOURCE_RELEASE_ID,
            "categorical_types": {
                "dataset_gate": ["pass", "source_dependent", "metadata_pending", "skipped", "dropped"],
                "evidence_state": ["supported", "indeterminate", "untestable", "not_applicable", "source_dependent", "skipped", "dropped"],
            },
            "display_contract": {
                "zero_is_negative_only_with_adequate_negative_rule": True,
                "missing_renders_as": "untestable_or_indeterminate",
                "moran_i_semantics": "spatial_organization_not_disease_direction_or_percentage_strength",
                "biological_and_technical_counts_are_separate": True,
            },
            "legacy_compatibility": {
                "field": "s5_spatial",
                "retained_for_releases": 1,
                "excluded_from_rankings": True,
                "excluded_from_active_layer_counts": True,
                "excluded_from_evidence_strength_bars": True,
                "excluded_from_default_UI_decisions": True,
            },
        },
    )
    semantics = output / "spatial_semantics.json"
    manifest.append({
        "release_id": RESOURCE_RELEASE_ID,
        "product": semantics.name,
        "source_relative_path": "generated_from_spatial_resource_contract",
        "bytes": semantics.stat().st_size,
        "sha256": sha256_file(semantics),
    })
    write_tsv(
        output / "product_manifest.tsv",
        ("release_id", "product", "source_relative_path", "bytes", "sha256"),
        manifest,
    )
    write_tsv(
        output / "READY",
        ("release_id", "status", "n_products", "product_manifest_sha256", "canonical_portal_written"),
        [{
            "release_id": RESOURCE_RELEASE_ID,
            "status": "candidate_portal_products_ready_for_adjudication",
            "n_products": len(manifest),
            "product_manifest_sha256": sha256_file(output / "product_manifest.tsv"),
            "canonical_portal_written": "FALSE",
        }],
    )
    print(f"WROTE {len(manifest)} candidate portal products to {output}")


if __name__ == "__main__":
    main()
