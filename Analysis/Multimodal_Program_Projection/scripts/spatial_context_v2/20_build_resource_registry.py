#!/usr/bin/env python3
"""Seal the manifest-driven spatial dataset registry for a candidate release."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from spatial_resource_lib import (
    PROGRAM_RELEASE_ID,
    RESOURCE_RELEASE_ID,
    SpatialResourceError,
    load_dataset_registry,
    sha256_file,
    write_json,
    write_tsv,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--candidate-root", type=Path, required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    project = args.project_root.resolve()
    config = Path(__file__).with_name("spatial_resource_datasets.tsv")
    rows = load_dataset_registry(config)
    output = args.candidate_root.resolve() / "registry"
    registry_json = output / "spatial_dataset_registry.json"

    payload = {
        "release_id": RESOURCE_RELEASE_ID,
        "program_release_id": PROGRAM_RELEASE_ID,
        "registry_contract_sha256": sha256_file(config),
        "datasets": [],
    }
    sources = []
    for row in rows:
        item = dict(row)
        item["coverage_required"] = row["coverage_required"].upper() == "TRUE"
        relative = row["gene_axis_source"].strip()
        if relative:
            source = project / relative
            if not source.is_file():
                raise SpatialResourceError(f"registered gene-axis source is missing: {source}")
            source_record = {
                "dataset_id": row["dataset_id"],
                "assay_id": row["assay_id"],
                "relative_path": relative,
                "bytes": source.stat().st_size,
                "sha256": sha256_file(source),
                "source_role": "assay_gene_axis_and_detection_source",
            }
            sources.append(source_record)
            item["gene_axis_source_sha256"] = source_record["sha256"]
            item["gene_axis_source_bytes"] = source_record["bytes"]
        else:
            item["gene_axis_source_sha256"] = ""
            item["gene_axis_source_bytes"] = None
        payload["datasets"].append(item)

    if args.check_only:
        observed = json.loads(registry_json.read_text(encoding="utf-8"))
        if observed != payload:
            raise SpatialResourceError("sealed spatial_dataset_registry.json drift")
        print(f"PASS registry check: {len(rows)} dataset/analysis decisions")
        return
    if registry_json.exists():
        raise SpatialResourceError(f"immutable candidate registry already exists: {registry_json}")
    output.mkdir(parents=True, exist_ok=True)
    write_json(registry_json, payload)
    write_tsv(
        output / "source_files.tsv",
        ("dataset_id", "assay_id", "relative_path", "bytes", "sha256", "source_role"),
        sources,
    )
    write_tsv(
        output / "READY",
        ("release_id", "status", "n_registry_rows", "n_coverage_assays", "registry_contract_sha256", "registry_json_sha256"),
        [{
            "release_id": RESOURCE_RELEASE_ID,
            "status": "sealed_candidate_registry_no_canonical_promotion",
            "n_registry_rows": len(rows),
            "n_coverage_assays": sum(row["coverage_required"].upper() == "TRUE" for row in rows),
            "registry_contract_sha256": sha256_file(config),
            "registry_json_sha256": sha256_file(registry_json),
        }],
    )
    print(f"WROTE {registry_json} ({len(rows)} decisions)")


if __name__ == "__main__":
    main()
