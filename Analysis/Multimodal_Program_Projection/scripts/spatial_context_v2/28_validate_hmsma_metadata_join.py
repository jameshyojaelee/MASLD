#!/usr/bin/env python3
"""Fail-closed preflight for an author-supplied HMSMA metadata join."""

from __future__ import annotations

import argparse
from pathlib import Path

from spatial_resource_lib import (
    SpatialResourceError,
    h5ad_column_values,
    load_dataset_registry,
    read_tsv,
    validate_hmsma_metadata_join,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--metadata-tsv", type=Path, required=True)
    parser.add_argument(
        "--inference-kind",
        choices=("donor", "clinical", "histology", "metabolite"),
        default="donor",
    )
    args = parser.parse_args()

    project = args.project_root.resolve()
    script_root = Path(__file__).resolve().parent
    registry = load_dataset_registry(script_root / "spatial_resource_datasets.tsv")
    hmsma = next(row for row in registry if row["dataset_id"] == "HRA007511_HMSMA")
    h5ad_path = Path(hmsma["gene_axis_source"])
    if not h5ad_path.is_absolute():
        h5ad_path = project / h5ad_path

    try:
        import h5py
    except ImportError as exc:
        raise SpatialResourceError("h5py is required to enumerate HMSMA arrays") from exc
    with h5py.File(h5ad_path, "r") as handle:
        if "sample_id" not in handle["obs"]:
            raise SpatialResourceError("HMSMA H5AD lacks registered obs/sample_id")
        expected_arrays = sorted(set(h5ad_column_values(handle["obs/sample_id"])))
    if len(expected_arrays) != 35:
        raise SpatialResourceError(f"expected 35 HMSMA arrays, found {len(expected_arrays)}")

    _, rows = read_tsv(args.metadata_tsv)
    validate_hmsma_metadata_join(rows, expected_arrays, args.inference_kind)
    print(
        "PASS HMSMA authoritative metadata join "
        f"arrays={len(expected_arrays)} inference_kind={args.inference_kind}"
    )


if __name__ == "__main__":
    main()
