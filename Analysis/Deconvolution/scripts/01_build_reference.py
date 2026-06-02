#!/usr/bin/env python3
"""Build a merged reference AnnData from Cell Ranger filtered matrices."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List

import anndata as ad
import scanpy as sc


def read_manifest(manifest_path: Path, species: str) -> List[Path]:
    filtered_h5 = []
    with manifest_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"species", "filtered_h5", "cellranger_done"}
        if not required.issubset(reader.fieldnames or []):
            missing = required - set(reader.fieldnames or [])
            raise ValueError(f"Manifest missing columns: {missing}")
        for row in reader:
            if row["species"] != species:
                continue
            if row["cellranger_done"] != "1":
                continue
            filtered_h5.append(Path(row["filtered_h5"]))
    return filtered_h5


def annotate_qc(adata: ad.AnnData) -> None:
    mt_mask = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["mt"] = mt_mask
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="run_manifest.tsv path")
    parser.add_argument("--species", required=True, choices=["mouse", "human"])
    parser.add_argument("--output", required=True, help="Output .h5ad path")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-cells", type=int, default=3)
    parser.add_argument("--max-genes", type=int, default=6000)
    parser.add_argument("--max-mito", type=float, default=20.0)
    parser.add_argument("--no-filter", action="store_true")
    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    inputs = read_manifest(manifest_path, args.species)
    if not inputs:
        raise SystemExit(f"No completed {args.species} runs found in {manifest_path}")

    adatas = []
    for h5 in inputs:
        if not h5.exists():
            raise SystemExit(f"Missing file: {h5}")
        adata = sc.read_10x_h5(str(h5))
        if adata.n_obs == 0:
            print(f"Warning: {h5} contains no cells. Skipping.")
            continue
        adata.var_names_make_unique()
        run_id = h5.parent.parent.name  # .../<run>/outs/filtered_feature_bc_matrix.h5
        adata.obs["sample"] = run_id
        adata.obs["run_id"] = run_id
        adata.obs["species"] = args.species
        # Preserve raw counts for scANVI
        adata.layers["counts"] = adata.X.copy()
        annotate_qc(adata)
        if not args.no_filter:
            sc.pp.filter_cells(adata, min_genes=args.min_genes)
            sc.pp.filter_genes(adata, min_cells=args.min_cells)
            if args.max_genes > 0:
                adata = adata[adata.obs.n_genes_by_counts < args.max_genes, :]
            if args.max_mito >= 0:
                adata = adata[adata.obs.pct_counts_mt < args.max_mito, :]
        # Skip samples that are now empty after filtering
        if adata.n_obs == 0:
            print(f"Warning: {run_id} has 0 cells after filtering. Skipping.")
            continue
        adatas.append(adata)

    merged = ad.concat(
        adatas,
        join="outer",
        label="sample",
        keys=[a.obs["sample"].iloc[0] for a in adatas],
    )
    merged.obs_names_make_unique()
    merged.var_names_make_unique()
    merged.write_h5ad(str(output_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
