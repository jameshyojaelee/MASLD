#!/usr/bin/env python3
"""
02_build_anndata.py — Convert SpaceRanger outputs to AnnData, QC, normalize.

Loads all SpaceRanger outputs, applies spot/gene QC, normalizes expression,
selects HVGs, and saves merged checkpoint. GPU-accelerated where available.

SLURM: --partition=cpu --cpus=8 --mem=64G --time=4:00:00
"""

import pathlib
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, SPATIAL_ROOT, RESULTS_DIR, load_config, load_dataset_config,
    init_spatial_gpu, get_processor, save_checkpoint, check_checkpoint,
    print_header, print_step,
)


def load_visium_sample(spaceranger_outs: pathlib.Path, sample_id: str) -> ad.AnnData:
    """Load one SpaceRanger output into AnnData with spatial coords."""
    adata = sc.read_visium(spaceranger_outs, count_file="filtered_feature_bc_matrix.h5")
    adata.obs["sample_id"] = sample_id
    adata.var_names_make_unique()
    # Store raw counts in layer (cell2location needs raw counts)
    adata.layers["counts"] = adata.X.copy()
    return adata


def spatial_qc(adata: ad.AnnData, params: dict, sample_id: str) -> ad.AnnData:
    """QC filtering with spatial autocorrelation check."""
    import squidpy as sq

    # Calculate QC metrics
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)

    # Spatial autocorrelation of QC metrics (detect tissue damage)
    try:
        sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)
        # Use a highly expressed gene as proxy for QC spatial autocorrelation
        top_gene = adata.var_names[np.array(adata.X.sum(axis=0)).flatten().argmax()]
        sq.gr.spatial_autocorr(adata, mode="moran", genes=[top_gene], n_perms=50)
    except Exception as e:
        print(f"    WARNING: Spatial autocorr failed for {sample_id}: {e}")

    # Filter spots
    n_before = adata.n_obs
    mask = (
        (adata.obs["n_genes_by_counts"] >= params["min_genes_per_spot"])
        & (adata.obs["total_counts"] >= params["min_counts_per_spot"])
        & (adata.obs["pct_counts_mt"] <= params["max_pct_mt"])
    )
    adata = adata[mask].copy()

    # Filter genes
    sc.pp.filter_genes(adata, min_cells=params["min_spots_per_gene"])

    pct = adata.n_obs / max(n_before, 1) * 100
    print(f"    QC: {n_before} → {adata.n_obs} spots ({pct:.1f}% retained), "
          f"{adata.n_vars} genes")
    return adata


def main():
    parser = argparse.ArgumentParser(description="Build AnnData from SpaceRanger outputs")
    parser.add_argument("--use-gpu", action="store_true", default=True,
                        help="Use GPU for normalization (default: True)")
    parser.add_argument("--force", action="store_true",
                        help="Reprocess even if checkpoint exists")
    parser.add_argument("--dataset", type=str, default=None,
                        help="Process specific dataset(s), comma-separated "
                             "(e.g., GSE192741,HRA007511_HMSMA). Default: all")
    parser.add_argument("--output-suffix", type=str, default="",
                        help="Suffix for output checkpoint files "
                             "(e.g., '_hmsma' → merged_spatial_hmsma.h5ad)")
    args = parser.parse_args()

    # Determine output names
    suffix = args.output_suffix
    checkpoint_raw = f"merged_spatial{suffix}_raw.h5ad"
    checkpoint_norm = f"merged_spatial{suffix}.h5ad"

    print_header("02: Build AnnData + QC + Normalize")
    if args.dataset:
        print(f"  Dataset filter: {args.dataset}")
    if suffix:
        print(f"  Output suffix: {suffix}")

    # Check checkpoint
    if not args.force and check_checkpoint(checkpoint_norm):
        print(f"  Checkpoint exists ({checkpoint_norm}). Use --force to reprocess.")
        return

    # Initialize GPU if requested
    use_gpu = args.use_gpu
    if use_gpu:
        try:
            init_spatial_gpu()
            print("  GPU initialized successfully")
        except Exception as e:
            print(f"  GPU init failed ({e}), falling back to CPU")
            use_gpu = False

    config = load_config()
    datasets = load_dataset_config()
    qc_params = config["qc"]
    norm_params = config["normalization"]

    # Load all samples
    adatas = []
    qc_summary = []

    # F008 fix: the canonical unsuffixed `merged_spatial.h5ad` must be the
    # GSE192741 DISCOVERY cohort only. Merging the validation-only Vu cohort
    # (FFPE CytAssist) into the same object pooled discovery + validation across
    # a fresh-frozen-vs-FFPE batch axis and collapsed the discovery/validation
    # split the manuscript claims. Default now = GSE192741 only; pass
    # `--dataset all` to explicitly merge every dataset, or name datasets.
    DEFAULT_PRIMARY_DATASET = "GSE192741"
    selected_datasets = None
    if args.dataset:
        if args.dataset.strip().lower() == "all":
            selected_datasets = None  # explicit opt-in to merge ALL datasets
        else:
            selected_datasets = set(d.strip() for d in args.dataset.split(","))
    elif not suffix:
        # unsuffixed canonical build → discovery cohort only
        selected_datasets = {DEFAULT_PRIMARY_DATASET}
        print(f"  (default) canonical object = {DEFAULT_PRIMARY_DATASET} only "
              f"(pass --dataset all to merge every dataset)")

    for dataset_name, ds_config in datasets["datasets"].items():
        # Skip datasets not in filter
        if selected_datasets and dataset_name not in selected_datasets:
            continue

        sr_base = RESULTS_DIR / "spaceranger" / dataset_name
        if not sr_base.exists():
            print(f"\n  Skipping {dataset_name}: no SpaceRanger output")
            continue

        print(f"\n  Loading {dataset_name}...")
        # Get per-sample metadata if available
        sample_meta = ds_config.get("sample_metadata", {})

        sample_dirs = sorted([
            d for d in sr_base.iterdir()
            if d.is_dir() and (d / "outs" / "filtered_feature_bc_matrix.h5").exists()
        ])

        for i, sample_dir in enumerate(sample_dirs):
            sample_id = sample_dir.name
            meta = sample_meta.get(sample_id, {})
            species = meta.get("species", ds_config.get("species", "unknown"))
            condition = meta.get("condition", ds_config.get("condition", "unknown"))

            # Skip mouse samples for human analysis pipeline
            if species == "mouse":
                print_step(f"Skipping {sample_id} (mouse)", i + 1, len(sample_dirs))
                continue

            print_step(f"Loading {sample_id} ({species}, {condition})", i + 1, len(sample_dirs))

            try:
                adata = load_visium_sample(sample_dir / "outs", sample_id)
                adata.obs["dataset"] = dataset_name
                adata.obs["species"] = species
                adata.obs["condition"] = condition
                adata.obs["individual"] = meta.get("individual", sample_id)
                n_before = adata.n_obs

                adata = spatial_qc(adata, qc_params, sample_id)

                qc_summary.append({
                    "sample_id": sample_id,
                    "dataset": dataset_name,
                    "species": species,
                    "condition": condition,
                    "individual": meta.get("individual", sample_id),
                    "n_spots_raw": n_before,
                    "n_spots_qc": adata.n_obs,
                    "n_genes": adata.n_vars,
                    "pct_retained": adata.n_obs / max(n_before, 1) * 100,
                    "median_genes": adata.obs["n_genes_by_counts"].median(),
                    "median_counts": adata.obs["total_counts"].median(),
                    "median_pct_mt": adata.obs["pct_counts_mt"].median(),
                })

                adatas.append(adata)

            except Exception as e:
                print(f"    ERROR loading {sample_id}: {e}")
                continue

    if not adatas:
        print("\n  ERROR: No samples loaded successfully!")
        sys.exit(1)

    # Save QC summary
    qc_df = pd.DataFrame(qc_summary)
    qc_file = RESULTS_DIR / "qc" / f"spot_qc_summary{suffix}.csv"
    qc_file.parent.mkdir(parents=True, exist_ok=True)
    qc_df.to_csv(qc_file, index=False)
    print(f"\n  Total: {len(adatas)} samples, {sum(a.n_obs for a in adatas)} spots")

    # Merge all samples
    print("\n  Merging samples...")
    merged = ad.concat(adatas, join="inner", merge="same")
    merged.obs_names_make_unique()
    print(f"  Merged: {merged.n_obs} spots × {merged.n_vars} genes")

    # Save raw checkpoint (before normalization)
    save_checkpoint(merged, checkpoint_raw)

    # Normalize
    print("\n  Normalizing...")
    proc = get_processor(use_gpu)
    proc.pp.normalize_total(merged, target_sum=norm_params["target_sum"])

    # HVG selection BEFORE log1p (seurat_v3 requires count-scale data)
    print(f"  Selecting {norm_params['n_top_genes']} HVGs...")
    sc.pp.highly_variable_genes(
        merged,
        n_top_genes=norm_params["n_top_genes"],
        flavor=norm_params["hvg_flavor"],
    )

    proc.pp.log1p(merged)
    n_hvg = merged.var["highly_variable"].sum()
    print(f"  HVGs selected: {n_hvg}")

    # Save normalized checkpoint
    save_checkpoint(merged, checkpoint_norm)

    # Summary stats
    print(f"\n  === Summary ===")
    print(f"  Samples: {merged.obs['sample_id'].nunique()}")
    print(f"  Datasets: {merged.obs['dataset'].nunique()}")
    print(f"  Spots: {merged.n_obs}")
    print(f"  Genes: {merged.n_vars}")
    print(f"  HVGs: {n_hvg}")
    for ds in merged.obs["dataset"].unique():
        n = (merged.obs["dataset"] == ds).sum()
        print(f"    {ds}: {n} spots")

    print_header("02: Complete")


if __name__ == "__main__":
    main()
