#!/usr/bin/env bash
#SBATCH --job-name=validate_load
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/validate_load_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/validate_load_%j.err
set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

python3 << 'PYEOF'
import scanpy as sc
import squidpy as sq
import pathlib

SR_BASE = pathlib.Path("Analysis/Spatial/results/spaceranger/GSE192741")

# Test loading all 15 samples
samples = sorted([d.name for d in SR_BASE.iterdir() if d.is_dir()])
print(f"Found {len(samples)} samples\n")

total_spots = 0
for sample in samples:
    outs = SR_BASE / sample / "outs"
    try:
        adata = sc.read_visium(outs, count_file="filtered_feature_bc_matrix.h5")
        adata.var_names_make_unique()
        n_spots = adata.n_obs
        n_genes = adata.n_vars
        total_spots += n_spots
        has_spatial = "spatial" in adata.obsm
        has_img = "spatial" in adata.uns and len(adata.uns["spatial"]) > 0
        print(f"  {sample}: {n_spots:,} spots × {n_genes:,} genes, "
              f"spatial={has_spatial}, img={has_img}")
    except Exception as e:
        print(f"  {sample}: FAILED ({e})")

print(f"\nTotal: {total_spots:,} spots across {len(samples)} samples")

# Test squidpy spatial neighbors on first sample
print("\nTesting squidpy on first sample...")
adata = sc.read_visium(SR_BASE / samples[0] / "outs")
adata.var_names_make_unique()

# QC metrics
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
print(f"  Median genes/spot: {adata.obs['n_genes_by_counts'].median():.0f}")
print(f"  Median counts/spot: {adata.obs['total_counts'].median():.0f}")
print(f"  Median %MT: {adata.obs['pct_counts_mt'].median():.1f}%")

# Spatial neighbors
sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)
print(f"  Spatial neighbors: OK (connectivity shape: {adata.obsp['spatial_connectivities'].shape})")

# Quick Moran's I on total_counts
sq.gr.spatial_autocorr(adata, mode="moran", genes=["total_counts"], n_perms=50)
morans_i = adata.uns["moranI"].loc["total_counts", "I"]
print(f"  Moran's I (total_counts): {morans_i:.3f}")

print("\nVALIDATION COMPLETE - Pipeline ready for GSE192741!")
PYEOF
