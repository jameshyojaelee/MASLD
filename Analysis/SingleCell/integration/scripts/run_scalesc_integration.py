#!/usr/bin/env python3
"""Run ScaleSC integration on prepared h5ad files.

Usage:
  python run_scalesc_integration.py --species human [--resolution 0.5] [--n_hvgs 4000]
  python run_scalesc_integration.py --species mouse [--resolution 0.3]

Expects h5ad files in:
  {BASE}/Analysis/SingleCell/integration/input_h5ad/{species}/
"""

import argparse
import logging
import os
import sys
import time
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INTEGRATION_DIR = BASE / "Analysis/SingleCell/integration"


def init_gpu():
    """Initialize GPU with PyTorch first, then RMM pool."""
    import torch
    if not torch.cuda.is_available():
        log.error("No GPU available")
        sys.exit(1)

    device = torch.cuda.current_device()
    gpu_name = torch.cuda.get_device_name(device)
    gpu_mem = torch.cuda.get_device_properties(device).total_memory // (1024**2)
    log.info(f"GPU: {gpu_name}, {gpu_mem} MiB")

    # Init CUBLAS via a small matmul
    _ = torch.randn(2, 2, device="cuda") @ torch.randn(2, 2, device="cuda")
    torch.cuda.synchronize()

    # RMM pool (30% of GPU memory)
    try:
        import rmm
        from rmm.allocators.cupy import rmm_cupy_allocator
        import cupy as cp

        pool_size = int(gpu_mem * 0.8) * 1024**2
        rmm.reinitialize(pool_allocator=True, initial_pool_size=pool_size)
        cp.cuda.set_allocator(rmm_cupy_allocator)
        log.info(f"RMM pool: {pool_size // (1024**3):.1f} GB")
    except Exception as e:
        log.warning(f"RMM init failed: {e}")

    return gpu_name, gpu_mem


def run_integration(species: str, resolution: float, n_hvgs: int, max_cell_batch: int):
    """Run ScaleSC integration pipeline."""
    import scalesc as ssc

    input_dir = INTEGRATION_DIR / "input_h5ad" / species
    output_dir = INTEGRATION_DIR / "output" / species
    output_dir.mkdir(parents=True, exist_ok=True)

    h5ad_files = sorted(input_dir.glob("*.h5ad"))
    log.info(f"Input: {len(h5ad_files)} h5ad files in {input_dir}")

    if not h5ad_files:
        log.error(f"No h5ad files found in {input_dir}")
        sys.exit(1)

    # --- ScaleSC Pipeline ---
    t0 = time.time()
    log.info("=" * 60)
    log.info(f"ScaleSC Integration: {species}")
    log.info("=" * 60)

    pipeline = ssc.ScaleSC(
        data_dir=str(input_dir),
        max_cell_batch=max_cell_batch,
        preload_on_cpu=True,
        preload_on_gpu=True,
        output_dir=str(output_dir),
    )
    log.info(f"Init + data load: {time.time() - t0:.1f}s")

    t1 = time.time()
    pipeline.calculate_qc_metrics()
    log.info(f"QC metrics: {time.time() - t1:.1f}s")

    t1 = time.time()
    pipeline.filter_genes(min_count=3)
    pipeline.filter_cells(min_count=200, max_count=10000)
    log.info(f"Filtering: {time.time() - t1:.1f}s")

    t1 = time.time()
    pipeline.highly_variable_genes(n_top_genes=n_hvgs)
    log.info(f"HVG selection: {time.time() - t1:.1f}s | {n_hvgs} HVGs")

    t1 = time.time()
    pipeline.normalize_log1p()
    log.info(f"Normalize + log1p: {time.time() - t1:.1f}s")

    t1 = time.time()
    pipeline.pca(n_components=50)
    log.info(f"PCA: {time.time() - t1:.1f}s")

    t1 = time.time()
    pipeline.harmony(sample_col_name="sample", max_iter_harmony=30)
    log.info(f"Harmony batch correction: {time.time() - t1:.1f}s")

    t1 = time.time()
    pipeline.neighbors(n_neighbors=30, n_pcs=50, use_rep="X_pca_harmony")
    log.info(f"Neighbors: {time.time() - t1:.1f}s")

    t1 = time.time()
    pipeline.leiden(resolution=resolution)
    log.info(f"Leiden clustering (res={resolution}): {time.time() - t1:.1f}s")

    t1 = time.time()
    pipeline.umap()
    log.info(f"UMAP: {time.time() - t1:.1f}s")

    # --- Save pre-annotation result ---
    pre_name = f"scalesc_{species}_integrated"
    pipeline.save(data_name=pre_name)
    log.info(f"Saved pre-annotation: {output_dir}/{pre_name}.h5ad")

    # --- Cluster merging + marker identification ---
    t1 = time.time()
    try:
        pipeline.find_cluster_and_merge()
        log.info(f"Cluster merge + markers: {time.time() - t1:.1f}s")
    except Exception as e:
        log.warning(f"Cluster merge failed (non-fatal): {e}")

    # --- Save final annotated result ---
    final_name = f"scalesc_{species}_annotated"
    pipeline.save(data_name=final_name)

    total_time = time.time() - t0
    log.info("=" * 60)
    log.info(f"ScaleSC Integration Complete: {total_time:.1f}s")
    log.info(f"  Species: {species}")
    log.info(f"  Samples: {len(h5ad_files)}")
    log.info(f"  Output: {output_dir}/{final_name}.h5ad")
    log.info("=" * 60)


def main():
    parser = argparse.ArgumentParser(description="ScaleSC single-cell integration")
    parser.add_argument("--species", required=True, choices=["human", "mouse"])
    parser.add_argument("--resolution", type=float, default=0.5)
    parser.add_argument("--n-hvgs", type=int, default=4000)
    parser.add_argument("--max-cell-batch", type=int, default=200000,
                        help="Max cells per GPU batch (default: 200K)")
    args = parser.parse_args()

    gpu_name, gpu_mem = init_gpu()

    run_integration(
        species=args.species,
        resolution=args.resolution,
        n_hvgs=args.n_hvgs,
        max_cell_batch=args.max_cell_batch,
    )


if __name__ == "__main__":
    main()
