#!/usr/bin/env python3
"""Validate rapids_singlecell installation on GPU node.

Tests: GPU detection, RAPIDS imports, basic workflow on synthetic data,
and benchmark vs scanpy on the same data.
"""

import time
import sys

print("=" * 60)
print("  rapids_singlecell Installation Validation")
print("=" * 60)

# --- Step 1: GPU detection ---
print("\n--- GPU Detection ---")
try:
    import cupy as cp
    print(f"CuPy version: {cp.__version__}")
    print(f"CUDA available: {cp.cuda.is_available()}")
    n_gpus = cp.cuda.runtime.getDeviceCount()
    print(f"GPU count: {n_gpus}")
    for i in range(n_gpus):
        props = cp.cuda.runtime.getDeviceProperties(i)
        name = props["name"].decode() if isinstance(props["name"], bytes) else props["name"]
        mem_gb = props["totalGlobalMem"] / 1e9
        print(f"  GPU {i}: {name} ({mem_gb:.1f} GB)")
except Exception as e:
    print(f"FAILED: {e}")
    sys.exit(1)

# --- Step 2: RAPIDS imports ---
print("\n--- RAPIDS Imports ---")
try:
    import rapids_singlecell as rsc
    print(f"rapids_singlecell v{rsc.__version__}")
    import cuml
    print(f"cuML v{cuml.__version__}")
    import scanpy as sc
    print(f"scanpy v{sc.__version__}")
    import anndata as ad
    print(f"anndata v{ad.__version__}")
    import rmm
    print(f"RMM v{rmm.__version__}")
    print("All imports OK")
except Exception as e:
    print(f"FAILED: {e}")
    sys.exit(1)

# --- Step 3: RMM pool allocator ---
print("\n--- RMM Pool Allocator ---")
try:
    from rmm.allocators.cupy import rmm_cupy_allocator
    rmm.reinitialize(pool_allocator=True, initial_pool_size=2**30)  # 1 GB
    cp.cuda.set_allocator(rmm_cupy_allocator)
    print("RMM pool allocator initialized (1 GB)")
except Exception as e:
    print(f"WARNING: Pool allocator failed, using default: {e}")

# --- Step 4: Synthetic data workflow ---
print("\n--- Synthetic Data Workflow (50K cells x 2000 genes) ---")
import numpy as np
from scipy.sparse import random as sp_random

np.random.seed(42)
n_cells = 50_000
n_genes = 2_000

print(f"Creating synthetic AnnData ({n_cells:,} cells x {n_genes:,} genes)...")
X = sp_random(n_cells, n_genes, density=0.1, format="csr", dtype=np.float32)
X.data = np.abs(X.data) * 100  # Make counts-like
adata = ad.AnnData(X=X)
adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
adata.obs_names = [f"Cell_{i}" for i in range(n_cells)]

# GPU workflow
print("\n  GPU workflow (rapids_singlecell):")
adata_gpu = adata.copy()
t0 = time.time()

rsc.get.anndata_to_GPU(adata_gpu)
rsc.pp.filter_genes(adata_gpu, min_cells=3)
rsc.pp.normalize_total(adata_gpu, target_sum=1e4)
rsc.pp.log1p(adata_gpu)
rsc.pp.highly_variable_genes(adata_gpu, n_top_genes=1000, flavor="seurat_v3")
adata_gpu = adata_gpu[:, adata_gpu.var["highly_variable"]].copy()
rsc.pp.scale(adata_gpu, max_value=10)
rsc.pp.pca(adata_gpu, n_comps=50)
rsc.pp.neighbors(adata_gpu)
rsc.tl.umap(adata_gpu)
rsc.tl.leiden(adata_gpu, resolution=0.5)

t_gpu = time.time() - t0
n_clusters_gpu = adata_gpu.obs["leiden"].nunique()
print(f"    Time: {t_gpu:.2f}s")
print(f"    Clusters: {n_clusters_gpu}")
print(f"    UMAP shape: {adata_gpu.obsm['X_umap'].shape}")

# CPU workflow
print("\n  CPU workflow (scanpy):")
adata_cpu = adata.copy()
t0 = time.time()

sc.pp.filter_genes(adata_cpu, min_cells=3)
sc.pp.normalize_total(adata_cpu, target_sum=1e4)
sc.pp.log1p(adata_cpu)
sc.pp.highly_variable_genes(adata_cpu, n_top_genes=1000, flavor="seurat_v3")
adata_cpu = adata_cpu[:, adata_cpu.var["highly_variable"]].copy()
sc.pp.scale(adata_cpu, max_value=10)
sc.pp.pca(adata_cpu, n_comps=50)
sc.pp.neighbors(adata_cpu)
sc.tl.umap(adata_cpu)
sc.tl.leiden(adata_cpu, resolution=0.5)

t_cpu = time.time() - t0
n_clusters_cpu = adata_cpu.obs["leiden"].nunique()
print(f"    Time: {t_cpu:.2f}s")
print(f"    Clusters: {n_clusters_cpu}")

# --- Step 5: Summary ---
print("\n" + "=" * 60)
print("  BENCHMARK SUMMARY")
print("=" * 60)
print(f"  GPU time:  {t_gpu:.2f}s")
print(f"  CPU time:  {t_cpu:.2f}s")
print(f"  Speedup:   {t_cpu / t_gpu:.1f}x")
print(f"  GPU memory used: {cp.cuda.Device(0).mem_info[1] - cp.cuda.Device(0).mem_info[0]:.0f} bytes "
      f"({(cp.cuda.Device(0).mem_info[1] - cp.cuda.Device(0).mem_info[0]) / 1e9:.2f} GB)")
print()
print("  VALIDATION: PASSED")
print("=" * 60)
