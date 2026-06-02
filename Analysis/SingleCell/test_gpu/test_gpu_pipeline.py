#!/usr/bin/env python3
"""End-to-end GPU pipeline test using existing Liver_Atlas samples.

Tests three scenarios:
  1. Per-sample scanpy_secondary.py (GPU vs CPU comparison)
  2. Multi-sample integration preprocessing (GPU vs CPU)
  3. Mini scVI training on GPU-preprocessed data

All outputs go to test_gpu/ — never overwrites production results.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time
from pathlib import Path

import numpy as np
import scanpy as sc

# Import GPU utilities
SCRIPT_DIR = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))
from gpu_utils import init_gpu, get_processor, to_gpu, from_gpu

# --- Seed pinning (added 2026-04-22 per T0.8) -------------------------------
import random
os.environ['PYTHONHASHSEED'] = '42'
random.seed(42)
np.random.seed(42)
import torch
torch.manual_seed(42)
torch.cuda.manual_seed_all(42)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
try:
    import scvi
    scvi.settings.seed = 42
except Exception:
    pass
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

# --- Test samples (GSE192740 human, from Liver_Atlas) ---
# Mix of large (~27M, snRNA nuclei) and small (~360K, sorted) samples
PROJECT_ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CELLRANGER_DIR = PROJECT_ROOT / "Liver_Atlas" / "cellranger" / "human"

TEST_SAMPLES = [
    "SRR17375012",  # ~27M, nuclei (large, ~10K+ cells)
    "SRR17375014",  # ~27M, nuclei (large)
    "SRR17375013",  # ~360K, sorted (small, ~500 cells)
]


def find_h5(sample_id: str) -> Path:
    p = CELLRANGER_DIR / sample_id / "outs" / "filtered_feature_bc_matrix.h5"
    if not p.exists():
        raise FileNotFoundError(f"Missing: {p}")
    return p


def load_samples(sample_ids: list[str]) -> sc.AnnData:
    """Load and concatenate multiple 10x samples."""
    adatas = []
    for sid in sample_ids:
        h5_path = find_h5(sid)
        log.info(f"Loading {sid} from {h5_path}")
        adata = sc.read_10x_h5(str(h5_path))
        adata.var_names_make_unique()
        adata.obs["sample"] = sid
        adata.obs["dataset"] = "GSE192740"
        log.info(f"  {sid}: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        adatas.append(adata)

    combined = sc.concat(adatas, join="outer", fill_value=0)
    combined.obs_names_make_unique()
    log.info(f"Combined: {combined.n_obs:,} cells × {combined.n_vars:,} genes")
    return combined


# ─── Test 1: Per-sample scanpy_secondary (GPU vs CPU) ────────────────────────

def test_per_sample(output_dir: Path, use_gpu: bool) -> dict:
    """Run scanpy_secondary-style pipeline on a single sample."""
    mode = "GPU" if use_gpu else "CPU"
    log.info(f"\n{'='*60}")
    log.info(f"TEST 1: Per-sample pipeline ({mode})")
    log.info(f"{'='*60}")

    pp, tl = get_processor(use_gpu)
    sample_id = TEST_SAMPLES[0]  # Large sample for meaningful benchmark
    adata = sc.read_10x_h5(str(find_h5(sample_id)))
    adata.var_names_make_unique()
    adata.obs["sample"] = sample_id

    t0 = time.time()

    # QC
    mt_mask = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["mt"] = mt_mask
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs.n_genes_by_counts < 6000, :].copy()
    adata = adata[adata.obs.pct_counts_mt < 20, :].copy()
    log.info(f"  After QC: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    if adata.n_obs < 50:
        log.warning(f"  Too few cells ({adata.n_obs}), skipping")
        return {"mode": mode, "cells": adata.n_obs, "time": 0, "status": "skipped"}

    # GPU transfer
    if use_gpu:
        to_gpu(adata)

    # Preprocessing
    pp.normalize_total(adata, target_sum=1e4)
    pp.log1p(adata)

    n_top = min(2000, adata.n_vars - 1)
    pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="cell_ranger")

    if use_gpu:
        from_gpu(adata)
    adata = adata[:, adata.var.highly_variable].copy()
    if use_gpu:
        to_gpu(adata)

    pp.scale(adata, max_value=10)
    pp.pca(adata)
    pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    tl.umap(adata, min_dist=0.3)
    tl.leiden(adata, resolution=0.5)

    if use_gpu:
        from_gpu(adata)

    elapsed = time.time() - t0

    # Save
    out_file = output_dir / f"test1_per_sample_{mode.lower()}.h5ad"
    adata.write_h5ad(str(out_file))
    log.info(f"  Saved: {out_file}")
    log.info(f"  Time: {elapsed:.2f}s | Cells: {adata.n_obs:,} | Clusters: {adata.obs['leiden'].nunique()}")

    return {
        "mode": mode,
        "cells": adata.n_obs,
        "genes": adata.n_vars,
        "clusters": int(adata.obs["leiden"].nunique()),
        "time": elapsed,
        "status": "ok",
    }


# ─── Test 2: Multi-sample integration preprocessing (GPU vs CPU) ─────────────

def test_integration_preprocess(output_dir: Path, use_gpu: bool) -> dict:
    """Run multi-sample preprocessing (the GPU-heavy part of scvi_integration)."""
    mode = "GPU" if use_gpu else "CPU"
    log.info(f"\n{'='*60}")
    log.info(f"TEST 2: Multi-sample integration preprocessing ({mode})")
    log.info(f"{'='*60}")

    pp, tl = get_processor(use_gpu)
    adata = load_samples(TEST_SAMPLES)

    t0 = time.time()

    # QC
    mt_mask = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["mt"] = mt_mask
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs.pct_counts_mt < 20, :].copy()
    log.info(f"  After QC: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # Store raw counts
    adata.layers["counts"] = adata.X.copy()

    # HVG selection on CPU (seurat_v3 with layer="counts" needs scipy sparse)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    n_top = min(4000, adata.n_vars - 1)
    try:
        sc.pp.highly_variable_genes(
            adata, n_top_genes=n_top, flavor="seurat_v3",
            batch_key="dataset", layer="counts",
        )
    except Exception as e:
        log.warning(f"  seurat_v3 HVG failed ({e}), trying cell_ranger flavor")
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="cell_ranger")

    n_hvg = adata.var.highly_variable.sum()
    log.info(f"  HVG selected: {n_hvg:,}")
    adata_hvg = adata[:, adata.var.highly_variable].copy()

    # GPU for compute-heavy operations (PCA, neighbors, UMAP, leiden)
    if use_gpu:
        to_gpu(adata_hvg)

    pp.scale(adata_hvg, max_value=10)
    pp.pca(adata_hvg)
    pp.neighbors(adata_hvg, n_neighbors=30, n_pcs=30)
    tl.umap(adata_hvg)
    tl.leiden(adata_hvg, resolution=1.0)

    if use_gpu:
        from_gpu(adata_hvg)

    elapsed = time.time() - t0

    # Save
    out_file = output_dir / f"test2_integration_preprocess_{mode.lower()}.h5ad"
    adata_hvg.write_h5ad(str(out_file))
    log.info(f"  Saved: {out_file}")
    log.info(f"  Time: {elapsed:.2f}s | Cells: {adata_hvg.n_obs:,} | HVG: {n_hvg:,} | Clusters: {adata_hvg.obs['leiden'].nunique()}")

    return {
        "mode": mode,
        "cells": adata_hvg.n_obs,
        "hvg": int(n_hvg),
        "clusters": int(adata_hvg.obs["leiden"].nunique()),
        "time": elapsed,
        "status": "ok",
    }


# ─── Test 3: Mini scVI training ──────────────────────────────────────────────

def test_mini_scvi(output_dir: Path, use_gpu: bool) -> dict:
    """Quick scVI training (5 epochs) to verify GPU tensor operations."""
    mode = "GPU" if use_gpu else "CPU"
    log.info(f"\n{'='*60}")
    log.info(f"TEST 3: Mini scVI training ({mode})")
    log.info(f"{'='*60}")

    try:
        import scvi
    except ImportError:
        log.warning("  scvi-tools not installed, skipping")
        return {"mode": mode, "status": "skipped", "time": 0}

    pp, tl = get_processor(use_gpu)
    adata = load_samples(TEST_SAMPLES[:2])  # Use 2 samples

    # Minimal QC
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    mt_mask = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["mt"] = mt_mask
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    adata = adata[adata.obs.pct_counts_mt < 20, :].copy()
    log.info(f"  After QC: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # HVG (on CPU — scVI needs counts layer)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    n_top = min(2000, adata.n_vars - 1)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top)
    adata = adata[:, adata.var.highly_variable].copy()
    log.info(f"  HVG subset: {adata.n_obs:,} × {adata.n_vars:,}")

    t0 = time.time()

    # scVI setup + train (PyTorch handles its own GPU)
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample")
    model = scvi.model.SCVI(adata, n_latent=10, n_hidden=64, n_layers=1)

    import torch
    use_cuda = use_gpu and torch.cuda.is_available()
    log.info(f"  scVI training on {'GPU (CUDA)' if use_cuda else 'CPU'}")
    model.train(max_epochs=5, accelerator="gpu" if use_cuda else "cpu")

    # Get latent representation (returns numpy on CPU)
    latent = model.get_latent_representation()
    adata.obsm["X_scVI"] = latent
    log.info(f"  Latent shape: {latent.shape}")

    # Free PyTorch GPU cache before rapids_singlecell ops
    if use_cuda:
        torch.cuda.empty_cache()

    # Post-scVI clustering (GPU-accelerated via rapids_singlecell)
    if use_gpu:
        to_gpu(adata)
    pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)
    tl.umap(adata)
    tl.leiden(adata, resolution=0.5)
    if use_gpu:
        from_gpu(adata)

    elapsed = time.time() - t0

    # Save
    out_file = output_dir / f"test3_mini_scvi_{mode.lower()}.h5ad"
    model_dir = output_dir / f"test3_scvi_model_{mode.lower()}"
    adata.write_h5ad(str(out_file))
    model.save(str(model_dir), overwrite=True)
    log.info(f"  Saved: {out_file}")
    log.info(f"  Model: {model_dir}")
    log.info(f"  Time: {elapsed:.2f}s | Clusters: {adata.obs['leiden'].nunique()}")

    return {
        "mode": mode,
        "cells": adata.n_obs,
        "latent_dim": latent.shape[1],
        "clusters": int(adata.obs["leiden"].nunique()),
        "time": elapsed,
        "status": "ok",
    }


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default=str(PROJECT_ROOT / "RNA-seq" / "Analysis" / "SingleCell" / "test_gpu"),
        help="Output directory for test results",
    )
    parser.add_argument(
        "--tests", nargs="+", type=int, default=[1, 2, 3],
        choices=[1, 2, 3],
        help="Which tests to run (default: all)",
    )
    parser.add_argument(
        "--gpu-only", action="store_true",
        help="Skip CPU comparison runs (faster)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Verify test samples exist
    for sid in TEST_SAMPLES:
        find_h5(sid)
    log.info(f"All {len(TEST_SAMPLES)} test samples found")

    # Initialize GPU
    use_gpu = init_gpu()
    if not use_gpu:
        log.error("GPU initialization failed — cannot run GPU tests")
        sys.exit(1)
    log.info("GPU initialized successfully")

    results = []

    test_funcs = {
        1: ("Per-sample pipeline", test_per_sample),
        2: ("Integration preprocessing", test_integration_preprocess),
        3: ("Mini scVI training", test_mini_scvi),
    }

    for test_id in args.tests:
        name, func = test_funcs[test_id]

        # GPU run
        log.info(f"\n>>> Running Test {test_id} ({name}) — GPU")
        try:
            gpu_result = func(output_dir, use_gpu=True)
            gpu_result["test"] = name
            results.append(gpu_result)
        except Exception as e:
            log.error(f"  GPU test {test_id} FAILED: {e}", exc_info=True)
            results.append({"test": name, "mode": "GPU", "status": f"FAILED: {e}", "time": 0})

        # CPU comparison (unless --gpu-only)
        if not args.gpu_only:
            log.info(f"\n>>> Running Test {test_id} ({name}) — CPU")
            try:
                cpu_result = func(output_dir, use_gpu=False)
                cpu_result["test"] = name
                results.append(cpu_result)
            except Exception as e:
                log.error(f"  CPU test {test_id} FAILED: {e}", exc_info=True)
                results.append({"test": name, "mode": "CPU", "status": f"FAILED: {e}", "time": 0})

    # ─── Summary ──────────────────────────────────────────────────────────

    log.info(f"\n{'='*70}")
    log.info("RESULTS SUMMARY")
    log.info(f"{'='*70}")
    log.info(f"{'Test':<35} {'Mode':<5} {'Cells':>8} {'Time':>8} {'Speedup':>8} {'Status'}")
    log.info(f"{'-'*35} {'-'*5} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")

    # Group by test name for speedup calc
    by_test = {}
    for r in results:
        by_test.setdefault(r["test"], {})[r["mode"]] = r

    for test_name, modes in by_test.items():
        for mode in ["GPU", "CPU"]:
            if mode not in modes:
                continue
            r = modes[mode]
            cells = r.get("cells", "—")
            t = r.get("time", 0)
            status = r.get("status", "—")

            speedup = "—"
            if mode == "GPU" and "CPU" in modes:
                cpu_t = modes["CPU"].get("time", 0)
                if cpu_t > 0 and t > 0:
                    speedup = f"{cpu_t / t:.1f}x"

            cells_str = f"{cells:,}" if isinstance(cells, int) else str(cells)
            log.info(f"{test_name:<35} {mode:<5} {cells_str:>8} {t:>7.1f}s {speedup:>8} {status}")

    # Save results
    import json
    results_file = output_dir / "test_results.json"
    with open(results_file, "w") as f:
        json.dump(results, f, indent=2, default=str)
    log.info(f"\nResults saved: {results_file}")
    log.info(f"Output dir: {output_dir}")

    # Check for failures
    failures = [r for r in results if "FAILED" in str(r.get("status", ""))]
    if failures:
        log.error(f"\n{len(failures)} test(s) FAILED!")
        return 1
    log.info("\nAll tests PASSED!")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
