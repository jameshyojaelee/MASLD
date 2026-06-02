"""GPU acceleration utilities for single-cell analysis.

Provides a unified interface that conditionally uses rapids_singlecell (GPU)
or scanpy (CPU) depending on availability and user preference.

Usage:
    from gpu_utils import init_gpu, get_processor

    use_gpu = init_gpu()  # Returns True if GPU initialized
    pp, tl = get_processor(use_gpu)  # Returns (pp_module, tl_module)

    # Then use pp/tl instead of sc.pp/sc.tl:
    pp.filter_genes(adata, min_cells=3)
    pp.normalize_total(adata)
    tl.umap(adata)
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)


def init_gpu(pool_size: int | None = None, pool_fraction: float = 0.3) -> bool:
    """Initialize GPU with RMM pool allocator.

    Initializes PyTorch CUDA context first (if available) to avoid
    CUBLAS_STATUS_NOT_INITIALIZED conflicts with the RMM memory pool.

    Args:
        pool_size: Initial pool size in bytes. If None, uses pool_fraction of GPU memory.
        pool_fraction: Fraction of GPU memory for RMM pool (default 0.3 = 30%).

    Returns:
        True if GPU initialized successfully, False otherwise.
    """
    try:
        import cupy as cp
        import rmm
        from rmm.allocators.cupy import rmm_cupy_allocator

        n_gpus = cp.cuda.runtime.getDeviceCount()
        if n_gpus == 0:
            logger.warning("No GPUs detected")
            return False

        props = cp.cuda.runtime.getDeviceProperties(0)
        name = props["name"].decode() if isinstance(props["name"], bytes) else props["name"]
        total_mem = props["totalGlobalMem"]
        logger.info(f"GPU: {name} ({total_mem / 1e9:.1f} GB)")

        # Initialize PyTorch CUDA context BEFORE RMM to prevent cuBLAS conflicts
        try:
            import torch
            if torch.cuda.is_available():
                torch.cuda.init()
                torch.zeros(1, device="cuda")  # Force cuBLAS init
                logger.info("PyTorch CUDA context initialized (pre-RMM)")
        except ImportError:
            pass

        if pool_size is None:
            pool_size = int(total_mem * pool_fraction)

        rmm.reinitialize(pool_allocator=True, initial_pool_size=pool_size)
        cp.cuda.set_allocator(rmm_cupy_allocator)
        logger.info(f"RMM pool allocator initialized ({pool_size / 1e9:.1f} GB)")
        return True

    except ImportError:
        logger.info("RAPIDS/CuPy not installed, using CPU")
        return False
    except Exception as e:
        logger.warning(f"GPU initialization failed: {e}")
        return False


def get_processor(use_gpu: bool = False):
    """Return (pp_module, tl_module) for either GPU or CPU processing.

    Args:
        use_gpu: If True, return rapids_singlecell modules.

    Returns:
        Tuple of (preprocessing_module, tools_module).
        For GPU: (rsc.pp, rsc.tl)
        For CPU: (sc.pp, sc.tl)
    """
    if use_gpu:
        import rapids_singlecell as rsc
        return rsc.pp, rsc.tl
    else:
        import scanpy as sc
        return sc.pp, sc.tl


def to_gpu(adata: AnnData) -> AnnData:
    """Transfer AnnData to GPU memory."""
    import rapids_singlecell as rsc
    rsc.get.anndata_to_GPU(adata)
    return adata


def from_gpu(adata: AnnData) -> AnnData:
    """Transfer AnnData back to CPU memory."""
    import rapids_singlecell as rsc
    rsc.get.anndata_to_CPU(adata)
    return adata
