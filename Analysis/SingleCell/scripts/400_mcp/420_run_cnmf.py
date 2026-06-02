#!/usr/bin/env python
"""
420_run_cnmf.py — unified driver for cNMF prepare / factorize / consensus stages.

Usage:
  # 1. Prepare (single job, once per target AnnData)
  python 420_run_cnmf.py prepare --input atlas_cnmf_global.h5ad --name global --k 10,12,14,16,18,20,22,25,28,32,36,40 --n-iter 200 --num-hvg 2500

  # 2. Factorize (run as SLURM array; each task processes one worker-index)
  python 420_run_cnmf.py factorize --name global --worker-index $SLURM_ARRAY_TASK_ID --total-workers 32

  # 3. Consensus (after all factorize workers finish)
  python 420_run_cnmf.py consensus --name global --k 10,12,14,16,18,20,22,25,28,32,36,40 --density-threshold 0.03
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
CNMF_ROOT = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/cnmf_runs"
CNMF_ROOT.mkdir(parents=True, exist_ok=True)
INPUT_DIR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs"


def _parse_klist(arg: str) -> list[int]:
    # Accept comma OR colon separated (colons are safe through sbatch --export which
    # itself uses commas). Whitespace also acts as a separator.
    import re
    return [int(x) for x in re.split(r"[,\s:]+", arg.strip()) if x]


def cmd_prepare(args: argparse.Namespace) -> None:
    from cnmf import cNMF
    import scanpy as sc
    import numpy as np
    import scipy.sparse as sp

    input_path = INPUT_DIR / args.input

    # cnmf.prepare internally normalizes by dividing by library size, assigning back
    # into adata.X. If X is int32/int64 the ufunc-divide fails:
    #   "Cannot cast ufunc 'divide' output from dtype('float64') to dtype('int32')"
    # Pre-cast to float32 if needed and write a sibling .f32.h5ad (do not mutate input).
    adata = sc.read(str(input_path))
    x = adata.X
    target_dtype = None
    if sp.issparse(x):
        if x.dtype.kind in "iu":  # integer sparse
            target_dtype = np.float32
    else:
        if x.dtype.kind in "iu":
            target_dtype = np.float32
    if target_dtype is not None:
        print(f"[420.prepare] casting adata.X from {x.dtype} to {target_dtype.__name__} "
              f"(cnmf.prepare requires float for library-size normalization)")
        if sp.issparse(x):
            adata.X = x.astype(target_dtype)
        else:
            adata.X = x.astype(target_dtype, copy=False)
        f32_path = input_path.with_suffix(".f32.h5ad")
        if not f32_path.exists() or f32_path.stat().st_mtime < input_path.stat().st_mtime:
            print(f"[420.prepare] writing float-cast h5ad -> {f32_path.name}")
            adata.write_h5ad(str(f32_path), compression="gzip")
        effective_input = f32_path
    else:
        effective_input = input_path

    cobj = cNMF(output_dir=str(CNMF_ROOT), name=args.name)
    cobj.prepare(
        counts_fn=str(effective_input),
        components=_parse_klist(args.k),
        n_iter=args.n_iter,
        seed=args.seed,
        num_highvar_genes=args.num_hvg,
        densify=False,
    )
    print(f"[420.prepare] Done. Outputs in: {CNMF_ROOT / args.name}")


def cmd_factorize(args: argparse.Namespace) -> None:
    from cnmf import cNMF
    cobj = cNMF(output_dir=str(CNMF_ROOT), name=args.name)
    cobj.factorize(worker_i=args.worker_index, total_workers=args.total_workers)
    print(f"[420.factorize] Worker {args.worker_index}/{args.total_workers} done.")


def cmd_consensus(args: argparse.Namespace) -> None:
    from cnmf import cNMF
    import scanpy as sc
    import numpy as np
    cobj = cNMF(output_dir=str(CNMF_ROOT), name=args.name)
    # Only combine the requested k values (skip missing iteration files for cancelled k)
    cobj.combine(components=_parse_klist(args.k), skip_missing_files=True)
    # Pre-load norm_counts and cast X to float64 (sklearn requires H.dtype == X.dtype)
    norm_counts = sc.read(cobj.paths['normalized_counts'])
    if norm_counts.X.dtype != np.float64:
        print(f"[420.consensus] casting norm_counts.X from {norm_counts.X.dtype} to float64")
        norm_counts.X = norm_counts.X.astype(np.float64)
    k_list = _parse_klist(args.k)
    for k in k_list:
        print(f"[420.consensus] k={k}")
        try:
            cobj.consensus(k=k, density_threshold=args.density_threshold, show_clustering=False,
                           close_clustergram_fig=True, norm_counts=norm_counts)
        except Exception as e:
            print(f"[420.consensus] k={k} FAILED: {e}")
    # Stability plot
    try:
        cobj.k_selection_plot(close_fig=True)
    except Exception as e:
        print(f"[420.consensus] k-selection plot failed: {e}")
    print("[420.consensus] DONE.")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="stage", required=True)

    p = sub.add_parser("prepare")
    p.add_argument("--input", required=True, help="Input h5ad filename under inputs/")
    p.add_argument("--name", required=True)
    p.add_argument("--k", required=True, help="Comma-separated list of k values")
    p.add_argument("--n-iter", type=int, default=200)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--num-hvg", type=int, default=2500)

    f = sub.add_parser("factorize")
    f.add_argument("--name", required=True)
    f.add_argument("--worker-index", type=int, required=True)
    f.add_argument("--total-workers", type=int, required=True)

    c = sub.add_parser("consensus")
    c.add_argument("--name", required=True)
    c.add_argument("--k", required=True)
    c.add_argument("--density-threshold", type=float, default=0.03)

    args = ap.parse_args(argv)
    if args.stage == "prepare":
        cmd_prepare(args)
    elif args.stage == "factorize":
        cmd_factorize(args)
    elif args.stage == "consensus":
        cmd_consensus(args)


if __name__ == "__main__":
    main()
