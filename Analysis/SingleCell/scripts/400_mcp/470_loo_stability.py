#!/usr/bin/env python
"""
470_loo_stability.py — Dataset leave-one-out cNMF stability.

For each of the 7 scRNA datasets:
  1. Re-run cNMF at the chosen k excluding that dataset's cells
  2. Compute gene-spectra Pearson r vs full-fit programs
  3. Take max-r assignment (Hungarian or greedy) → program recovery Jaccard of top-100 genes

Outputs:
  results_gpu_v2/mcp/benchmarks/loo_stability.tsv
  results_gpu_v2/mcp/benchmarks/loo_recovery_heatmap.tsv

This script is a DRIVER — it launches 7 sub-jobs to SLURM (one per dataset) using the standard
cnmf prepare/factorize/consensus flow, then aggregates. Each sub-run writes to
cnmf_runs/global_loo_{dataset}/.
"""
from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
SCRIPT_DIR = PROJECT_ROOT / "Analysis/SingleCell/scripts/400_mcp"
INPUTS = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs"
CNMF_ROOT = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/cnmf_runs"
BENCH = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/benchmarks"
BENCH.mkdir(parents=True, exist_ok=True)


def build_loo_inputs():
    atlas = ad.read_h5ad(INPUTS / "atlas_cnmf_global.h5ad")
    datasets = sorted(atlas.obs["dataset"].astype(str).unique())
    made = []
    for d in datasets:
        sub = atlas[atlas.obs["dataset"] != d].copy()
        out = INPUTS / f"atlas_cnmf_global_loo_{d}.h5ad"
        if not out.exists():
            sub.write_h5ad(out, compression="gzip")
        made.append((d, sub.shape[0]))
        print(f"[470] LOO {d}: kept {sub.shape[0]} cells -> {out.name}")
    return datasets


def submit_loo_jobs(datasets, k):
    jobs = []
    for d in datasets:
        # Prepare
        cmd_prep = [
            "sbatch",
            "--export", f"CNMF_INPUT=atlas_cnmf_global_loo_{d}.h5ad,CNMF_NAME=global_loo_{d},CNMF_KLIST={k},CNMF_NITER=150,CNMF_NHVG=2500",
            str(SCRIPT_DIR / "420_cnmf_prepare.sbatch"),
        ]
        r = subprocess.run(cmd_prep, capture_output=True, text=True)
        prep_jid = r.stdout.strip().split()[-1]
        # Factorize
        cmd_fac = [
            "sbatch",
            f"--dependency=afterok:{prep_jid}",
            "--export", f"CNMF_NAME=global_loo_{d},TOTAL_WORKERS=32",
            str(SCRIPT_DIR / "420_cnmf_factorize.sbatch"),
        ]
        r2 = subprocess.run(cmd_fac, capture_output=True, text=True)
        fac_jid = r2.stdout.strip().split()[-1]
        # Consensus
        cmd_cons = [
            "sbatch",
            f"--dependency=afterok:{fac_jid}",
            "--export", f"CNMF_NAME=global_loo_{d},CNMF_KLIST={k},CNMF_DTHRESH=0.03",
            str(SCRIPT_DIR / "420_cnmf_consensus.sbatch"),
        ]
        r3 = subprocess.run(cmd_cons, capture_output=True, text=True)
        jobs.append((d, prep_jid, fac_jid, r3.stdout.strip().split()[-1]))
        print(f"[470] LOO {d}: jobs {prep_jid}, {fac_jid}, {r3.stdout.strip().split()[-1]}")
    pd.DataFrame(jobs, columns=["dataset", "prepare_jid", "factorize_jid", "consensus_jid"]).to_csv(
        BENCH / "loo_job_ids.tsv", sep="\t", index=False,
    )


def aggregate(datasets, k, dthresh=0.03):
    """Compare LOO programs to full-fit programs, compute Jaccard top-100."""
    dt_tag = f"{dthresh:.2f}".replace(".", "_")
    full_f = CNMF_ROOT / "global" / "global" / f"global.gene_spectra_score.k_{k}.dt_{dt_tag}.txt"
    if not full_f.exists():
        full_f = CNMF_ROOT / "global" / f"global.gene_spectra_score.k_{k}.dt_{dt_tag}.txt"
    full = pd.read_csv(full_f, sep="\t", index_col=0)
    top_full = {p: set(full.loc[p].nlargest(100).index) for p in full.index}

    rows = []
    for d in datasets:
        loo_f = CNMF_ROOT / f"global_loo_{d}" / f"global_loo_{d}" / f"global_loo_{d}.gene_spectra_score.k_{k}.dt_{dt_tag}.txt"
        if not loo_f.exists():
            loo_f = CNMF_ROOT / f"global_loo_{d}" / f"global_loo_{d}.gene_spectra_score.k_{k}.dt_{dt_tag}.txt"
        if not loo_f.exists():
            print(f"[470] LOO {d} k={k} missing — skip")
            continue
        loo = pd.read_csv(loo_f, sep="\t", index_col=0)
        top_loo = {p: set(loo.loc[p].nlargest(100).index) for p in loo.index}
        # Greedy Hungarian: for each full program, find best LOO match
        for pf, sf in top_full.items():
            best_j = 0
            best_lp = None
            for pl, sl in top_loo.items():
                j = len(sf & sl) / max(1, len(sf | sl))
                if j > best_j:
                    best_j, best_lp = j, pl
            rows.append({"dataset_left_out": d, "full_program": pf,
                         "loo_best_match": best_lp, "jaccard_top100": best_j})
    pd.DataFrame(rows).to_csv(BENCH / f"loo_stability_k{k}.tsv", sep="\t", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--stage", choices=["setup", "submit", "aggregate"], default="aggregate")
    ap.add_argument(
        "--datasets",
        type=str,
        default="",
        help="Comma-separated dataset list to restrict to. Empty = all 7.",
    )
    ap.add_argument(
        "--exclude-dataset",
        type=str,
        default="",
        help="Comma-separated dataset list to exclude from the run (e.g. already done).",
    )
    args = ap.parse_args()
    all_ds = sorted(pd.read_csv(INPUTS / "donor_metadata.tsv", sep="\t")["dataset"].astype(str).unique())
    if args.datasets:
        keep = [d.strip() for d in args.datasets.split(",") if d.strip()]
        datasets = [d for d in all_ds if d in keep]
    else:
        datasets = list(all_ds)
    if args.exclude_dataset:
        drop = {d.strip() for d in args.exclude_dataset.split(",") if d.strip()}
        datasets = [d for d in datasets if d not in drop]
    if args.stage == "setup":
        build_loo_inputs()
    elif args.stage == "submit":
        submit_loo_jobs(datasets, args.k)
    elif args.stage == "aggregate":
        aggregate(datasets, args.k)
