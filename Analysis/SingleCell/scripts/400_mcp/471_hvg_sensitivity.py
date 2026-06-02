#!/usr/bin/env python
"""
471_hvg_sensitivity.py — cNMF stability across HVG counts (1k/2k/4k/8k).

Driver script. Submits sbatch jobs with --export CNMF_NAME=global_hvg_{n} CNMF_NHVG={n}.
After completion, aggregates program recovery Jaccard via the same function as 470_loo_stability.

Usage:
  python 471_hvg_sensitivity.py --stage submit --k 20
  python 471_hvg_sensitivity.py --stage aggregate --k 20
"""
from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
SCRIPT_DIR = PROJECT_ROOT / "Analysis/SingleCell/scripts/400_mcp"
CNMF_ROOT = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/cnmf_runs"
BENCH = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/benchmarks"
BENCH.mkdir(parents=True, exist_ok=True)

HVG_VALUES = [1000, 2000, 4000, 8000]


def submit(k: int) -> None:
    jobs = []
    for n in HVG_VALUES:
        name = f"global_hvg_{n}"
        prep = subprocess.run([
            "sbatch",
            "--export", f"CNMF_INPUT=atlas_cnmf_global.h5ad,CNMF_NAME={name},CNMF_KLIST={k},CNMF_NITER=150,CNMF_NHVG={n}",
            str(SCRIPT_DIR / "420_cnmf_prepare.sbatch"),
        ], capture_output=True, text=True)
        prep_j = prep.stdout.strip().split()[-1]
        fac = subprocess.run([
            "sbatch", f"--dependency=afterok:{prep_j}",
            "--export", f"CNMF_NAME={name},TOTAL_WORKERS=32",
            str(SCRIPT_DIR / "420_cnmf_factorize.sbatch"),
        ], capture_output=True, text=True)
        fac_j = fac.stdout.strip().split()[-1]
        cons = subprocess.run([
            "sbatch", f"--dependency=afterok:{fac_j}",
            "--export", f"CNMF_NAME={name},CNMF_KLIST={k},CNMF_DTHRESH=0.03",
            str(SCRIPT_DIR / "420_cnmf_consensus.sbatch"),
        ], capture_output=True, text=True)
        jobs.append((name, prep_j, fac_j, cons.stdout.strip().split()[-1]))
    pd.DataFrame(jobs, columns=["name", "prepare_jid", "factorize_jid", "consensus_jid"]).to_csv(
        BENCH / "hvg_job_ids.tsv", sep="\t", index=False)
    print("[471] submitted:", jobs)


def aggregate(k: int, dthresh: float = 0.03) -> None:
    dt_tag = f"{dthresh:.2f}".replace(".", "_")
    full_f = CNMF_ROOT / "global" / "global" / f"global.gene_spectra_score.k_{k}.dt_{dt_tag}.txt"
    if not full_f.exists():
        full_f = CNMF_ROOT / "global" / f"global.gene_spectra_score.k_{k}.dt_{dt_tag}.txt"
    full = pd.read_csv(full_f, sep="\t", index_col=0)
    top_full = {p: set(full.loc[p].nlargest(100).index) for p in full.index}

    rows = []
    for n in HVG_VALUES:
        name = f"global_hvg_{n}"
        fn = CNMF_ROOT / name / name / f"{name}.gene_spectra_score.k_{k}.dt_{dt_tag}.txt"
        if not fn.exists():
            fn = CNMF_ROOT / name / f"{name}.gene_spectra_score.k_{k}.dt_{dt_tag}.txt"
        if not fn.exists():
            continue
        alt = pd.read_csv(fn, sep="\t", index_col=0)
        top_alt = {p: set(alt.loc[p].nlargest(100).index) for p in alt.index}
        for pf, sf in top_full.items():
            best = max(
                (len(sf & sa) / max(1, len(sf | sa)) for sa in top_alt.values()),
                default=0,
            )
            rows.append({"n_hvg": n, "full_program": pf, "jaccard_top100_best": best})
    pd.DataFrame(rows).to_csv(BENCH / f"hvg_sensitivity_k{k}.tsv", sep="\t", index=False)
    print(f"[471] wrote hvg_sensitivity_k{k}.tsv")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["submit", "aggregate"], required=True)
    ap.add_argument("--k", type=int, required=True)
    args = ap.parse_args()
    if args.stage == "submit":
        submit(args.k)
    else:
        aggregate(args.k)
