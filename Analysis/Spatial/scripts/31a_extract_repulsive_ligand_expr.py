#!/usr/bin/env python3
"""
31a_extract_repulsive_ligand_expr.py — per-spot log-normalized ligand expression.

Companion extractor for 31_spatial_immune_exclusion.R (F121). The immune-
exclusion repulsive-ligand test needs spot-level expression of IDO1/TGFB1/
VEGFA/CXCL12 on the LOG1P-CPM scale, but the deconvolved h5ad `.X` is raw
integer counts (cell2location overwrites .X with the counts layer) and no R
script in this pipeline reads h5ad. Following the pipeline convention (Python
extracts -> CSV -> R consumes), this writes a small per-spot matrix:

    immune_exclusion/repulsive_ligand_spot_expr.csv
      columns: spot, sample_id, condition, <one column per FOUND ligand>

Missing ligands (e.g. IDO1 is filtered out of GSE192741's var) are skipped and
reported, so 31.R applies BH only over genes that actually exist.

Env: spatial (needs scanpy via spatial_stats.ensure_lognorm). Pure data-prep;
does not run any pipeline / snakemake / sbatch.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import RESULTS_DIR, load_deconvolved_adata
from spatial_stats import ensure_lognorm, _get_expr

# Repulsive / immune-excluding ligands tested in 31.R step 5.
REPULSIVE_LIGANDS = ["IDO1", "TGFB1", "VEGFA", "CXCL12"]
DONOR_COL = "sample_id"


def main():
    adata = load_deconvolved_adata()
    print(f"  Loaded {adata.n_obs} spots, {adata.obs[DONOR_COL].nunique()} donors")

    # F038-style: read expression from a log1p-CPM layer, never raw .X counts.
    lognorm_layer = ensure_lognorm(adata)

    present = [g for g in REPULSIVE_LIGANDS if g in set(adata.var_names)]
    missing = [g for g in REPULSIVE_LIGANDS if g not in set(adata.var_names)]
    if missing:
        print(f"  NOTE: ligands absent from var (skipped): {missing}")
    if not present:
        raise SystemExit("No repulsive ligands present in the dataset var_names.")

    out = pd.DataFrame({
        "spot": adata.obs_names.astype(str),
        "sample_id": adata.obs[DONOR_COL].astype(str).values,
        "condition": adata.obs["condition"].astype(str).values,
    })
    for g in present:
        out[g] = _get_expr(adata, g, layer=lognorm_layer)

    outdir = RESULTS_DIR / "immune_exclusion"
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "repulsive_ligand_spot_expr.csv"
    out.to_csv(outpath, index=False)
    print(f"  Wrote {outpath} ({len(out)} spots, ligands: {present})")


if __name__ == "__main__":
    main()
