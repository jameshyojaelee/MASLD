#!/usr/bin/env python3
"""
16c_run_commot.py — Run COMMOT optimal transport cell-cell communication.

Runs spatial communication inference per condition (Healthy/Steatotic) using
the CellChat LR database filtered in 16b. Uses optimal transport with
distance threshold and regularization parameters from config.

SLURM: --partition=cpu --cpus=16 --mem=128G --time=12:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, load_deconvolved_adata,
    save_checkpoint, save_csv, print_header, print_step,
)


def load_filtered_lr_database():
    """Load filtered LR database from 16b."""
    path = RESULTS_DIR / "commot" / "lr_database_cellchat_filtered.csv"
    if not path.exists():
        print(f"  ERROR: {path} not found. Run 16b_prepare_commot.py first.")
        sys.exit(1)
    df = pd.read_csv(path, index_col=0)
    print(f"  Loaded filtered LR database: {len(df)} pairs")
    return df


def run_commot_per_condition(adata, condition, commot_config, lr_db):
    """Run COMMOT spatial communication for one condition."""
    import commot as ct

    adata_sub = adata[adata.obs["condition"] == condition].copy()
    n_spots = adata_sub.n_obs
    print(f"    Spots: {n_spots}")

    if n_spots < 100:
        print(f"    WARNING: Only {n_spots} spots for {condition}, skipping")
        return None

    # Ensure spatial coordinates exist
    if "spatial" not in adata_sub.obsm:
        print("    ERROR: No spatial coordinates in obsm['spatial']")
        return None

    # Ensure normalized expression (COMMOT expects log-normalized)
    if adata_sub.X.max() > 50:
        print("    Normalizing expression...")
        sc.pp.normalize_total(adata_sub, target_sum=1e4)
        sc.pp.log1p(adata_sub)

    # Run COMMOT spatial communication
    dis_thr = commot_config.get("dis_thr", 500)
    cot_eps_p = commot_config.get("cot_eps_p", 0.1)
    cot_rho = commot_config.get("cot_rho", 10)

    print(f"    Running COMMOT (dis_thr={dis_thr}, eps_p={cot_eps_p}, rho={cot_rho})...")
    ct.tl.spatial_communication(
        adata_sub,
        database_name="CellChat",
        df_ligrec=lr_db,
        dis_thr=dis_thr,
        cot_eps_p=cot_eps_p,
        cot_rho=cot_rho,
        heteromeric=True,
        pathway_sum=True,
    )

    # Report results
    commot_keys = [k for k in adata_sub.obsp.keys() if "commot" in k.lower()]
    uns_keys = [k for k in adata_sub.uns.keys() if "commot" in k.lower()]
    print(f"    COMMOT obsp keys: {len(commot_keys)}")
    print(f"    COMMOT uns keys: {len(uns_keys)}")

    return adata_sub


def main():
    print_header("16c: Run COMMOT Spatial Communication")

    config = load_config()
    commot_config = config["commot"]
    output_dir = RESULTS_DIR / "commot"
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Load data ──
    print_step("Loading data", 1, 3)
    adata = load_deconvolved_adata()
    print(f"  Loaded: {adata.n_obs} spots × {adata.n_vars} genes")

    lr_db = load_filtered_lr_database()

    # ── Run per condition ──
    print_step("Running COMMOT per condition", 2, 3)
    conditions = sorted(adata.obs["condition"].unique().tolist())
    print(f"  Conditions: {conditions}")

    results = {}
    for i, condition in enumerate(conditions):
        print(f"\n  [{i+1}/{len(conditions)}] Condition: {condition}")
        adata_sub = run_commot_per_condition(adata, condition, commot_config, lr_db)
        if adata_sub is not None:
            results[condition] = adata_sub

    # ── Save results ──
    print_step("Saving results", 3, 3)
    for condition, adata_sub in results.items():
        cond_label = condition.replace(" ", "_").replace("/", "_")
        save_checkpoint(adata_sub, f"adata_commot_{cond_label}.h5ad", subdir="commot")

    # ── Summary ──
    print(f"\n  Summary:")
    print(f"    Conditions processed: {len(results)}/{len(conditions)}")
    for condition, adata_sub in results.items():
        commot_keys = [k for k in adata_sub.obsp.keys() if "commot" in k.lower()]
        print(f"    {condition}: {adata_sub.n_obs} spots, "
              f"{len(commot_keys)} communication matrices")

    print_header("16c: Complete")


if __name__ == "__main__":
    main()
