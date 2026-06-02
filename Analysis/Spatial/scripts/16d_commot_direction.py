#!/usr/bin/env python3
"""
16d_commot_direction.py — Compute COMMOT signaling direction vector fields.

For each condition and key MASLD pathway, computes directional sender/receiver
vector fields per spot using COMMOT's communication_direction.

SLURM: --partition=cpu --cpus=8 --mem=64G --time=4:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config,
    save_csv, print_header, print_step,
)
from spatial_stats import get_commot_pathways, match_key_pathways


def load_commot_adata(condition):
    """Load COMMOT results h5ad for a condition."""
    cond_label = condition.replace(" ", "_").replace("/", "_")
    path = RESULTS_DIR / "commot" / f"adata_commot_{cond_label}.h5ad"
    if not path.exists():
        print(f"  WARNING: {path} not found")
        return None
    adata = sc.read_h5ad(path)
    print(f"  Loaded {condition}: {adata.n_obs} spots")
    return adata


def get_available_pathways(adata):
    """Extract pathway names available in COMMOT results.

    Fixed (F081/F082): the real pathway-level matrices live in ``adata.obsp``
    as ``commot-CellChat-<PATHWAY>`` (no 2nd '-'), NOT in obsm. obsm only holds
    the global ``commot-CellChat-sum-sender/-receiver`` aggregate, so the old
    obsm/'sum' parse returned ['sum'] and never the MASLD pathways. Delegates to
    the shared canonical-pathway discovery (excludes LR pairs + 'total'/'sum').
    """
    return get_commot_pathways(adata, db="CellChat")


def compute_direction_for_pathway(adata, pathway, condition):
    """Compute communication direction for a specific pathway."""
    import commot as ct

    try:
        ct.tl.communication_direction(
            adata,
            database_name="CellChat",
            pathway_name=pathway,
            k=5,
        )
    except Exception as e:
        print(f"      WARNING: communication_direction failed for "
              f"{pathway}/{condition}: {e}")
        return None

    # Extract sender/receiver direction vectors from obsm
    sender_key = f"commot_sender_vf-CellChat-{pathway}"
    receiver_key = f"commot_receiver_vf-CellChat-{pathway}"

    records = []
    spatial_coords = adata.obsm["spatial"]

    for i, barcode in enumerate(adata.obs_names):
        record = {
            "barcode": barcode,
            "condition": condition,
            "pathway": pathway,
            "x": spatial_coords[i, 0],
            "y": spatial_coords[i, 1],
        }

        # Extract sender vectors if available
        if sender_key in adata.obsm:
            vf = adata.obsm[sender_key]
            record["sender_vx"] = vf[i, 0]
            record["sender_vy"] = vf[i, 1]
            record["sender_magnitude"] = np.sqrt(vf[i, 0]**2 + vf[i, 1]**2)

        # Extract receiver vectors if available
        if receiver_key in adata.obsm:
            vf = adata.obsm[receiver_key]
            record["receiver_vx"] = vf[i, 0]
            record["receiver_vy"] = vf[i, 1]
            record["receiver_magnitude"] = np.sqrt(vf[i, 0]**2 + vf[i, 1]**2)

        records.append(record)

    return pd.DataFrame(records)


def main():
    print_header("16d: COMMOT Communication Direction")

    config = load_config()
    commot_config = config["commot"]
    key_pathways = commot_config.get("key_pathways", [])
    output_dir = RESULTS_DIR / "commot"
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Load COMMOT results per condition ──
    print_step("Loading COMMOT results", 1, 2)
    conditions = ["Healthy", "Steatotic"]
    adata_dict = {}
    for condition in conditions:
        adata = load_commot_adata(condition)
        if adata is not None:
            adata_dict[condition] = adata

    if not adata_dict:
        print("  ERROR: No COMMOT results found. Run 16c_run_commot.py first.")
        sys.exit(1)

    # ── Compute direction per pathway per condition ──
    print_step("Computing direction vector fields", 2, 2)
    all_direction_dfs = []

    for condition, adata in adata_dict.items():
        available = get_available_pathways(adata)
        print(f"\n  {condition}: {len(available)} pathways available")

        # Filter to key pathways that are available (case-insensitive +
        # WNT-family match; F082: config 'TGFB'/'WNT' vs CellChat 'TGFb'/'ncWNT')
        pathways_to_process = match_key_pathways(key_pathways, available)
        if not pathways_to_process:
            # Fall back to all available pathways
            pathways_to_process = available[:10]  # Cap at 10
            print(f"    No key pathways matched config, using top {len(pathways_to_process)} available")

        for i, pathway in enumerate(pathways_to_process):
            print(f"    [{i+1}/{len(pathways_to_process)}] {pathway}...")
            direction_df = compute_direction_for_pathway(adata, pathway, condition)
            if direction_df is not None and len(direction_df) > 0:
                all_direction_dfs.append(direction_df)

                # Save per pathway/condition
                cond_label = condition.replace(" ", "_").replace("/", "_")
                save_csv(
                    direction_df,
                    f"direction_{pathway}_{cond_label}.csv",
                    subdir="commot",
                )

                # Report summary statistics
                if "sender_magnitude" in direction_df.columns:
                    mean_s = direction_df["sender_magnitude"].mean()
                    print(f"      Mean sender magnitude: {mean_s:.4f}")
                if "receiver_magnitude" in direction_df.columns:
                    mean_r = direction_df["receiver_magnitude"].mean()
                    print(f"      Mean receiver magnitude: {mean_r:.4f}")

    # ── Summary ──
    if all_direction_dfs:
        combined = pd.concat(all_direction_dfs, ignore_index=True)
        save_csv(combined, "direction_all_pathways.csv", subdir="commot")

        print(f"\n  Summary:")
        print(f"    Total direction records: {len(combined)}")
        print(f"    Conditions: {combined['condition'].nunique()}")
        print(f"    Pathways: {combined['pathway'].nunique()}")
        for pw in combined["pathway"].unique():
            pw_df = combined[combined["pathway"] == pw]
            n_conds = pw_df["condition"].nunique()
            print(f"      {pw}: {n_conds} conditions, {len(pw_df)} records")

    print_header("16d: Complete")


if __name__ == "__main__":
    main()
