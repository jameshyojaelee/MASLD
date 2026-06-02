#!/usr/bin/env python3
"""
07_run_da_test.py — Standalone DA test on SnapATAC2 processed h5ad.

Fixes non-unique obs_names and runs diff_test per cell type.
"""
import snapatac2 as snap
import anndata as ad
import numpy as np
import pandas as pd
import os
import sys
import argparse
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="snapatac2_processed.h5ad")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    log.info("Loading %s", args.input)
    adata = ad.read_h5ad(args.input)
    log.info("Shape: %d x %d", adata.n_obs, adata.n_vars)
    log.info("obs columns: %s", list(adata.obs.columns))

    # Fix non-unique obs_names (in-memory only; save separately to avoid race)
    if not adata.obs_names.is_unique:
        n_dup = adata.n_obs - adata.obs_names.nunique()
        log.info("Fixing %d duplicate obs_names (in-memory)", n_dup)
        adata.obs_names_make_unique()

    # Save fixed version to separate file (avoid corrupting original during read)
    fixed_path = args.input.replace(".h5ad", "_fixed.h5ad")
    adata.write(fixed_path)
    log.info("Saved fixed h5ad: %s", fixed_path)

    # Identify groups
    conditions = adata.obs["condition"].unique().tolist()
    log.info("Conditions: %s", conditions)

    masld_labels = [c for c in conditions
                    if str(c).upper() in ["MASLD", "NASH", "NAFLD", "NAFL", "MASH", "MASL"]]
    normal_labels = [c for c in conditions
                     if str(c).upper() in ["NORMAL", "HEALTHY", "CONTROL"]]

    if not masld_labels or not normal_labels:
        log.error("Cannot identify groups: %s", conditions)
        sys.exit(1)

    log.info("MASLD: %s, Normal: %s", masld_labels, normal_labels)

    cell_types = adata.obs["cell_type"].unique().tolist()
    all_results = []

    for ct in cell_types:
        ct_mask = adata.obs["cell_type"] == ct
        masld_mask = ct_mask & adata.obs["condition"].isin(masld_labels)
        normal_mask = ct_mask & adata.obs["condition"].isin(normal_labels)

        n_masld = int(masld_mask.sum())
        n_normal = int(normal_mask.sum())

        if n_masld < 5 or n_normal < 5:
            log.info("  %s: skipping (MASLD=%d, Normal=%d)", ct, n_masld, n_normal)
            continue

        log.info("  %s: MASLD=%d, Normal=%d", ct, n_masld, n_normal)

        masld_cells = adata.obs_names[masld_mask].tolist()
        normal_cells = adata.obs_names[normal_mask].tolist()

        try:
            result = snap.tl.diff_test(
                adata,
                cell_group1=masld_cells,
                cell_group2=normal_cells,
            )
            if hasattr(result, 'to_pandas'):
                df = result.to_pandas()
            else:
                df = result
            df["cell_type"] = ct
            all_results.append(df)

            # Count significant
            for col in ["adjusted p-value", "adjusted_p_value", "padj", "fdr"]:
                if col in df.columns:
                    n_sig = int((df[col] < 0.05).sum())
                    log.info("  %s: %d DA regions (padj<0.05)", ct, n_sig)
                    break
        except Exception as e:
            log.error("  %s: FAILED: %s", ct, str(e))

    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        out_path = os.path.join(args.output_dir, "scatac_da_results.csv")
        combined.to_csv(out_path, index=False)
        log.info("Saved %d DA results -> %s", len(combined), out_path)
    else:
        log.warning("No DA results generated")


if __name__ == "__main__":
    main()
