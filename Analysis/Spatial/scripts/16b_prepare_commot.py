#!/usr/bin/env python3
"""
16b_prepare_commot.py — Prepare COMMOT ligand-receptor database and annotations.

Loads the CellChat human LR database via COMMOT, filters by expression
frequency, and annotates MASLD-relevant pathway families.

SLURM: --partition=cpu --cpus=8 --mem=64G --time=2:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, load_deconvolved_adata,
    save_csv, print_header, print_step,
)


def load_cellchat_database(signaling_type="Secreted Signaling"):
    """Load CellChat human LR database via COMMOT."""
    import commot as ct

    print(f"  Loading CellChat database (type={signaling_type})...")
    df_lr = ct.pp.ligand_receptor_database(
        species="human", signaling_type=signaling_type,
    )
    print(f"  Raw LR pairs: {len(df_lr)}")
    return df_lr


def compute_expression_frequency(adata, genes):
    """Compute fraction of spots expressing each gene (>0 counts)."""
    gene_mask = adata.var_names.isin(genes)
    if hasattr(adata.X, "toarray"):
        expr_matrix = adata.X[:, gene_mask].toarray()
    else:
        expr_matrix = np.asarray(adata.X[:, gene_mask])

    freq = (expr_matrix > 0).mean(axis=0)
    freq_dict = dict(zip(adata.var_names[gene_mask], freq))
    return freq_dict


def _get_lr_columns(df_lr):
    """Identify ligand, receptor, pathway columns in CellChat DB.

    CellChat DB from COMMOT has numeric column names: 0=ligand, 1=receptor,
    2=pathway, 3=signaling_type. Some versions may use string names.
    """
    # Try named columns first
    lig_col = next((c for c in df_lr.columns if isinstance(c, str) and "ligand" in c.lower()), None)
    rec_col = next((c for c in df_lr.columns if isinstance(c, str) and "receptor" in c.lower()), None)
    pw_col = next((c for c in df_lr.columns if isinstance(c, str) and "pathway" in c.lower()), None)

    # Fall back to positional (CellChat default: 0=lig, 1=rec, 2=pathway, 3=type)
    if lig_col is None:
        cols = list(df_lr.columns)
        lig_col = cols[0]
        rec_col = cols[1]
        pw_col = cols[2] if len(cols) > 2 else None

    return lig_col, rec_col, pw_col


def filter_lr_by_expression(df_lr, freq_dict, min_freq=0.05):
    """Filter LR pairs requiring both ligand and receptor expressed in >=min_freq spots."""
    lig_col, rec_col, _ = _get_lr_columns(df_lr)

    keep = []
    for idx, row in df_lr.iterrows():
        ligand = row[lig_col]
        receptor = row[rec_col]

        # Handle multi-subunit receptors (e.g., "TGFBR1_TGFBR2")
        lig_genes = str(ligand).split("_") if "_" in str(ligand) else [str(ligand)]
        rec_genes = str(receptor).split("_") if "_" in str(receptor) else [str(receptor)]

        lig_ok = all(freq_dict.get(g, 0) >= min_freq for g in lig_genes)
        rec_ok = all(freq_dict.get(g, 0) >= min_freq for g in rec_genes)

        if lig_ok and rec_ok:
            keep.append(idx)

    df_filtered = df_lr.loc[keep].copy()
    print(f"  After expression filter (>={min_freq*100:.0f}% spots): "
          f"{len(df_filtered)}/{len(df_lr)} LR pairs")
    return df_filtered


def annotate_masld_pathways(df_lr, key_pathways):
    """Classify each LR pair into MASLD-relevant pathway families."""
    lig_col, rec_col, pathway_col = _get_lr_columns(df_lr)

    if pathway_col is None:
        print("  WARNING: No pathway column found in LR database, "
              "annotating from gene names")

        pathway_assignments = []
        for _, row in df_lr.iterrows():
            pair_str = f"{row[lig_col]}_{row[rec_col]}".upper()
            assigned = "Other"
            for kp in key_pathways:
                if kp.upper() in pair_str:
                    assigned = kp
                    break
            pathway_assignments.append(assigned)
        df_lr = df_lr.copy()
        df_lr["masld_pathway"] = pathway_assignments
    else:
        # Match pathway column against key pathways
        pathway_assignments = []
        for _, row in df_lr.iterrows():
            pw = str(row[pathway_col]).upper()
            assigned = "Other"
            for kp in key_pathways:
                if kp.upper() in pw:
                    assigned = kp
                    break
            pathway_assignments.append(assigned)
        df_lr = df_lr.copy()
        df_lr["masld_pathway"] = pathway_assignments

    # Summary
    for kp in key_pathways:
        n = (df_lr["masld_pathway"] == kp).sum()
        if n > 0:
            print(f"    {kp}: {n} LR pairs")
    n_other = (df_lr["masld_pathway"] == "Other").sum()
    print(f"    Other: {n_other} LR pairs")

    return df_lr


def main():
    print_header("16b: Prepare COMMOT LR Database")

    config = load_config()
    commot_config = config["commot"]
    output_dir = RESULTS_DIR / "commot"
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Load CellChat database ──
    print_step("Loading CellChat LR database", 1, 4)
    signaling_type = commot_config.get("signaling_type", "Secreted")
    # COMMOT expects full signaling type string
    if signaling_type == "Secreted":
        signaling_type = "Secreted Signaling"
    df_lr = load_cellchat_database(signaling_type=signaling_type)

    # ── Step 2: Load spatial data and compute expression frequency ──
    print_step("Computing gene expression frequency", 2, 4)
    adata = load_deconvolved_adata()
    print(f"  Loaded: {adata.n_obs} spots × {adata.n_vars} genes")

    # Collect all genes mentioned in LR database
    lig_col, rec_col, _ = _get_lr_columns(df_lr)

    all_lr_genes = set()
    for _, row in df_lr.iterrows():
        for g in str(row[lig_col]).split("_"):
            all_lr_genes.add(g)
        for g in str(row[rec_col]).split("_"):
            all_lr_genes.add(g)

    present_genes = all_lr_genes & set(adata.var_names)
    print(f"  LR genes in database: {len(all_lr_genes)}, "
          f"present in adata: {len(present_genes)}")

    freq_dict = compute_expression_frequency(adata, present_genes)

    # ── Step 3: Filter by expression ──
    print_step("Filtering LR pairs by expression", 3, 4)
    min_freq = 0.05  # 5% of spots
    df_filtered = filter_lr_by_expression(df_lr, freq_dict, min_freq=min_freq)

    # ── Step 4: Annotate MASLD pathways ──
    print_step("Annotating MASLD-relevant pathways", 4, 4)
    key_pathways = commot_config.get("key_pathways", [])
    df_annotated = annotate_masld_pathways(df_filtered, key_pathways)

    # ── Save outputs ──
    save_csv(df_annotated, "lr_database_cellchat_filtered.csv", subdir="commot")

    # Save pathway annotations separately
    pathway_summary = df_annotated.groupby("masld_pathway").size().reset_index(name="n_pairs")
    pathway_summary = pathway_summary.sort_values("n_pairs", ascending=False)
    save_csv(pathway_summary, "masld_pathway_annotations.csv", subdir="commot")

    # ── Summary ──
    print(f"\n  Summary:")
    print(f"    Raw LR pairs: {len(df_lr)}")
    print(f"    After expression filter: {len(df_filtered)}")
    print(f"    Key MASLD pathways annotated: {len(key_pathways)}")
    print(f"    LR pairs in key pathways: "
          f"{(df_annotated['masld_pathway'] != 'Other').sum()}")

    print_header("16b: Complete")


if __name__ == "__main__":
    main()
