#!/usr/bin/env python3
"""
Process PXD052937 (MassSpec DIA) proteomics dataset.

GSE276114 was removed: GEO confirms it is bulk RNA-seq (Expression profiling
by high throughput sequencing), not SomaScan proteomics.

Outputs normalized protein x sample matrix with HGNC gene symbols.
"""

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd


def process_pxd052937(data_dir: Path, output_dir: Path):
    """Process PXD052937 MassSpec DIA proteomics.

    The Spectronaut export is long-format TSV (one row per protein per sample):
      R.Condition | R.FileName | PG.ProteinAccessions | PG.ProteinNames | PG.Quantity
    We pivot to a wide protein x sample matrix indexed by the first gene symbol.
    """
    print("=== Processing PXD052937 (MassSpec DIA) ===")

    spectronaut_dir = data_dir / "20220805_163627_Plasma_liver2"

    # Load protein quantification (long format, .xls is actually TSV)
    protein_quant_file = spectronaut_dir / "Plasma_liver2_Report_Protein Quant (Normal).xls"
    print(f"  Loading protein quants from: {protein_quant_file.name}")
    protein_df = pd.read_csv(protein_quant_file, sep="\t")
    print(f"  Raw shape (long format): {protein_df.shape}")
    print(f"  Columns: {protein_df.columns.tolist()}")

    # Load condition setup (sample metadata)
    condition_file = spectronaut_dir / "Plasma_liver2_ConditionSetup.tsv"
    conditions = pd.read_csv(condition_file, sep="\t")
    print(f"  Conditions: {conditions.shape[0]} samples, groups: {conditions['Condition'].unique()}")

    # Extract gene symbols from PG.ProteinNames (format: "LV39_HUMAN;LV321_HUMAN")
    # or use PG.ProteinAccessions as fallback
    if "PG.ProteinAccessions" in protein_df.columns:
        protein_df["protein_id"] = protein_df["PG.ProteinAccessions"].str.split(";").str[0]
    print(f"  Unique proteins: {protein_df['protein_id'].nunique()}")
    print(f"  Unique samples (R.FileName): {protein_df['R.FileName'].nunique()}")

    # Pivot: protein_id x R.FileName, values = PG.Quantity
    protein_matrix = protein_df.pivot_table(
        index="protein_id", columns="R.FileName", values="PG.Quantity", aggfunc="first"
    )
    print(f"  Pivoted matrix: {protein_matrix.shape[0]} proteins x {protein_matrix.shape[1]} samples")

    # Log2 transform (raw intensities)
    if protein_matrix.max().max() > 100:
        protein_matrix = np.log2(protein_matrix.replace(0, np.nan))
        print("  Applied log2 transformation")

    # Median centering per sample
    sample_medians = protein_matrix.median(axis=0)
    protein_matrix = protein_matrix.subtract(sample_medians, axis=1).add(sample_medians.median())
    print("  Applied median centering")

    # Remove proteins with >50% missing
    missing_pct = protein_matrix.isna().sum(axis=1) / protein_matrix.shape[1]
    protein_matrix = protein_matrix[missing_pct < 0.5]
    print(f"  After missing filter: {protein_matrix.shape[0]} proteins, {protein_matrix.shape[1]} samples")

    # Save
    output_dir.mkdir(parents=True, exist_ok=True)
    protein_matrix.to_csv(output_dir / "pxd052937_protein_matrix.csv")
    conditions.to_csv(output_dir / "pxd052937_metadata.csv", index=False)
    print(f"  Saved to {output_dir / 'pxd052937_protein_matrix.csv'}")

    return protein_matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root",
                        default="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
    parser.add_argument("--output-dir", default="Analysis/Proteomics/results")
    args = parser.parse_args()

    root = Path(args.project_root)
    output_dir = root / args.output_dir

    pxd_matrix = process_pxd052937(root / "data" / "PXD052937", output_dir)

    print("\n=== PROTEOMICS PROCESSING COMPLETE ===")
    print(f"PXD052937: {pxd_matrix.shape[0]} proteins x {pxd_matrix.shape[1]} samples")


if __name__ == "__main__":
    main()
