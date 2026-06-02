#!/usr/bin/env python3
"""
gene_overlap_share.py - Visualization of Unique vs Shared Gene Contributions

Generates a stacked bar plot showing how many genes in each dataset are
unique to that dataset vs shared with at least one other dataset.

Includes:
1. 7 Transcriptomic Datasets (MCD, Govaere, Hoang)
2. GWAS Genes (from Closest_genes.csv)

Filters:
- Transcriptomics: padj < 0.1, log2FC > 0.8, TPM > 1.0
- GWAS: No filter (All genes mapped to Ensembl)
"""

import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Definitions
ROOT = Path(__file__).resolve().parents[3] # Go up 3 levels: scripts/ -> cutoff_analysis/ -> RNA-seq/ -> Cas13/
DATA_DIR = ROOT / "streamlit_deg_explorer" / "data"

INHOUSE_MCD_FILES = {
    "MCD Week pooled (combined)": "mcd_week_pooled_combined.tsv.gz",
}

EXTERNAL_MCD_FILES = {
    "GSE156918 (external MCD)": "other_mcd_gse156918.tsv.gz",
    "GSE205974 (external MCD)": "other_mcd_gse205974.tsv.gz",
}

PATIENT_FILES = {
    "GSE130970": {
        "nas_1plus": ROOT / "RNA-seq/patient_RNAseq/results/GSE130970/deseq2_strict/differential_expression.csv",
        "fibrosis": "gse130970_fibrosis.csv.gz",
    },
    "GSE135251": {
        "nas_1plus": ROOT / "RNA-seq/patient_RNAseq/results/GSE135251/deseq2_strict/differential_expression.csv",
        "fibrosis": "gse135251_fibrosis.csv.gz",
    },
}

GWAS_FILE = "Closest_genes.csv"
ORTHOLOG_FILENAME = "mouse_human_orthologs.tsv.gz"
ORTHOLOG_PATH = DATA_DIR / ORTHOLOG_FILENAME

# TPM values sourced from bundled data files.
# NOTE: DEG results use NAS1+ vs strict control from nas_threshold_sensitivity/cumulative_nas/
TPM_FILES = {
    "GSE130970": DATA_DIR / "gse130970_nas_high.csv.gz",
    "GSE135251": DATA_DIR / "gse135251_nas_high.csv.gz"
}

def strip_version(gene_id: str) -> str:
    return str(gene_id).split(".")[0]

def load_df(path, tpm_path=None):
    sep = "\t" if str(path).endswith(".tsv.gz") else ","
    df = pd.read_csv(path, sep=sep)
    cols = {c.lower(): c for c in df.columns}
    rename = {}
    if "gene" in cols: rename[cols["gene"]] = "gene_id"
    elif "gene_id" in cols: rename[cols["gene_id"]] = "gene_id"
    
    if "gene_symbol" in cols: rename[cols["gene_symbol"]] = "gene_symbol"
    elif "symbol" in cols: rename[cols["symbol"]] = "gene_symbol"
    
    if "log2foldchange" in cols: rename[cols["log2foldchange"]] = "log2FoldChange"
    if "padj" in cols: rename[cols["padj"]] = "padj"
    
    if "tpm_mean" in cols: rename[cols["tpm_mean"]] = "tpm_mean"
    elif "tpm" in cols: rename[cols["tpm"]] = "tpm_mean"
    
    df = df.rename(columns=rename)
    
    # Merge TPM if missing
    if "tpm_mean" not in df.columns and tpm_path and os.path.exists(tpm_path):
        print(f"  Loading TPMs from {tpm_path}")
        tpm_sep = "\t" if str(tpm_path).endswith(".tsv.gz") else ","
        tpm_df = pd.read_csv(tpm_path, sep=tpm_sep)
        t_cols = {c.lower(): c for c in tpm_df.columns}
        t_gene = t_cols.get("gene_id") or t_cols.get("gene")
        t_val = t_cols.get("tpm_mean") or t_cols.get("tpm")
        
        if t_gene and t_val:
            # Map based on stripped version
            tmap = dict(zip(tpm_df[t_gene].astype(str).map(strip_version), tpm_df[t_val]))
            df["tpm_mean"] = df["gene_id"].astype(str).map(strip_version).map(tmap)
            
    return df

def strip_version(gene_id: str) -> str:
    return str(gene_id).split(".")[0]

def load_ortholog_map(path: Path):
    df = pd.read_csv(path, sep="\t")
    cols = {c.lower(): c for c in df.columns}
    mouse_col = cols.get("mouse_ensembl_gene_id") or cols.get("ensembl_gene_id") or cols.get("mouse_gene_id")
    human_col = cols.get("human_ensembl_gene_id") or cols.get("hsapiens_homolog_ensembl_gene")
    ortho_col = cols.get("orthology_type") or cols.get("hsapiens_homolog_orthology_type")
    
    df = df.rename(columns={mouse_col: "mouse_ensembl_gene_id", human_col: "human_ensembl_gene_id"})
    
    # Filter One2One
    if ortho_col:
        # print("  Filtering for one-to-one orthologs explicitly...")
        # df = df[df[ortho_col] == "ortholog_one2one"]
        pass

    df["mouse_ensembl_gene_id"] = df["mouse_ensembl_gene_id"].map(strip_version)
    df["human_ensembl_gene_id"] = df["human_ensembl_gene_id"].map(strip_version)
    
    # Map Mouse -> Human and Human -> Mouse
    m2h = {}
    h2m = {}
    for row in df.itertuples():
        m2h.setdefault(row.mouse_ensembl_gene_id, set()).add(row.human_ensembl_gene_id)
        h2m.setdefault(row.human_ensembl_gene_id, set()).add(row.mouse_ensembl_gene_id)
    return m2h, h2m

def main():
    print("Loading datasets...")
    datasets = {}
    
    # Load Transcriptomics
    datasets["Cas13 MCD mouse"] = {
        "df": load_df(DATA_DIR / INHOUSE_MCD_FILES["MCD Week pooled (combined)"]), 
        "species": "mouse"
    }
    datasets["Paquette MCD mouse"] = {
        "df": load_df(DATA_DIR / EXTERNAL_MCD_FILES["GSE156918 (external MCD)"]), 
        "species": "mouse"
    }
    datasets["Yue MCD mouse"] = {
        "df": load_df(DATA_DIR / EXTERNAL_MCD_FILES["GSE205974 (external MCD)"]), 
        "species": "mouse"
    }
    # Patients
    datasets["Hoang NAS1+"] = {
        "df": load_df(DATA_DIR / PATIENT_FILES["GSE130970"]["nas_1plus"], TPM_FILES["GSE130970"]), 
        "species": "human"
    }
    datasets["Hoang Fibrosis"] = {
        "df": load_df(DATA_DIR / PATIENT_FILES["GSE130970"]["fibrosis"], TPM_FILES["GSE130970"]), 
        "species": "human"
    }
    datasets["Govaere NAS1+"] = {
        "df": load_df(DATA_DIR / PATIENT_FILES["GSE135251"]["nas_1plus"], TPM_FILES["GSE135251"]), 
        "species": "human"
    }
    datasets["Govaere Fibrosis"] = {
        "df": load_df(DATA_DIR / PATIENT_FILES["GSE135251"]["fibrosis"], TPM_FILES["GSE135251"]), 
        "species": "human"
    }
    
    # Load GWAS
    gwas_df = pd.read_csv(DATA_DIR / GWAS_FILE, header=None, names=["gene_symbol"])
    datasets["GWAS Genes"] = {"df": gwas_df, "species": "human", "is_gwas": True}
    
    # Build Symbol -> Ensembl Map from Patient Data
    print("Building Symbol -> Ensembl Map...")
    symbol_to_id = {}
    for label, info in datasets.items():
        if info["species"] == "human" and not info.get("is_gwas"):
            df = info["df"]
            if "gene_symbol" in df.columns and "gene_id" in df.columns:
                temp_map = pd.Series(
                    df["gene_id"].astype(str).map(strip_version).values,
                    index=df["gene_symbol"]
                ).to_dict()
                symbol_to_id.update(temp_map)
    
    print(f"  Mapped {len(symbol_to_id)} symbols to Ensembl IDs.")
    
    # Load Orthologs
    m2h, h2m = load_ortholog_map(ORTHOLOG_PATH)
    print(f"Ortholog Map One-to-One Sizes: m2h={len(m2h)}, h2m={len(h2m)}")
    
    # Process Gene Sets (Apply Filters and Map to Human)
    padj_cut = 0.1
    lfc_cut = 0.8
    tpm_cut = 1.0 # Added filter
    
    final_sets = {}
    
    for label, info in datasets.items():
        df = info["df"]
        
        # Filter
        genes = set()
        if info.get("is_gwas"):
            # Map GWAS Symbols to IDs
            raw_symbols = set(df["gene_symbol"].dropna().astype(str))
            for sym in raw_symbols:
                if sym in symbol_to_id:
                    genes.add(symbol_to_id[sym])
                else:
                    # Try uppercase or simple match just in case
                    pass
            print(f"  GWAS: {len(genes)} IDs mapped from {len(raw_symbols)} symbols")
        else:
            # Transcriptomics Filtering
            mask = (df["padj"] < padj_cut) & (df["log2FoldChange"] > lfc_cut)
            
            # Apply TPM filter if available
            if "tpm_mean" in df.columns:
                mask &= (df["tpm_mean"] >= tpm_cut)
                
            genes = set(df.loc[mask, "gene_id"].astype(str).map(strip_version))
            
        # Map to Human if Mouse, or Filter Human if unchecked
        mapped = set()
        if info["species"] == "mouse":
            for g in genes:
                targets = m2h.get(g)
                if targets:
                    mapped.update(targets)
        else:
            # Human: Filter to only those with mouse orthologs (Strict Match)
            for g in genes:
                 if g in h2m:
                     mapped.add(g)
            
        final_sets[label] = mapped
        print(f"  {label}: {len(mapped)} genes (mapped/filtered)")
        
    # Calculate Unique vs Shared
    results = []
    
    all_genes_union = set()
    degs_only_union = set()
    
    for label, s in final_sets.items():
        all_genes_union.update(s)
        if label != "GWAS Genes":
            degs_only_union.update(s)
            
    print(f"\nTotal Dedup Genes (including GWAS): {len(all_genes_union)}")
    print(f"Total Dedup DEGs (excluding GWAS): {len(degs_only_union)}")
        
    for label, genes in final_sets.items():
        other_genes = set()
        for other_label, other_set in final_sets.items():
            if other_label != label:
                other_genes.update(other_set)
        
        shared = genes.intersection(other_genes)
        unique = genes - other_genes
        
        results.append({
            "Dataset": label,
            "unique_only": len(unique),
            "shared_any": len(shared),
            "total": len(genes)
        })
        
    df_res = pd.DataFrame(results).sort_values("total", ascending=True)
    print("\nResults:")
    print(df_res)
    
    # Plotting
    print("Generating plot...")
    
    # Colors: Magenta (Unique), Orange (Shared) - matching user screenshot style
    COLOR_UNIQUE = "#c2185b" # Dark Magenta
    COLOR_SHARED = "#ffb74d" # Light Orange
    # Set context and style
    sns.set_context("paper") 
    sns.set_style("whitegrid", {'axes.grid': False})

    # CRITICAL: Configure fonts for Illustrator (Type 42 = TrueType)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    # Font Family (Helvetica/Arial)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'sans-serif']

    # Strict Font Sizes (Journal Standards)
    plt.rcParams['font.size'] = 14           # Default text size
    plt.rcParams['axes.titlesize'] = 16      # Title size
    plt.rcParams['axes.labelsize'] = 14      # X and Y label size
    plt.rcParams['xtick.labelsize'] = 12     # X tick size
    plt.rcParams['ytick.labelsize'] = 12     # Y tick size
    plt.rcParams['legend.fontsize'] = 12     # Legend size
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    y_pos = np.arange(len(df_res))
    
    p1 = ax.barh(y_pos, df_res["unique_only"], color=COLOR_UNIQUE, label="Unique Only", edgecolor='white')
    p2 = ax.barh(y_pos, df_res["shared_any"], left=df_res["unique_only"], color=COLOR_SHARED, label="Shared (Any)", edgecolor='white')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df_res["Dataset"])
    ax.set_xlabel("Number of Genes", fontsize=16, fontweight='bold')
    ax.set_title(f"Unique vs Shared Gene Contributions\n(Shared = overlaps with ANY other dataset | TPM > {tpm_cut})", fontweight='bold', fontsize=18)
    
    ax.legend(loc='lower right')
    
    # Annotate with totals - REMOVED per user request
    # for i, (total, unique, shared) in enumerate(zip(df_res["total"], df_res["unique_only"], df_res["shared_any"])):
    #     ax.text(total + max(df_res["total"])*0.01, i, str(total), va='center', fontweight='bold', color='#333333')
        
    plt.tight_layout()
    
    plot_dir = Path(__file__).parent.parent / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    
    plot_path = plot_dir / "gene_overlap_share.png"
    plt.savefig(plot_path, dpi=300)
    print(f"Saved {plot_path}")

if __name__ == "__main__":
    main()
