
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
BIOTYPE_PATH = DATA_DIR / "ensembl_gene_biotypes.tsv.gz"

INHOUSE_MCD_FILES = {
    "MCD Week pooled (combined)": "mcd_week_pooled_combined.tsv.gz",
}

EXTERNAL_MCD_FILES = {
    "GSE156918 (external MCD)": "other_mcd_gse156918.tsv.gz",
    "GSE205974 (external MCD)": "other_mcd_gse205974.tsv.gz",
}

def load_df(path):
    sep = "\t" if path.name.endswith(".tsv.gz") else ","
    df = pd.read_csv(path, sep=sep)
    # Normalize columns
    cols = {c.lower(): c for c in df.columns}
    rename = {}
    if "gene" in cols: rename[cols["gene"]] = "gene_id"
    elif "gene_id" in cols: rename[cols["gene_id"]] = "gene_id"
    
    if "log2foldchange" in cols: rename[cols["log2foldchange"]] = "log2FoldChange"
    if "padj" in cols: rename[cols["padj"]] = "padj"
    
    if "tpm_mean" in cols: rename[cols["tpm_mean"]] = "tpm_mean"
    elif "tpm" in cols: rename[cols["tpm"]] = "tpm_mean"
    
    df = df.rename(columns=rename)
    return df

def strip_version(gene_id: str) -> str:
    return str(gene_id).split(".")[0]

def load_biotype_lncrna_ids(path: Path) -> set:
    print(f"Loading biotypes from {path}...")
    df = pd.read_csv(path, sep="\t")
    # Filter for Mouse lncRNAs
    # Note: species column is likely "mouse" or "human"
    mask = (df["species"] == "mouse") & (df["gene_biotype"] == "lncRNA")
    lncrnas = set(df[mask]["ensembl_gene_id"].astype(str).map(strip_version))
    print(f"Found {len(lncrnas)} mouse lncRNA IDs.")
    return lncrnas

PADJ_CUTOFF = 0.1

def run_analysis(datasets, valid_lncrnas, log2fc_cutoff, tpm_range):
    print(f"\n--- Running Analysis for log2FC > {log2fc_cutoff} ---")
    results = []

    for tpm_cut in tpm_range:
        # Identify Upregulated Sets filtered by lncRNA
        current_sets = {}
        for label, info in datasets.items():
            df = info["df"]
            # Filter DEGs
            mask = (df["padj"] < PADJ_CUTOFF) & (df["log2FoldChange"] > log2fc_cutoff)
            if "tpm_mean" in df.columns:
                mask &= (df["tpm_mean"] >= tpm_cut)
            
            # Get IDs and strip version
            genes = df.loc[mask, "gene_id"].astype(str).map(strip_version)
            
            # Filter for lncRNA membership
            lncrna_genes = set([g for g in genes if g in valid_lncrnas])
            
            current_sets[label] = lncrna_genes
            
        # Union of Mouse lncRNAs (Integrative pool)
        dedup_union = set()
        step_dataset_counts = {}
        
        for label, genes in current_sets.items():
            step_dataset_counts[label] = len(genes)
            dedup_union.update(genes)
            
        union_size = len(dedup_union)
        
        # Record stats
        row = {
            "tpm_cutoff": tpm_cut,
            "union_size": union_size,
        }
        row.update(step_dataset_counts)
        results.append(row)

    res_df = pd.DataFrame(results)
    
    # Colors
    COLOR_MAP = {
        "Cas13 MCD mouse": "#FFCD69",
        "Paquette MCD mouse": "#2A8C7D",
        "Yue MCD mouse": "#ED5565"
    }

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
    
    # Plotting
    print("Generating plots...")
    
    # 1. Main Line Plot
    plt.figure(figsize=(8, 6))
    
    # Plot Union
    plt.plot(res_df["tpm_cutoff"], res_df["union_size"], label="Total Unique lncRNAs (Union)", 
             linewidth=4, color="#333333", marker="o", zorder=10)
    
    # Plot individual datasets
    for col in res_df.columns:
        if col not in ["tpm_cutoff", "union_size"]:
            color = COLOR_MAP.get(col, "#888888")
            plt.plot(res_df["tpm_cutoff"], res_df[col], label=col, 
                     alpha=0.8, linewidth=2, linestyle="-", color=color)
            
    plt.xlabel("Global TPM Cutoff", fontweight='bold', fontsize=16)
    plt.ylabel("Number of lncRNAs", fontweight='bold', fontsize=16)
    plt.title(f"Sensitivity of lncRNA Count to TPM Cutoff\n(padj < {PADJ_CUTOFF}, log2FC > {log2fc_cutoff})", fontsize=18)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
    plt.tight_layout()
    
    fname_suffix = f"lfc{log2fc_cutoff}"
    plot_dir = Path(__file__).parent.parent / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    
    plot_path = plot_dir / f"lncrna_sensitivity_plot_{fname_suffix}.png"
    plt.savefig(plot_path, dpi=300)
    print(f"Saved {plot_path}")
    plt.close()

    # 2. Violin Plot (TPM Distribution of lncRNAs)
    plt.figure(figsize=(8, 6))
    all_vp_data = []
    
    for label, info in datasets.items():
        df = info["df"]
        mask = (df["padj"] < PADJ_CUTOFF) & (df["log2FoldChange"] > log2fc_cutoff)
        sig_df = df[mask].copy()
        
        # Filter for lncRNA membership
        sig_df["clean_id"] = sig_df["gene_id"].astype(str).map(strip_version)
        lncrna_df = sig_df[sig_df["clean_id"].isin(valid_lncrnas)]
        
        if "tpm_mean" in lncrna_df.columns:
            vals = lncrna_df["tpm_mean"].dropna()
            vals = vals[vals > 0]
            if len(vals) > 5000: vals = vals.sample(5000)
            
            temp = pd.DataFrame({
                "TPM": vals, 
                "Dataset": label, 
                "Species": "Mouse" # All are mouse
            })
            all_vp_data.append(temp)
            
    vp_df = pd.concat(all_vp_data)
    vp_df["log2(TPM+1)"] = np.log2(vp_df["TPM"] + 1)
    
    order = [
        "Cas13 MCD mouse", 
        "Paquette MCD mouse", 
        "Yue MCD mouse"
    ]
    order = [l for l in order if l in vp_df["Dataset"].unique()]
    
    sns.violinplot(
        data=vp_df, 
        x="Dataset", 
        y="log2(TPM+1)", 
        hue="Dataset",
        palette=COLOR_MAP,
        order=order,
        inner="quartile", 
        cut=0, 
        linewidth=1,
        density_norm="width",
        dodge=False,
        legend=False
    )
    plt.title(f"TPM Distribution (log2FC > {log2fc_cutoff})", fontsize=18, fontweight='bold')
    plt.xlabel("")
    plt.ylabel("log2(TPM + 1)", fontweight='bold', fontsize=16)
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    
    vp_path = plot_dir / f"lncrna_distribution_violin_{fname_suffix}.png"
    plt.savefig(vp_path, dpi=300)
    print(f"Saved {vp_path}")
    plt.close()

def main():
    print("Loading data...")
    # Load Datasets
    datasets = {}
    
    # In-house MCD - Only Pooled
    if "MCD Week pooled (combined)" in INHOUSE_MCD_FILES:
        fname = INHOUSE_MCD_FILES["MCD Week pooled (combined)"]
        datasets["Cas13 MCD mouse"] = {"df": load_df(DATA_DIR / fname), "species": "mouse"}
        
    # External MCD
    datasets["Paquette MCD mouse"] = {"df": load_df(DATA_DIR / EXTERNAL_MCD_FILES["GSE156918 (external MCD)"]), "species": "mouse"}
    datasets["Yue MCD mouse"] = {"df": load_df(DATA_DIR / EXTERNAL_MCD_FILES["GSE205974 (external MCD)"]), "species": "mouse"}

    print(f"Loaded {len(datasets)} datasets.")

    # Load lncRNA allowlist
    valid_lncrnas = load_biotype_lncrna_ids(BIOTYPE_PATH)
    
    # Analysis Configuration
    log2fc_cutoffs = [0.5, 0.8, 1.0, 2.0]
    tpm_range = np.arange(0, 10.25, 0.25) # 0 to 10 step 0.25
    
    for lfc in log2fc_cutoffs:
        run_analysis(datasets, valid_lncrnas, lfc, tpm_range)

if __name__ == "__main__":
    main()
