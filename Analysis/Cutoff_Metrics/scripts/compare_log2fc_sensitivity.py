
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
ORTHOLOG_PATH = DATA_DIR / "mouse_human_orthologs.tsv.gz"

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
        # "nas_low": "gse130970_nas_low.csv.gz",
        "fibrosis": "gse130970_fibrosis.csv.gz",
    },
    "GSE135251": {
        "nas_1plus": ROOT / "RNA-seq/patient_RNAseq/results/GSE135251/deseq2_strict/differential_expression.csv",
        # "nas_low": "gse135251_nas_low.csv.gz",
        "fibrosis": "gse135251_fibrosis.csv.gz",
    },
}

# TPM values sourced from bundled data files.
# NOTE: DEG results use NAS1+ vs strict control from nas_threshold_sensitivity/cumulative_nas/
TPM_FILES = {
    "GSE130970": DATA_DIR / "gse130970_nas_high.csv.gz",
    "GSE135251": DATA_DIR / "gse135251_nas_high.csv.gz"
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
# Strict Font Sizes (Journal Standards)
plt.rcParams['font.size'] = 14           # Default text size
plt.rcParams['axes.titlesize'] = 16      # Title size
plt.rcParams['axes.labelsize'] = 14      # X and Y label size
plt.rcParams['xtick.labelsize'] = 12     # X tick size
plt.rcParams['ytick.labelsize'] = 12     # Y tick size
plt.rcParams['legend.fontsize'] = 12     # Legend size

def strip_version(gene_id: str) -> str:
    return str(gene_id).split(".")[0]

def load_df(path, tpm_path=None):
    sep = "\t" if str(path).endswith(".tsv.gz") else ","
    df = pd.read_csv(path, sep=sep)
    cols = {c.lower(): c for c in df.columns}
    rename = {}
    if "gene" in cols: rename[cols["gene"]] = "gene_id"
    elif "gene_id" in cols: rename[cols["gene_id"]] = "gene_id"
    if "log2foldchange" in cols: rename[cols["log2foldchange"]] = "log2FoldChange"
    if "padj" in cols: rename[cols["padj"]] = "padj"
    if "tpm_mean" in cols: rename[cols["tpm_mean"]] = "tpm_mean"
    elif "tpm" in cols: rename[cols["tpm"]] = "tpm_mean"
    
    df = df.rename(columns=rename)
    
    # Merge TPM if missing and source provided
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

def load_ortholog_map(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    cols = {c.lower(): c for c in df.columns}
    mouse_col = cols.get("mouse_ensembl_gene_id") or cols.get("ensembl_gene_id") or cols.get("mouse_gene_id")
    human_col = cols.get("human_ensembl_gene_id") or cols.get("hsapiens_homolog_ensembl_gene")
    df = df.rename(columns={mouse_col: "mouse_ensembl_gene_id", human_col: "human_ensembl_gene_id"})
    df["mouse_ensembl_gene_id"] = df["mouse_ensembl_gene_id"].map(strip_version)
    df["human_ensembl_gene_id"] = df["human_ensembl_gene_id"].map(strip_version)
    return df

def build_maps(df):
    m2h = {}
    h2m = {}
    for row in df.itertuples():
        m2h.setdefault(row.mouse_ensembl_gene_id, set()).add(row.human_ensembl_gene_id)
        h2m.setdefault(row.human_ensembl_gene_id, set()).add(row.mouse_ensembl_gene_id)
    return m2h, h2m

def load_biotype_lncrna_ids(path: Path) -> set:
    df = pd.read_csv(path, sep="\t")
    mask = (df["species"] == "mouse") & (df["gene_biotype"] == "lncRNA")
    return set(df[mask]["ensembl_gene_id"].astype(str).map(strip_version))

PADJ_CUTOFF = 0.1

# Custom Palette (Cyan -> Pink -> Magenta -> Purple)
# 0.5: Cyan (#43afb4)
# 0.8: Light Pink (#f06ca9)
# 1.0: Deep Magenta (#b83674)
# 2.0: Dark Purple (#703d74)
CUSTOM_PALETTE = {
    0.5: "#43afb4",
    0.8: "#f06ca9",
    1.0: "#b83674",
    2.0: "#703d74"
}

def run_pcg_comparison(datasets, m2h, h2m, start, end, step):
    cutoffs = [0.5, 0.8, 1.0, 2.0]
    tpm_range = np.arange(start, end + step, step)
    
    comparisons = []
    
    print("\n[PCG] Calculating unions...")
    for lfc in cutoffs:
        print(f"  Processing log2FC > {lfc}")
        for tpm in tpm_range:
            current_sets = {}
            for label, info in datasets.items():
                df = info["df"]
                mask = (df["padj"] < PADJ_CUTOFF) & (df["log2FoldChange"] > lfc)
                if "tpm_mean" in df.columns:
                    mask &= (df["tpm_mean"] >= tpm)
                
                genes = set(df.loc[mask, "gene_id"].astype(str).map(strip_version))
                
                # Map to Human (Strict Orthology)
                mapped = set()
                if info["species"] == "human":
                     for g in genes:
                         if g in h2m: mapped.add(g)
                elif info["species"] == "mouse":
                    for g in genes:
                        if g in m2h: mapped.update(m2h[g])
                
                current_sets[label] = mapped

            # Union
            dedup_union = set()
            for s in current_sets.values():
                dedup_union.update(s)
            
            comparisons.append({
                "cutoff_val": lfc,
                "tpm": tpm,
                "count": len(dedup_union)
            })

    df_res = pd.DataFrame(comparisons)
    
    # Plotting
    # Plotting
    plt.figure(figsize=(8, 6))
    sns.lineplot(data=df_res, x="tpm", y="count", hue="cutoff_val", 
                 palette=CUSTOM_PALETTE, linewidth=2.5, marker="o")
    
    plt.xlabel("Global TPM Cutoff", fontweight='bold', fontsize=16)
    plt.ylabel("Total Deduplicated DEGs (Union)", fontweight='bold', fontsize=16)
    plt.title(f"PCG Sensitivity: Log2FC Comparison\n(padj < {PADJ_CUTOFF}, Strict Orthology)", fontsize=18)
    plt.legend(title="Log2FC Cutoff")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    plot_dir = Path(__file__).parent.parent / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    
    out_path = plot_dir / "pcg_log2fc_comparison.png"
    plt.savefig(out_path, dpi=300)
    print(f"Saved {out_path}")


def run_lncrna_comparison(datasets, valid_lncrnas, start, end, step):
    # Only use mouse datasets for lncRNA
    mouse_datasets = {k: v for k, v in datasets.items() if v["species"] == "mouse"}
    
    cutoffs = [0.5, 0.8, 1.0, 2.0]
    tpm_range = np.arange(start, end + step, step)
    
    comparisons = []
    
    print("\n[lncRNA] Calculating unions...")
    for lfc in cutoffs:
        print(f"  Processing log2FC > {lfc}")
        for tpm in tpm_range:
            dedup_union = set()
            
            for label, info in mouse_datasets.items():
                df = info["df"]
                mask = (df["padj"] < PADJ_CUTOFF) & (df["log2FoldChange"] > lfc)
                if "tpm_mean" in df.columns:
                    mask &= (df["tpm_mean"] >= tpm)
                
                genes = df.loc[mask, "gene_id"].astype(str).map(strip_version)
                # Filter lncRNA
                lnc_genes = [g for g in genes if g in valid_lncrnas]
                dedup_union.update(lnc_genes)
            
            comparisons.append({
                "cutoff_val": lfc,
                "tpm": tpm,
                "count": len(dedup_union)
            })

    df_res = pd.DataFrame(comparisons)
    
    # Plotting
    # Plotting
    plt.figure(figsize=(8, 6))
    sns.lineplot(data=df_res, x="tpm", y="count", hue="cutoff_val", 
                 palette=CUSTOM_PALETTE, linewidth=2.5, marker="o")
    
    plt.xlabel("Global TPM Cutoff", fontweight='bold', fontsize=16)
    plt.ylabel("Total Unique lncRNAs (Union)", fontweight='bold', fontsize=16)
    plt.title(f"lncRNA Sensitivity: Log2FC Comparison\n(padj < {PADJ_CUTOFF}, Mouse Datasets Only)", fontsize=18)
    plt.legend(title="Log2FC Cutoff")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    plot_dir = Path(__file__).parent.parent / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    
    out_path = plot_dir / "lncrna_log2fc_comparison.png"
    plt.savefig(out_path, dpi=300)
    print(f"Saved {out_path}")

def main():
    print("Loading data...")
    datasets = {}
    
    # Load all (some unused for lncRNA but cheap to load)
    if "MCD Week pooled (combined)" in INHOUSE_MCD_FILES:
        fname = INHOUSE_MCD_FILES["MCD Week pooled (combined)"]
        datasets["Cas13 MCD mouse"] = {"df": load_df(DATA_DIR / fname), "species": "mouse"}
        
    for label, fname in EXTERNAL_MCD_FILES.items():
        datasets[f"MCD (external) | {label}"] = {"df": load_df(DATA_DIR / fname), "species": "mouse"}
        
    for ds_name, files in PATIENT_FILES.items():
        tpath = TPM_FILES.get(ds_name)
        if "nas_1plus" in files:
            datasets[f"Patient | {ds_name} | nas_1plus"] = {"df": load_df(DATA_DIR / files["nas_1plus"], tpath), "species": "human"}
        if "fibrosis" in files:
            datasets[f"Patient | {ds_name} | fibrosis"] = {"df": load_df(DATA_DIR / files["fibrosis"], tpath), "species": "human"}

    print(f"Loaded {len(datasets)} datasets.")

    # PCG Requirements
    ortho_df = load_ortholog_map(ORTHOLOG_PATH)
    m2h, h2m = build_maps(ortho_df)
    
    # lncRNA Requirements
    valid_lncrnas = load_biotype_lncrna_ids(BIOTYPE_PATH)
    
    # Run Comparisons
    # PCG: Step 0.5
    run_pcg_comparison(datasets, m2h, h2m, 0, 10, 0.5)
    
    # lncRNA: Step 0.25
    run_lncrna_comparison(datasets, valid_lncrnas, 0, 10, 0.25)

if __name__ == "__main__":
    main()
