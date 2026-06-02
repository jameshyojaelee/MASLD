import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Definitions copied/adapted from app.py to avoid import side-effects
ROOT = Path(__file__).resolve().parents[3] # Go up 3 levels: scripts/ -> cutoff_analysis/ -> RNA-seq/ -> Cas13/
DATA_DIR = ROOT / "streamlit_deg_explorer" / "data"

INHOUSE_MCD_FILES = {
    "MCD Week 1": "mcd_week1.tsv.gz",
    "MCD Week 2": "mcd_week2.tsv.gz",
    "MCD Week 3": "mcd_week3.tsv.gz",
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
        # "fibrosis": "gse130970_fibrosis.csv.gz", # Assuming fibrosis specific analysis not run/needed in strict pass or same file? 
        # Actually user asked for NAS1+ vs Strict. Fibrosis specific comparison might be redundant or removed.
        # Let's verify if fibrosis specific comparison exists or if we should just use the main file.
        # For now, pointing nas_1plus to the main result.
        # If 'fibrosis' key implies a different comparison (e.g. F >= 1 vs F0), we haven't explicitely run that as 'fibrosis' named output. 
        # But the strict analysis captures 'Disease' vs 'Control'.
        # I will start by pointing nas_1plus. If fibrosis key is used, I should check if I need to run it or if it's legacy.
        # The script iterates over keys.
    },
    "GSE135251": {
        "nas_1plus": ROOT / "RNA-seq/patient_RNAseq/results/GSE135251/deseq2_strict/differential_expression.csv",
    },
}

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

def load_ortholog_map(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    cols = {c.lower(): c for c in df.columns}
    mouse_col = cols.get("mouse_ensembl_gene_id") or cols.get("ensembl_gene_id") or cols.get("mouse_gene_id")
    human_col = cols.get("human_ensembl_gene_id") or cols.get("hsapiens_homolog_ensembl_gene")
    ortho_col = cols.get("orthology_type") or cols.get("hsapiens_homolog_orthology_type")
    
    df = df.rename(columns={
        mouse_col: "mouse_ensembl_gene_id",
        human_col: "human_ensembl_gene_id",
        ortho_col: "orthology_type",
    })
    df["mouse_ensembl_gene_id"] = df["mouse_ensembl_gene_id"].map(strip_version)
    df["human_ensembl_gene_id"] = df["human_ensembl_gene_id"].map(strip_version)
    return df

def build_maps(df, one2one_only):
    if one2one_only:
        df = df[df["orthology_type"] == "ortholog_one2one"]
    m2h = {}
    h2m = {}
    for row in df.itertuples():
        m2h.setdefault(row.mouse_ensembl_gene_id, set()).add(row.human_ensembl_gene_id)
        h2m.setdefault(row.human_ensembl_gene_id, set()).add(row.mouse_ensembl_gene_id)
    return m2h, h2m

# Global Params
PADJ_CUTOFF = 0.1

def run_analysis(datasets, m2h, h2m, log2fc_cutoff, tpm_range):
    print(f"\n--- Running Analysis for log2FC > {log2fc_cutoff} ---")
    include_unmapped_human = False
    include_unmapped_mouse = False
    
    results = []

    for tpm_cut in tpm_range:
        # 1. Identify Upregulated Sets
        current_sets = {}
        for label, info in datasets.items():
            df = info["df"]
            # Filter
            mask = (df["padj"] < PADJ_CUTOFF) & (df["log2FoldChange"] > log2fc_cutoff)
            if "tpm_mean" in df.columns:
                mask &= (df["tpm_mean"] >= tpm_cut)
            
            genes = set(df.loc[mask, "gene_id"].astype(str).map(strip_version))
            current_sets[label] = genes
            
        # 2. Canonicalize/Deduplicate
        dedup_union = set()
        
        step_dataset_counts = {}
        
        for label, genes in current_sets.items():
            species = datasets[label]["species"]
            mapped_set = set()
            
            if species == "human":
                for g in genes:
                    if include_unmapped_human:
                        mapped_set.add(g)
                    elif g in h2m: # Has a mouse match
                        mapped_set.add(g)
                        
            elif species == "mouse":
                for g in genes:
                    targets = m2h.get(g)
                    if targets:
                        mapped_set.update(targets)
                    elif include_unmapped_mouse:
                        mapped_set.add(f"MOUSE:{g}")
            
            # current_sets[label] = mapped_set # Not needed unless debugging
            step_dataset_counts[label] = len(mapped_set)
            dedup_union.update(mapped_set)
            
        union_size = len(dedup_union)
        
        # Record stats
        row = {
            "tpm_cutoff": tpm_cut,
            "union_size": union_size,
        }
        row.update(step_dataset_counts)
        results.append(row)

    res_df = pd.DataFrame(results)
    
    # --- Plotting ---
    
    # Helper to color mapping (For Line Plot - Distinguish shades)
    # Helper to color mapping (Matched to plot_tpm_distributions_publication.R)
    # UPDATED: Magenta, Pink, Purple theme
    
    # Mouse (Purples)
    # Cas13 -> #9C27B0 (Purple)
    # Paquette -> #AB47BC (Lighter Purple)
    # Yue -> #7B1FA2 (Dark Purple)
    
    # Human (Pinks/Red-Violets)
    # Hoang NAS1+ -> #E91E63 (Pink)
    # Hoang Fibrosis -> #C2185B (Dark Pink)
    # Govaere NAS1+ -> #FF4081 (Light Pink)
    # Govaere Fibrosis -> #F50057 (Strong Pink)

    COLOR_MAP = {
        "Cas13 MCD mouse": "#9C27B0",
        "Paquette MCD mouse": "#AB47BC",
        "Yue MCD mouse": "#7B1FA2",
        "Hoang NAS1+": "#E91E63",
        "Hoang Fibrosis": "#C2185B",
        "Govaere NAS1+": "#FF4081",
        "Govaere Fibrosis": "#F50057"
    }

    def get_line_color(label):
        # Clean label if needed (remove 'Human', etc if passed full label map key)
        # The keys here match the values in LABEL_MAP below
        
        # If label is one of the keys in LABEL_MAP, map it first
        short_label = LABEL_MAP.get(label, label)
        return COLOR_MAP.get(short_label, "#888888")

    # Label Mapping
    LABEL_MAP = {
        "Patient | GSE135251 | nas_1plus": "Govaere NAS1+",
        "Patient | GSE130970 | nas_1plus": "Hoang NAS1+",
        "Patient | GSE135251 | fibrosis": "Govaere Fibrosis",
        "Patient | GSE130970 | fibrosis": "Hoang Fibrosis",
        "MCD (in-house) | MCD Week pooled (combined)": "Cas13 MCD mouse",
        "MCD (external) | GSE156918 (external MCD)": "Paquette MCD mouse",
        "MCD (external) | GSE205974 (external MCD)": "Yue MCD mouse"
    }
    
    # Define Palette for Species (used in violin plot)
    # Using Generic colors if specific dataset not passed, but violin uses dataset hue
    # Actually violin plot logic below uses SPECIES_PALETTE for 'hue="Species"'
    # But current violin plot code uses 'hue="Species"'? 
    # Let's check violin plot code. It sets hue="Dataset".
    
    # We will modify the violin plot call later in this file.
    # For now, let's keep SPECIES_PALETTE compatible just in case, but we intend to change the plot.
    SPECIES_PALETTE = {
        "Human": "#E91E63", # Pink
        "Mouse": "#9C27B0"  # Purple
    }
    
    def get_species(label):
        if "Govaere" in label or "Hoang" in label or "(Human)" in label: return "Human"
        if "mouse" in label.lower(): return "Mouse"
        return "Unknown"

    print("Generating plots...")
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
    plt.rcParams.update({'axes.grid': True, 'grid.alpha': 0.3})
    
    # 1. Main Line Plot
    # 1. Main Line Plot
    plt.figure(figsize=(8, 6))
    
    # Plot Union
    plt.plot(res_df["tpm_cutoff"], res_df["union_size"], label="Total Deduplicated DEGs (Union)", 
             linewidth=4, color="#333333", marker="o", zorder=10)
    
    # Plot individual datasets
    for col in res_df.columns:
        if col not in ["tpm_cutoff", "union_size"]:
            display_label = LABEL_MAP.get(col, col)
            color = get_line_color(col)
            plt.plot(res_df["tpm_cutoff"], res_df[col], label=display_label, 
                     alpha=0.8, linewidth=2, linestyle="-", color=color)
            
    plt.xlabel("Global TPM Cutoff", fontweight='bold', fontsize=16)
    plt.ylabel("Number of Genes", fontweight='bold', fontsize=16)
    plt.title(f"Sensitivity of DEG Count to TPM Cutoff\n(padj < {PADJ_CUTOFF}, log2FC > {log2fc_cutoff}, Strict Orthology)", fontsize=18)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
    plt.tight_layout()
    
    fname_suffix = f"lfc{log2fc_cutoff}"
    plot_dir = Path(__file__).parent.parent / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    plot_path = plot_dir / f"tpm_sensitivity_plot_{fname_suffix}.pdf"
    plt.savefig(plot_path, dpi=300, transparent=True)
    print(f"Saved {plot_path}")
    plt.close() # Close to free memory

    # 2. Violin Plot
    plt.figure(figsize=(10, 6))
    # Collect all TPM values
    all_vp_data = []
    for label, info in datasets.items():
        df = info["df"]
        mask = (df["padj"] < PADJ_CUTOFF) & (df["log2FoldChange"] > log2fc_cutoff)
        sig_df = df[mask]
        
        if "tpm_mean" in sig_df.columns:
            vals = sig_df["tpm_mean"].dropna()
            vals = vals[vals > 0]
            if len(vals) > 5000: vals = vals.sample(5000)
            
            display_label = LABEL_MAP.get(label, label)
            species = get_species(display_label)
            
            temp = pd.DataFrame({
                "TPM": vals, 
                "Dataset": display_label, 
                "Species": species
            })
            all_vp_data.append(temp)
            
    vp_df = pd.concat(all_vp_data)
    vp_df["log2(TPM+1)"] = np.log2(vp_df["TPM"] + 1)
    
    # Sort datasets
    order = [
        "Govaere NAS1+", "Govaere Fibrosis",
        "Hoang NAS1+", "Hoang Fibrosis",
        "Cas13 MCD mouse", "Paquette MCD mouse", "Yue MCD mouse"
    ]
    # Filter to only those present
    order = [l for l in order if l in vp_df["Dataset"].unique()]
    
    sns.violinplot(
        data=vp_df, 
        x="Dataset", 
        y="log2(TPM+1)", 
        hue="Dataset",
        palette=COLOR_MAP,
        order=order,
        inner="quartile", 
        linewidth=1,
        density_norm="width",
        dodge=False,
        legend=False
    )
    plt.title(f"TPM Distribution (padj < {PADJ_CUTOFF}, log2FC > {log2fc_cutoff})", fontsize=18, fontweight='bold')
    plt.xlabel("")
    plt.ylabel("log2(TPM + 1)", fontweight='bold', fontsize=16)
    plt.xticks(rotation=25, ha="right")
    # plt.legend(title="Species", frameon=False) # Redundant with x-axis labels
    plt.tight_layout()
    
    vp_path = plot_dir / f"tpm_distribution_violin_{fname_suffix}.pdf"
    plt.savefig(vp_path, dpi=300, transparent=True)
    print(f"Saved {vp_path}")
    plt.close()

def main():
    print("Loading data...")
    # Load Datasets
    datasets = {}
    
    # In-house MCD - Only Pooled
    if "MCD Week pooled (combined)" in INHOUSE_MCD_FILES:
        fname = INHOUSE_MCD_FILES["MCD Week pooled (combined)"]
        datasets["MCD (in-house) | MCD Week pooled (combined)"] = {"df": load_df(DATA_DIR / fname), "species": "mouse"}
        
    # External MCD - All (GSE156918, GSE205974)
    for label, fname in EXTERNAL_MCD_FILES.items():
        datasets[f"MCD (external) | {label}"] = {"df": load_df(DATA_DIR / fname), "species": "mouse"}
        
    # Patients - NAS 1+ AND Fibrosis
    for ds_name, files in PATIENT_FILES.items():
        tpath = TPM_FILES.get(ds_name)
        # NAS 1+
        if "nas_1plus" in files:
            fname = files["nas_1plus"]
            label = f"Patient | {ds_name} | nas_1plus"
            datasets[label] = {"df": load_df(DATA_DIR / fname, tpath), "species": "human"}
        # Fibrosis
        if "fibrosis" in files:
            fname = files["fibrosis"]
            label = f"Patient | {ds_name} | fibrosis"
            datasets[label] = {"df": load_df(DATA_DIR / fname, tpath), "species": "human"}

    print(f"Loaded {len(datasets)} datasets.")

    # Load Orthologs
    print("Loading ortholog map...")
    ortho_df = load_ortholog_map(ORTHOLOG_PATH)
    m2h, h2m = build_maps(ortho_df, one2one_only=False)
    
    # Analysis Configuration
    log2fc_cutoffs = [0.5, 0.8, 1.0, 2.0]
    tpm_range = np.arange(0, 10.5, 0.5) # 0 to 10 step 0.5
    
    for lfc in log2fc_cutoffs:
        run_analysis(datasets, m2h, h2m, lfc, tpm_range)

if __name__ == "__main__":
    main()
