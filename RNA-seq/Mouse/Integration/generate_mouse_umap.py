
# generate_mouse_umap.py
# Goal: Integrate Mouse RNA-seq datasets using Scanpy and generate UMAP plots.
# Supports GPU acceleration via USE_GPU environment variable.

import os
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import gzip
import csv
import logging

# GPU acceleration support
_GPU_UTILS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "..", "Analysis", "SingleCell", "scripts")
sys.path.insert(0, _GPU_UTILS)
from gpu_utils import init_gpu, get_processor, to_gpu, from_gpu

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

# GPU setup — set USE_GPU=1 environment variable to enable
USE_GPU = False
if os.environ.get("USE_GPU", "0") == "1":
    USE_GPU = init_gpu()
    if USE_GPU:
        print("GPU mode: ON (rapids_singlecell)")
    else:
        print("GPU mode: FAILED — falling back to CPU (scanpy)")
pp, tl = get_processor(USE_GPU)

# Set Global Plotting Settings for High Quality
sc.set_figure_params(dpi=300, dpi_save=300, format='pdf', vector_friendly=True, transparent=True)
plt.rcParams['pdf.fonttype'] = 42 # TrueType fonts

# Define Paths
BASE_DIR = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
INHOUSE_COUNTS = os.path.join(BASE_DIR, "Mouse/InHouse_MCD/counts/featurecounts/gene_counts.txt")
OTHER_COUNTS = os.path.join(BASE_DIR, "Mouse/Public_Diet_Models/counts/featurecounts/gene_counts.txt")
PUBLIC_MCD_COUNTS = os.path.join(BASE_DIR, "Mouse/Public_MCD/counts/featurecounts/gene_counts.txt")

INHOUSE_META = os.path.join(BASE_DIR, "Mouse/InHouse_MCD/metadata/samples.tsv")
OTHER_META = os.path.join(BASE_DIR, "Mouse/Public_Diet_Models/metadata/samples.tsv")
PUBLIC_MCD_META = os.path.join(BASE_DIR, "Mouse/Public_MCD/metadata/samples.tsv")
OUT_DIR = os.path.join(BASE_DIR, "Mouse/Integration")

# Ensure output directory exists
os.makedirs(OUT_DIR, exist_ok=True)

# Set Scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150, frameon=False, vector_friendly=True, format='pdf')
sc.settings.figdir = OUT_DIR

# ==============================================================================
# 1. Load and Clean Counts
# ==============================================================================

print("Loading count matrices...")

def load_featurecounts(path):
    # Determine skip rows
    with open(path) as f:
        first_line = f.readline()
    skip = 1 if not first_line.startswith("Geneid") else 0
    
    df = pd.read_csv(path, sep="\t", skiprows=skip)
    # Filter columns: keep Geneid and samples (cols > 6)
    meta_cols = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
    cols = df.columns
    sample_cols = [c for c in cols if c not in meta_cols]
    
    # Extract counts
    counts = df[sample_cols].copy()
    counts.index = df["Geneid"]
    
    return counts

cnts_inhouse = load_featurecounts(INHOUSE_COUNTS)
cnts_other = load_featurecounts(OTHER_COUNTS)
cnts_mcd = load_featurecounts(PUBLIC_MCD_COUNTS)

# Clean Column Names
# In-house: remove path, keep basename (e.g., 1CF)
cnts_inhouse.columns = [os.path.basename(c).split(".Aligned")[0] for c in cnts_inhouse.columns]

# Other: remove path, keep basename (e.g., SRR...)
cnts_other.columns = [os.path.basename(c).split(".Aligned")[0] for c in cnts_other.columns]
cnts_mcd.columns = [os.path.basename(c).split(".Aligned")[0] for c in cnts_mcd.columns]

print(f"In-house samples: {len(cnts_inhouse.columns)}")
print(f"Other samples: {len(cnts_other.columns)}")
print(f"Public MCD samples: {len(cnts_mcd.columns)}")

# Merge Counts (Inner Join on Genes)
counts_merged = pd.merge(cnts_inhouse, cnts_other, left_index=True, right_index=True, how="inner")
counts_merged = pd.merge(counts_merged, cnts_mcd, left_index=True, right_index=True, how="inner")
print(f"Merged dimensions: {counts_merged.shape}")

# Create AnnData
adata = sc.AnnData(X=counts_merged.T)

# ==============================================================================
# 2. Load and Prepare Metadata
# ==============================================================================

print("Loading metadata...")

meta_inhouse = pd.read_csv(INHOUSE_META, sep="\t")
meta_other = pd.read_csv(OTHER_META, sep="\t")
meta_mcd = pd.read_csv(PUBLIC_MCD_META, sep="\t")

# Map In-house
# Use 'week' column to create finer labels for MCD
# Condition = diet_week if diet == MCD, else Control
def map_inhouse_cond(row):
    diet = str(row['diet'])
    if diet == 'MCD':
        wk = str(row['week'])
        return f"MCD ({wk}wk)"
    return diet

meta_inhouse["Condition"] = meta_inhouse.apply(map_inhouse_cond, axis=1)
meta_inhouse["Dataset"] = "In-house_MCD"
meta_inhouse = meta_inhouse.set_index("sample_id")

# Map Other
meta_other["Condition"] = meta_other["condition"]
meta_other["Dataset"] = meta_other["dataset"]
meta_other = meta_other.set_index("sample_id")

# Map Public MCD
# Using 'diet' column (MCD vs Control)
meta_mcd["Condition"] = meta_mcd["diet"]
# Refine External MCD Timepoints based on GSM ID
# GSE156918 (GSM474...) -> 6 weeks
# GSE205974 (GSM623...) -> 4 weeks
meta_mcd = meta_mcd.set_index("sample_id") # Set index first for row.name access

def map_external_mcd_cond(row):
    cond = row['Condition']
    sid = str(row.name) # sample_id index
    if cond == 'MCD':
        if sid.startswith('GSM474'):
            return 'MCD (6wk)'
        elif sid.startswith('GSM623'):
            return 'MCD (4wk)'
    return cond

meta_mcd["Condition"] = meta_mcd.apply(map_external_mcd_cond, axis=1)
meta_mcd["Dataset"] = "External MCD"

# Combine
metadata_combined = pd.concat([
    meta_inhouse[["Condition", "Dataset"]],
    meta_other[["Condition", "Dataset"]],
    meta_mcd[["Condition", "Dataset"]]
])

# Annotate AnnData
common_samples = adata.obs_names.intersection(metadata_combined.index)
print(f"Matching samples: {len(common_samples)}")

adata = adata[common_samples]
adata.obs = metadata_combined.loc[common_samples]

# ==============================================================================
# 2b. Filter & Rename (User Request Refinements)
# ==============================================================================

print("Applying dataset exclusions and label simplifications...")

# 1. Exclude Datasets/Conditions
# Remove GSE263273 completely
adata = adata[adata.obs['Dataset'] != 'GSE263273']

# Remove GSE225616 (GAN) completely
adata = adata[adata.obs['Dataset'] != 'GSE225616']

# Remove High_Fat_Diet_DDD86481 from GSE224069
adata = adata[~((adata.obs['Dataset'] == 'GSE224069') & (adata.obs['Condition'] == 'High_Fat_Diet_DDD86481'))]

# 2. Rename Datasets
# 2. Rename Datasets
# Splitting GSE162876 (CDAHFD/FPC) by timepoint for better integration
# Identify 7wk vs 20wk samples based on Condition (which was merged from 'condition' column)
# Note: Metadata merging step above assigned 'Condition' from 'condition' column in samples.tsv
# We need to use the ORIGINAL 'condition' values if possible, or infer from current 'Condition' if mapped.
# Actually, 'Condition' in adata.obs is currently the raw condition from samples.tsv (before normalization step 3).
mask_162876 = adata.obs['Dataset'] == 'GSE162876'
adata.obs.loc[mask_162876 & adata.obs['Condition'].str.contains('7wk', case=False, na=False), 'Dataset'] = 'GSE162876_7wk'
adata.obs.loc[mask_162876 & adata.obs['Condition'].str.contains('20wk', case=False, na=False), 'Dataset'] = 'GSE162876_20wk'

# Splitting GSE159911 (LIDPAD) by timepoint using external metadata
print("Splitting GSE159911 (LIDPAD) by timepoint...")
gse159911_dir = os.path.join(BASE_DIR, "Mouse/Public_Diet_Models/GSE159911/metadata")
run_table_path = os.path.join(gse159911_dir, "GSE159911_SRARunTable.csv")
series_matrix_path = os.path.join(gse159911_dir, "GSE159911_series_matrix.txt.gz")

if os.path.exists(run_table_path) and os.path.exists(series_matrix_path):
    try:
        # 1. Map SRR (Run) -> GSM (Sample Name)
        srr_to_gsm = {}
        with open(run_table_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if 'Run' in row and 'Sample Name' in row:
                    srr_to_gsm[row['Run']] = row['Sample Name']
                    
        # 2. Map GSM -> Timepoint (from Series Matrix)
        gsm_to_time = {}
        with gzip.open(series_matrix_path, 'rt') as f:
            lines = f.readlines()
            
        geo_accessions = []
        time_characteristics = []
        
        # Parse matrix header
        for line in lines:
            line = line.strip()
            if line.startswith('!Sample_geo_accession'):
                # First col is label, rest are samples
                geo_accessions = line.split('\t')[1:]
                # Remove quotes if present
                geo_accessions = [g.strip('"') for g in geo_accessions]
            
            if line.startswith('!Sample_characteristics_ch1') and 'diet duration (weeks)' in line:
                vals = line.split('\t')[1:]
                time_characteristics = vals
                    
        # Build map
        if geo_accessions and time_characteristics:
             for i, gsm in enumerate(geo_accessions):
                 if i < len(time_characteristics):
                     desc = time_characteristics[i].strip('"')
                     if 'diet duration (weeks):' in desc:
                         try:
                             val = desc.split(':')[-1].strip()
                             val = val.split()[0]
                             gsm_to_time[gsm] = f'{val}wk'
                         except:
                             pass
        
        # 3. Update adata
        mask_lidpad = adata.obs['Dataset'].str.contains('GSE159911')
        count_updated = 0
        for idx in adata.obs[mask_lidpad].index:
            # Index is typical sample ID
            srr = idx
            if srr in srr_to_gsm:
                gsm = srr_to_gsm[srr]
                if gsm in gsm_to_time:
                    tp = gsm_to_time[gsm]
                    # Directly set to Condition to preserve Batch ID (Dataset) for RegressOut
                    adata.obs.loc[idx, 'Condition'] = f'LIDPAD ({tp})'
                    count_updated += 1

                    
        print(f"Updated {count_updated} LIDPAD samples with timepoints.")
        
    except Exception as e:
        print(f"Failed to refine LIDPAD metadata: {e}")
else:
    print("Warning: GSE159911 metadata files not found. Skipping LIDPAD refinement.")

ds_map = {
    'In-house_MCD': 'In-House MCD',
    'GSE159911': 'LIDPAD', # Fallback
    'GSE159911_3wk': 'LIDPAD (3wk)',
    'GSE159911_6wk': 'LIDPAD (6wk)',
    'GSE159911_12wk': 'LIDPAD (12wk)',
    'GSE162876_7wk': 'CDAHFD (7wk)',
    'GSE162876_20wk': 'FPC (20wk)',
    'GSE274914': 'HFD #1',
    'GSE224069': 'HFD #2'
}
# Ensure string type for reliable replacement
adata.obs['Dataset'] = adata.obs['Dataset'].astype(str).replace(ds_map)

# 3. Rename Conditions
cond_map = {
    '7wks_LFD_Ctrl': 'Control',
    '7wks_CDAHFD': 'CDAHFD (7wk)',
    '20wks_LFD_Ctrl': 'Control',
    '20wks_FPC': 'FPC (20wk)',
    '7w_LFD': 'Control',
    '7w_HFD': 'HFD (7wk)',
    '52w_LFD': 'Control',
    '52w_HFD': 'HFD (52wk)',
    '52w_HFD': 'HFD (52wk)',
    'Chow_Diet_Vehicle': 'Control',
    'High_Fat_Diet_Vehicle': 'HFD (19wk)'
}
adata.obs['Condition'] = adata.obs['Condition'].astype(str).replace(cond_map)

# 4. Map External MCD to explicit label (Already handled in Step 2)
# adata.obs.loc[(adata.obs['Dataset'] == 'External MCD') & (adata.obs['Condition'] == 'MCD'), 'Condition'] = 'MCD (External)'

# 5. Reorder Categories to force Control first
# 5. Reorder Categories to force Temporal Order
# Explicitly define the desired order for the legend
desired_order = [
    'Control',
    # LIDPAD Timecourse
    'LIDPAD (1wk)', 'LIDPAD (3wk)', 'LIDPAD (4wk)', 'LIDPAD (8wk)', 
    'LIDPAD (12wk)', 'LIDPAD (16wk)', 'LIDPAD (32wk)', 'LIDPAD (48wk)',
    # MCD Timecourse
    'MCD (1wk)', 'MCD (2wk)', 'MCD (3wk)', 'MCD (4wk)', 'MCD (6wk)', 'MCD (8wk)', 'MCD (External)',
    # HFD Timecourse
    'HFD (7wk)', 'HFD (19wk)', 'HFD (52wk)', 'HFD #1', 'HFD #2',
    # Other Models
    'CDAHFD (7wk)',
    'FPC (20wk)'
]

# Get actual unique categories present in the data
present_conditions = adata.obs['Condition'].unique().tolist()
# Filter desired_order to only include present conditions to avoid processing errors
ordered_present = [c for c in desired_order if c in present_conditions]
# Add any missing conditions at the end (just in case)
missing = [c for c in present_conditions if c not in ordered_present]
final_order = ordered_present + sorted(missing)

adata.obs['Condition'] = pd.Categorical(adata.obs['Condition'], categories=final_order, ordered=True)

print(f"Post-filtering samples: {adata.shape[0]}")
print(f"Unique Datasets: {adata.obs['Dataset'].unique()}")
print(f"Unique Conditions: {adata.obs['Condition'].unique()}")

# ==============================================================================
# 3. Preprocessing & Normalization
# ==============================================================================

print("Preprocessing...")

# Transfer to GPU for accelerated preprocessing
if USE_GPU:
    print("Transferring to GPU...")
    to_gpu(adata)

# Filter genes
pp.filter_genes(adata, min_cells=3)

# Normalize (CPM for bulk UMAP)
pp.normalize_total(adata, target_sum=1e6)
pp.log1p(adata)

# Highly Variable Genes (Batch Aware)
print("Selecting HVGs per batch...")
pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key='Dataset')

if USE_GPU:
    from_gpu(adata)
adata.raw = adata  # Store full normalized data
adata = adata[:, adata.var.highly_variable]

print(f"HVGs: {adata.shape[1]}")

# Regress out Dataset effect (always CPU — regress_out uses statsmodels)
print("Regressing out Dataset effect (Strict Correction)...")
sc.pp.regress_out(adata, ['Dataset'])

# Scale + PCA on GPU if available
if USE_GPU:
    to_gpu(adata)

pp.scale(adata, max_value=10)
pp.pca(adata)

if USE_GPU:
    from_gpu(adata)

# QC/Outlier Check: Remove small disconnected clusters
# DISABLED per user request (Step 1307)
# print("Checking for outliers...")
# sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')
# sc.tl.leiden(adata, resolution=0.5, key_added='leiden_qc')
# cluster_counts = adata.obs['leiden_qc'].value_counts()
# print("Cluster sizes:")
# print(cluster_counts)
# 
# # Identify small clusters (< 10 samples) as outliers
# small_clusters = cluster_counts[cluster_counts < 10].index
# if len(small_clusters) > 0:
#     # print(f"Removing {len(small_clusters)} outlier clusters: {small_clusters.tolist()}")
#     # adata = adata[~adata.obs['leiden_qc'].isin(small_clusters)].copy()
#     # print(f"New shape: {adata.shape}")
#     print(f"Skipping removal of {len(small_clusters)} outlier clusters (User Requested).")
# else:
#     print("No obvious outlier clusters found.")

# ==============================================================================
# 4. Batch Correction
# ==============================================================================

print("Batch Correction...")
# Use Harmony with stricter diversity penalty if needed
try:
    import scanpy.external as sce
    print("Using Harmony integration...")
    # theta=2 default. Increase to force mixing if persistent separation.
    # theta=3 for Mouse Severe Integration (keep consistent with v2 tuning, regression adds power)
    sce.pp.harmony_integrate(adata, 'Dataset', theta=3, max_iter_harmony=20)
    use_rep = 'X_pca_harmony'
except (ImportError, AttributeError):
    print("Harmony not found. Using Regress Out / Combat logic...")
    # sc.pp.combat(adata, key='Dataset') # ComBat often robust
    # Or regress out
    sc.pp.regress_out(adata, ['Dataset'])
    sc.tl.pca(adata, svd_solver='arpack') # Re-run PCA on regressed data
    use_rep = 'X_pca'

# ==============================================================================
# 5. UMAP and Plotting
# ==============================================================================

print("Running UMAP...")
# GPU-accelerated neighbors + UMAP
if USE_GPU:
    to_gpu(adata)

pp.neighbors(adata, n_neighbors=50, n_pcs=40, use_rep=use_rep)
tl.umap(adata, min_dist=0.5)

if USE_GPU:
    from_gpu(adata)

# Verification: Print counts per dataset
print("Final Dataset Counts:")
print(adata.obs['Dataset'].value_counts())

print("Plotting...")
# Randomize sample order for reproducible z-ordering
np.random.seed(42)
idx = np.arange(adata.shape[0])
np.random.shuffle(idx)
adata = adata[idx].copy()

# Plot UMAP
# ------------------------------------------------------------------------------

# Define Custom Palette
color_dict = {
    'Control': '#d3d3d3',  # Light Grey for all controls
    
    # LIDPAD Timecourse (Greens: Light -> Dark)
    'LIDPAD (1wk)': '#E5F5E0',
    'LIDPAD (3wk)': '#C7E9C0',
    'LIDPAD (4wk)': '#A1D99B', 
    'LIDPAD (8wk)': '#74C476',
    'LIDPAD (12wk)': '#41AB5D',
    'LIDPAD (16wk)': '#238B45',
    'LIDPAD (32wk)': '#006D2C',
    'LIDPAD (48wk)': '#00441B',

    # HFD Timecourse (Purples: Light -> Dark)
    'HFD (7wk)': '#BCBDDC',
    'HFD (19wk)': '#9E9AC8',
    'HFD (52wk)': '#756BB1',
    'HFD #1': '#BCBDDC', # Fallback if needed
    'HFD #2': '#756BB1', # Fallback

    # Other Models (Reds/Oranges)
    'MCD': '#CB181D',         # Fallback
    'MCD (External)': '#CB181D',
    'MCD (1wk)': '#FCBBA1',
    'MCD (2wk)': '#FB6A4A',
    'MCD (3wk)': '#CB181D', 
    'MCD (4wk)': '#67000D',   # Very Dark Red
    'MCD (6wk)': '#800026',   # Darkest Red
    'MCD (8wk)': '#49000a',   # Ultra Dark Red
    'In-House MCD': '#EF3B2C',
    'CDAHFD (7wk)': '#008B8B', # Dark Cyan (Distinct from Reds)
    'FPC (20wk)': '#8B4513'    # Saddle Brown (Distinct from Reds)
}

# 1. Plot Condition with Custom Palette
print("Plotting Condition UMAP...")
sc.pl.umap(adata, color='Condition', 
           palette=color_dict,
           save="_mouse_integrated_condition.pdf", 
           title="Integrated Mouse Atlas: Disease Condition",
           frameon=False,
           size=30,           
           legend_loc='right margin',
           legend_fontsize=8)

# 2. Plot Dataset with Default Palette (Distinct batches)
print("Plotting Dataset UMAP...")
sc.pl.umap(adata, color='Dataset', 
           save="_mouse_integrated_dataset.pdf", 
           title="Integrated Mouse Atlas: Dataset Source",
           frameon=False,
           size=30,           
           legend_loc='right margin',
           legend_fontsize=8)

print(f"Done! Saved plots to {OUT_DIR}/umap_mouse_integrated.pdf")

print("Generating Trajectory Analysis...")

# 6. Diffusion Pseudotime (DPT)
# -----------------------------
print("Calculating DPT...")
# Find root: Control samples
#Condition is 'Control'
healthy_indices = np.where(adata.obs['Condition'] == 'Control')[0]

if len(healthy_indices) > 0:
    # Pick the one with the lowest pseudotime? No DPT yet.
    # Simple approach: Pick the first one. 
    root_idx = healthy_indices[0] 
    print(f"Set root to cell index {root_idx} (Control)")
    
    adata.uns['iroot'] = root_idx
    
    # Run Diffusion Map first
    sc.tl.diffmap(adata)
    
    # Run DPT
    sc.tl.dpt(adata)
    
    # Check for infinite values (disconnected)
    n_inf = np.isinf(adata.obs['dpt_pseudotime']).sum()
    if n_inf > 0:
        print(f"WARNING: {n_inf} samples have infinite pseudotime (disconnected components).")
        print("Try increasing n_neighbors further if this persists.")
    
    # Plot DPT
    print("Plotting Pseudotime...")
    sc.pl.umap(adata, color=['dpt_pseudotime'], color_map='viridis', 
               save="_mouse_integrated_pseudotime.pdf",
               frameon=False, size=30, legend_loc='right margin')
    
    # 7. Grid-based Arrow Flow
    # ------------------------
    print("generating Arrow Flow...")
    
    def plot_arrows(adata_obj, basis='umap', grid_size=20):
        import numpy as np
        from scipy.interpolate import griddata, NearestNDInterpolator
        from scipy.ndimage import gaussian_filter
        
        # Get coordinates
        coords = adata_obj.obsm[f'X_{basis}']
        x, y = coords[:, 0], coords[:, 1]
        
        # Get pseudotime and handle non-finite values
        pt = adata_obj.obs['dpt_pseudotime'].values
        valid = np.isfinite(pt)
        x, y, pt = x[valid], y[valid], pt[valid]
        
        # Create grid
        x_min, x_max = x.min(), x.max()
        y_min, y_max = y.min(), y.max()
        
        # Expand slightly
        margin_x = (x_max - x_min) * 0.05
        margin_y = (y_max - y_min) * 0.05
        
        x_grid = np.linspace(x_min - margin_x, x_max + margin_x, grid_size)
        y_grid = np.linspace(y_min - margin_y, y_max + margin_y, grid_size)
        
        grid_x, grid_y = np.meshgrid(x_grid, y_grid)
        
        # Interpolate DPT values onto grid
        # Use linear interpolation first
        grid_dpt = griddata((x, y), pt, (grid_x, grid_y), method='linear')
        
        # Handle NaNs (outside convex hull) using nearest
        mask = np.isnan(grid_dpt)
        if np.any(mask):
             interp = NearestNDInterpolator(list(zip(x,y)), pt)
             grid_dpt_filled = interp(grid_x, grid_y)
             # Smooth the surface
             grid_dpt_smooth = gaussian_filter(grid_dpt_filled, sigma=1)
        else:
             grid_dpt_smooth = gaussian_filter(grid_dpt, sigma=1)
             
        # ROI Mask to prevent arrows in empty space?
        # Ideally we only show arrows near data.
        # But for now, showing global field is okay, or we can mask by density.
             
        # Calculate gradient
        # dy, dx = np.gradient(grid_dpt_smooth) 
        # Note: Gradient points in direction of increase (Healthy -> Disease)
        dy, dx = np.gradient(grid_dpt_smooth)
        
        # Plot
        plt.figure(figsize=(8,8))
        # Plot points
        plt.scatter(x, y, c='lightgrey', s=20, alpha=0.5)
        
        # Plot arrows
        plt.quiver(grid_x, grid_y, dx, dy, color='black', scale=None) # Auto scale
        plt.title('Disease Trajectory (DPT Gradient)')
        plt.axis('off')
        plt.savefig(os.path.join(OUT_DIR, "mouse_integrated_arrows.pdf"))
        plt.close()

    try:
        plot_arrows(adata)
    except Exception as e:
        print(f"Arrow plotting failed: {e}")

else:
    print("WARNING: No 'Control' samples found! Skipping DPT.")
