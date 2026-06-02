"""
spatial_utils.py — Shared utilities for the spatial transcriptomics pipeline.

Provides configuration loading, checkpoint management, data loading helpers,
and GPU initialization following the project-wide patterns.
"""

import pathlib
import sys
import yaml
import numpy as np
import pandas as pd

# Lazy imports: scanpy/anndata may not be available in all envs (e.g. gsmap)
try:
    import scanpy as sc
    import anndata as ad
except ImportError:
    sc = None
    ad = None

# ── Constants ────────────────────────────────────────────────────────────────
# cell2location prefixes cell type names with this string in obsm columns
C2L_PREFIX = "q05cell_abundance_w_sf_means_per_cluster_mu_fg_"

# ── Paths ────────────────────────────────────────────────────────────────────
PROJECT_ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SPATIAL_ROOT = PROJECT_ROOT / "Analysis/Spatial"
CONFIG_DIR = SPATIAL_ROOT / "config"
RESULTS_DIR = SPATIAL_ROOT / "results"
LOGS_DIR = SPATIAL_ROOT / "logs"


# ── Configuration ────────────────────────────────────────────────────────────

def load_config(config_name: str = "spatial_params.yaml") -> dict:
    """Load spatial pipeline YAML configuration."""
    with open(CONFIG_DIR / config_name) as f:
        return yaml.safe_load(f)


def load_dataset_config(config_name: str = "spatial_datasets.yaml") -> dict:
    """Load spatial dataset registry."""
    with open(CONFIG_DIR / config_name) as f:
        return yaml.safe_load(f)


def resolve_path(relative_path: str) -> pathlib.Path:
    """Resolve a project-relative path to absolute."""
    return PROJECT_ROOT / relative_path


# ── GPU Initialization ───────────────────────────────────────────────────────

def init_spatial_gpu():
    """Initialize GPU for cell2location / PyTorch-based models.

    Ensures CUDA is available for PyTorch. The spatial env does not use
    RAPIDS/RMM, so this is simpler than the SingleCell pipeline GPU init.
    """
    import torch
    if torch.cuda.is_available():
        device = torch.cuda.get_device_name(0)
        print(f"  GPU initialized: {device}")
    else:
        print("  WARNING: No GPU available, models will run on CPU")


def get_processor(use_gpu: bool = False):
    """Return scanpy module for preprocessing.

    The spatial pipeline uses scanpy (CPU) for preprocessing — the data
    is small enough (~200K spots) that GPU acceleration is unnecessary.
    """
    return sc


# ── Checkpoint Management ────────────────────────────────────────────────────

def checkpoint_path(name: str, subdir: str = "preprocessed") -> pathlib.Path:
    """Get path for a checkpoint file."""
    return RESULTS_DIR / subdir / name


def check_checkpoint(name: str, subdir: str = "preprocessed") -> bool:
    """Check if a checkpoint file exists (for resume-safe execution)."""
    path = checkpoint_path(name, subdir)
    if path.exists():
        size_gb = path.stat().st_size / 1e9
        print(f"  Checkpoint exists: {path} ({size_gb:.1f} GB)")
        return True
    return False


def save_checkpoint(adata: ad.AnnData, name: str, subdir: str = "preprocessed"):
    """Save AnnData checkpoint with size reporting."""
    outdir = RESULTS_DIR / subdir
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / name
    adata.write_h5ad(path)
    size_gb = path.stat().st_size / 1e9
    print(f"  Checkpoint saved: {path} ({size_gb:.1f} GB)")


def save_csv(df: pd.DataFrame, name: str, subdir: str = "preprocessed"):
    """Save DataFrame checkpoint as CSV."""
    outdir = RESULTS_DIR / subdir
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / name
    df.to_csv(path, index=True)
    print(f"  CSV saved: {path} ({len(df)} rows)")


# ── Data Loading ─────────────────────────────────────────────────────────────

def load_spatial_adata(name: str = "merged_spatial.h5ad",
                       subdir: str = "preprocessed") -> ad.AnnData:
    """Load preprocessed spatial AnnData from checkpoint."""
    path = checkpoint_path(name, subdir)
    print(f"Loading spatial data: {path}")
    return sc.read_h5ad(path)


def load_deconvolved_adata(name: str = "spatial_deconvolved.h5ad") -> ad.AnnData:
    """Load cell2location-deconvolved spatial AnnData."""
    path = checkpoint_path(name, "cell2location/spatial_model")
    print(f"Loading deconvolved data: {path}")
    return sc.read_h5ad(path)


def load_dream_degs(padj_thresh: float = 0.1, lfc_thresh: float = 0.0) -> pd.DataFrame:
    """Load dream mega-analysis results, optionally filtered.

    Adds 'symbol' column by mapping Ensembl IDs via the multi-evidence atlas.
    """
    config = load_config()
    path = resolve_path(config["paths"]["dream_results"])
    df = pd.read_csv(path)
    # Support both column naming conventions
    padj_col = "padj" if "padj" in df.columns else "adj.P.Val"
    if padj_thresh < 1.0 or lfc_thresh > 0:
        df = df[(df[padj_col] < padj_thresh) & (df["logFC"].abs() > lfc_thresh)]

    # Add symbol column if missing
    if "symbol" not in df.columns and "gene" in df.columns:
        ens_to_sym = build_ensembl_to_symbol_map()
        df["symbol"] = df["gene"].map(lambda g: ens_to_sym.get(g) or ens_to_sym.get(g.split(".")[0] if isinstance(g, str) else g))

    print(f"Loaded {len(df)} dream DEGs (padj<{padj_thresh}, |LFC|>{lfc_thresh})")
    return df


def build_ensembl_to_symbol_map() -> dict:
    """Build Ensembl ID → gene symbol mapping from multi-evidence atlas.

    Returns dict like {ENSG00000000003: "TSPAN6", ENSG00000000003.17: "TSPAN6", ...}.
    Falls back to dream_results gene column (versioned Ensembl IDs without symbols).
    """
    config = load_config()
    # Primary source: multi-evidence atlas has ensembl_id + human_symbol
    atlas_path = resolve_path(config["paths"]["multi_evidence_atlas"])
    mapping = {}
    if atlas_path.exists():
        df = pd.read_csv(atlas_path, usecols=["ensembl_id", "human_symbol"])
        df = df.dropna(subset=["ensembl_id", "human_symbol"])
        mapping = dict(zip(df["ensembl_id"], df["human_symbol"]))
        # Add version-stripped fallback
        for ensembl_id, symbol in list(mapping.items()):
            if isinstance(ensembl_id, str):
                stripped = ensembl_id.split(".")[0]
                if stripped not in mapping:
                    mapping[stripped] = symbol

    # Secondary: dream_results gene column → add versioned IDs if atlas missed them
    dream_path = resolve_path(config["paths"]["dream_results"])
    if dream_path.exists() and len(mapping) > 0:
        dream = pd.read_csv(dream_path, usecols=["gene"])
        for gene_id in dream["gene"]:
            stripped = gene_id.split(".")[0] if isinstance(gene_id, str) else gene_id
            if stripped in mapping and gene_id not in mapping:
                mapping[gene_id] = mapping[stripped]

    print(f"  Ensembl→symbol map: {len(mapping)} entries")
    return mapping


def load_deconv_scores() -> pd.DataFrame:
    """Load bulk deconvolution attribution scores."""
    config = load_config()
    path = resolve_path(config["paths"]["deconv_scores"])
    return pd.read_csv(path)


def load_multi_evidence_atlas() -> pd.DataFrame:
    """Load the multi-evidence target atlas."""
    config = load_config()
    path = resolve_path(config["paths"]["multi_evidence_atlas"])
    return pd.read_csv(path)


def load_conserved() -> list:
    """Load Conserved gene list from cross-species concordance atlas."""
    atlas_path = PROJECT_ROOT / "Analysis/Cross_Species_Concordance/results/unified_atlas"
    candidates = list(atlas_path.glob("cross_species_concordance_atlas*.csv"))
    if not candidates:
        # Fallback: try to get from multi-evidence atlas
        atlas = load_multi_evidence_atlas()
        if "is_conserved" in atlas.columns:
            sym_col = "human_symbol" if "human_symbol" in atlas.columns else "symbol"
            return atlas[atlas["is_conserved"] == True][sym_col].tolist()
        print("WARNING: Conserved gene list not found")
        return []
    df = pd.read_csv(candidates[0])
    col = "human_symbol" if "human_symbol" in df.columns else "symbol"
    cat_col = "primary_category" if "primary_category" in df.columns else None
    if cat_col:
        return df[df[cat_col] == "Conserved"][col].tolist()
    return df[col].tolist()


# ── Cell Type Harmonization ──────────────────────────────────────────────────
# Maps scATAC cell type names → cell2location cluster names (for MultiSP/ONTraC)
ATAC_TO_C2L = {
    "Hepatocyte": "Hepatocytes",
    "Endothelial": "Endothelial cells",
    "LSEC": "Endothelial cells",
    "Stellate_Cell": "Fibroblasts",
    "Kupffer_Cell": "Macrophages",
    "Macrophage": "Macrophages",
    "Cholangiocyte": "Cholangiocytes",
    "B_Cell": "B cells",
    "Plasma_Cell": "Plasma cells",
    "NK_T_Cell": "T cells",
    "T_Cell": "T cells",
    "NK_Cell": "T cells",
}


def harmonize_cell_types(cell_types, source="atac"):
    """Harmonize cell type names between modalities.

    Parameters
    ----------
    cell_types : list or pd.Series
        Cell type labels to harmonize.
    source : str
        Source modality ('atac' for scATAC → cell2location mapping).

    Returns
    -------
    list
        Harmonized cell type names.
    """
    mapping = ATAC_TO_C2L if source == "atac" else {}
    if isinstance(cell_types, pd.Series):
        return cell_types.map(lambda x: mapping.get(x, x)).tolist()
    return [mapping.get(ct, ct) for ct in cell_types]


def strip_c2l_prefix(name):
    """Strip cell2location column prefix from a cell type name."""
    if isinstance(name, str) and name.startswith(C2L_PREFIX):
        return name[len(C2L_PREFIX):]
    return name


# ── Enhanced Results Loaders ────────────────────────────────────────────────

def load_coloc_results(pp4_threshold=0.5):
    """Load COLOC results from all GWAS, filtered by PP.H4 threshold."""
    coloc_dir = PROJECT_ROOT / "RNA-seq/results/causal_inference"
    gwas_dirs = {
        "UKBB_ALT": coloc_dir / "broadaway_ukbb",
        "UKBB_AST": coloc_dir / "broadaway_ukbb_ast",
        "UKBB_GGT": coloc_dir / "broadaway_ukbb_ggt",
        "PDFF": coloc_dir / "broadaway_pdff",
    }
    all_coloc = []
    for gwas_name, gwas_dir in gwas_dirs.items():
        coloc_path = gwas_dir / "coloc_results.csv"
        if not coloc_path.exists():
            continue
        df = pd.read_csv(coloc_path)
        pp4_col = [c for c in df.columns if "PP.H4" in c or "pp4" in c.lower()]
        if pp4_col:
            df = df[df[pp4_col[0]] > pp4_threshold]
        df["gwas"] = gwas_name
        all_coloc.append(df)
    if not all_coloc:
        return pd.DataFrame()
    return pd.concat(all_coloc, ignore_index=True)


def load_drug_targets():
    """Load multi-layer drug targets."""
    config = load_config()
    path = resolve_path(config["paths"]["drug_targets"])
    if path.exists():
        return pd.read_csv(path)
    return pd.DataFrame()


def load_zonation_scores():
    """Load per-spot zonation scores from 04a results."""
    path = RESULTS_DIR / "zonation" / "spot_zonation_scores.csv"
    if path.exists():
        return pd.read_csv(path, index_col=0)
    return pd.DataFrame()


def fisher_test(set_a, set_b, universe):
    """2x2 Fisher's exact test for enrichment of set_a in set_b."""
    from scipy.stats import fisher_exact
    a = len(set_a & set_b)
    b = len(set_a - set_b)
    c = len(set_b - set_a)
    d = len(universe - set_a - set_b)
    if min(a + b, c + d, a + c, b + d) == 0:
        return np.nan, np.nan, 0
    odds, pval = fisher_exact([[a, b], [c, d]])
    return odds, pval, a


# ── SpaceRanger Helpers ──────────────────────────────────────────────────────

def parse_spaceranger_metrics(outs_dir: pathlib.Path) -> dict:
    """Parse SpaceRanger metrics_summary.csv into a dict."""
    csv_path = outs_dir / "metrics_summary.csv"
    if not csv_path.exists():
        return {"status": "MISSING", "path": str(outs_dir)}
    df = pd.read_csv(csv_path)
    metrics = {}
    for col in df.columns:
        val = df[col].iloc[0]
        if isinstance(val, str) and "%" in val:
            metrics[col] = float(val.replace("%", "").replace(",", ""))
        elif isinstance(val, str) and "," in val:
            metrics[col] = int(val.replace(",", ""))
        else:
            metrics[col] = val
    return metrics


def assess_spaceranger_qc(metrics: dict, thresholds: dict) -> tuple:
    """Assess SpaceRanger QC metrics against thresholds.

    Returns: (pass_qc: bool, failures: list[str])
    """
    checks = {
        "Reads Mapped Confidently to Genome": ("min_reads_mapped_pct", ">"),
        "Median Genes per Spot": ("min_median_genes", ">"),
        "Number of Spots Under Tissue": ("min_spots_under_tissue", ">"),
        "Sequencing Saturation": ("min_saturation", ">"),
    }
    failures = []
    for metric_name, (threshold_key, direction) in checks.items():
        if metric_name in metrics and threshold_key in thresholds:
            val = metrics[metric_name]
            thresh = thresholds[threshold_key]
            if direction == ">" and val < thresh:
                failures.append(f"{metric_name}={val} < {thresh}")
    return len(failures) == 0, failures


# ── Printing / Logging ───────────────────────────────────────────────────────

def print_header(title: str, width: int = 70):
    """Print a formatted section header."""
    print(f"\n{'=' * width}")
    print(f"  {title}")
    print(f"{'=' * width}\n")


def print_step(step: str, n: int = None, total: int = None):
    """Print a step progress indicator."""
    if n is not None and total is not None:
        print(f"  [{n}/{total}] {step}")
    else:
        print(f"  → {step}")
