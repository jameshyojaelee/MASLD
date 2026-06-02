#!/usr/bin/env python3
"""
100_bulkformer_embeddings.py
Extract per-sample embeddings from the pre-trained BulkFormer model (147M).

BulkFormer: 147M-parameter GNN+Performer pre-trained on >500K bulk RNA-seq
profiles. Combines graph convolutions (gene-gene interactions from TCGA PPI)
with Performer attention for global expression dependencies.

Pipeline:
  1. Clone BulkFormer repo + download pretrained weights + data from Zenodo
  2. Install required dependencies (performer-pytorch, torch-geometric, gdown)
  3. Load our expression data from merged_dge.rds (full gene set, raw counts)
  4. Convert raw counts -> log1p(TPM) using BulkFormer's gene_length_df
  5. Align to BulkFormer's 20,010-gene vocabulary (missing genes -> -10)
  6. Run BulkFormer encoder -> per-sample 643-dim embeddings (mean aggregation)
  7. Save embeddings as CSV

Fallback: If BulkFormer inference fails (dependency/GPU/gene mismatch), fall
back to the existing VAE 64-dim embeddings and log the error.

Input:
  - results/integration/merged_dge.rds (full counts matrix)
  - results/staging_classifier/prepared_data.h5 (sample IDs, gene names)
  - results/staging_classifier/embeddings_all_samples.csv (fallback VAE)

Output:
  - results/staging_classifier/bulkformer_embeddings.csv (samples x embed_dim)
  - results/staging_classifier/bulkformer_gene_mapping.csv (gene mapping log)

SLURM: gpu partition, 1xL40S, 8 CPUs, 64GB RAM, 48h
Env:   micromamba activate rapids_singlecell

Usage:
  sbatch --job-name=stg100_bulkformer \
         --partition=gpu --gres=gpu:1 \
         --cpus-per-task=8 --mem=64G --time=48:00:00 \
         --output=logs/100_bulkformer_%j.out \
         --error=logs/100_bulkformer_%j.err \
         --wrap="bash -c 'eval \"\\$(micromamba shell hook --shell bash)\" && \\
                 micromamba activate rapids_singlecell && \\
                 cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \\
                 python 100_bulkformer_embeddings.py'"
"""

import os
import sys
import time
import subprocess
import warnings
import logging
import shutil
from pathlib import Path
from collections import OrderedDict

import numpy as np
import pandas as pd
import h5py

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
INT = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUTDIR = os.path.join(INT, "results/staging_classifier")
LOGDIR = os.path.join(OUTDIR, "logs")
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(LOGDIR, exist_ok=True)

# BulkFormer will be cloned here
BULKFORMER_DIR = os.path.join(BASE, "tools/BulkFormer")
BULKFORMER_REPO = "https://github.com/KangBoming/BulkFormer.git"

# Zenodo data URLs (required model data files)
ZENODO_BASE = "https://zenodo.org/records/15744294/files"
ZENODO_DATA_FILES = {
    "G_tcga.pt": f"{ZENODO_BASE}/G_tcga.pt?download=1",
    "G_tcga_weight.pt": f"{ZENODO_BASE}/G_tcga_weight.pt?download=1",
    "esm2_feature_concat.pt": f"{ZENODO_BASE}/esm2_feature_concat.pt?download=1",
    "interested_gene_list.pt": f"{ZENODO_BASE}/interested_gene_list.pt?download=1",
    "bulkformer_gene_info.csv": f"{ZENODO_BASE}/bulkformer_gene_info.csv?download=1",
}

# Google Drive ID for the 147M model
MODEL_GDRIVE_ID = "1UtqN_vCh3669Fs-GU5CTE7F7UnuQCAzN"
MODEL_FILENAME = "bulkformer_147M.pt"

# Inputs
MERGED_DGE_RDS = os.path.join(INT, "results/integration/merged_dge.rds")
H5_PATH = os.path.join(OUTDIR, "prepared_data.h5")
FALLBACK_EMBEDDINGS = os.path.join(OUTDIR, "embeddings_all_samples.csv")

# Outputs
OUTPUT_EMBEDDINGS = os.path.join(OUTDIR, "bulkformer_embeddings.csv")
OUTPUT_GENE_MAP = os.path.join(OUTDIR, "bulkformer_gene_mapping.csv")

# Rscript path
RSCRIPT = "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(os.path.join(LOGDIR, "100_bulkformer.log")),
    ],
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------
def run_cmd(cmd, desc="", check=True, timeout=3600):
    """Run a shell command with logging."""
    log.info(f"CMD [{desc}]: {cmd}")
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, timeout=timeout
    )
    if result.stdout.strip():
        log.info(f"  stdout: {result.stdout.strip()[:500]}")
    if result.stderr.strip():
        log.warning(f"  stderr: {result.stderr.strip()[:500]}")
    if check and result.returncode != 0:
        raise RuntimeError(
            f"Command failed (rc={result.returncode}): {cmd}\n{result.stderr[:1000]}"
        )
    return result


def download_file(url, dest_path, desc=""):
    """Download a file via wget (curl fallback)."""
    if os.path.exists(dest_path):
        log.info(f"  Already exists: {dest_path}")
        return
    log.info(f"  Downloading {desc}: {url} -> {dest_path}")
    try:
        run_cmd(
            f'wget -q --no-check-certificate -O "{dest_path}" "{url}"',
            desc=f"wget {desc}",
            timeout=7200,
        )
    except RuntimeError:
        run_cmd(
            f'curl -sL -o "{dest_path}" "{url}"',
            desc=f"curl {desc}",
            timeout=7200,
        )


def use_fallback(reason):
    """Copy existing VAE embeddings as fallback output."""
    log.warning(f"FALLBACK: {reason}")
    log.warning(f"Using existing VAE embeddings from {FALLBACK_EMBEDDINGS}")
    if not os.path.exists(FALLBACK_EMBEDDINGS):
        log.error("Fallback VAE embeddings not found either. Cannot proceed.")
        sys.exit(1)
    df = pd.read_csv(FALLBACK_EMBEDDINGS)
    # Rename columns to indicate fallback
    rename = {c: c.replace("z", "bf") for c in df.columns if c.startswith("z")}
    df = df.rename(columns=rename)
    df.to_csv(OUTPUT_EMBEDDINGS, index=False)
    log.info(f"Fallback embeddings saved: {OUTPUT_EMBEDDINGS} ({df.shape})")
    log.info(f"NOTE: These are VAE embeddings (64-dim), not BulkFormer (643-dim)")
    return df


# ============================================================================
# STEP 0: Environment setup — install missing dependencies
# ============================================================================
def install_dependencies():
    """Install BulkFormer's required packages into the current environment."""
    import importlib

    pip_exe = sys.executable.replace("python", "pip")
    if not os.path.exists(pip_exe):
        pip_exe = f"{sys.executable} -m pip"

    deps = {
        "performer_pytorch": "performer-pytorch",
        "torch_geometric": "torch-geometric",
        "gdown": "gdown",
    }

    for module_name, pip_name in deps.items():
        try:
            importlib.import_module(module_name)
            log.info(f"  {module_name}: already installed")
        except ImportError:
            log.info(f"  Installing {pip_name}...")
            if pip_name == "torch-geometric":
                # torch-geometric needs special install for compatibility
                run_cmd(
                    f"{sys.executable} -m pip install --no-deps torch-geometric",
                    desc=f"pip install {pip_name}",
                    timeout=600,
                )
                # Also install pyg-lib, torch-sparse, torch-scatter, torch-cluster
                import torch
                torch_ver = torch.__version__.split("+")[0]
                cuda_tag = torch.__version__.split("+")[1] if "+" in torch.__version__ else "cpu"
                pyg_url = f"https://data.pyg.org/whl/torch-{torch_ver}+{cuda_tag}.html"
                for pyg_dep in ["torch-sparse", "torch-scatter", "torch-cluster"]:
                    try:
                        run_cmd(
                            f"{sys.executable} -m pip install --no-deps {pyg_dep} -f {pyg_url}",
                            desc=f"pip install {pyg_dep}",
                            timeout=600,
                        )
                    except RuntimeError:
                        log.warning(f"  Could not install {pyg_dep} from PyG wheels, trying plain pip")
                        try:
                            run_cmd(
                                f"{sys.executable} -m pip install --no-deps {pyg_dep}",
                                desc=f"pip install {pyg_dep} (fallback)",
                                timeout=600,
                            )
                        except RuntimeError:
                            log.warning(f"  {pyg_dep} installation failed — will try to proceed")
            else:
                run_cmd(
                    f"{sys.executable} -m pip install --no-deps {pip_name}",
                    desc=f"pip install {pip_name}",
                    timeout=600,
                )


# ============================================================================
# STEP 1: Clone BulkFormer and download data
# ============================================================================
def setup_bulkformer():
    """Clone repo, download model weights and data files from Zenodo."""
    bf_root = Path(BULKFORMER_DIR)
    bf_inner = bf_root / "BulkFormer"

    # Clone if needed
    if not (bf_inner / "utils" / "BulkFormer.py").exists():
        log.info("Cloning BulkFormer repository...")
        bf_root.mkdir(parents=True, exist_ok=True)
        if bf_inner.exists():
            shutil.rmtree(bf_inner)
        run_cmd(
            f"git clone {BULKFORMER_REPO} {bf_inner}",
            desc="git clone BulkFormer",
            timeout=300,
        )
    else:
        log.info("BulkFormer repo already cloned.")

    # Download Zenodo data files -> BulkFormer/data/
    data_dir = bf_inner / "data"
    data_dir.mkdir(exist_ok=True)
    for fname, url in ZENODO_DATA_FILES.items():
        download_file(url, str(data_dir / fname), desc=fname)

    # Download model weights from Google Drive -> BulkFormer/model/
    model_dir = bf_inner / "model"
    model_dir.mkdir(exist_ok=True)
    model_path = model_dir / MODEL_FILENAME
    if not model_path.exists():
        log.info(f"Downloading {MODEL_FILENAME} from Google Drive...")
        try:
            import gdown
            gdown.download(
                id=MODEL_GDRIVE_ID,
                output=str(model_path),
                quiet=False,
            )
        except Exception as e:
            log.warning(f"gdown failed: {e}")
            # Try direct URL fallback
            gdrive_url = f"https://drive.google.com/uc?export=download&id={MODEL_GDRIVE_ID}"
            try:
                download_file(gdrive_url, str(model_path), desc="model weights (direct)")
            except Exception as e2:
                log.warning(f"Direct download also failed: {e2}")
                # Try with gdown CLI
                try:
                    run_cmd(
                        f"{sys.executable} -m gdown '{MODEL_GDRIVE_ID}' -O '{model_path}' --fuzzy",
                        desc="gdown CLI",
                        timeout=3600,
                    )
                except RuntimeError:
                    raise RuntimeError(
                        f"Cannot download model weights. Please manually download "
                        f"bulkformer_147M.pt from Google Drive ID {MODEL_GDRIVE_ID} "
                        f"and place it at {model_path}"
                    )
    else:
        log.info(f"Model weights already exist: {model_path}")

    return bf_inner


# ============================================================================
# STEP 2: Load our expression data (full gene set from merged_dge.rds)
# ============================================================================
def load_expression_data():
    """Load raw counts from merged_dge.rds and sample IDs from H5."""
    log.info("Loading expression data...")

    # Get sample IDs from H5
    with h5py.File(H5_PATH, "r") as h5:
        sample_ids = np.array(
            [s.decode() if isinstance(s, bytes) else s for s in h5["sample_ids"][:]]
        )
    log.info(f"  H5 sample IDs: {len(sample_ids)}")

    # Convert merged_dge.rds to CSV via Rscript
    counts_csv = os.path.join(OUTDIR, "_bulkformer_raw_counts.csv")
    genes_csv = os.path.join(OUTDIR, "_bulkformer_gene_info.csv")

    if not os.path.exists(counts_csv) or not os.path.exists(genes_csv):
        log.info("  Converting merged_dge.rds -> CSV via Rscript...")
        r_code = f"""
library(edgeR)
dge <- readRDS("{MERGED_DGE_RDS}")
# Write raw counts (genes x samples)
counts_df <- as.data.frame(dge$counts)
counts_df$gene_id <- rownames(dge$counts)
write.csv(counts_df, "{counts_csv}", row.names=FALSE)
# Write gene info
gene_info <- data.frame(
    gene_id = rownames(dge$genes),
    gene_symbol = dge$genes$gene_name,
    gene_length = dge$genes$gene_length,
    stringsAsFactors = FALSE
)
write.csv(gene_info, "{genes_csv}", row.names=FALSE)
cat("Counts:", nrow(dge$counts), "genes x", ncol(dge$counts), "samples\\n")
cat("Gene info:", nrow(gene_info), "rows\\n")
"""
        r_script_path = os.path.join(OUTDIR, "_bulkformer_export.R")
        with open(r_script_path, "w") as f:
            f.write(r_code)
        run_cmd(
            f"{RSCRIPT} {r_script_path}",
            desc="export counts from RDS",
            timeout=600,
        )

    # Load counts
    log.info("  Loading counts CSV...")
    counts_df = pd.read_csv(counts_csv)
    gene_id_col = counts_df.pop("gene_id")

    # Load gene info
    gene_info = pd.read_csv(genes_csv)
    gene_info = gene_info.set_index("gene_id")

    # Ensure column order matches H5 sample IDs
    # counts_df columns are sample IDs
    available_samples = [s for s in sample_ids if s in counts_df.columns]
    if len(available_samples) < len(sample_ids):
        log.warning(
            f"  {len(sample_ids) - len(available_samples)} samples in H5 not found in DGE"
        )
    counts_df = counts_df[available_samples]
    counts_df.index = gene_id_col.values

    log.info(f"  Raw counts: {counts_df.shape[0]} genes x {counts_df.shape[1]} samples")
    log.info(f"  Gene info: {gene_info.shape[0]} genes")

    return counts_df, gene_info, np.array(available_samples)


# ============================================================================
# STEP 3: Preprocess for BulkFormer (log1p-TPM + gene alignment)
# ============================================================================
def preprocess_for_bulkformer(counts_df, gene_info, bf_dir):
    """Convert raw counts to log1p(TPM), align to BulkFormer gene vocabulary."""
    log.info("Preprocessing for BulkFormer...")

    # Load BulkFormer gene info (vocabulary)
    bf_gene_info_path = bf_dir / "data" / "bulkformer_gene_info.csv"
    bf_gene_info = pd.read_csv(bf_gene_info_path)
    bf_gene_list = bf_gene_info["ensg_id"].tolist()
    log.info(f"  BulkFormer vocabulary: {len(bf_gene_list)} genes")
    log.info(f"  First 3 BulkFormer genes: {bf_gene_list[:3]}")

    # Load gene length dictionary from BulkFormer
    bf_gene_length_path = bf_dir / "data" / "gene_length_df.csv"
    if bf_gene_length_path.exists():
        bf_gene_length = pd.read_csv(bf_gene_length_path)
        gene_length_dict = bf_gene_length.set_index("ensg_id")["length"].to_dict()
        log.info(f"  BulkFormer gene lengths: {len(gene_length_dict)} entries")
    else:
        log.warning("  gene_length_df.csv not found, using our gene lengths")
        gene_length_dict = {}

    # Create gene ID mapping: our versioned Ensembl IDs -> BulkFormer's IDs
    # Our genes: ENSG00000136709.13 ; BulkFormer may use ENSG00000136709 (no version)
    our_genes = counts_df.index.tolist()

    # Strip version from our gene IDs
    our_gene_base = {g.split(".")[0]: g for g in our_genes}

    # Check BulkFormer gene format
    bf_has_version = any("." in g for g in bf_gene_list[:100])
    if bf_has_version:
        log.info("  BulkFormer uses versioned Ensembl IDs")
        # Direct matching with version
        gene_map = {g: g for g in our_genes if g in set(bf_gene_list)}
        # Also try base ID matching for misses
        bf_base = {g.split(".")[0]: g for g in bf_gene_list}
        for our_g in our_genes:
            if our_g not in gene_map:
                base = our_g.split(".")[0]
                if base in bf_base:
                    gene_map[our_g] = bf_base[base]
    else:
        log.info("  BulkFormer uses unversioned Ensembl IDs")
        bf_gene_set = set(bf_gene_list)
        gene_map = {}
        for our_g in our_genes:
            base = our_g.split(".")[0]
            if base in bf_gene_set:
                gene_map[our_g] = base

    log.info(f"  Mapped {len(gene_map)} / {len(our_genes)} genes to BulkFormer vocab")
    log.info(
        f"  Missing from BulkFormer: {len(bf_gene_list) - len(gene_map)} genes "
        f"will be filled with -10"
    )

    # Save gene mapping for reference
    map_df = pd.DataFrame(
        [
            {"our_gene_id": k, "bulkformer_gene_id": v, "matched": True}
            for k, v in gene_map.items()
        ]
    )
    unmatched = [g for g in our_genes if g not in gene_map]
    map_df = pd.concat(
        [
            map_df,
            pd.DataFrame(
                [
                    {"our_gene_id": g, "bulkformer_gene_id": "", "matched": False}
                    for g in unmatched[:100]  # cap for file size
                ]
            ),
        ],
        ignore_index=True,
    )
    map_df.to_csv(OUTPUT_GENE_MAP, index=False)
    log.info(f"  Gene mapping saved: {OUTPUT_GENE_MAP}")

    # BulkFormer normalization: raw counts -> log1p(TPM)
    # TPM = (count / gene_length_kb) / sum(count / gene_length_kb) * 1e6
    log.info("  Computing log1p(TPM)...")

    # Use BulkFormer gene lengths where available, fallback to our lengths
    counts_matrix = counts_df.values.astype(np.float64)  # (genes, samples)

    # Get gene lengths
    gene_lengths = np.ones(counts_matrix.shape[0]) * 1000  # default 1kb
    for i, g in enumerate(our_genes):
        # Try BulkFormer length dict first (may use base ID)
        base_g = g.split(".")[0]
        if g in gene_length_dict:
            gene_lengths[i] = gene_length_dict[g]
        elif base_g in gene_length_dict:
            gene_lengths[i] = gene_length_dict[base_g]
        elif g in gene_info.index and "gene_length" in gene_info.columns:
            gl = gene_info.loc[g, "gene_length"]
            if pd.notna(gl) and gl > 0:
                gene_lengths[i] = gl

    gene_lengths_kb = gene_lengths / 1000.0  # convert to kb

    # TPM calculation (per-sample)
    rate = counts_matrix / gene_lengths_kb[:, np.newaxis]  # (genes, samples)
    rate_sum = rate.sum(axis=0, keepdims=True)  # (1, samples)
    rate_sum[rate_sum == 0] = 1e-6
    tpm = rate / rate_sum * 1e6
    log_tpm = np.log1p(tpm)  # log1p(TPM)

    # Create DataFrame with our gene IDs -> BulkFormer gene IDs
    # Transpose: (genes, samples) -> (samples, genes) for main_gene_selection
    log_tpm_df = pd.DataFrame(
        log_tpm.T,
        index=counts_df.columns,  # sample IDs
        columns=[gene_map.get(g, g.split(".")[0]) for g in our_genes],
    )

    # Remove duplicate columns (multiple versions mapping to same base)
    log_tpm_df = log_tpm_df.loc[:, ~log_tpm_df.columns.duplicated(keep="first")]
    log.info(f"  log1p(TPM) matrix: {log_tpm_df.shape}")

    # Align to BulkFormer gene vocabulary (main_gene_selection)
    bf_genes_present = set(log_tpm_df.columns)
    to_fill = [g for g in bf_gene_list if g not in bf_genes_present]
    log.info(f"  Genes to fill with -10 (missing): {len(to_fill)}")
    log.info(
        f"  Gene coverage: {len(bf_gene_list) - len(to_fill)}/{len(bf_gene_list)} "
        f"({100*(len(bf_gene_list) - len(to_fill))/len(bf_gene_list):.1f}%)"
    )

    # Pad missing genes with -10 (BulkFormer's mask token)
    if to_fill:
        padding = pd.DataFrame(
            np.full((log_tpm_df.shape[0], len(to_fill)), -10.0),
            columns=to_fill,
            index=log_tpm_df.index,
        )
        log_tpm_df = pd.concat([log_tpm_df, padding], axis=1)

    # Reorder to BulkFormer vocabulary order
    input_df = log_tpm_df[bf_gene_list]
    log.info(f"  Final input matrix: {input_df.shape}")

    # Compute global features BulkFormer uses
    mask_rate = len(to_fill) / len(bf_gene_list)
    log.info(f"  Gene missing rate: {mask_rate:.4f}")

    return input_df, mask_rate


# ============================================================================
# STEP 4: Run BulkFormer inference
# ============================================================================
def run_bulkformer_inference(input_df, bf_dir, mask_rate):
    """Load BulkFormer model and extract sample-level embeddings."""
    import torch
    from torch_geometric.typing import SparseTensor

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log.info(f"Running BulkFormer inference on {device}...")

    # Add BulkFormer to Python path
    bf_code_dir = str(bf_dir)
    if bf_code_dir not in sys.path:
        sys.path.insert(0, bf_code_dir)

    from utils.BulkFormer import BulkFormer
    from model.config import model_params

    # Load graph structure
    data_dir = bf_dir / "data"
    model_dir = bf_dir / "model"

    log.info("  Loading graph structure...")
    graph_raw = torch.load(
        str(data_dir / "G_tcga.pt"), map_location="cpu", weights_only=False
    )
    weights = torch.load(
        str(data_dir / "G_tcga_weight.pt"), map_location="cpu", weights_only=False
    )
    graph = SparseTensor(
        row=graph_raw[1], col=graph_raw[0], value=weights
    ).t().to(device)

    # Load ESM2 gene embeddings
    log.info("  Loading ESM2 gene embeddings...")
    gene_emb = torch.load(
        str(data_dir / "esm2_feature_concat.pt"), map_location="cpu", weights_only=False
    )

    # Configure model
    model_params["graph"] = graph
    model_params["gene_emb"] = gene_emb

    log.info(f"  Model config: dim={model_params['dim']}, "
             f"gene_length={model_params['gene_length']}, "
             f"p_repeat={model_params['p_repeat']}")

    # Build model
    model = BulkFormer(**model_params).to(device)

    # Load pretrained weights
    model_path = model_dir / MODEL_FILENAME
    log.info(f"  Loading pretrained weights: {model_path}")
    ckpt = torch.load(str(model_path), map_location=device, weights_only=False)
    new_state_dict = OrderedDict()
    for key, value in ckpt.items():
        new_key = key[7:] if key.startswith("module.") else key
        new_state_dict[new_key] = value
    model.load_state_dict(new_state_dict)
    model.eval()

    n_params = sum(p.numel() for p in model.parameters())
    log.info(f"  Model loaded: {n_params/1e6:.1f}M parameters")

    # Load interested gene indices (for filtering)
    interested_gene_idx = None
    interested_path = data_dir / "interested_gene_list.pt"
    if interested_path.exists():
        interested_gene_idx = torch.load(
            str(interested_path), map_location="cpu", weights_only=False
        )
        log.info(f"  Interested gene indices loaded: {len(interested_gene_idx)} genes")

    # Run inference in batches
    expr_array = input_df.values.astype(np.float32)
    n_samples = expr_array.shape[0]
    batch_size = 8  # conservative for 147M model on L40S
    all_embeddings = []

    log.info(f"  Running inference: {n_samples} samples, batch_size={batch_size}")

    # Set mask_prob based on actual missing gene rate
    mask_prob = mask_rate

    with torch.no_grad(), torch.amp.autocast("cuda", enabled=(device.type == "cuda")):
        for start in range(0, n_samples, batch_size):
            end = min(start + batch_size, n_samples)
            batch = torch.tensor(
                expr_array[start:end], dtype=torch.float32, device=device
            )

            # Forward pass — get gene-level embeddings
            gene_emb_out = model(batch, mask_prob=mask_prob, output_expr=False)
            # gene_emb_out shape: (batch, n_genes, hidden_dim)
            gene_emb_np = gene_emb_out.detach().cpu().numpy()

            # Filter to interested genes if available
            if interested_gene_idx is not None:
                gene_emb_np = gene_emb_np[:, interested_gene_idx, :]

            # Aggregate: mean across genes -> sample-level embedding
            sample_emb = np.mean(gene_emb_np, axis=1)  # (batch, hidden_dim)
            all_embeddings.append(sample_emb)

            del gene_emb_out, batch
            if device.type == "cuda":
                torch.cuda.empty_cache()

            if (start // batch_size + 1) % 20 == 0 or end == n_samples:
                log.info(f"  Processed {end}/{n_samples} samples")

    embeddings = np.vstack(all_embeddings)
    log.info(f"  Embeddings shape: {embeddings.shape}")
    return embeddings


# ============================================================================
# Main
# ============================================================================
def main():
    print("=" * 70)
    print("100: BulkFormer Embedding Extraction")
    print("=" * 70)
    t0 = time.time()

    # Validate inputs
    if not os.path.exists(H5_PATH):
        log.error(f"H5 not found: {H5_PATH}")
        sys.exit(1)
    if not os.path.exists(MERGED_DGE_RDS):
        log.error(f"merged_dge.rds not found: {MERGED_DGE_RDS}")
        sys.exit(1)

    import torch
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log.info(f"Device: {device}")
    if device.type == "cuda":
        log.info(f"  GPU: {torch.cuda.get_device_name(0)}")
        log.info(f"  VRAM: {torch.cuda.get_device_properties(0).total_mem / 1e9:.1f} GB")

    # -----------------------------------------------------------------------
    # Step 0: Install dependencies
    # -----------------------------------------------------------------------
    log.info("Step 0: Installing dependencies...")
    try:
        install_dependencies()
    except Exception as e:
        use_fallback(f"Dependency installation failed: {e}")
        return

    # -----------------------------------------------------------------------
    # Step 1: Setup BulkFormer (clone + download)
    # -----------------------------------------------------------------------
    log.info("Step 1: Setting up BulkFormer...")
    try:
        bf_dir = setup_bulkformer()
    except Exception as e:
        use_fallback(f"BulkFormer setup failed: {e}")
        return

    # -----------------------------------------------------------------------
    # Step 2: Load expression data
    # -----------------------------------------------------------------------
    log.info("Step 2: Loading expression data...")
    try:
        counts_df, gene_info, sample_ids = load_expression_data()
    except Exception as e:
        use_fallback(f"Expression data loading failed: {e}")
        return

    # -----------------------------------------------------------------------
    # Step 3: Preprocess for BulkFormer
    # -----------------------------------------------------------------------
    log.info("Step 3: Preprocessing for BulkFormer...")
    try:
        input_df, mask_rate = preprocess_for_bulkformer(counts_df, gene_info, bf_dir)
    except Exception as e:
        use_fallback(f"Preprocessing failed: {e}")
        return

    # -----------------------------------------------------------------------
    # Step 4: Run BulkFormer inference
    # -----------------------------------------------------------------------
    log.info("Step 4: Running BulkFormer inference...")
    try:
        embeddings = run_bulkformer_inference(input_df, bf_dir, mask_rate)
    except Exception as e:
        use_fallback(f"BulkFormer inference failed: {e}")
        return

    # -----------------------------------------------------------------------
    # Step 5: Save embeddings
    # -----------------------------------------------------------------------
    log.info("Step 5: Saving embeddings...")
    embed_cols = [f"bf{i}" for i in range(embeddings.shape[1])]
    embed_df = pd.DataFrame(embeddings, columns=embed_cols)
    embed_df.insert(0, "sample_id", sample_ids)
    embed_df.to_csv(OUTPUT_EMBEDDINGS, index=False)

    log.info(f"  BulkFormer embeddings saved: {OUTPUT_EMBEDDINGS}")
    log.info(f"  Shape: {embed_df.shape[0]} samples x {embeddings.shape[1]} dimensions")

    # Summary stats
    log.info(f"  Embedding stats: mean={embeddings.mean():.4f}, "
             f"std={embeddings.std():.4f}, "
             f"min={embeddings.min():.4f}, max={embeddings.max():.4f}")

    # Clean up temp files
    for tmp in ["_bulkformer_raw_counts.csv", "_bulkformer_gene_info.csv",
                "_bulkformer_export.R"]:
        tmp_path = os.path.join(OUTDIR, tmp)
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    elapsed = time.time() - t0
    log.info(f"\nDone in {elapsed/60:.1f} minutes.")
    log.info(f"Output: {OUTPUT_EMBEDDINGS}")


if __name__ == "__main__":
    main()
