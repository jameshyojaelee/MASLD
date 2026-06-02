#!/usr/bin/env python3
"""
03_chromvar_motifs.py — chromVAR TF motif enrichment analysis (Module 2c)

Computes per-cell TF deviation scores from SnapATAC2 peak matrix using the
chromVAR algorithm, then tests for differential TF activity between MASLD
conditions per cell type.

Pipeline position:
  01_cellranger_arc.sh → 02_snapatac2_processing.py → **03_chromvar_motifs.py**

Inputs:
  - SnapATAC2 processed AnnData with peak matrix and cell type labels

Outputs:
  - results/chromvar/chromvar_tf_activity.csv    Per-celltype TF stats
  - results/chromvar/chromvar_deviations.h5ad     Per-cell TF deviation scores
  - results/chromvar/tf_activity_heatmap.pdf      Top variable TF heatmap

Environment:
  micromamba activate atac_env

Usage:
  cd Analysis/ATAC/Human_Multiome
  python scripts/03_chromvar_motifs.py \\
      --input results/snapatac2/snapatac2_processed.h5ad \\
      --output-dir results/chromvar \\
      --genome hg38

Author: MASLD-Atlas pipeline
"""

from __future__ import annotations

import argparse
import gc
import logging
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

# Donor-level pseudobulk helpers (fixes per-cell pseudoreplication; F024/F028).
# utils_pseudobulk.py lives in this same scripts/ dir.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils_pseudobulk import aggregate_cells_to_donors, donor_groupwise_test

warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Donor-level differential-TF configuration (pseudoreplication fix)
# ---------------------------------------------------------------------------
# The differential test runs at DONOR resolution, not per-cell. Donor labels
# come from the curated bridge; the per-donor chromVAR means (if already
# exported by 03b) are preferred over re-aggregating the per-cell h5ad.
_ATAC_DIR = Path(__file__).resolve().parent.parent  # Analysis/ATAC/Human_Multiome
DONOR_META_TSV = _ATAC_DIR / "metadata" / "donor_metadata_curated.tsv"
# Existing per-donor chromVAR means, in preference order (first that exists wins)
PER_DONOR_CANDIDATES = [
    _ATAC_DIR / "results" / "chromvar_v2" / "chromvar_per_donor.tsv.gz",
    _ATAC_DIR / "results" / "chromvar" / "chromvar_per_donor.tsv.gz",
]
# Condition -> donor group. MASL/MASH = disease (reported as "MASLD"); NORMAL = healthy.
DISEASE_CONDITIONS = ("MASL", "MASH")
HEALTHY_CONDITIONS = ("NORMAL",)
GROUP_DISEASE = "MASLD"   # output naming preserved from the legacy per-cell table
GROUP_HEALTHY = "Normal"
MIN_DONORS_PER_GROUP = 2  # donor-level floor for a defined Mann-Whitney U


# ---------------------------------------------------------------------------
# Dependency checks
# ---------------------------------------------------------------------------
def _check_imports():
    """Verify that required packages are importable and report versions."""
    required = {
        "anndata": "anndata",
        "scanpy": "scanpy",
    }
    optional = {
        "snapatac2": "snapatac2",
        "gimmemotifs": "gimmemotifs",
        "pyjaspar": "pyjaspar",
    }
    for label, mod_name in required.items():
        try:
            mod = __import__(mod_name)
            log.info(f"  {label}: {getattr(mod, '__version__', 'found')}")
        except ImportError:
            log.error(f"Required package '{label}' not installed. "
                      "Run: micromamba activate atac_env")
            sys.exit(1)

    available = {}
    for label, mod_name in optional.items():
        try:
            mod = __import__(mod_name)
            log.info(f"  {label}: {getattr(mod, '__version__', 'found')}")
            available[label] = True
        except ImportError:
            log.info(f"  {label}: not installed (optional)")
            available[label] = False
    return available


# ---------------------------------------------------------------------------
# Motif scanning helpers
# ---------------------------------------------------------------------------
def _get_peak_sequences(adata, genome: str = "hg38"):
    """Extract genomic sequences for peaks from a reference genome.

    Uses pysam or pyfaidx to pull sequences from a local FASTA.  Falls back
    to a simple coordinate representation if no FASTA is available.

    Returns:
        dict mapping peak name (chr:start-end) to DNA sequence string.
    """
    ref_paths = {
        "hg38": "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/"
                "refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa",
    }
    fasta_path = ref_paths.get(genome)

    peak_names = adata.var_names.tolist()
    sequences: dict[str, str] = {}

    if fasta_path and os.path.exists(fasta_path):
        try:
            import pysam
            fa = pysam.FastaFile(fasta_path)
            for pname in peak_names:
                try:
                    chrom, rest = pname.split(":")
                    start, end = rest.split("-")
                    seq = fa.fetch(chrom, int(start), int(end)).upper()
                    sequences[pname] = seq
                except Exception:
                    sequences[pname] = ""
            fa.close()
            log.info(f"Extracted sequences for {len(sequences):,} peaks from {fasta_path}")
            return sequences
        except ImportError:
            log.warning("pysam not installed; trying pyfaidx")
        except Exception as exc:
            log.warning(f"pysam FASTA extraction failed: {exc}")

        try:
            from pyfaidx import Fasta
            fa = Fasta(fasta_path)
            for pname in peak_names:
                try:
                    chrom, rest = pname.split(":")
                    start, end = rest.split("-")
                    sequences[pname] = str(fa[chrom][int(start):int(end)]).upper()
                except Exception:
                    sequences[pname] = ""
            log.info(f"Extracted sequences for {len(sequences):,} peaks via pyfaidx")
            return sequences
        except ImportError:
            log.warning("pyfaidx not installed")
        except Exception as exc:
            log.warning(f"pyfaidx extraction failed: {exc}")

    log.warning("No genome FASTA accessible; motif scanning will use "
                "SnapATAC2 built-in methods or gimmemotifs defaults")
    return None


def _scan_motifs_snapatac2(adata, genome: str):
    """SnapATAC2 motif_enrichment does set-level enrichment, not per-cell
    chromVAR deviations. Skip this and use pyjaspar/gimmemotifs instead."""
    log.info("SnapATAC2 motif_enrichment is region-level (not per-cell); skipping")
    return None


def _scan_motifs_gimmemotifs(adata, peak_sequences: dict | None):
    """Scan peaks for TF binding motifs using gimmemotifs."""
    from gimmemotifs.motif import read_motifs
    from gimmemotifs.scanner import Scanner

    log.info("Loading JASPAR 2024 vertebrate motifs via gimmemotifs...")

    # gimmemotifs ships with JASPAR databases; use vertebrate core
    try:
        from gimmemotifs.motif import default_motifs
        motifs = default_motifs()
        log.info(f"Loaded {len(motifs)} default motifs from gimmemotifs")
    except Exception:
        # Fallback: load JASPAR PFM file if gimmemotifs default fails
        jaspar_pfm = os.path.expanduser(
            "~/.local/share/gimmemotifs/motif_databases/JASPAR2024_vertebrates.pfm"
        )
        if os.path.exists(jaspar_pfm):
            motifs = read_motifs(jaspar_pfm)
            log.info(f"Loaded {len(motifs)} JASPAR 2024 vertebrate motifs")
        else:
            log.error("Cannot locate JASPAR motif database. Install gimmemotifs "
                      "databases: `gimme motif2factors --database JASPAR2024`")
            return None, None

    if peak_sequences is None:
        log.error("No peak sequences available for motif scanning")
        return None, None

    # Write peak sequences to temp FASTA for scanner
    import tempfile
    tmp_fasta = tempfile.NamedTemporaryFile(
        suffix=".fa", delete=False, mode="w"
    )
    for pname, seq in peak_sequences.items():
        if seq:
            tmp_fasta.write(f">{pname}\n{seq}\n")
    tmp_fasta.close()

    log.info(f"Scanning {len(peak_sequences):,} peaks against {len(motifs)} motifs...")
    scanner = Scanner()
    scanner.set_motifs(motifs)
    scanner.set_genome(tmp_fasta.name)

    peak_list = [p for p in adata.var_names if peak_sequences.get(p, "")]
    n_peaks = len(peak_list)
    n_motifs = len(motifs)

    # Build binary peak x motif match matrix
    match_matrix = np.zeros((n_peaks, n_motifs), dtype=np.float32)
    peak_to_idx = {p: i for i, p in enumerate(peak_list)}

    results = scanner.scan(tmp_fasta.name, nreport=1, score_threshold=5.0)
    for motif_idx, motif_matches in enumerate(results):
        for match in motif_matches:
            seq_name = match.seq_id if hasattr(match, "seq_id") else str(match[0])
            if seq_name in peak_to_idx:
                match_matrix[peak_to_idx[seq_name], motif_idx] = 1.0

    os.unlink(tmp_fasta.name)

    motif_names = [m.id if hasattr(m, "id") else str(m) for m in motifs]
    log.info(f"Motif scanning complete: {int(match_matrix.sum()):,} peak-motif matches")

    # Filter motifs with too few peak matches — rare motifs produce unstable
    # deviations because expected counts approach zero, causing division overflow
    MIN_PEAKS_PER_MOTIF = 50
    peaks_per_motif = match_matrix.sum(axis=0).astype(int)
    keep_mask = peaks_per_motif >= MIN_PEAKS_PER_MOTIF
    n_dropped = int((~keep_mask).sum())
    if n_dropped > 0:
        match_matrix = match_matrix[:, keep_mask]
        motif_names = [m for m, k in zip(motif_names, keep_mask) if k]
        log.info(f"  Filtered {n_dropped} motifs with <{MIN_PEAKS_PER_MOTIF} peak matches "
                 f"({len(motif_names)} motifs retained)")

    return match_matrix, motif_names


# ---------------------------------------------------------------------------
# chromVAR deviation scoring
# ---------------------------------------------------------------------------
def _compute_gc_content(sequences: dict[str, str]) -> dict[str, float]:
    """Compute GC fraction for each peak."""
    gc = {}
    for name, seq in sequences.items():
        if not seq:
            gc[name] = 0.5
            continue
        s = seq.upper()
        n = len(s)
        gc[name] = (s.count("G") + s.count("C")) / n if n > 0 else 0.5
    return gc


def _compute_chromvar_deviations(
    adata,
    match_matrix: np.ndarray,
    motif_names: list[str],
    peak_sequences: dict[str, str] | None,
    n_bg: int = 50,
    seed: int = 42,
):
    """Compute chromVAR-style deviation z-scores.

    Implements the core chromVAR algorithm:
    1. For each peak, compute expected accessibility from mean accessibility
    2. Match background peaks by GC content and mean accessibility
    3. For each motif, compute deviation = (observed - expected) / expected_sd

    Args:
        adata: AnnData with peak count matrix in .X
        match_matrix: (n_peaks, n_motifs) binary motif match matrix
        motif_names: list of motif IDs
        peak_sequences: dict of peak name -> DNA sequence (for GC matching)
        n_bg: number of background iterations
        seed: random seed

    Returns:
        AnnData with cells x motifs deviation scores
    """
    import anndata as ad

    rng = np.random.RandomState(seed)
    n_cells, n_peaks = adata.shape
    n_motifs = match_matrix.shape[1]

    log.info(f"Computing chromVAR deviations: {n_cells:,} cells x {n_peaks:,} peaks x {n_motifs} motifs")

    # Get count matrix as sparse
    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.astype(np.float64)

    # Per-cell total counts and per-peak mean accessibility
    cell_totals = np.array(X.sum(axis=1)).ravel()  # (n_cells,)
    peak_means = np.array(X.mean(axis=0)).ravel()   # (n_peaks,)
    peak_means_norm = peak_means / peak_means.sum()  # fraction of total signal per peak

    # GC content for background matching
    if peak_sequences is not None:
        gc_vals = _compute_gc_content(peak_sequences)
        peak_gc = np.array([gc_vals.get(p, 0.5) for p in adata.var_names[:n_peaks]])
    else:
        peak_gc = np.full(n_peaks, 0.5)

    # Bin peaks by GC content + mean accessibility for background matching
    n_gc_bins = 10
    n_acc_bins = 10
    gc_bins = np.digitize(peak_gc, np.linspace(0, 1, n_gc_bins + 1)) - 1
    gc_bins = np.clip(gc_bins, 0, n_gc_bins - 1)

    log_acc = np.log1p(peak_means)
    acc_quantiles = np.linspace(0, 100, n_acc_bins + 1)
    acc_edges = np.percentile(log_acc, acc_quantiles)
    acc_bins = np.digitize(log_acc, acc_edges) - 1
    acc_bins = np.clip(acc_bins, 0, n_acc_bins - 1)

    combined_bins = gc_bins * n_acc_bins + acc_bins

    # Build background peak sets: for each peak, sample peers from same bin
    bin_members: dict[int, list[int]] = {}
    for i in range(n_peaks):
        b = combined_bins[i]
        bin_members.setdefault(b, []).append(i)

    # =========================================================================
    # Vectorized deviation computation using sparse matrix multiplication
    # Key insight: observed = X @ match_matrix computes all motif sums at once
    # =========================================================================
    log.info("Computing deviation z-scores (vectorized)...")

    # Convert match_matrix to sparse for efficient multiplication
    match_sparse = sp.csc_matrix(match_matrix.astype(np.float64))

    # Foreground: observed counts per cell per motif (single sparse matmul)
    log.info("  Computing observed counts (X @ match_matrix)...")
    observed_all = X @ match_sparse  # (n_cells, n_motifs) sparse
    observed_all = np.asarray(observed_all.todense())

    # Expected: cell_totals_i * sum(peak_means_norm_j * match_j_m) for each motif m
    motif_fracs = peak_means_norm @ match_matrix  # (n_motifs,)
    expected_all = np.outer(cell_totals, motif_fracs)  # (n_cells, n_motifs)

    # Foreground deviation = (observed - expected) / max(expected, eps)
    # Use eps=1e-4 (not 1e-10) to prevent overflow from rare motifs
    # with near-zero expected counts
    DEVIATION_EPS = 1e-4
    fg_deviation_all = (observed_all - expected_all) / np.maximum(expected_all, DEVIATION_EPS)

    # Background: for each background iteration, permute peaks within GC/acc bins
    # and compute background deviations
    log.info(f"  Computing {n_bg} background iterations...")
    bg_mean_all = np.zeros((n_cells, n_motifs), dtype=np.float64)
    bg_sq_all = np.zeros((n_cells, n_motifs), dtype=np.float64)

    # Pre-build permutation lookup: for each peak, candidates in same bin
    bin_candidates = np.empty(n_peaks, dtype=object)
    for i in range(n_peaks):
        bin_candidates[i] = bin_members[combined_bins[i]]

    for bg_i in range(n_bg):
        if bg_i % 10 == 0:
            log.info(f"    Background iteration {bg_i + 1}/{n_bg}")

        # Permute peaks within bins (vectorized)
        perm = np.array([rng.choice(bin_candidates[i]) for i in range(n_peaks)])

        # Background match matrix: permuted rows of original
        bg_match = match_matrix[perm, :]  # (n_peaks, n_motifs)
        bg_match_sparse = sp.csc_matrix(bg_match.astype(np.float64))

        # Background observed
        bg_observed = X @ bg_match_sparse
        bg_observed = np.asarray(bg_observed.todense())

        # Background expected
        bg_fracs = peak_means_norm @ bg_match  # (n_motifs,)
        bg_expected = np.outer(cell_totals, bg_fracs)

        # Background deviation (same eps as foreground)
        bg_dev = (bg_observed - bg_expected) / np.maximum(bg_expected, DEVIATION_EPS)

        # Online mean/variance (Welford's algorithm)
        bg_mean_all += bg_dev
        bg_sq_all += bg_dev ** 2

    bg_mean_all /= n_bg
    bg_var_all = bg_sq_all / n_bg - bg_mean_all ** 2
    bg_sd_all = np.sqrt(np.maximum(bg_var_all, 1e-20))

    # Z-score — clip extreme values and ensure standard float32
    deviations = (fg_deviation_all - bg_mean_all) / np.maximum(bg_sd_all, 1e-6)
    # Clip to [-50, 50] to prevent downstream overflow in mean computation
    deviations = np.clip(deviations, -50, 50).astype(np.float32)
    log.info("  Deviation z-scores computed (clipped to [-50, 50])")

    # Build output AnnData
    dev_adata = ad.AnnData(
        X=deviations,
        obs=adata.obs.copy(),
        var=pd.DataFrame({"motif_name": motif_names}, index=motif_names),
    )
    log.info(f"Deviation matrix computed: {dev_adata.shape}")
    return dev_adata


# ---------------------------------------------------------------------------
# Differential TF activity
# ---------------------------------------------------------------------------
def _load_donor_to_group() -> dict[str, str]:
    """Map donor_id (D##) -> donor group {GROUP_DISEASE, GROUP_HEALTHY}.

    Reads the curated bridge donor_metadata_curated.tsv. The per-donor chromVAR
    TSV and the deviations h5ad both key donors by donor_id (D##), which matches
    the metadata `donor_id` column (NOT `donor_id_atac` = MM_##; see
    09_chromvar_propagation.py lines 247-249).
    """
    if not DONOR_META_TSV.exists():
        log.error(f"Donor metadata not found: {DONOR_META_TSV}")
        return {}
    meta = pd.read_csv(DONOR_META_TSV, sep="\t")
    if "donor_id" not in meta.columns or "condition" not in meta.columns:
        log.error(f"Donor metadata missing donor_id/condition columns. "
                  f"Have: {list(meta.columns)}")
        return {}
    d2g: dict[str, str] = {}
    for _, r in meta.iterrows():
        cond = str(r["condition"])
        if cond in DISEASE_CONDITIONS:
            d2g[str(r["donor_id"])] = GROUP_DISEASE
        elif cond in HEALTHY_CONDITIONS:
            d2g[str(r["donor_id"])] = GROUP_HEALTHY
    n_dis = sum(v == GROUP_DISEASE for v in d2g.values())
    n_hea = sum(v == GROUP_HEALTHY for v in d2g.values())
    log.info(f"  Donor groups from {DONOR_META_TSV.name}: "
             f"{n_dis} {GROUP_DISEASE} (MASL+MASH), {n_hea} {GROUP_HEALTHY} (NORMAL)")
    return d2g


def _rename_donor_test_columns(res: pd.DataFrame, ct: str) -> pd.DataFrame:
    """Rename donor_groupwise_test output to the legacy per-cell table schema.

    Legacy columns (preserved): tf_name, cell_type, mean_deviation_masld,
    mean_deviation_normal, logFC_deviation, pvalue, padj, n_masld, n_normal.
    Donor-level additions: n_donors, n_donors_masld, n_donors_normal, u_stat.
    The helper returns mean_<GROUP_DISEASE>/mean_<GROUP_HEALTHY>; map those to
    the *_masld / *_normal names and surface mean_diff as logFC_deviation.
    """
    out = res.rename(
        columns={
            "feature": "tf_name",
            f"mean_{GROUP_DISEASE}": "mean_deviation_masld",
            f"mean_{GROUP_HEALTHY}": "mean_deviation_normal",
            "mean_diff": "logFC_deviation",
            "n_donors_a": "n_donors_masld",
            "n_donors_b": "n_donors_normal",
        }
    )
    out["cell_type"] = ct
    out["n_donors"] = out["n_donors_masld"] + out["n_donors_normal"]
    # Legacy n_masld/n_normal columns now carry DONOR counts (the test unit)
    out["n_masld"] = out["n_donors_masld"]
    out["n_normal"] = out["n_donors_normal"]
    for c in ("mean_deviation_masld", "mean_deviation_normal", "logFC_deviation"):
        out[c] = out[c].round(6)
    return out


def _differential_tf_activity_from_per_donor(per_donor_path: Path,
                                             donor_to_group: dict[str, str]):
    """Donor-level differential TF activity from an existing per-donor TSV.

    The TSV (03b output) is long-form: donor_id_atac | cell_type | TF |
    mean_deviation | ... where donor_id_atac holds D## ids. For each cell type
    we pivot to a (donor x TF) matrix and run donor_groupwise_test.
    """
    log.info(f"  Source: per-donor means {per_donor_path}")
    pd_df = pd.read_csv(per_donor_path, sep="\t")
    donor_col = "donor_id_atac" if "donor_id_atac" in pd_df.columns else (
        "donor_id" if "donor_id" in pd_df.columns else None
    )
    tf_col = "TF" if "TF" in pd_df.columns else ("motif_name" if "motif_name" in pd_df.columns else None)
    if donor_col is None or tf_col is None or "cell_type" not in pd_df.columns \
            or "mean_deviation" not in pd_df.columns:
        log.error(f"  Per-donor TSV missing required columns. Have: {list(pd_df.columns)}")
        return pd.DataFrame()

    pd_df[donor_col] = pd_df[donor_col].astype(str)
    frames = []
    for ct, g in pd_df.groupby("cell_type", observed=True):
        # (donor x TF) matrix of mean deviations
        mat = g.pivot_table(index=donor_col, columns=tf_col,
                            values="mean_deviation", aggfunc="mean")
        donors = mat.index.astype(str).tolist()
        res = donor_groupwise_test(
            mat.values, donors, donor_to_group,
            GROUP_DISEASE, GROUP_HEALTHY,
            feature_names=mat.columns.tolist(),
            min_per_group=MIN_DONORS_PER_GROUP,
        )
        n_a, n_b = int(res["n_donors_a"].iloc[0]), int(res["n_donors_b"].iloc[0])
        if n_a < MIN_DONORS_PER_GROUP or n_b < MIN_DONORS_PER_GROUP:
            log.info(f"  Skipping {ct}: {GROUP_DISEASE}={n_a}, {GROUP_HEALTHY}={n_b} "
                     f"donors (need >={MIN_DONORS_PER_GROUP})")
            continue
        frames.append(_rename_donor_test_columns(res, ct))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def _differential_tf_activity_from_h5ad(dev_adata, donor_to_group: dict[str, str],
                                        cell_type_col: str, donor_col: str = "donor_id"):
    """Donor-level differential TF activity by aggregating the per-cell h5ad.

    Fallback path when no per-donor TSV exists. Aggregates the EXISTING per-cell
    deviation z-scores to per-donor means (does NOT recompute chromVAR) and runs
    donor_groupwise_test per cell type.
    """
    log.info("  Source: aggregating per-cell deviations h5ad to donor means")
    if donor_col not in dev_adata.obs.columns:
        log.error(f"  Column '{donor_col}' not in obs; cannot aggregate to donors. "
                  f"Available: {list(dev_adata.obs.columns)}")
        return pd.DataFrame()
    if cell_type_col not in dev_adata.obs.columns:
        log.error(f"  Column '{cell_type_col}' not in obs. "
                  f"Available: {list(dev_adata.obs.columns)}")
        return pd.DataFrame()

    motif_names = dev_adata.var_names.tolist()
    donors_obs = dev_adata.obs[donor_col].astype(str).values
    ct_obs = dev_adata.obs[cell_type_col].astype(str).values
    X = dev_adata.X
    frames = []
    for ct in pd.unique(ct_obs):
        ct_mask = ct_obs == ct
        X_ct = X[ct_mask, :]
        donors_ct = donors_obs[ct_mask]
        # Aggregate the existing per-cell z-scores to per-donor means
        donors, M, _ = aggregate_cells_to_donors(X_ct, donors_ct, agg="mean")
        if not donors:
            continue
        res = donor_groupwise_test(
            M, donors, donor_to_group,
            GROUP_DISEASE, GROUP_HEALTHY,
            feature_names=motif_names,
            min_per_group=MIN_DONORS_PER_GROUP,
        )
        n_a, n_b = int(res["n_donors_a"].iloc[0]), int(res["n_donors_b"].iloc[0])
        if n_a < MIN_DONORS_PER_GROUP or n_b < MIN_DONORS_PER_GROUP:
            log.info(f"  Skipping {ct}: {GROUP_DISEASE}={n_a}, {GROUP_HEALTHY}={n_b} "
                     f"donors (need >={MIN_DONORS_PER_GROUP})")
            continue
        frames.append(_rename_donor_test_columns(res, ct))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def _differential_tf_activity(
    dev_adata,
    cell_type_col: str = "cell_type",
    condition_col: str = "condition",
):
    """Test differential TF deviation scores between MASLD and Normal per cell type.

    DONOR-LEVEL Mann-Whitney U with BH correction (the donor, not the cell, is
    the experimental unit; fixes pseudoreplication F024/F028). Prefers existing
    per-donor chromVAR means (03b export); otherwise aggregates the per-cell
    deviation h5ad to donor means. chromVAR deviations are never recomputed here.

    Returns:
        DataFrame with columns: tf_name, cell_type, mean_deviation_masld,
        mean_deviation_normal, logFC_deviation, pvalue, padj, n_masld, n_normal
        (now donor counts), plus n_donors, n_donors_masld, n_donors_normal, u_stat.
    """
    from statsmodels.stats.multitest import multipletests

    log.info("Computing DONOR-LEVEL differential TF activity per cell type...")

    donor_to_group = _load_donor_to_group()
    if not donor_to_group:
        log.error("No donor->group mapping available; cannot run donor-level test")
        return pd.DataFrame()

    # Prefer an existing per-donor means table; fall back to aggregating the h5ad.
    per_donor_path = next((p for p in PER_DONOR_CANDIDATES if p.exists()), None)
    if per_donor_path is not None:
        df = _differential_tf_activity_from_per_donor(per_donor_path, donor_to_group)
    else:
        log.info("  No per-donor TSV found; aggregating per-cell deviations h5ad")
        df = _differential_tf_activity_from_h5ad(
            dev_adata, donor_to_group, cell_type_col=cell_type_col
        )

    if df.empty:
        log.warning("No differential TF activity results produced")
        return df

    # BH correction across all (TF x cell_type) donor-level tests (NaN-safe)
    pvals = df["pvalue"].to_numpy(dtype=float)
    padj = np.full(pvals.shape, np.nan)
    ok = ~np.isnan(pvals)
    if ok.any():
        _, padj_ok, _, _ = multipletests(pvals[ok], method="fdr_bh")
        padj[ok] = padj_ok
    df["padj"] = padj

    # Preserve legacy column order, then append donor-level diagnostics
    legacy_cols = [
        "tf_name", "cell_type", "mean_deviation_masld", "mean_deviation_normal",
        "pvalue", "logFC_deviation", "n_masld", "n_normal", "padj",
    ]
    extra_cols = ["n_donors", "n_donors_masld", "n_donors_normal", "u_stat"]
    df = df[[c for c in legacy_cols if c in df.columns]
            + [c for c in extra_cols if c in df.columns]]

    df = df.sort_values("padj", na_position="last").reset_index(drop=True)

    n_sig = int((df["padj"] < 0.05).sum())
    n_tested = int(df["pvalue"].notna().sum())
    log.info(f"Differential TF activity (donor-level): {len(df):,} rows, "
             f"{n_tested:,} with a defined p-value, {n_sig:,} significant (padj<0.05)")
    return df


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------
def _plot_tf_heatmap(dev_adata, output_path: str, n_top: int = 50,
                     cell_type_col: str = "cell_type"):
    """Plot heatmap of top variable TFs across cell types."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    log.info(f"Generating TF activity heatmap (top {n_top} variable TFs)...")

    if cell_type_col not in dev_adata.obs.columns:
        log.warning(f"Column '{cell_type_col}' not in obs; skipping heatmap")
        return

    # Mean deviation per cell type
    cell_types = dev_adata.obs[cell_type_col].dropna().unique()
    motif_names = dev_adata.var_names.tolist()

    mean_dev = np.zeros((len(cell_types), len(motif_names)))
    for i, ct in enumerate(cell_types):
        mask = dev_adata.obs[cell_type_col] == ct
        if mask.sum() > 0:
            mean_dev[i, :] = np.mean(dev_adata.X[mask.values, :], axis=0)

    mean_df = pd.DataFrame(mean_dev, index=cell_types, columns=motif_names)

    # Select top variable motifs
    var_per_motif = mean_df.var(axis=0)
    top_motifs = var_per_motif.nlargest(n_top).index.tolist()
    plot_df = mean_df[top_motifs]

    # Plot
    fig_height = max(6, len(cell_types) * 0.4 + 2)
    fig_width = max(10, n_top * 0.25 + 3)

    g = sns.clustermap(
        plot_df.T,
        cmap="RdBu_r",
        center=0,
        figsize=(fig_width, fig_height),
        xticklabels=True,
        yticklabels=True,
        dendrogram_ratio=(0.1, 0.15),
        cbar_kws={"label": "Mean deviation z-score"},
        row_cluster=True,
        col_cluster=True,
    )
    g.ax_heatmap.set_xlabel("Cell type")
    g.ax_heatmap.set_ylabel("TF motif")
    g.fig.suptitle("TF Motif Activity Across Cell Types (chromVAR)", y=1.02)
    plt.tight_layout()
    g.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info(f"Heatmap saved: {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="chromVAR TF motif enrichment analysis (Module 2c)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="SnapATAC2 processed AnnData (.h5ad) with peak matrix",
    )
    parser.add_argument(
        "--output-dir", "-o",
        default="results/chromvar",
        help="Output directory for results",
    )
    parser.add_argument(
        "--genome", "-g",
        default="hg38",
        help="Genome assembly for sequence extraction",
    )
    parser.add_argument(
        "--cell-type-col",
        default="cell_type",
        help="Column in obs with cell type labels",
    )
    parser.add_argument(
        "--condition-col",
        default="condition",
        help="Column in obs with disease condition labels",
    )
    parser.add_argument(
        "--n-background",
        type=int,
        default=50,
        help="Number of background iterations for deviation scoring",
    )
    parser.add_argument(
        "--n-top-heatmap",
        type=int,
        default=50,
        help="Number of top variable TFs to show in heatmap",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for background sampling",
    )
    args = parser.parse_args()

    t_start = time.time()
    log.info("=" * 70)
    log.info("chromVAR TF Motif Enrichment Analysis (Module 2c)")
    log.info("=" * 70)

    # Check dependencies
    log.info("Checking dependencies...")
    available = _check_imports()

    import anndata as ad

    # Create output directory
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------------
    # Step 1: Load SnapATAC2 AnnData
    # -----------------------------------------------------------------------
    log.info(f"Loading SnapATAC2 AnnData: {args.input}")
    if not os.path.exists(args.input):
        log.error(f"Input file not found: {args.input}")
        sys.exit(1)

    adata = ad.read_h5ad(args.input)
    log.info(f"  Shape: {adata.shape[0]:,} cells x {adata.shape[1]:,} peaks")
    log.info(f"  obs columns: {list(adata.obs.columns)}")

    if args.cell_type_col in adata.obs.columns:
        cts = adata.obs[args.cell_type_col].value_counts()
        log.info(f"  Cell types ({len(cts)}):")
        for ct, n in cts.items():
            log.info(f"    {ct}: {n:,}")

    if args.condition_col in adata.obs.columns:
        conds = adata.obs[args.condition_col].value_counts()
        log.info(f"  Conditions: {dict(conds)}")

    # -----------------------------------------------------------------------
    # Step 1b: Filter to accessible tiles/peaks for tractable motif scanning
    # -----------------------------------------------------------------------
    log.info("-" * 50)
    log.info("Step 1b: Filtering to top accessible regions")

    max_features = 50000  # Cap at 50K most-accessible tiles for tractable chromVAR
    if adata.shape[1] > max_features:
        log.info(f"  Input has {adata.shape[1]:,} features; filtering to top {max_features:,}")
        X = adata.X
        if sp.issparse(X):
            feature_sums = np.array(X.sum(axis=0)).ravel()
        else:
            feature_sums = np.asarray(X.sum(axis=0)).ravel()
        top_idx = np.argsort(feature_sums)[-max_features:]
        top_idx = np.sort(top_idx)
        adata = adata[:, top_idx].copy()
        log.info(f"  Filtered: {adata.shape[0]:,} cells x {adata.shape[1]:,} features")
    else:
        log.info(f"  Features ({adata.shape[1]:,}) below threshold; using all")

    gc.collect()

    # -----------------------------------------------------------------------
    # Step 2: Motif scanning
    # -----------------------------------------------------------------------
    log.info("-" * 50)
    log.info("Step 2: TF motif scanning")
    match_matrix = None
    motif_names = None

    # Strategy 1: gimmemotifs (if installed)
    if match_matrix is None and available.get("gimmemotifs", False):
        log.info("Using gimmemotifs for motif scanning...")
        peak_sequences = _get_peak_sequences(adata, args.genome)
        match_matrix, motif_names = _scan_motifs_gimmemotifs(adata, peak_sequences)

    # Strategy 2: pyjaspar + Bio.motifs (vectorized PWM scoring)
    if match_matrix is None and available.get("pyjaspar", False):
        log.info("Using pyjaspar + Bio.motifs for motif scanning...")
        peak_sequences = _get_peak_sequences(adata, args.genome)
        if peak_sequences is not None:
            try:
                from pyjaspar import jaspardb
                from Bio.motifs import matrix as bio_matrix

                jdb = jaspardb()
                motifs = jdb.fetch_motifs(
                    collection="CORE",
                    tax_group=["vertebrates"],
                )
                log.info(f"Loaded {len(motifs)} JASPAR vertebrate motifs")

                peak_list = list(adata.var_names)
                n_peaks = len(peak_list)
                n_motifs_j = len(motifs)
                match_mat = np.zeros((n_peaks, n_motifs_j), dtype=np.float32)
                motif_names = []

                # Build Bio.motifs PSSM objects for efficient scoring
                pssms = []
                for motif in motifs:
                    motif_names.append(motif.name)
                    try:
                        pwm = motif.counts.normalize(pseudocounts=0.5)
                        pssm = pwm.log_odds()
                        pssms.append(pssm)
                    except Exception:
                        pssms.append(None)

                # Scan peaks in batches
                report_interval = max(1, n_peaks // 20)
                for p_idx, pname in enumerate(peak_list):
                    if p_idx % report_interval == 0:
                        log.info(f"  Scanning peak {p_idx:,}/{n_peaks:,} "
                                 f"({100 * p_idx / n_peaks:.0f}%)")
                    seq = peak_sequences.get(pname, "")
                    if len(seq) < 6:
                        continue
                    for m_idx, pssm in enumerate(pssms):
                        if pssm is None:
                            continue
                        motif_len = pssm.length
                        if len(seq) < motif_len:
                            continue
                        try:
                            scores = pssm.calculate(seq)
                            if isinstance(scores, (list, np.ndarray)):
                                max_score = max(scores) if len(scores) > 0 else -999
                            else:
                                max_score = float(scores)
                            # Stricter threshold: ~80% of max possible score
                            # Targets ~5-15% match density (vs 44% with 0.5)
                            max_possible = sum(max(pssm[nt][i] for nt in "ACGT")
                                               for i in range(motif_len))
                            if max_score > max_possible * 0.8:
                                match_mat[p_idx, m_idx] = 1.0
                        except Exception:
                            continue

                match_matrix = match_mat
                n_matches = int(match_mat.sum())
                log.info(f"pyjaspar scanning complete: {n_matches:,} peak-motif matches "
                         f"({n_matches / (n_peaks * n_motifs_j) * 100:.2f}% density)")
            except Exception as exc:
                log.warning(f"pyjaspar motif scanning failed: {exc}")

    if match_matrix is None:
        log.error(
            "No motif scanning backend available. Install one of:\n"
            "  - gimmemotifs >= 0.18\n"
            "  - pyjaspar >= 3.0 + biopython\n"
            "Aborting."
        )
        sys.exit(1)

    # -----------------------------------------------------------------------
    # Step 3: Compute chromVAR deviation scores
    # -----------------------------------------------------------------------
    log.info("-" * 50)
    log.info("Step 3: Computing chromVAR deviation scores")

    peak_sequences = _get_peak_sequences(adata, args.genome)
    dev_adata = _compute_chromvar_deviations(
        adata,
        match_matrix,
        motif_names,
        peak_sequences,
        n_bg=args.n_background,
            seed=args.seed,
        )

    gc.collect()

    # -----------------------------------------------------------------------
    # Step 4: Save per-cell deviation AnnData
    # -----------------------------------------------------------------------
    log.info("-" * 50)
    log.info("Step 4: Saving per-cell TF deviations")

    dev_h5ad_path = out_dir / "chromvar_deviations.h5ad"
    # Ensure X is standard numpy float32 (not GPU TF32) before saving
    if sp.issparse(dev_adata.X):
        dev_adata.X = np.asarray(dev_adata.X.todense(), dtype=np.float32)
    elif not isinstance(dev_adata.X, np.ndarray) or dev_adata.X.dtype != np.float32:
        dev_adata.X = np.asarray(dev_adata.X, dtype=np.float32)
    dev_adata.write_h5ad(dev_h5ad_path)
    log.info(f"  Saved: {dev_h5ad_path} ({dev_adata.shape})")

    # -----------------------------------------------------------------------
    # Step 5: Differential TF activity per cell type
    # -----------------------------------------------------------------------
    log.info("-" * 50)
    log.info("Step 5: Differential TF activity (MASLD vs Normal)")

    tf_activity_df = _differential_tf_activity(
        dev_adata,
        cell_type_col=args.cell_type_col,
        condition_col=args.condition_col,
    )

    tf_csv_path = out_dir / "chromvar_tf_activity.csv"
    if not tf_activity_df.empty:
        tf_activity_df.to_csv(tf_csv_path, index=False)
        log.info(f"  Saved: {tf_csv_path} ({len(tf_activity_df):,} rows)")

        # Summary
        sig_df = tf_activity_df[tf_activity_df["padj"] < 0.05]
        if not sig_df.empty:
            top10 = sig_df.head(10)[["tf_name", "cell_type", "logFC_deviation", "padj"]]
            log.info("  Top 10 significant TF-celltype pairs:")
            for _, row in top10.iterrows():
                log.info(f"    {row['tf_name']:>20s} | {row['cell_type']:>20s} | "
                         f"logFC={row['logFC_deviation']:+.3f} | padj={row['padj']:.2e}")
    else:
        log.warning("  No differential TF activity results to save")

    # -----------------------------------------------------------------------
    # Step 6: TF activity heatmap
    # -----------------------------------------------------------------------
    log.info("-" * 50)
    log.info("Step 6: Generating TF activity heatmap")

    heatmap_path = out_dir / "tf_activity_heatmap.pdf"
    try:
        _plot_tf_heatmap(
            dev_adata,
            str(heatmap_path),
            n_top=args.n_top_heatmap,
            cell_type_col=args.cell_type_col,
        )
    except Exception as exc:
        log.warning(f"Heatmap generation failed: {exc}")

    # -----------------------------------------------------------------------
    # Done
    # -----------------------------------------------------------------------
    elapsed = time.time() - t_start
    log.info("=" * 70)
    log.info(f"chromVAR analysis COMPLETE in {elapsed / 60:.1f} min")
    log.info(f"  Deviations: {dev_h5ad_path}")
    log.info(f"  TF activity: {tf_csv_path}")
    log.info(f"  Heatmap:     {heatmap_path}")
    log.info("=" * 70)


if __name__ == "__main__":
    main()
