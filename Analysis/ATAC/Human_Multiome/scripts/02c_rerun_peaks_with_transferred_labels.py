#!/usr/bin/env python3
"""
02c_rerun_peaks_with_transferred_labels.py — Re-run MACS3 peak calling with
label-transferred cell types from 02b.

Reads the label-transferred AnnData, recreates the AnnDataSet from per-donor
h5ads (for fragment file access required by MACS3), and runs per-cell-type
peak calling using corrected labels. Compares old vs new peak sets.

Pipeline position:
  02_snapatac2_processing.py → 02b_label_transfer_from_rna.py
                             → **02c_rerun_peaks_with_transferred_labels.py**

Inputs:
  - results/label_transfer/snapatac2_label_transferred.h5ad (from 02b)
  - results/snapatac2/per_donor/*.h5ad (backed h5ads with fragment refs)

Outputs:
  - results/label_transfer/cell_type_peak_sets_v2/  Per-cell-type BED files
  - results/label_transfer/peak_comparison.csv       Old vs new peak set comparison
  - results/label_transfer/snapatac2_relabeled.h5ad  h5ad with cell_type = transferred

Environment:
  micromamba activate snapatac2

Usage:
  cd Analysis/ATAC/Human_Multiome
  python scripts/02c_rerun_peaks_with_transferred_labels.py \
      --label_transferred results/label_transfer/snapatac2_label_transferred.h5ad \
      --per_donor_dir results/snapatac2/per_donor \
      --old_peak_dir results/snapatac2/cell_type_peak_sets \
      --output_dir results/label_transfer
"""

from __future__ import annotations

import argparse
import gc
import logging
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore", category=FutureWarning)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser(
        description="Re-run MACS3 peak calling with label-transferred cell types"
    )
    p.add_argument(
        "--label_transferred",
        default="results/label_transfer/snapatac2_label_transferred.h5ad",
        help="Label-transferred AnnData from 02b",
    )
    p.add_argument(
        "--per_donor_dir",
        default="results/snapatac2/per_donor",
        help="Per-donor backed h5ads with fragment references",
    )
    p.add_argument(
        "--old_peak_dir",
        default="results/snapatac2/cell_type_peak_sets",
        help="Original per-cell-type peak BED directory",
    )
    p.add_argument(
        "--output_dir",
        default="results/label_transfer",
        help="Output directory",
    )
    return p.parse_args()


def load_and_relabel(label_transferred_path: str) -> sc.AnnData:
    """Load label-transferred h5ad, set cell_type to transferred labels.

    Filters out Low_confidence and Unassigned cells before peak calling.
    """
    log.info("Loading label-transferred AnnData: %s", label_transferred_path)
    adata = sc.read_h5ad(label_transferred_path)

    # Preserve original labels
    if "cell_type_gene_activity" not in adata.obs.columns:
        adata.obs["cell_type_gene_activity"] = adata.obs["cell_type"].copy()

    # Replace cell_type with transferred labels (as plain strings, not Categorical,
    # because AnnDataSet.obs cannot accept Categorical dtype).
    # Sanitize names: replace / and other special chars with _ for MACS3 file output.
    import re
    adata.obs["cell_type"] = (
        adata.obs["cell_type_transferred"]
        .astype(str)
        .apply(lambda x: re.sub(r'[/\\:*?"<>| ]', '_', x))
        .values
    )

    # Keep Low_confidence/Unassigned cells — MACS3 will create separate peak sets
    # for them which can be ignored. Removing them would cause cell count mismatch
    # with the AnnDataSet (which has all cells from per-donor h5ads).
    n_low = (adata.obs["cell_type"].isin(["Low_confidence", "Unassigned"])).sum()
    if n_low > 0:
        log.info(
            "%d cells are Low_confidence/Unassigned (kept for AnnDataSet alignment, "
            "peak sets will be ignored downstream)", n_low,
        )

    log.info("Relabeled cell types:\n%s", adata.obs["cell_type"].value_counts().to_string())
    return adata


def create_anndataset(per_donor_dir: str, adata: sc.AnnData):
    """Recreate AnnDataSet from per-donor backed h5ads for MACS3 fragment access."""
    import snapatac2 as snap

    h5ad_files = sorted(Path(per_donor_dir).glob("*.h5ad"))
    if not h5ad_files:
        raise FileNotFoundError(f"No h5ad files in {per_donor_dir}")

    log.info("Opening %d per-donor backed h5ads...", len(h5ad_files))
    # SnapATAC2 AnnDataSet expects list of (name_str, AnnData) tuples
    processed = []
    for f in h5ad_files:
        a = snap.read(str(f), backed="r")
        processed.append((f.stem, a))
        log.info("  %s: %d cells", f.stem, a.n_obs)

    # Create temporary AnnDataSet
    temp_h5ads = str(Path(per_donor_dir).parent / "temp_peak_recall.h5ads")
    dataset = snap.AnnDataSet(
        adatas=processed,
        filename=temp_h5ads,
        add_key="donor_id",
    )
    log.info("AnnDataSet: %d cells x %d features", dataset.n_obs, dataset.n_vars)

    return dataset, processed, temp_h5ads


def run_macs3(dataset, adata, output_dir: Path):
    """Run per-cell-type MACS3 peak calling with transferred labels.

    Uses snap.tl.macs3() but monkey-patches tempfile.TemporaryDirectory to
    persist MACS3 output files. The SnapATAC2 wrapper crashes at the final
    adata.uns write (HDF5 version mismatch in backed AnnDataSet), but by then
    MACS3 has already written narrowPeak files to the temp dir. We intercept
    those files before cleanup.
    """
    import tempfile

    import snapatac2 as snap

    peak_dir = output_dir / "cell_type_peak_sets_v2"
    peak_dir.mkdir(parents=True, exist_ok=True)

    # Transfer cell_type labels to AnnDataSet
    try:
        ct_values = adata.obs["cell_type"].astype(str).tolist()
        dataset.obs["cell_type"] = ct_values
        log.info("Set cell_type on AnnDataSet for per-cell-type peak calling")
    except Exception as e:
        log.error("Cannot set cell_type on AnnDataSet: %s", e)
        return {}

    # Persistent directory for MACS3 output (survives the uns write crash)
    macs3_tmpdir = str(output_dir / "macs3_tmp")
    os.makedirs(macs3_tmpdir, exist_ok=True)

    # Monkey-patch TemporaryDirectory so snap.tl.macs3 writes to our
    # persistent dir instead of a real temp dir that gets cleaned up on error.
    # This is necessary because snap.tl.macs3 creates a TemporaryDirectory
    # internally and deletes it before returning peak results.
    import threading
    _patch_lock = threading.Lock()
    OrigTmpDir = tempfile.TemporaryDirectory

    class PersistentTmpDir:
        """Drop-in replacement that writes to a persistent directory."""
        def __init__(self, *a, **kw):
            self.name = macs3_tmpdir

        def __enter__(self):
            return self.name

        def __exit__(self, *a):
            pass  # don't delete — we need the output files

        def cleanup(self):
            pass

    with _patch_lock:
        tempfile.TemporaryDirectory = PersistentTmpDir

    log.info("Running MACS3 (groupby=cell_type) with persistent tmpdir...")
    try:
        peak_result = snap.tl.macs3(
            dataset,
            groupby="cell_type",
            qvalue=0.05,
            shift=-100,
            extsize=200,
        )
        log.info("MACS3 completed without error")
    except RuntimeError as e:
        if "H5Gcreate2" in str(e) or "bad object header" in str(e):
            log.info(
                "Caught expected HDF5 write error (peaks computed successfully). "
                "Recovering from dataset.uns or tmpdir.",
            )
            peak_result = None
        else:
            tempfile.TemporaryDirectory = OrigTmpDir
            raise
    finally:
        with _patch_lock:
            tempfile.TemporaryDirectory = OrigTmpDir

    # --- Extract peaks from return value, dataset.uns, or narrowPeak files ---
    peak_counts = {}

    # Strategy 1: Use return value from snap.tl.macs3()
    peaks_data = peak_result
    log.info("peak_result type: %s, is None: %s", type(peak_result).__name__, peak_result is None)

    if peaks_data is None:
        # Strategy 2: Check dataset.uns with explicit key names
        for key_name in ["peaks", "macs3", "peak_annotation"]:
            try:
                val = dataset.uns[key_name]
                if val is not None:
                    peaks_data = val
                    log.info("Found peaks in dataset.uns['%s'] (type: %s)", key_name, type(val).__name__)
                    break
            except (KeyError, Exception) as e:
                log.debug("  dataset.uns['%s'] not found: %s", key_name, e)

    if peaks_data is None:
        # Strategy 2b: List all uns keys
        try:
            uns_keys = list(dataset.uns.keys())
            log.info("dataset.uns keys: %s", uns_keys)
            for key in uns_keys:
                try:
                    val = dataset.uns[key]
                    log.info("  uns['%s'] type=%s", key, type(val).__name__)
                    if isinstance(val, dict):
                        peaks_data = val
                        log.info("  Using uns['%s'] as peaks_data", key)
                        break
                except Exception as e2:
                    log.warning("  uns['%s'] access error: %s", key, e2)
        except Exception as e:
            log.warning("Cannot list dataset.uns keys: %s", e)

    if peaks_data is not None:
        log.info("Extracting peaks from return value/uns (type: %s)", type(peaks_data).__name__)
        if isinstance(peaks_data, dict):
            for group_name, peaks_obj in peaks_data.items():
                if hasattr(peaks_obj, "to_pandas"):
                    peaks_df = peaks_obj.to_pandas()
                elif isinstance(peaks_obj, pd.DataFrame):
                    peaks_df = peaks_obj
                else:
                    log.warning("  %s: unknown peak format %s", group_name, type(peaks_obj))
                    continue
                if len(peaks_df) > 0:
                    bed_path = peak_dir / f"{group_name}_peaks.bed"
                    # DEDUPE: macs3 groupby can emit coordinate-duplicate intervals
                    # (~24% for abundant cell types); writing them verbatim inflates
                    # every downstream peak-count and FDR denominator. Keep unique
                    # (chr,start,end) only.
                    bed3 = peaks_df.iloc[:, :3].drop_duplicates()
                    bed3.to_csv(bed_path, sep="\t", header=False, index=False)
                    peak_counts[group_name] = len(bed3)
                    log.info("  %s: %d peaks -> %s", group_name, len(peaks_df), bed_path)
        else:
            log.warning("peaks_data is not a dict: %s", type(peaks_data))

    # Strategy 3: Fall back to narrowPeak files in tmpdir
    if not peak_counts:
        import glob
        narrowpeaks = glob.glob(os.path.join(macs3_tmpdir, "**", "*.narrowPeak"), recursive=True)
        if not narrowpeaks:
            narrowpeaks = glob.glob(os.path.join(macs3_tmpdir, "*.narrowPeak"))
        log.info("Fallback: found %d narrowPeak files in %s", len(narrowpeaks), macs3_tmpdir)

        for np_file in sorted(narrowpeaks):
            basename = os.path.basename(np_file)
            ct_name = basename.replace("_peaks.narrowPeak", "").replace(".narrowPeak", "")
            n_peaks = sum(1 for _ in open(np_file))
            if n_peaks == 0:
                continue
            bed_path = peak_dir / f"{ct_name}_peaks.bed"
            _seen = set()
            with open(np_file) as fin, open(bed_path, "w") as fout:
                for line in fin:
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        key = (parts[0], parts[1], parts[2])
                        if key in _seen:          # dedupe coordinate-duplicate intervals
                            continue
                        _seen.add(key)
                        fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")
            peak_counts[ct_name] = len(_seen)
            log.info("  %s: %d peaks -> %s", ct_name, n_peaks, bed_path)

    # Remove Low_confidence and Unassigned peak sets (not biologically meaningful)
    # and flag cell types with <500 cells as underpowered
    for exclude_ct in ["Low_confidence", "Unassigned"]:
        bed = peak_dir / f"{exclude_ct}_peaks.bed"
        if bed.exists():
            bed.unlink()
            peak_counts.pop(exclude_ct, None)
            log.info("  Removed %s peak set (not biologically meaningful)", exclude_ct)

    for ct, n_peaks in list(peak_counts.items()):
        ct_cells = (adata.obs["cell_type"] == ct).sum()
        if ct_cells < 500:
            log.warning(
                "  %s: only %d cells — peak set (%d peaks) is underpowered",
                ct, ct_cells, n_peaks,
            )

    log.info("Peak calling complete: %d cell types with peaks", len(peak_counts))
    return peak_counts


def _load_peaks_as_intervals(bed_path):
    """Load BED file as list of (chrom, start, end) tuples with integer coords."""
    peaks = []
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                peaks.append((parts[0], int(parts[1]), int(parts[2])))
    return peaks


def _overlap_jaccard(peaks_a, peaks_b, min_overlap_frac=0.5):
    """Compute overlap-based Jaccard between two peak sets.

    Two peaks overlap if they share the same chromosome and have
    reciprocal overlap >= min_overlap_frac of the shorter peak.
    """
    if not peaks_a or not peaks_b:
        return 0.0, 0, 0, 0

    # Index peaks_b by chromosome for fast lookup
    from collections import defaultdict
    b_by_chrom = defaultdict(list)
    for chrom, start, end in peaks_b:
        b_by_chrom[chrom].append((start, end))
    # Sort by start for binary search
    for chrom in b_by_chrom:
        b_by_chrom[chrom].sort()

    n_a_matched = 0
    matched_b = set()

    for chrom, a_start, a_end in peaks_a:
        a_len = a_end - a_start
        for i, (b_start, b_end) in enumerate(b_by_chrom.get(chrom, [])):
            if b_start > a_end:
                break  # sorted, no more overlaps possible
            if b_end < a_start:
                continue
            # Compute overlap
            ovl = min(a_end, b_end) - max(a_start, b_start)
            min_len = min(a_len, b_end - b_start)
            if min_len > 0 and ovl / min_len >= min_overlap_frac:
                n_a_matched += 1
                matched_b.add((chrom, b_start, b_end))
                break  # count each A peak at most once

    intersection = n_a_matched  # ~= len(matched_b)
    union = len(peaks_a) + len(peaks_b) - intersection
    jaccard = intersection / union if union > 0 else 0.0
    return jaccard, intersection, n_a_matched, len(matched_b)


def compare_peak_sets(old_peak_dir: str, new_peak_dir: Path, output_dir: Path):
    """Compare old vs new per-cell-type peak sets using overlap-based Jaccard."""
    old_dir = Path(old_peak_dir)
    comparison = []

    if not old_dir.exists():
        log.warning("Old peak directory not found: %s", old_dir)
        return

    for new_bed in sorted(new_peak_dir.glob("*_peaks.bed")):
        ct_name = new_bed.stem.replace("_peaks", "")
        old_bed = old_dir / new_bed.name

        new_peaks = _load_peaks_as_intervals(new_bed)

        if old_bed.exists():
            old_peaks = _load_peaks_as_intervals(old_bed)
            jaccard, intersection, _, _ = _overlap_jaccard(old_peaks, new_peaks)

            comparison.append(
                {
                    "cell_type": ct_name,
                    "old_n_peaks": len(old_peaks),
                    "new_n_peaks": len(new_peaks),
                    "overlap_intersection": intersection,
                    "jaccard_overlap": jaccard,
                }
            )
            log.info(
                "  %s: old=%d, new=%d, Jaccard(overlap)=%.3f",
                ct_name, len(old_peaks), len(new_peaks), jaccard,
            )
        else:
            comparison.append(
                {
                    "cell_type": ct_name,
                    "old_n_peaks": 0,
                    "new_n_peaks": len(new_peaks),
                    "overlap_intersection": 0,
                    "jaccard_overlap": 0.0,
                }
            )
            log.info("  %s: NEW cell type (no old peaks), %d peaks", ct_name, len(new_peaks))

    # Check for old cell types that disappeared
    for old_bed in sorted(old_dir.glob("*_peaks.bed")):
        ct_name = old_bed.stem.replace("_peaks", "")
        new_bed = new_peak_dir / old_bed.name
        if not new_bed.exists():
            old_peaks_count = sum(1 for _ in open(old_bed))
            comparison.append(
                {
                    "cell_type": ct_name,
                    "old_n_peaks": old_peaks_count,
                    "new_n_peaks": 0,
                    "overlap_intersection": 0,
                    "jaccard_overlap": 0.0,
                }
            )
            log.info("  %s: REMOVED cell type, was %d peaks", ct_name, old_peaks_count)

    if comparison:
        df = pd.DataFrame(comparison)
        df.to_csv(output_dir / "peak_comparison.csv", index=False)
        log.info("\nPeak comparison saved to %s", output_dir / "peak_comparison.csv")

        mean_jaccard = df["jaccard_overlap"].mean()
        log.info("Mean overlap Jaccard across cell types: %.3f", mean_jaccard)


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # File logging
    fh = logging.FileHandler(output_dir / "peak_recalling.log")
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    log.addHandler(fh)

    log.info("=" * 60)
    log.info("Re-run MACS3 Peak Calling with Transferred Labels")
    log.info("=" * 60)

    # --- Load and relabel ---
    adata = load_and_relabel(args.label_transferred)

    # Save relabeled h5ad (cell_type = transferred) for downstream chromVAR/SCENIC+
    relabeled_path = output_dir / "snapatac2_relabeled.h5ad"
    adata.write_h5ad(str(relabeled_path))
    log.info("Saved relabeled h5ad: %s", relabeled_path)

    # --- Create AnnDataSet for fragment access ---
    dataset, processed, temp_h5ads = create_anndataset(args.per_donor_dir, adata)

    # --- Run MACS3 ---
    try:
        peak_counts = run_macs3(dataset, adata, output_dir)
    finally:
        # Clean up AnnDataSet
        dataset.close()
        try:
            Path(temp_h5ads).unlink(missing_ok=True)
        except Exception:
            pass
        for _, a in processed:
            try:
                if hasattr(a, "file") and a.file:
                    a.file.close()
            except Exception:
                pass

    # --- Compare old vs new peak sets ---
    log.info("\n=== Peak Set Comparison ===")
    new_peak_dir = output_dir / "cell_type_peak_sets_v2"
    compare_peak_sets(args.old_peak_dir, new_peak_dir, output_dir)

    log.info("\n" + "=" * 60)
    log.info("Peak re-calling COMPLETE")
    log.info("  New peaks: %s", new_peak_dir)
    log.info("  Relabeled h5ad: %s", relabeled_path)
    log.info("=" * 60)


if __name__ == "__main__":
    main()
