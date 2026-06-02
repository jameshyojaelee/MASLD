#!/usr/bin/env python3
"""Annotate DA peaks with nearest gene symbols for L8 integration.

Solves two issues:
  Module 2: DA results have peak coordinates (chr:start-end) but no gene symbols.
            Maps each peak to nearest gene TSS within a promoter window.
  Module 2c: chromVAR results are per-TF per-cell-type, not per-gene.
             Maps TFs to genes via Module 3 regulon targets.

Outputs annotated CSVs that the L8 integration script can consume.

Usage:
    python 08_annotate_peaks_for_l8.py \
        --da-results results/snapatac2/scatac_da_results.csv \
        --chromvar results/chromvar/chromvar_tf_activity.csv \
        --regulons scenic_plus/hepatocyte_regulons.csv \
        --gtf /path/to/gencode.v49.gtf.gz \
        --output-dir results/l8_annotated
"""

import argparse
import gzip
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)


def parse_tss_from_gtf(gtf_path):
    """Parse gene TSS positions from GENCODE GTF.

    Returns DataFrame: gene_name, chrom, tss, strand
    """
    log.info("Parsing TSS from GTF: %s", gtf_path)
    records = []
    opener = gzip.open if gtf_path.endswith(".gz") else open

    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            chrom = fields[0]
            if not chrom.startswith("chr"):
                continue

            strand = fields[6]
            start = int(fields[3])  # 1-based
            end = int(fields[4])

            attrs = fields[8]
            gene_name = None
            gene_type = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_name"):
                    gene_name = attr.split('"')[1] if '"' in attr else attr.split(" ")[1]
                elif attr.startswith("gene_type"):
                    gene_type = attr.split('"')[1] if '"' in attr else attr.split(" ")[1]

            if gene_name is None:
                continue

            tss = start if strand == "+" else end
            records.append({
                "gene_name": gene_name,
                "chrom": chrom,
                "tss": tss,
                "strand": strand,
                "gene_type": gene_type or "",
            })

    df = pd.DataFrame(records)
    # Keep first occurrence per gene (primary annotation)
    df = df.drop_duplicates(subset="gene_name", keep="first")
    log.info("  Parsed %d gene TSS positions", len(df))
    return df


def annotate_da_peaks(da_path, tss_df, promoter_window=5000):
    """Map DA peaks to nearest gene within promoter window.

    For each peak, find the nearest gene TSS. If within promoter_window,
    annotate the peak with the gene symbol.

    Returns annotated DataFrame with gene_symbol column.
    """
    log.info("Annotating DA peaks with nearest gene...")
    da = pd.read_csv(da_path)
    log.info("  DA results: %d rows, columns: %s", len(da), list(da.columns))

    # Parse peak coordinates
    peak_col = "feature name"
    if peak_col not in da.columns:
        for c in da.columns:
            if "feature" in c.lower() or "peak" in c.lower() or "region" in c.lower():
                peak_col = c
                break

    coords = da[peak_col].str.extract(r"(chr\w+):(\d+)-(\d+)")
    coords.columns = ["chrom", "start", "end"]
    coords["start"] = coords["start"].astype(int)
    coords["end"] = coords["end"].astype(int)
    coords["mid"] = (coords["start"] + coords["end"]) // 2

    da = pd.concat([da, coords], axis=1)

    # Build TSS lookup by chromosome
    chrom_tss = {}
    for chrom in tss_df["chrom"].unique():
        ct = tss_df[tss_df["chrom"] == chrom].sort_values("tss")
        chrom_tss[chrom] = ct[["gene_name", "tss"]].values  # (N, 2)

    # Map each peak to nearest gene
    gene_symbols = []
    distances = []

    for _, row in da.iterrows():
        chrom = row["chrom"]
        mid = row["mid"]

        if chrom not in chrom_tss or pd.isna(chrom):
            gene_symbols.append(None)
            distances.append(np.nan)
            continue

        genes_on_chrom = chrom_tss[chrom]
        tss_positions = genes_on_chrom[:, 1].astype(int)

        # Binary search for nearest TSS
        idx = np.searchsorted(tss_positions, mid)
        candidates = []
        if idx > 0:
            candidates.append(idx - 1)
        if idx < len(tss_positions):
            candidates.append(idx)

        best_gene = None
        best_dist = np.inf
        for c in candidates:
            dist = abs(int(tss_positions[c]) - mid)
            if dist < best_dist:
                best_dist = dist
                best_gene = genes_on_chrom[c, 0]

        if best_dist <= promoter_window:
            gene_symbols.append(best_gene)
            distances.append(best_dist)
        else:
            gene_symbols.append(best_gene)  # still annotate, let L8 decide
            distances.append(best_dist)

    da["gene_symbol"] = gene_symbols
    da["distance_to_tss"] = distances

    # For peaks with multiple genes at same distance, keep nearest
    annotated = da.dropna(subset=["gene_symbol"])
    log.info("  Annotated: %d / %d peaks mapped to genes", len(annotated), len(da))
    log.info("  Peaks within %dkb of TSS: %d",
             promoter_window // 1000,
             (annotated["distance_to_tss"] <= promoter_window).sum())

    # Rename columns to match L8 expectations
    col_map = {
        "log2(fold_change)": "logFC",
        "adjusted p-value": "padj",
        "p-value": "pvalue",
    }
    annotated = annotated.rename(columns=col_map)

    return annotated


def annotate_chromvar_with_genes(chromvar_path, regulon_path):
    """Map chromVAR TF results to genes via Module 3 regulon targets.

    For each (TF, cell_type) row in chromVAR, expand to gene-level rows
    using the TF's target genes from the SCENIC+ regulons.

    Returns DataFrame with gene_symbol column added.
    """
    log.info("Mapping chromVAR TFs to genes via regulon targets...")
    cv = pd.read_csv(chromvar_path)
    log.info("  chromVAR: %d rows", len(cv))

    regulons = pd.read_csv(regulon_path)
    log.info("  Regulons: %d rows", len(regulons))

    # Build TF -> target genes mapping
    tf_targets = {}
    for _, row in regulons.iterrows():
        tf = row.get("tf_name", "")
        target = row.get("target_gene", "")
        if tf and target:
            tf_targets.setdefault(tf, set()).add(target)

    log.info("  %d TFs with regulon targets", len(tf_targets))

    # Expand chromVAR rows: for each TF in chromVAR that has regulon targets,
    # create one row per target gene
    expanded = []
    for _, row in cv.iterrows():
        tf = row["tf_name"]
        # Normalize TF name (chromVAR may use different capitalization or composite names)
        tf_clean = tf.split("::")[0]  # Handle composite motifs like ARNT::HIF1A

        targets = tf_targets.get(tf_clean, set())
        if not targets:
            # Try case-insensitive match
            for reg_tf in tf_targets:
                if reg_tf.upper() == tf_clean.upper():
                    targets = tf_targets[reg_tf]
                    break

        if targets:
            for gene in targets:
                new_row = row.to_dict()
                new_row["target_gene"] = gene
                expanded.append(new_row)

    if not expanded:
        log.warning("  No TF-gene mappings found. ChromVAR cannot be mapped to genes.")
        return pd.DataFrame()

    result = pd.DataFrame(expanded)
    log.info("  Expanded to %d gene-level rows (%d unique genes)",
             len(result), result["target_gene"].nunique())
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Annotate DA peaks and chromVAR results with gene symbols for L8",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--da-results", required=True,
                        help="scATAC DA results CSV (peak-level)")
    parser.add_argument("--chromvar", required=True,
                        help="chromVAR TF activity CSV")
    parser.add_argument("--regulons", required=True,
                        help="Module 3 hepatocyte_regulons.csv")
    parser.add_argument("--gtf", required=True,
                        help="GENCODE GTF for TSS positions")
    parser.add_argument("--output-dir", default="results/l8_annotated",
                        help="Output directory for annotated CSVs")
    parser.add_argument("--promoter-window", type=int, default=5000,
                        help="Promoter window in bp for peak-gene mapping")
    args = parser.parse_args()

    t0 = time.time()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Parse TSS
    tss_df = parse_tss_from_gtf(args.gtf)

    # Module 2: Annotate DA peaks
    log.info("=" * 60)
    log.info("Module 2: Annotating DA peaks with gene symbols")
    da_annotated = annotate_da_peaks(args.da_results, tss_df,
                                     promoter_window=args.promoter_window)
    da_out = out_dir / "scatac_da_gene_annotated.csv"
    da_annotated.to_csv(da_out, index=False)
    log.info("  Saved: %s (%d rows)", da_out, len(da_annotated))

    # Cell type distribution
    if "cell_type" in da_annotated.columns:
        log.info("  Per cell type:")
        for ct, count in da_annotated["cell_type"].value_counts().items():
            log.info("    %s: %d peaks", ct, count)

    # Module 2c: Annotate chromVAR with gene symbols
    log.info("=" * 60)
    log.info("Module 2c: Mapping chromVAR TFs to genes via regulons")
    cv_annotated = annotate_chromvar_with_genes(args.chromvar, args.regulons)
    cv_out = out_dir / "chromvar_gene_annotated.csv"
    if not cv_annotated.empty:
        cv_annotated.to_csv(cv_out, index=False)
        log.info("  Saved: %s (%d rows)", cv_out, len(cv_annotated))
    else:
        log.warning("  No chromVAR gene annotations produced")

    log.info("=" * 60)
    log.info("Done in %.1fs", time.time() - t0)


if __name__ == "__main__":
    main()
