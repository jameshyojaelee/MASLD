#!/usr/bin/env python3
"""Per-cell-type scATAC peak landscape: promoter / genic / distal annotation.

For each cell-type peak BED file, classify each peak relative to GENCODE v49:
  - promoter: peak midpoint within +/- 2 kb of a gene TSS
  - genic   : peak overlaps a gene body (any exon/intron span) but is NOT a promoter
  - distal  : neither (intergenic)

Output:
    Analysis/ATAC/Human_Multiome/results/snapatac2/peak_landscape_per_ct.tsv
    columns: cell_type, n_peaks_total, n_promoter, n_genic, n_distal

Run via run_31_peak_landscape.sbatch (spatial env, 1 CPU, 8 GB, 1 h, cpu partition).
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
from collections import defaultdict
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

BASE = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
PEAK_DIR = BASE / "Analysis/ATAC/Human_Multiome/results/label_transfer/cell_type_peak_sets_v2"
OUT_TSV = BASE / "Analysis/ATAC/Human_Multiome/results/snapatac2/peak_landscape_per_ct.tsv"
GTF = Path(
    "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
)

PROMOTER_WIN = 2000  # +/- 2 kb of TSS


def parse_gencode(gtf_path: Path) -> tuple[dict, dict]:
    """Parse GTF and return (tss_by_chrom, gene_body_by_chrom).

    tss_by_chrom[chrom]      -> sorted np.array of TSS positions
    gene_body_by_chrom[chrom] -> sorted np.array of (start, end) pairs
    """
    log.info("Parsing GENCODE: %s", gtf_path)
    tss_by_chrom: dict[str, list[int]] = defaultdict(list)
    body_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)

    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    n_genes = 0
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
            # skip patch/hapl/scaff non-canonical contigs by requiring main contigs only
            # (still keep chrX/Y/M as they are main)
            start = int(fields[3])  # 1-based, inclusive
            end = int(fields[4])
            strand = fields[6]
            tss = start if strand == "+" else end
            tss_by_chrom[chrom].append(tss)
            body_by_chrom[chrom].append((start, end))
            n_genes += 1

    # Sort
    tss_sorted = {c: np.array(sorted(v), dtype=np.int64) for c, v in tss_by_chrom.items()}
    # For body: sort by start, then collapse to a single sorted-start array + end array
    body_sorted: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for c, vs in body_by_chrom.items():
        vs_sorted = sorted(vs, key=lambda x: x[0])
        starts = np.array([v[0] for v in vs_sorted], dtype=np.int64)
        ends = np.array([v[1] for v in vs_sorted], dtype=np.int64)
        # cumulative max of ends to enable interval-overlap via binary search
        max_ends = np.maximum.accumulate(ends)
        body_sorted[c] = (starts, ends, max_ends)

    log.info("  Parsed %d gene records across %d chromosomes", n_genes, len(tss_sorted))
    return tss_sorted, body_sorted


def classify_peaks(
    bed_path: Path,
    tss: dict,
    body: dict,
    promoter_win: int = PROMOTER_WIN,
) -> tuple[int, int, int, int]:
    """Classify peaks in a BED file as promoter / genic / distal.

    Returns (n_total, n_promoter, n_genic, n_distal).
    """
    log.info("Classifying %s", bed_path.name)
    df = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
        dtype={"chrom": str, "start": np.int64, "end": np.int64},
    )
    df["mid"] = (df["start"] + df["end"]) // 2
    n_total = len(df)

    n_prom = 0
    n_genic = 0
    n_distal = 0

    # Process per chromosome for cache locality
    for chrom, sub in df.groupby("chrom", sort=False):
        mids = sub["mid"].values
        starts = sub["start"].values
        ends = sub["end"].values

        # ── promoter: |mid - nearest TSS| <= promoter_win ──
        if chrom in tss:
            tss_arr = tss[chrom]
            idx = np.searchsorted(tss_arr, mids)
            # candidate neighbours: idx-1 and idx
            left = np.clip(idx - 1, 0, len(tss_arr) - 1)
            right = np.clip(idx, 0, len(tss_arr) - 1)
            d_left = np.abs(mids - tss_arr[left])
            d_right = np.abs(mids - tss_arr[right])
            nearest = np.minimum(d_left, d_right)
            is_prom = nearest <= promoter_win
        else:
            is_prom = np.zeros(len(sub), dtype=bool)

        # ── genic: peak interval overlaps a gene body ──
        if chrom in body:
            b_starts, b_ends, b_max_ends = body[chrom]
            # An interval (s,e) overlaps any [bs, be] iff bs <= e and be >= s.
            # Find candidate gene index range: first bs > e, all bodies with index < that
            # whose max_end >= s overlap. We test each peak independently for simplicity
            # (vectorised with searchsorted on b_starts).
            # b_starts is sorted ascending.
            right_bound = np.searchsorted(b_starts, ends, side="right")
            # iterate scalar loop only over peaks that didn't pass promoter test
            is_genic = np.zeros(len(sub), dtype=bool)
            for i in range(len(sub)):
                if is_prom[i]:
                    continue
                rb = right_bound[i]
                if rb == 0:
                    continue
                # any body with start <= end[i] and end >= start[i]?
                # walk back from rb-1 while b_starts >= start[i] to find one with end >= start[i],
                # but simpler: among first rb bodies, max end (b_max_ends[rb-1]) must >= start[i]
                if b_max_ends[rb - 1] >= starts[i]:
                    # scan candidate ends in [0, rb) for any >= start[i]
                    # do binary-ish scan: since we only need existence, iterate from rb-1 backward
                    # until first body with start < some threshold; cheap because gene density bounded
                    # Use vectorised any on this slice (typically small after promoter pre-filter)
                    if (b_ends[:rb] >= starts[i]).any():
                        is_genic[i] = True
        else:
            is_genic = np.zeros(len(sub), dtype=bool)

        is_distal = ~(is_prom | is_genic)

        n_prom += int(is_prom.sum())
        n_genic += int(is_genic.sum())
        n_distal += int(is_distal.sum())

    log.info(
        "  total=%d  promoter=%d  genic=%d  distal=%d",
        n_total,
        n_prom,
        n_genic,
        n_distal,
    )
    return n_total, n_prom, n_genic, n_distal


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peak-dir", default=str(PEAK_DIR), type=str)
    parser.add_argument("--gtf", default=str(GTF), type=str)
    parser.add_argument("--out", default=str(OUT_TSV), type=str)
    parser.add_argument("--promoter-win", default=PROMOTER_WIN, type=int)
    args = parser.parse_args()

    peak_dir = Path(args.peak_dir)
    gtf = Path(args.gtf)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    if not peak_dir.is_dir():
        log.error("Peak dir not found: %s", peak_dir)
        return 1
    if not gtf.exists():
        log.error("GTF not found: %s", gtf)
        return 1

    bed_files = sorted(peak_dir.glob("*_peaks.bed"))
    if not bed_files:
        log.error("No *_peaks.bed under %s", peak_dir)
        return 1
    log.info("Found %d cell-type BED files", len(bed_files))

    tss, body = parse_gencode(gtf)

    rows = []
    for bed in bed_files:
        ct = bed.stem.replace("_peaks", "")
        n_total, n_prom, n_genic, n_distal = classify_peaks(
            bed, tss, body, promoter_win=args.promoter_win
        )
        rows.append(
            {
                "cell_type": ct,
                "n_peaks_total": n_total,
                "n_promoter": n_prom,
                "n_genic": n_genic,
                "n_distal": n_distal,
            }
        )

    df = pd.DataFrame(rows).sort_values("n_peaks_total", ascending=False)
    df.to_csv(out, sep="\t", index=False)
    log.info("Wrote: %s", out)
    log.info("\n%s", df.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
