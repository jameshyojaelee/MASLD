#!/usr/bin/env python3
"""
Mark, for one chromosome, which reference-panel SNPs fall inside each of the
123 gene-program windows built by 01_build_gene_windows.py.

The output is a single "thin" annotation file per chromosome holding all 123
columns at once. Computing LD scores for 123 annotations in one pass over the
genotypes is far cheaper than 123 separate passes; the columns are split apart
afterwards by 03_split_ldscores.py.

SNP order is copied from the PLINK .bim file without reordering. Stratified LD
score regression assumes the annotation rows and the genotype rows are the same
SNPs in the same order, and it does not check.

Usage: 02_make_annot.py <chromosome> <bfile_prefix_dir> <out_dir>
"""

import gzip
import os
import sys
from collections import defaultdict

import numpy as np

REPO = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
BASE = os.path.join(REPO, "GWAS/ldsc/program_heritability_20260830")
BED = os.path.join(BASE, "gene_sets/annotation_windows.bed")
MANIFEST = os.path.join(BASE, "gene_sets/annotation_manifest.tsv")


def annotation_order():
    names = []
    with open(MANIFEST) as fh:
        fh.readline()
        for line in fh:
            names.append(line.split("\t")[0])
    return names


def load_windows(chrom_label):
    spans = defaultdict(list)
    with open(BED) as fh:
        for line in fh:
            chrom, lo, hi, name = line.rstrip("\n").split("\t")
            if chrom == chrom_label:
                spans[name].append((int(lo), int(hi)))
    return spans


def main():
    chrom = sys.argv[1]
    bfile_dir = sys.argv[2]
    out_dir = sys.argv[3]
    os.makedirs(out_dir, exist_ok=True)

    bim = os.path.join(bfile_dir, "1000G.EUR.QC.%s.bim" % chrom)
    pos = np.loadtxt(bim, usecols=3, dtype=np.int64)
    if not np.all(np.diff(pos) >= 0):
        sys.exit("%s is not sorted by base pair; refusing to build annotations" % bim)

    names = annotation_order()
    spans = load_windows("chr%s" % chrom)
    n_snp = pos.shape[0]
    mat = np.zeros((n_snp, len(names)), dtype=np.uint8)

    for j, name in enumerate(names):
        for lo, hi in spans.get(name, []):
            # BED is half-open and 0-based; .bim positions are 1-based.
            i = np.searchsorted(pos, lo + 1, side="left")
            k = np.searchsorted(pos, hi, side="right")
            mat[i:k, j] = 1

    # Write the ASCII table directly. Two bytes per value (digit + separator)
    # keeps this a single buffer copy instead of 780k formatted rows.
    ncol = len(names)
    out = np.full((n_snp, 2 * ncol), 9, dtype=np.uint8)  # tab
    out[:, 0::2] = np.array([48, 49], dtype=np.uint8)[mat]
    out[:, -1] = 10  # newline

    path = os.path.join(out_dir, "programs.%s.annot.gz" % chrom)
    with gzip.open(path, "wb", compresslevel=6) as fh:
        fh.write(("\t".join(names) + "\n").encode())
        fh.write(out.tobytes())

    counts = mat.sum(axis=0)
    summary = os.path.join(out_dir, "programs.%s.snp_counts.tsv" % chrom)
    with open(summary, "w") as fh:
        fh.write("annotation\tn_snps\tn_snps_total\n")
        for name, c in zip(names, counts):
            fh.write("%s\t%d\t%d\n" % (name, c, n_snp))
    print("chr%s: %d SNPs, %d annotations, min=%d median=%d max=%d SNPs per annotation"
          % (chrom, n_snp, ncol, counts.min(), int(np.median(counts)), counts.max()))


if __name__ == "__main__":
    main()
