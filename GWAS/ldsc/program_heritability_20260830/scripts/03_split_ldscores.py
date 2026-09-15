#!/usr/bin/env python3
"""
Split the 123-column LD score output into one file set per annotation.

LD scores are computed once for all 123 annotations together, but the
partitioned regression is run one program at a time (each program tested on its
own, conditional on the baselineLD model and on its own cell-type background).
That requires each annotation to sit in its own prefix, in the four-file layout
LDSC expects:

    <slug>.<chr>.annot.gz        the 0/1 membership column
    <slug>.<chr>.l2.ldscore.gz   the LD score column
    <slug>.<chr>.l2.M            SNP count
    <slug>.<chr>.l2.M_5_50       SNP count, common variants only

The .annot.gz files are needed because --overlap-annot re-reads them to work out
how much each pair of annotations overlaps before converting coefficients into
enrichment.

Usage: 03_split_ldscores.py <chromosome> <annot_dir> <ldscore_dir> <out_dir>
"""

import gzip
import os
import sys

import numpy as np
import pandas as pd

REPO = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
BASE = os.path.join(REPO, "GWAS/ldsc/program_heritability_20260830")
MANIFEST = os.path.join(BASE, "gene_sets/annotation_manifest.tsv")

ASCII = np.array([48, 49], dtype=np.uint8)  # '0', '1'


def read_manifest():
    df = pd.read_csv(MANIFEST, sep="\t", dtype=str)
    return list(df["annotation"]), list(df["file_slug"])


def main():
    chrom, annot_dir, ldscore_dir, out_dir = sys.argv[1:5]
    os.makedirs(out_dir, exist_ok=True)
    names, slugs = read_manifest()

    # --- annotation columns ---
    src = os.path.join(annot_dir, "programs.%s.annot.gz" % chrom)
    annot = pd.read_csv(src, sep="\t", dtype=np.uint8)
    if list(annot.columns) != names:
        sys.exit("annotation names in %s do not match the manifest" % src)
    n_snp = annot.shape[0]
    for name, slug in zip(names, slugs):
        col = annot[name].to_numpy()
        buf = np.full((n_snp, 2), 10, dtype=np.uint8)  # newline in column 1
        buf[:, 0] = ASCII[col]
        with gzip.open(os.path.join(out_dir, "%s.%s.annot.gz" % (slug, chrom)),
                       "wb", compresslevel=6) as fh:
            fh.write((name + "\n").encode())
            fh.write(buf.tobytes())
    del annot

    # --- LD score columns ---
    src = os.path.join(ldscore_dir, "programs.%s.l2.ldscore.gz" % chrom)
    ld = pd.read_csv(src, sep="\t")
    expected = ["CHR", "SNP", "BP"] + [n + "L2" for n in names]
    if list(ld.columns) != expected:
        sys.exit("unexpected column layout in %s" % src)
    for name, slug in zip(names, slugs):
        out = ld[["CHR", "SNP", "BP", name + "L2"]]
        out.to_csv(os.path.join(out_dir, "%s.%s.l2.ldscore.gz" % (slug, chrom)),
                   sep="\t", index=False, compression="gzip")
    del ld

    # --- SNP counts ---
    for suffix in ("l2.M", "l2.M_5_50"):
        src = os.path.join(ldscore_dir, "programs.%s.%s" % (chrom, suffix))
        with open(src) as fh:
            values = fh.read().split()
        if len(values) != len(names):
            sys.exit("%s holds %d values but there are %d annotations"
                     % (src, len(values), len(names)))
        for slug, v in zip(slugs, values):
            with open(os.path.join(out_dir, "%s.%s.%s" % (slug, chrom, suffix)), "w") as fh:
                fh.write(v + "\n")

    print("chr%s: %d SNPs split into %d annotation file sets under %s"
          % (chrom, n_snp, len(names), out_dir))


if __name__ == "__main__":
    main()
