#!/usr/bin/env python3
"""
How many genuinely independent tests are there among the 117 programs?

A gene belongs to exactly one program within a cell type, but many genes are
picked up in several cell types, and the 100 kb windows overlap further. The 117
tests run within one trait are therefore correlated, and reporting them as 117
independent tests overstates the multiple testing burden.

This computes the correlation matrix of the 117 program annotations across the
reference panel SNPs and summarises it two ways:

  Li and Ji (2005)        effective number of tests from the eigenvalue spectrum
  Cheverud and Nyholt     effective number from the variance of the eigenvalues

Neither replaces the Benjamini-Hochberg correction, which is already valid under
positive dependence. They are reported so a reader can see how far 117 is from
the number of independent questions actually being asked.

Correlations are accumulated one chromosome at a time so the full SNP by
annotation matrix is never held in memory.

Usage: 04_effective_tests.py <annot_dir> <out_tsv>
"""

import os
import sys

import numpy as np
import pandas as pd

REPO = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
BASE = os.path.join(REPO, "GWAS/ldsc/program_heritability_20260830")
MANIFEST = os.path.join(BASE, "gene_sets/annotation_manifest.tsv")


def main():
    annot_dir, out_tsv = sys.argv[1:3]
    man = pd.read_csv(MANIFEST, sep="\t", dtype=str)
    names = list(man["annotation"])
    prog_names = list(man.loc[man["kind"] == "program", "annotation"])
    m = len(prog_names)

    n = 0
    colsum = np.zeros(m)
    gram = np.zeros((m, m))
    for chrom in range(1, 23):
        path = os.path.join(annot_dir, "programs.%d.annot.gz" % chrom)
        block = pd.read_csv(path, sep="\t", dtype=np.uint8)
        if list(block.columns) != names:
            sys.exit("annotation names in %s do not match the manifest" % path)
        x = block[prog_names].to_numpy(dtype=np.float64)
        n += x.shape[0]
        colsum += x.sum(axis=0)
        gram += x.T @ x
        print("chr%d: %d SNPs" % (chrom, x.shape[0]))

    mean = colsum / n
    cov = gram / n - np.outer(mean, mean)
    sd = np.sqrt(np.clip(np.diag(cov), 0, None))
    denom = np.outer(sd, sd)
    corr = np.where(denom > 0, cov / np.where(denom > 0, denom, 1.0), 0.0)
    np.fill_diagonal(corr, 1.0)

    eig = np.sort(np.abs(np.linalg.eigvalsh(corr)))[::-1]
    li_ji = float(np.sum((eig >= 1).astype(float) + (eig - np.floor(eig))))
    var_eig = float(np.var(eig, ddof=1))
    cheverud = 1.0 + (m - 1.0) * (1.0 - var_eig / m)

    off = corr[np.triu_indices(m, 1)]
    with open(out_tsv, "w") as fh:
        fh.write("quantity\tvalue\n")
        fh.write("n_program_annotations\t%d\n" % m)
        fh.write("n_reference_snps\t%d\n" % n)
        fh.write("mean_pairwise_correlation\t%.4f\n" % off.mean())
        fh.write("median_pairwise_correlation\t%.4f\n" % float(np.median(off)))
        fh.write("max_pairwise_correlation\t%.4f\n" % off.max())
        fh.write("effective_tests_li_ji\t%.2f\n" % li_ji)
        fh.write("effective_tests_cheverud_nyholt\t%.2f\n" % cheverud)

    pd.DataFrame(corr, index=prog_names, columns=prog_names).to_csv(
        os.path.join(os.path.dirname(out_tsv), "program_annotation_correlation.tsv"),
        sep="\t", float_format="%.5f")
    print("effective tests: Li-Ji %.2f, Cheverud-Nyholt %.2f (of %d)"
          % (li_ji, cheverud, m))


if __name__ == "__main__":
    main()
