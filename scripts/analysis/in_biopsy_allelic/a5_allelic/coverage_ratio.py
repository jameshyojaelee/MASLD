#!/usr/bin/env python3
"""3'/5' coverage ratio of one library (PRESPEC_A RNA integrity proxy).

Genes with exactly one transcript in the GTF; the 500 with most reads in this
library's featureCounts table. For each, the first and last 20% of the transcript
(in transcript coordinates, strand-aware) are mapped back to genomic exon pieces;
`samtools bedcov` sums per-base depth over each piece. Per gene: mean depth over
the 3' segment / mean depth over the 5' segment. Reported: the median over genes.
"""
import argparse
import json
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

SAMTOOLS = "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/samtools"


def single_isoform_exons(gtf):
    tx_of_gene, exons, strand = defaultdict(set), defaultdict(list), {}
    for line in open(gtf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] != "exon":
            continue
        g = f[8].split('gene_id "', 1)[1].split('"', 1)[0]
        t = f[8].split('transcript_id "', 1)[1].split('"', 1)[0]
        tx_of_gene[g].add(t)
        exons[g].append((f[0], int(f[3]), int(f[4])))
        strand[g] = f[6]
    return {g: (sorted(exons[g], key=lambda e: e[1]), strand[g]) for g, t in tx_of_gene.items() if len(t) == 1}


def segment(exons, strand, frac, from_5prime):
    """Genomic pieces covering the first (5') or last (3') `frac` of the transcript."""
    order = exons if strand == "+" else exons[::-1]
    total = sum(e - s + 1 for _, s, e in exons)
    need = max(1, int(total * frac))
    if not from_5prime:
        order = order[::-1]
    pieces = []
    for c, s, e in order:
        length = e - s + 1
        take = min(length, need)
        at_start = (strand == "+") == from_5prime
        pieces.append((c, s, s + take - 1) if at_start else (c, e - take + 1, e))
        need -= take
        if need <= 0:
            break
    return pieces


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-genes", type=int, default=500)
    a = ap.parse_args()
    fc = pd.read_csv(a.counts, sep="\t", comment="#", index_col=0)
    reads = fc.iloc[:, -1]
    genes = single_isoform_exons(a.gtf)
    top = reads[reads.index.isin(genes.keys())].sort_values(ascending=False).index[: a.n_genes]
    rows = []
    for g in top:
        ex, st = genes[g]
        for end, pieces in (("5p", segment(ex, st, 0.2, True)), ("3p", segment(ex, st, 0.2, False))):
            for c, s, e in pieces:
                rows.append((c, s - 1, e, g, end))
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as fh:
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
        bed = fh.name
    res = subprocess.run([SAMTOOLS, "bedcov", bed, a.bam], capture_output=True, text=True, check=True).stdout
    Path(bed).unlink()
    cov = pd.DataFrame([l.split("\t") for l in res.strip().split("\n")], columns=["chrom", "start", "end", "gene", "end5or3", "sum"])
    cov[["start", "end", "sum"]] = cov[["start", "end", "sum"]].astype(float)
    cov["len"] = cov["end"] - cov["start"]
    per = cov.groupby(["gene", "end5or3"]).agg(s=("sum", "sum"), l=("len", "sum")).reset_index()
    per["mean_depth"] = per["s"] / per["l"]
    wide = per.pivot(index="gene", columns="end5or3", values="mean_depth")
    wide = wide[(wide["5p"] > 0) & (wide["3p"] > 0)]
    ratio = np.log2(wide["3p"] / wide["5p"])
    Path(a.out).write_text(json.dumps({"genes_used": int(len(ratio)),
                                       "median_log2_3p_over_5p": float(ratio.median()),
                                       "iqr_log2_3p_over_5p": float(ratio.quantile(0.75) - ratio.quantile(0.25))}, indent=2))


if __name__ == "__main__":
    main()
