"""Per-gene, stranded exon intervals of the tag genes (B-Dup, for 07_gene_dup.py).

01_regions.py wrote only the union of the tag-gene exons, with no gene label or strand.
07 assigns each read to one tag gene by its aligned blocks and the library's transcribed
strand, so it needs each gene's own exons and strand. This script rebuilds them from the
same GTF with the same rule as 01 (a gene's exons on its own chromosome, merged where
they overlap or touch) and checks the result against 01's outputs.

Inputs
  --tag-genes  tag_genes.tsv from 01 (gene_key, chrom, exon_bp, ...). Gene list only.
  --gtf        the GENCODE v49 primary-chromosome GTF that 01 used.
  --union-bed  tag_gene_exons.bed from 01 (merged union), for the check.
Output
  --out        BED6: chrom, start, end (0-based, half-open), gene_key, 0, strand.
Checks (the script stops if any fails)
  every tag gene found; per-gene merged exon bp equals tag_genes.tsv exon_bp;
  the union of all per-gene intervals equals tag_gene_exons.bed exactly.
Refuses to overwrite --out.
"""
import argparse
import os
import re
import sys
from collections import defaultdict

import pandas as pd

ATTR = re.compile(r'(\S+) "([^"]*)"')
CHROM_ORDER = {f"chr{i}": i for i in range(1, 23)} | {"chrX": 23, "chrY": 24, "chrM": 25}


def merge(intervals):
    """Merge 0-based half-open intervals that overlap or touch (same rule as 01)."""
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag-genes", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--union-bed", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if os.path.exists(args.out):
        sys.exit(f"refusing to overwrite {args.out}")

    tag = pd.read_csv(args.tag_genes, sep="\t", dtype={"gene_key": str, "chrom": str})
    keys = set(tag["gene_key"])
    gene_info, exons = {}, defaultdict(list)
    with open(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("gene", "exon"):
                continue
            key = dict(ATTR.findall(f[8]))["gene_id"].split(".")[0]
            if key not in keys:
                continue
            if f[2] == "gene":
                gene_info[key] = (f[0], f[6])
            else:
                exons[key].append((f[0], int(f[3]) - 1, int(f[4])))
    missing = sorted(keys - set(gene_info))
    if missing:
        sys.exit(f"{len(missing)} tag genes not in the GTF, e.g. {missing[:5]}")

    rows, union = [], defaultdict(list)
    expected_bp = dict(zip(tag["gene_key"], tag["exon_bp"]))
    for key in sorted(keys):
        chrom, strand = gene_info[key]
        if strand not in "+-":
            sys.exit(f"{key}: strand {strand!r}")
        gene_exons = merge([(s, e) for c, s, e in exons[key] if c == chrom])
        bp = sum(e - s for s, e in gene_exons)
        if bp != expected_bp[key]:
            sys.exit(f"{key}: exon bp {bp} != tag_genes.tsv {expected_bp[key]}")
        for s, e in gene_exons:
            rows.append((chrom, s, e, key, 0, strand))
            union[chrom].append((s, e))

    rebuilt = [(c, s, e) for c in sorted(union, key=lambda c: CHROM_ORDER.get(c, 99)) for s, e in merge(union[c])]
    with open(args.union_bed) as fh:
        given = [(c, int(s), int(e)) for c, s, e in (line.split("\t")[:3] for line in fh)]
    if rebuilt != given:
        sys.exit(f"union of per-gene exons ({len(rebuilt)} intervals) != {args.union_bed} ({len(given)})")

    rows.sort(key=lambda r: (CHROM_ORDER.get(r[0], 99), r[1], r[2], r[3]))
    with open(args.out + ".part", "w") as out:
        for r in rows:
            out.write("\t".join(map(str, r)) + "\n")
    os.replace(args.out + ".part", args.out)
    n_plus = sum(1 for k in keys if gene_info[k][1] == "+")
    print(f"genes {len(keys)} (+ {n_plus}, - {len(keys) - n_plus}); per-gene intervals {len(rows)}; "
          f"union intervals {len(rebuilt)} = {args.union_bed}")


if __name__ == "__main__":
    main()
