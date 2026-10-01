"""Tag-gene intervals and the library manifest for the duplicate-read summary (B-Dup).

Inputs
  --tags       t3/tags_final.tsv. Only the gene_key column is read; tag positions and
               alleles are not used.
  --gtf        GENCODE v49 primary-chromosome GTF (same coordinates as the 706-contig
               BAMs on chr1-22).
  --bam-root   a1-full .../bai/<cohort>/<run>.bam (symlinks to the STAR BAMs, each with
               a .bam.bai beside it).
  --crosswalk  identity/sample_to_individual.tsv (run, cohort, individual_id).
Outputs (in --out)
  tag_genes.tsv        one row per tag gene: gene_key, gene_id, gene_name, chrom, start,
                       end (1-based, inclusive), n_exon_intervals, exon_bp.
  tag_gene_exons.bed   union of all annotated exons of the tag genes, merged
                       (0-based, half-open). Reads are counted here.
  tag_gene_bodies.bed  gene start..end spans, merged. Reads are selected here before
                       duplicate marking so both mates of a spliced pair are present.
  libraries.tsv        task, cohort, run, individual_id, bam, star_log (one row per library).
Rules
  Genes are matched on the unversioned Ensembl id. A tag gene absent from the GTF, a
  library without a BAM, index or STAR Log.final.out stops the script.
"""
import argparse
import os
import re
import sys
from collections import defaultdict

import pandas as pd

CHROM_ORDER = {f"chr{i}": i for i in range(1, 23)} | {"chrX": 23, "chrY": 24, "chrM": 25}
ATTR = re.compile(r'(\S+) "([^"]*)"')


def merge(intervals):
    """Merge 0-based half-open (start, end) intervals that overlap or touch."""
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def merged_bp(by_chrom):
    return sum(end - start for chrom in by_chrom for start, end in merge(by_chrom[chrom]))


def read_gtf(path, keys):
    genes, exons = {}, defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("gene", "exon"):
                continue
            attrs = dict(ATTR.findall(f[8]))
            key = attrs["gene_id"].split(".")[0]
            if key not in keys:
                continue
            if f[2] == "gene":
                genes[key] = (attrs["gene_id"], attrs.get("gene_name", ""), f[0], int(f[3]), int(f[4]))
            else:
                exons[key].append((f[0], int(f[3]) - 1, int(f[4])))
    return genes, exons


def write_bed(by_chrom, path):
    with open(path, "w") as out:
        for chrom in sorted(by_chrom, key=lambda c: CHROM_ORDER.get(c, 99)):
            for start, end in merge(by_chrom[chrom]):
                out.write(f"{chrom}\t{start}\t{end}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tags", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--bam-root", required=True)
    ap.add_argument("--crosswalk", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    keys = set(pd.read_csv(args.tags, sep="\t", usecols=["gene_key"])["gene_key"])
    genes, exons = read_gtf(args.gtf, keys)
    missing = sorted(keys - set(genes))
    if missing:
        sys.exit(f"{len(missing)} tag genes not in the GTF, e.g. {missing[:5]}")

    exon_by_chrom, body_by_chrom, rows = defaultdict(list), defaultdict(list), []
    for key in sorted(keys):
        gene_id, name, chrom, start, end = genes[key]
        gene_exons = merge([(s, e) for c, s, e in exons[key] if c == chrom])
        exon_by_chrom[chrom].extend(map(tuple, gene_exons))
        body_by_chrom[chrom].append((start - 1, end))
        rows.append((key, gene_id, name, chrom, start, end, len(gene_exons), sum(e - s for s, e in gene_exons)))
    pd.DataFrame(rows, columns=["gene_key", "gene_id", "gene_name", "chrom", "start", "end",
                                "n_exon_intervals", "exon_bp"]).to_csv(
        os.path.join(args.out, "tag_genes.tsv"), sep="\t", index=False)
    write_bed(exon_by_chrom, os.path.join(args.out, "tag_gene_exons.bed"))
    write_bed(body_by_chrom, os.path.join(args.out, "tag_gene_bodies.bed"))

    # keep_default_na=False: pandas would otherwise turn empty or "None" ids into NaN
    xw = pd.read_csv(args.crosswalk, sep="\t", dtype=str, keep_default_na=False)
    xw = xw[["cohort", "run", "individual_id"]].sort_values(["cohort", "run"]).reset_index(drop=True)
    xw["bam"] = [os.path.join(args.bam_root, c, f"{r}.bam") for c, r in zip(xw["cohort"], xw["run"])]
    logs, problems = [], []
    for bam in xw["bam"]:
        target = os.path.realpath(bam)
        log = target.replace(".Aligned.sortedByCoord.out.bam", ".Log.final.out")
        logs.append(log)
        for need in (bam, bam + ".bai", log):
            if not os.path.isfile(need):
                problems.append(need)
    if problems:
        sys.exit(f"{len(problems)} missing BAM/index/log files, e.g. {problems[:5]}")
    xw["star_log"] = logs
    xw.insert(0, "task", range(1, len(xw) + 1))
    xw.to_csv(os.path.join(args.out, "libraries.tsv"), sep="\t", index=False)

    print(f"tag genes: {len(keys)}; exon bp (merged): {merged_bp(exon_by_chrom)}; "
          f"gene-body bp (merged): {merged_bp(body_by_chrom)}; libraries: {len(xw)}")


if __name__ == "__main__":
    main()
