"""Independent re-derivation of the B-Dup duplicate counts for a few libraries.

Reads the same gene-body subset that 03_dup_library.sbatch builds (primary, MAPQ 255,
reads overlapping tag_gene_bodies.bed, each read once), but only the columns QNAME,
FLAG, RNAME, POS and CIGAR (no sequence, so no allele is ever seen), and recomputes:
  single-end   duplicates = reads - distinct (chrom, unclipped 5' position, strand)
  paired-end   duplicate reads = 2 * (pairs - distinct unordered pairs of mate ends),
               each mate end = (chrom, unclipped 5' position, strand). Reads whose mate
               is outside the subset are counted (mate_absent_reads) but not re-derived
               here. markdup treats them as single reads and compares each one's end
               with the other single reads and with both ends of every complete pair;
               R-DATA reproduced markdup's DUPLICATE SINGLE exactly with that rule
               (SRR9883072 7,442; SRR21623070 10,848). This script does not check
               exon-level duplicate counts; R-DATA bounded them independently and the
               step-3 values fall inside the bounds (SRR17442829, SRR9883072, SRR21623070).
  exon reads   reads whose aligned span (POS .. reference end, introns included, as
               samtools view -L uses it) overlaps a merged tag-gene exon.
Then prints these next to samtools markdup's DUPLICATE SINGLE / DUPLICATE PAIR and the
exon_reads from 03's per-library row.

Inputs : --libraries libraries.tsv; --regions dir with tag_gene_bodies.bed and
         tag_gene_exons.bed; --results dir with per_library/<run>.tsv and
         markdup_stats/<run>.txt; --runs run accessions to check.
Output : --out TSV (counts only); refuses to overwrite.
"""
import argparse
import bisect
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

import pandas as pd

SAMTOOLS = "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/samtools"
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
REF_OPS = set("MDN=X")


def read_bed(path):
    by_chrom = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            chrom, start, end = line.split("\t")[:3]
            by_chrom[chrom].append((int(start), int(end)))
    return {c: (sorted(s for s, _ in iv), sorted(e for _, e in iv)) for c, iv in by_chrom.items()}


def overlaps(exons, chrom, start, end):
    """True if 0-based half-open [start, end) overlaps a merged interval on chrom."""
    if chrom not in exons:
        return False
    starts, ends = exons[chrom]
    i = bisect.bisect_left(starts, end) - 1  # last interval starting before end
    return i >= 0 and ends[i] > start


def read_end(pos, cigar):
    """Leading soft clip, reference span length, trailing soft clip."""
    ops = CIGAR.findall(cigar)
    lead = int(ops[0][0]) if ops[0][1] == "S" else 0
    trail = int(ops[-1][0]) if ops[-1][1] == "S" else 0
    span = sum(int(n) for n, op in ops if op in REF_OPS)
    return lead, span, trail


def check(bam, bodies, exons):
    cmd = [SAMTOOLS, "view", "-M", "-L", bodies, "-q", "255", "-F", "0x904", bam]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    singles, mates = Counter(), {}
    pair_keys = Counter()
    n_reads = n_exon = n_single = 0
    for line in proc.stdout:
        qname, flag, chrom, pos, _, cigar = line.split("\t", 6)[:6]
        flag, pos = int(flag), int(pos)
        lead, span, trail = read_end(pos, cigar)
        start0, end0 = pos - 1, pos - 1 + span
        n_reads += 1
        n_exon += overlaps(exons, chrom, start0, end0)
        reverse = bool(flag & 16)
        end_key = (chrom, end0 + trail, "-") if reverse else (chrom, start0 - lead, "+")
        if not flag & 1:
            n_single += 1
            singles[end_key] += 1
            continue
        other = mates.pop(qname, None)
        if other is None:
            mates[qname] = end_key
        else:
            pair_keys[tuple(sorted((end_key, other)))] += 1
    if proc.wait() != 0:
        sys.exit(f"samtools failed on {bam}")
    n_pairs = sum(pair_keys.values())
    return {
        "reads": n_reads,
        "exon_reads_rederived": n_exon,
        "single_reads": n_single,
        "single_dup_rederived": n_single - len(singles),
        "pairs": n_pairs,
        "pair_dup_reads_rederived": 2 * (n_pairs - len(pair_keys)),
        "mate_absent_reads": len(mates),
    }


def read_markdup(path):
    values = {}
    with open(path) as fh:
        for line in fh:
            key, _, value = line.partition(":")
            if key in ("EXAMINED", "DUPLICATE SINGLE", "DUPLICATE PAIR"):
                values[key] = int(value)
    return values


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--libraries", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--runs", required=True, nargs="+")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if os.path.exists(args.out):
        sys.exit(f"refusing to overwrite {args.out}")

    libs = pd.read_csv(args.libraries, sep="\t", dtype=str, keep_default_na=False).set_index("run")
    exons = read_bed(f"{args.regions}/tag_gene_exons.bed")
    rows = []
    for run in args.runs:
        got = check(libs.loc[run, "bam"], f"{args.regions}/tag_gene_bodies.bed", exons)
        md = read_markdup(f"{args.results}/markdup_stats/{run}.txt")
        row03 = pd.read_csv(f"{args.results}/per_library/{run}.tsv", sep="\t").iloc[0]
        rows.append({"run": run, "cohort": libs.loc[run, "cohort"], **got,
                     "markdup_examined": md["EXAMINED"],
                     "markdup_dup_single": md["DUPLICATE SINGLE"],
                     "markdup_dup_pair": md["DUPLICATE PAIR"],
                     "exon_reads_03": int(row03["exon_reads"])})
    table = pd.DataFrame(rows)
    table.to_csv(args.out, sep="\t", index=False)
    print(table.T.to_string())


if __name__ == "__main__":
    main()
