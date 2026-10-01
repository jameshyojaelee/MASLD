"""Per-library x per-tag-gene duplicate counts under one read-end rule (B-Dup step 7).

Why: the step-3 column exon_dup_fraction pools all tag-exon reads, so most of its weight
sits in the most-read genes (positional saturation, which depends on expression and
depth), and it uses the pair rule for paired-end and the one-end rule for single-end
libraries. This step writes the per-gene counts from which a covariate defined the same
way for all 1,466 libraries, and matched for expression and depth, is computed (08).

Reads (the same filters as step 3)
  samtools view -M -L tag_gene_bodies.bed -q 255 -F 0x904 <bam> | cut -f 1-6
  i.e. primary, uniquely mapped (STAR MAPQ 255), not unmapped/secondary/supplementary,
  overlapping a tag-gene body, each read once. Only QNAME, FLAG, RNAME, POS, MAPQ and
  CIGAR reach this script; the read sequence and base qualities are never read.
One-end reads
  Read 1 (FLAG 0x40) of paired-end libraries; every read of single-end libraries.
  Read 2 is used only to form pair keys.
End key
  (RNAME, unclipped 5' position, strand): POS minus a leading soft clip for a forward
  read; reference end plus a trailing soft clip for a reverse read. This is the key
  samtools markdup uses for one read end (05_verify_markdup.py reproduced markdup's
  single-end duplicate counts exactly with it).
Gene assignment
  A one-end read is assigned to gene g when its aligned blocks (CIGAR M, = and X
  segments; introns, deletions and soft clips are not blocks) overlap the merged exons
  of g and of no other tag gene on the library's transcribed strand. The library strand
  is called from the reads: among one-end reads whose blocks overlap the exons of exactly
  one tag gene (either strand), the share aligned antisense to that gene is
  >= 0.9 -> "antisense" (read 1 opposite to the gene), <= 0.1 -> "sense",
  otherwise "unstranded" (strand ignored; a read is assigned if it overlaps one gene).
Per gene g (written to per_gene/<run>.tsv.gz)
  reads1          assigned one-end reads
  distinct1       distinct end keys among them; duplicates = reads1 - distinct1
  pairs           (paired-end) assigned read-1s whose mate is in the subset
  distinct_pairs  distinct unordered pairs of end keys (read 1, read 2) among those pairs;
                  pair duplicates = pairs - distinct_pairs (markdup's pair rule)
  mate_absent1    (paired-end) assigned read-1s whose mate is not in the subset
  reads1_sub, distinct1_sub
                  the same one-end counts on a random subsample of exactly N_SUB = 50,000
                  assigned one-end reads per library (numpy default_rng(20260929), without
                  replacement); NA when a library has fewer than 50,000 assigned reads.
Per library (per_library/<run>.tsv)
  layout, strand call and its counts, read totals, the header check with --no-PG
  (program names in @PG, and whether any STAR command line in @PG/@CO contains
  --bamRemoveDuplicatesType), and span_* check columns: the one-end rule over one-end
  reads whose span (POS to reference end, introns included, as samtools -L uses it)
  overlaps the merged tag exons on either strand, with no gene assignment. span_* only
  reproduces the R-DATA reviewer's recount; it is not a covariate.
Counts only; no allele, genotype or stage information is read or written.
"""
import argparse
import bisect
import gzip
import os
import re
import subprocess
import sys
from array import array
from collections import defaultdict

import numpy as np

SAMTOOLS = "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/samtools"
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
N_SUB = 50_000
SEED = 20260929
POS_OFFSET = 1 << 20  # keeps unclipped 5' positions left of base 1 non-negative
STRAND_CUT = 0.9


def read_gene_exons(path):
    """Disjoint segments per chromosome, each with the frozenset of genes whose exons cover it."""
    genes, index, by_chrom = [], {}, defaultdict(list)
    with open(path) as fh:
        for line in fh:
            chrom, start, end, key, _, strand = line.rstrip("\n").split("\t")
            if key not in index:
                index[key] = len(genes)
                genes.append((key, strand))
            by_chrom[chrom].append((int(start), int(end), index[key]))
    segments = {}
    for chrom, intervals in by_chrom.items():
        events = defaultdict(list)
        for start, end, g in intervals:
            events[start].append((g, 1))
            events[end].append((g, -1))
        points = sorted(events)
        starts, ends, sets, active = [], [], [], defaultdict(int)
        for k, point in enumerate(points[:-1]):
            for g, d in events[point]:
                active[g] += d
            covering = frozenset(g for g, n in active.items() if n > 0)
            if covering:
                starts.append(point)
                ends.append(points[k + 1])
                sets.append(covering)
        segments[chrom] = (starts, ends, sets)
    return genes, segments


def unique_or_code(candidates):
    """Gene index if exactly one candidate, -1 if none, -2 if several."""
    if len(candidates) == 1:
        return next(iter(candidates))
    return -1 if not candidates else -2


def header_check(bam):
    text = subprocess.run([SAMTOOLS, "view", "-H", "--no-PG", bam], check=True,
                          capture_output=True, text=True).stdout
    programs, dedup = set(), 0
    for line in text.splitlines():
        if line.startswith("@PG"):
            programs.update(f[3:] for f in line.split("\t") if f.startswith("PN:"))
        if line.startswith(("@PG", "@CO")) and "bamRemoveDuplicates" in line:
            dedup = 1
    return ",".join(sorted(programs)) or "none", dedup


def stream(bam, bodies, genes, segments):
    view = subprocess.Popen([SAMTOOLS, "view", "-M", "-L", bodies, "-q", "255", "-F", "0x904", bam],
                            stdout=subprocess.PIPE)
    cut = subprocess.Popen(["cut", "-f", "1-6"], stdin=view.stdout, stdout=subprocess.PIPE, text=True)
    view.stdout.close()
    strand_of = [s for _, s in genes]
    chrom_id = {}
    key1, g_sense, g_anti, g_any = array("q"), array("i"), array("i"), array("i")
    span_hit, mate = array("b"), array("q")
    pending = {}  # qname -> -(index + 1) for a stored read 1, None for an unstored read 1, key >= 0 for a read 2
    n = {"reads_streamed": 0, "one_end_reads": 0, "r2_reads": 0, "paired_flag_reads": 0,
         "strand_informative": 0, "strand_antisense": 0}
    empty = ((), (), ())
    for line in cut.stdout:
        qname, flag, chrom, pos, _, cigar = line.rstrip("\n").split("\t")
        flag, pos = int(flag), int(pos)
        n["reads_streamed"] += 1
        ops = CIGAR.findall(cigar)
        lead = int(ops[0][0]) if ops[0][1] == "S" else 0
        trail = int(ops[-1][0]) if ops[-1][1] == "S" else 0
        ref, blocks = pos - 1, []
        for length, op in ops:
            length = int(length)
            if op in "M=X":
                blocks.append((ref, ref + length))
                ref += length
            elif op in "DN":
                ref += length
        reverse = flag & 16
        pos5 = ref + trail if reverse else pos - 1 - lead
        cid = chrom_id.setdefault(chrom, len(chrom_id))
        key = (cid << 36) | ((pos5 + POS_OFFSET) << 1) | (1 if reverse else 0)
        paired = flag & 1
        if paired:
            n["paired_flag_reads"] += 1
        if paired and not flag & 0x40:  # read 2: pair key only
            n["r2_reads"] += 1
            other = pending.pop(qname, KeyError)
            if other is KeyError:
                pending[qname] = key
            elif other is not None:
                mate[-other - 1] = key
            continue
        n["one_end_reads"] += 1
        starts, ends, sets = segments.get(chrom, empty)
        hit = set()
        for bs, be in blocks:
            i = bisect.bisect_right(ends, bs)
            while i < len(starts) and starts[i] < be:
                hit |= sets[i]
                i += 1
        i = bisect.bisect_right(ends, pos - 1)
        on_span = i < len(starts) and starts[i] < ref
        if not hit and not on_span:
            if paired:
                if pending.pop(qname, KeyError) is KeyError:
                    pending[qname] = None
            continue
        read_strand = "-" if reverse else "+"
        if len(hit) == 1:
            n["strand_informative"] += 1
            n["strand_antisense"] += strand_of[next(iter(hit))] != read_strand
        key1.append(key)
        g_any.append(unique_or_code(hit))
        g_sense.append(unique_or_code({g for g in hit if strand_of[g] == read_strand}))
        g_anti.append(unique_or_code({g for g in hit if strand_of[g] != read_strand}))
        span_hit.append(1 if on_span else 0)
        if not paired:
            mate.append(-1)
            continue
        other = pending.pop(qname, KeyError)
        if other is KeyError:
            mate.append(-2)  # mate not seen yet; stays -2 if it never comes
            pending[qname] = -len(key1)
        else:
            mate.append(other)
    cut.wait()
    if view.wait() != 0 or cut.returncode != 0:
        sys.exit(f"samtools view | cut failed on {bam}")
    n["r2_without_read1"] = sum(1 for v in pending.values() if v is not None and v >= 0)
    arrays = {name: np.frombuffer(a, dtype=dt) if len(a) else np.zeros(0, dtype=dt)
              for name, a, dt in (("key1", key1, np.int64), ("g_sense", g_sense, np.int32),
                                  ("g_anti", g_anti, np.int32), ("g_any", g_any, np.int32),
                                  ("span", span_hit, np.int8), ("mate", mate, np.int64))}
    return n, arrays


def distinct_per_gene(genes_idx, *keys, n_genes):
    """Number of distinct (gene, key...) combinations per gene."""
    if len(genes_idx) == 0:
        return np.zeros(n_genes, dtype=np.int64)
    order = np.lexsort(tuple(reversed(keys)) + (genes_idx,))
    g = genes_idx[order]
    new = np.ones(len(g), dtype=bool)
    new[1:] = g[1:] != g[:-1]
    for k in keys:
        ks = k[order]
        new[1:] |= ks[1:] != ks[:-1]
    return np.bincount(g[new], minlength=n_genes)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--regions", required=True, help="dir with tag_gene_bodies.bed, tag_gene_exons_by_gene.bed")
    ap.add_argument("--out", required=True, help="restricted output dir")
    args = ap.parse_args()

    row_path = os.path.join(args.out, "per_library", f"{args.run}.tsv")
    gene_path = os.path.join(args.out, "per_gene", f"{args.run}.tsv.gz")
    if os.path.exists(row_path):
        print(f"{args.run} already done")
        return
    os.makedirs(os.path.dirname(row_path), exist_ok=True)
    os.makedirs(os.path.dirname(gene_path), exist_ok=True)

    genes, segments = read_gene_exons(os.path.join(args.regions, "tag_gene_exons_by_gene.bed"))
    n_genes = len(genes)
    programs, dedup_option = header_check(args.bam)
    n, a = stream(args.bam, os.path.join(args.regions, "tag_gene_bodies.bed"), genes, segments)

    layout = "paired" if n["paired_flag_reads"] else "single"
    anti_frac = n["strand_antisense"] / n["strand_informative"] if n["strand_informative"] else float("nan")
    if anti_frac >= STRAND_CUT:
        strand_call, g = "antisense", a["g_anti"]
    elif anti_frac <= 1 - STRAND_CUT:
        strand_call, g = "sense", a["g_sense"]
    else:
        strand_call, g = "unstranded", a["g_any"]

    assigned = g >= 0
    ga, k1, mate = g[assigned], a["key1"][assigned], a["mate"][assigned]
    reads1 = np.bincount(ga, minlength=n_genes)
    distinct1 = distinct_per_gene(ga, k1, n_genes=n_genes)
    if layout == "paired":
        has_mate = mate >= 0
        lo, hi = np.minimum(k1, mate)[has_mate], np.maximum(k1, mate)[has_mate]
        pairs = np.bincount(ga[has_mate], minlength=n_genes)
        distinct_pairs = distinct_per_gene(ga[has_mate], lo, hi, n_genes=n_genes)
        mate_absent1 = np.bincount(ga[mate == -2], minlength=n_genes)
    if len(ga) >= N_SUB:
        pick = np.random.default_rng(SEED).choice(len(ga), size=N_SUB, replace=False)
        reads1_sub = np.bincount(ga[pick], minlength=n_genes)
        distinct1_sub = distinct_per_gene(ga[pick], k1[pick], n_genes=n_genes)

    def col(values, i):
        return "NA" if values is None else str(int(values[i]))

    with gzip.open(gene_path + ".part", "wt") as out:
        out.write("gene_key\treads1\tdistinct1\tpairs\tdistinct_pairs\tmate_absent1\treads1_sub\tdistinct1_sub\n")
        paired_cols = (pairs, distinct_pairs, mate_absent1) if layout == "paired" else (None, None, None)
        sub_cols = (reads1_sub, distinct1_sub) if len(ga) >= N_SUB else (None, None)
        for i, (key, _) in enumerate(genes):
            out.write("\t".join([key, col(reads1, i), col(distinct1, i)]
                                + [col(c, i) for c in paired_cols] + [col(c, i) for c in sub_cols]) + "\n")
    os.replace(gene_path + ".part", gene_path)

    span = a["span"] == 1
    span_keys = a["key1"][span]
    row = {
        "cohort": args.cohort, "run": args.run, "layout": layout,
        "pg_programs_noPG": programs, "star_dedup_option_in_header": dedup_option,
        **n,
        "strand_antisense_fraction": f"{anti_frac:.6f}", "strand_call": strand_call,
        "one_end_assigned": int(assigned.sum()),
        "one_end_multi_gene": int((g == -2).sum()),
        "one_end_sub_drawn": N_SUB if len(ga) >= N_SUB else 0,
        "span_one_end_reads": int(span.sum()),
        "span_one_end_distinct": int(len(np.unique(span_keys))),
        "seed": SEED,
    }
    with open(row_path + ".part", "w") as out:
        out.write("\t".join(row) + "\n" + "\t".join(str(v) for v in row.values()) + "\n")
    os.replace(row_path + ".part", row_path)
    print(row)


if __name__ == "__main__":
    main()
