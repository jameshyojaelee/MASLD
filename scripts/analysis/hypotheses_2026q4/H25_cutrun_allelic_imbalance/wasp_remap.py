#!/usr/bin/env python3
"""WASP-equivalent mapping-bias filter for paired-end reads (H25).

make:   read a name-sorted BAM of proper pairs (MAPQ >= 30 already applied), find pairs whose
        bases at listed SNVs equal REF or ALT, write every allele-swapped version of each such
        pair as FASTQ (all 2^k combinations except the original; pairs with k > MAX_K are
        dropped), and record each original pair's expected mate starts.
check:  read the realigned versions (Bowtie2 output, same settings) and keep an original pair
        only if every version realigned as a proper pair, MAPQ >= 30 on both mates, to the
        same contig and the same start of both mates. Writes the kept original pairs
        (pairs with no swappable SNV are kept without a test) to a BAM.
"""
import argparse
import bisect
import json
from collections import defaultdict

import pysam

MAX_K = 6
MIN_MAPQ = 30
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def load_sites(path):
    pos = defaultdict(list)
    allele = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci, pi, ri, ai = (hdr.index(c) for c in ("chrom", "pos", "ref", "alt"))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            p0 = int(f[pi]) - 1
            pos[f[ci]].append(p0)
            allele[(f[ci], p0)] = (f[ri], f[ai])
    for c in pos:
        pos[c].sort()
    return pos, allele


def pairs(bam):
    """Yield (r1, r2) from a name-sorted BAM; groups that are not exactly one R1 + one R2 are skipped."""
    group = []
    name = None
    for r in bam.fetch(until_eof=True):
        if r.query_name != name:
            if len(group) == 2:
                yield group
            group = []
            name = r.query_name
        group.append(r)
    if len(group) == 2:
        yield group


def order(group):
    a, b = group
    if a.is_read1 and b.is_read2:
        return a, b
    if b.is_read1 and a.is_read2:
        return b, a
    return None


def snv_hits(read, pos, allele):
    """Return {ref_pos0: (qpos, base)} for listed SNVs where the read base is REF or ALT."""
    chrom = read.reference_name
    plist = pos.get(chrom)
    if not plist:
        return {}, 0
    lo = bisect.bisect_left(plist, read.reference_start)
    hi = bisect.bisect_left(plist, read.reference_end)
    if lo == hi:
        return {}, 0
    want = set(plist[lo:hi])
    seq = read.query_sequence
    out = {}
    n_over = 0
    for qp, rp in read.get_aligned_pairs(matches_only=True):
        if rp in want:
            n_over += 1
            b = seq[qp]
            if b in allele[(chrom, rp)]:
                out[rp] = (qp, b)
    return out, n_over


def make(a):
    pos, allele = load_sites(a.sites)
    bam = pysam.AlignmentFile(a.bam, "rb")
    f1 = open(a.out_prefix + ".remap_1.fq", "w")
    f2 = open(a.out_prefix + ".remap_2.fq", "w")
    exp = open(a.out_prefix + ".expected.tsv", "w")
    noswap = open(a.out_prefix + ".noswap_names.txt", "w")
    st = defaultdict(int)
    for g in pairs(bam):
        o = order(g)
        if o is None:
            st["bad_pair"] += 1
            continue
        r1, r2 = o
        if r1.reference_name != r2.reference_name:
            st["diff_contig"] += 1
            continue
        st["pairs"] += 1
        h1, n1 = snv_hits(r1, pos, allele)
        h2, n2 = snv_hits(r2, pos, allele)
        snvs = sorted(set(h1) | set(h2))
        if not snvs:
            if n1 or n2:
                st["overlap_other_base_only"] += 1
            noswap.write(r1.query_name + "\n")
            st["noswap"] += 1
            continue
        k = len(snvs)
        if k > MAX_K:
            st["drop_k_gt_max"] += 1
            continue
        s1, s2 = r1.query_sequence, r2.query_sequence
        q1 = pysam.qualities_to_qualitystring(r1.query_qualities)
        q2 = pysam.qualities_to_qualitystring(r2.query_qualities)
        nver = 0
        chrom = r1.reference_name
        for bits in range(2 ** k):
            l1, l2 = list(s1), list(s2)
            for j, rp in enumerate(snvs):
                base = allele[(chrom, rp)][(bits >> j) & 1]
                if rp in h1:
                    l1[h1[rp][0]] = base
                if rp in h2:
                    l2[h2[rp][0]] = base
            n1s, n2s = "".join(l1), "".join(l2)
            if n1s == s1 and n2s == s2:
                continue
            o1 = revcomp(n1s) if r1.is_reverse else n1s
            o2 = revcomp(n2s) if r2.is_reverse else n2s
            oq1 = q1[::-1] if r1.is_reverse else q1
            oq2 = q2[::-1] if r2.is_reverse else q2
            vname = f"{r1.query_name}.v{nver}"
            f1.write(f"@{vname}\n{o1}\n+\n{oq1}\n")
            f2.write(f"@{vname}\n{o2}\n+\n{oq2}\n")
            nver += 1
        exp.write(f"{r1.query_name}\t{nver}\t{chrom}\t{r1.reference_start}\t{r2.reference_start}\n")
        st["test"] += 1
        st["versions"] += nver
    for fh in (f1, f2, exp, noswap):
        fh.close()
    json.dump(dict(st), open(a.out_prefix + ".make_stats.json", "w"), indent=1)
    print(dict(st))


def check(a):
    expected = {}
    with open(a.out_prefix + ".expected.tsv") as fh:
        for line in fh:
            n, nv, c, s1, s2 = line.rstrip("\n").split("\t")
            expected[n] = [int(nv), c, int(s1), int(s2), 0]
    rb = pysam.AlignmentFile(a.remapped, "r" if a.remapped.endswith(".sam") else "rb")
    for g in pairs(rb):
        o = order(g)
        if o is None:
            continue
        r1, r2 = o
        orig = r1.query_name.rsplit(".v", 1)[0]
        e = expected.get(orig)
        if e is None:
            continue
        ok = (r1.is_proper_pair and r2.is_proper_pair and not r1.is_secondary and not r2.is_secondary
              and not r1.is_supplementary and not r2.is_supplementary
              and r1.mapping_quality >= MIN_MAPQ and r2.mapping_quality >= MIN_MAPQ
              and r1.reference_name == e[1] and r2.reference_name == e[1]
              and r1.reference_start == e[2] and r2.reference_start == e[3])
        if ok:
            e[4] += 1
    keep = set(n.strip() for n in open(a.out_prefix + ".noswap_names.txt"))
    n_pass = n_fail = 0
    for n, e in expected.items():
        if e[0] > 0 and e[4] == e[0]:
            keep.add(n)
            n_pass += 1
        else:
            n_fail += 1
    src = pysam.AlignmentFile(a.bam, "rb")
    out = pysam.AlignmentFile(a.out_bam, "wb", template=src)
    n_w = 0
    for g in pairs(src):
        if g[0].query_name in keep:
            for r in g:
                out.write(r)
            n_w += 1
    out.close()
    st = dict(tested=len(expected), tested_pass=n_pass, tested_fail=n_fail, kept_pairs=n_w)
    json.dump(st, open(a.out_prefix + ".check_stats.json", "w"), indent=1)
    print(st)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["make", "check"])
    ap.add_argument("--bam", required=True, help="name-sorted filtered BAM")
    ap.add_argument("--sites", help="sites TSV (chrom, pos 1-based, ref, alt)")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--remapped", help="check: realigned versions (SAM or BAM)")
    ap.add_argument("--out-bam", help="check: kept original pairs")
    a = ap.parse_args()
    make(a) if a.mode == "make" else check(a)


if __name__ == "__main__":
    main()
