#!/usr/bin/env python3
"""H25: random duplicate removal and fragment-level allele counts at listed SNVs.

Input is the coordinate-sorted, indexed BAM of pairs kept by wasp_remap.py check.
Duplicates: fragments with the same contig, R1 start, R2 start and R1 strand; one is kept,
chosen by a uniform draw from numpy default_rng(seed) in BAM order (cannot depend on allele).
Counting: one count per kept fragment; base quality >= MIN_BQ; if both mates cover the site
with passing bases that disagree, the fragment is not counted at that site.
Writes per-site ref, alt, other counts. Computes no imbalance statistic.
"""
import argparse
import json
from collections import defaultdict

import numpy as np
import pysam

MIN_BQ = 20


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--seed", type=int, default=20260923)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    bam = pysam.AlignmentFile(a.bam, "rb")

    best = {}
    n_frag = 0
    for r in bam.fetch(until_eof=True):
        if not r.is_read1 or r.is_unmapped:
            continue
        n_frag += 1
        key = (r.reference_id, r.reference_start, r.next_reference_start, r.is_reverse)
        u = rng.random()
        cur = best.get(key)
        if cur is None or u < cur[0]:
            best[key] = (u, r.query_name)
    keep = {v[1] for v in best.values()}
    del best

    sites = []
    with open(a.sites) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            sites.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))

    with open(a.out, "w") as out:
        out.write("sample\tchrom\tpos\tref\talt\tcls\tref_n\talt_n\tother_n\tdiscordant_mates\n")
        for s in sites:
            chrom, p0 = s["chrom"], int(s["pos"]) - 1
            frag = defaultdict(list)
            for r in bam.fetch(chrom, p0, p0 + 1):
                if r.query_name not in keep:
                    continue
                for qp, rp in r.get_aligned_pairs(matches_only=True):
                    if rp == p0:
                        if r.query_qualities[qp] >= MIN_BQ:
                            frag[r.query_name].append(r.query_sequence[qp])
                        break
            nr = na = no = nd = 0
            for bases in frag.values():
                if not bases:
                    continue
                if len(set(bases)) > 1:
                    nd += 1
                    continue
                b = bases[0]
                if b == s["ref"]:
                    nr += 1
                elif b == s["alt"]:
                    na += 1
                else:
                    no += 1
            out.write(f"{a.sample}\t{chrom}\t{s['pos']}\t{s['ref']}\t{s['alt']}\t{s.get('cls', 'NA')}\t{nr}\t{na}\t{no}\t{nd}\n")
    json.dump(dict(sample=a.sample, fragments_in=n_frag, fragments_after_dedup=len(keep), seed=a.seed),
              open(a.out + ".dedup.json", "w"), indent=1)


if __name__ == "__main__":
    main()
