#!/usr/bin/env python3
"""Independent check of variants.tsv for one library: index the BAM in a private
directory, then count alleles with pysam's pileup engine (not the CIGAR walk used by
extract_reads_one_library.py).

Usage: indep_variant_check.py BAM INDEX_DIR VARIANTS_TSV OUT_TSV
"""
import csv
import os
import sys
from collections import Counter

import pysam

bam_path, idx_dir, var_path, out_path = sys.argv[1:5]
if os.path.exists(out_path):
    sys.exit("refusing to overwrite " + out_path)
os.makedirs(idx_dir, exist_ok=True)
idx = os.path.join(idx_dir, os.path.basename(bam_path) + ".bai")
if not os.path.exists(idx):
    pysam.index("-b", bam_path, idx)
af = pysam.AlignmentFile(bam_path, "rb", index_filename=idx)

with open(var_path) as fh:
    variants = list(csv.DictReader(fh, delimiter="\t"))

BASES = "ACGT"
rows = []


def count_locus(chrom, pos1, accept_mq3, ins_anchor=False):
    """Return Counter of alleles over distinct fragments, plus alignment-level MAPQ census."""
    c = Counter()
    seen = set()
    census = Counter()
    for col in af.pileup(chrom, pos1 - 1, pos1, truncate=True, stepper="nofilter", min_base_quality=0,
                         max_depth=10_000_000, ignore_overlaps=False, ignore_orphans=False):
        for pr in col.pileups:
            a = pr.alignment
            census[(a.mapping_quality, bool(a.flag & 0x100))] += 1
            if a.flag & 0x700:
                continue
            mq = a.mapping_quality
            if not (mq == 255 or (accept_mq3 and mq == 3)):
                continue
            if a.query_name in seen:
                continue
            if ins_anchor:
                if pr.is_del or pr.is_refskip or pr.query_position is None:
                    continue
                seen.add(a.query_name)
                c["ins" if pr.indel > 0 else ("del_after" if pr.indel < 0 else "noins")] += 1
                continue
            if pr.is_del:
                seen.add(a.query_name)
                c["del"] += 1
                continue
            if pr.is_refskip or pr.query_position is None:
                continue
            if a.query_qualities[pr.query_position] < 20:
                continue
            seen.add(a.query_name)
            c[a.query_sequence[pr.query_position].upper()] += 1
    return c, census


for v in variants:
    chrom, pos = v["chrom"], int(v["pos"])
    if v["kind"] == "snv":
        multi = v["alt_contig"] != "NA"
        c, census = count_locus(chrom, pos, multi)
        tot = Counter(c)
        if multi:
            c2, census2 = count_locus(v["alt_contig"], int(v["alt_pos"]), True)
            tot.update(c2)
            alt_ct = sum(c2[b] for b in BASES)
        else:
            alt_ct = 0
        ref, alt = v["ref"], v["alt"]
        n_ref, n_alt = tot[ref], tot[alt]
        n_oth = sum(tot[b] for b in BASES if b not in (ref, alt))
        rows.append([v["rsid"], n_ref, n_alt, n_oth, tot["del"], n_ref + n_alt + n_oth,
                     alt_ct, ";".join(f"mq{k[0]}{'s' if k[1] else ''}={n}" for k, n in sorted(census.items()))])
    else:
        # insertion: compare query span with reference span between flank bases (CIGAR-length based)
        L, R = pos - 3, pos + 4
        seen, c, census = set(), Counter(), Counter()
        for a in af.fetch(chrom, pos - 6, pos + 6):
            census[(a.mapping_quality, bool(a.flag & 0x100))] += 1
            if a.flag & 0x700 or a.mapping_quality != 255 or a.query_name in seen:
                continue
            d = {r_: q_ for q_, r_ in a.get_aligned_pairs() if r_ is not None and q_ is not None}
            if (L - 1) not in d or (R - 1) not in d:
                continue
            seen.add(a.query_name)
            il = (d[R - 1] - d[L - 1]) - (R - L)
            c["ins1" if il == 1 else ("noindel" if il == 0 else "other_indel")] += 1
        rows.append([v["rsid"], c["noindel"], c["ins1"], c["other_indel"], 0, c["noindel"] + c["ins1"] + c["other_indel"], 0,
                     ";".join(f"mq{k[0]}{'s' if k[1] else ''}={n}" for k, n in sorted(census.items()))])

with open(out_path, "w") as fo:
    fo.write("library\trsid\tn_ref\tn_alt\tn_other\tn_del\tdepth\tdepth_on_alt_contig\tmapq_census_at_locus\n")
    lib = os.path.basename(bam_path).split(".")[0]
    for r in rows:
        fo.write(lib + "\t" + "\t".join(map(str, r)) + "\n")
print("done", lib)
