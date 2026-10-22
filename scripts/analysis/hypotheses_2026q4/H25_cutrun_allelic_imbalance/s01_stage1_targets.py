#!/usr/bin/env python3
"""H25 stage 1: build the target-site table and the k-mer bait FASTA.

Targets = tier-1/2 credible-set SNVs (lifted to hg38) inside source disease-gained or
disease-lost H3K27ac regions. REF/ALT come from the 1kGP high-coverage panel when the
panel SNV alleles equal the credible-set pair, else from the hg38 FASTA base.
Writes targets.tsv and bait_windows.fa (+/- FLANK bp, REF and ALT versions).
Reads no allele counts.
"""
import argparse
import subprocess
from pathlib import Path

import pysam

FLANK = 75


def load_bed(path):
    out = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        out.append(f)
    return out


def in_any(chrom, pos1, intervals):
    # intervals: list of (chrom, start0, end, *rest); pos1 is 1-based
    hits = [iv for iv in intervals if iv[0] == chrom and int(iv[1]) < pos1 <= int(iv[2])]
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cs-bed", required=True, help="lifted hg38 credible-set SNVs, name = id|pip|tier")
    ap.add_argument("--diff-bed", required=True)
    ap.add_argument("--blacklist", required=True)
    ap.add_argument("--kgp-dir", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--bcftools", required=True)
    ap.add_argument("--max-tier", type=int, default=2)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=False)

    diff = load_bed(a.diff_bed)
    black = load_bed(a.blacklist)
    fa = pysam.FastaFile(a.fasta)
    rows = []
    for f in load_bed(a.cs_bed):
        chrom, s0, e = f[0], int(f[1]), int(f[2])
        cid, pip, tier = f[3].split("|")
        if int(tier) > a.max_tier:
            continue
        hits = in_any(chrom, e, diff)
        if not hits:
            continue
        classes = sorted({h[4] for h in hits})
        if len(classes) != 1:
            raise SystemExit(f"site in both gained and lost regions: {cid}")
        a1, a2 = cid.split(":")[2:4]
        pos = e
        fasta_base = fa.fetch(chrom, pos - 1, pos).upper()
        vcf = Path(a.kgp_dir) / f"1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        q = subprocess.run([a.bcftools, "query", "-r", f"{chrom}:{pos}", "-f", "%POS\t%REF\t%ALT\t%AF\n", str(vcf)],
                           capture_output=True, text=True, check=True).stdout
        ref = alt = None
        af = "NA"
        src = "none"
        for line in q.splitlines():
            p, r, al, afv = line.split("\t")
            if int(p) == pos and len(r) == 1 and len(al) == 1 and {r, al} == {a1, a2}:
                ref, alt, af, src = r, al, afv, "1kGP"
        if ref is None and fasta_base in (a1, a2):
            ref = fasta_base
            alt = a2 if a1 == fasta_base else a1
            src = "fasta"
        if ref is None:
            continue
        if ref != fasta_base:
            raise SystemExit(f"REF mismatch with FASTA at {chrom}:{pos}")
        blk = bool(in_any(chrom, pos, black))
        region = hits[0][3]
        rows.append(dict(chrom=chrom, pos=pos, ref=ref, alt=alt, cls="CS_" + classes[0], cs_id=cid, pip=pip,
                         tier=tier, af_1kgp=af, allele_source=src, region=region, blacklist=int(blk)))
    rows.sort(key=lambda r: (r["chrom"], r["pos"]))
    cols = ["chrom", "pos", "ref", "alt", "cls", "cs_id", "pip", "tier", "af_1kgp", "allele_source", "region", "blacklist"]
    with open(out / "targets.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(out / "bait_windows.fa", "w") as fh:
        for r in rows:
            seq = fa.fetch(r["chrom"], r["pos"] - 1 - FLANK, r["pos"] + FLANK).upper()
            assert seq[FLANK] == r["ref"]
            fh.write(f">{r['chrom']}_{r['pos']}_ref\n{seq}\n")
            fh.write(f">{r['chrom']}_{r['pos']}_alt\n{seq[:FLANK] + r['alt'] + seq[FLANK + 1:]}\n")
    n_by = {}
    for r in rows:
        n_by[(r["cls"], r["blacklist"])] = n_by.get((r["cls"], r["blacklist"]), 0) + 1
    print("targets", len(rows), n_by)


if __name__ == "__main__":
    main()
