#!/usr/bin/env python3
"""Model A pivot, step T3: allele-specific mappability filter for tag SNVs.

For each unique tag SNV: L-bp windows (L = --read-length, default 100) lying inside one
exon of its gene and containing the SNV are written with the REF and with the ALT base;
all reads are aligned with STAR to --index (unspliced reads expected). Window starts are
taken per exon on a grid of step max(1, (hi - lo) // 20) from the lowest start, and only
the 20 lowest starts are kept, so windows with the SNV near the read start are rarely
simulated (they are not evenly spaced over read offsets). T3' round 2 uses
t3prime/scripts/t3prime_sim.py in the obs-20260926T155651Z execution instead.
A read "returns" if its primary record is unique (NH:i:1) and lies on its own chromosome
at its own start (SAM POS = window start, so a read whose first base is soft-clipped does
not return; a read whose last bases are soft-clipped does return, although bcftools
mpileup does not count a clipped SNV base). On an index with ALT/patch/scaffold contigs, a
read that also aligns equally well to such a contig has NH >= 2, and one that aligns
better there is placed there; neither returns. A tag is kept if both alleles have >= 2
windows, return rates differ by <= 0.01, and both are >= 0.5. Writes
tags_mappability.tsv and tags_final.tsv.
"""
import argparse
import subprocess
from pathlib import Path

import pandas as pd

MAXW = 20
STAR = "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/STAR"
SAMTOOLS = "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/samtools"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tags-dir", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--read-length", type=int, default=100, help="simulated read (mate) length")
    a = ap.parse_args()
    L = a.read_length
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    tags = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(Path(a.tags_dir).glob("chr*.tags.tsv"))], ignore_index=True)
    tags["gene_key"] = tags["gene_id"].str.split(".").str[0]
    snv = tags.drop_duplicates(["chrom", "pos", "ref", "alt", "gene_key"])
    ex = {}
    for line in open(a.gtf):
        if line.startswith("#"):
            continue
        f = line.split("\t", 9)
        if f[2] != "exon":
            continue
        g = f[8].split('gene_id "', 1)[1].split('"', 1)[0].split(".")[0]
        ex.setdefault((f[0], g), set()).add((int(f[3]), int(f[4])))
    windows = []
    for _, t in snv.iterrows():
        starts = set()
        for s, e in ex.get((t["chrom"], t["gene_key"]), ()):
            if s <= t["pos"] <= e and e - s + 1 >= L:
                lo, hi = max(s, t["pos"] - L + 1), min(t["pos"], e - L + 1)
                if lo <= hi:
                    step = max(1, (hi - lo) // MAXW)
                    starts.update(range(lo, hi + 1, step))
        for w in sorted(starts)[:MAXW]:
            windows.append((t["chrom"], t["pos"], t["ref"], t["alt"], w))
    reg = out / "windows.txt"
    reg.write_text("".join(f"{c}:{w}-{w + L - 1}\n" for c, p, r, x, w in windows))
    seqs = subprocess.run([SAMTOOLS, "faidx", a.fasta, "-r", str(reg)], capture_output=True, text=True, check=True).stdout
    recs, name, buf = {}, None, []
    for line in seqs.splitlines():
        if line.startswith(">"):
            if name:
                recs[name] = "".join(buf)
            name, buf = line[1:], []
        else:
            buf.append(line.strip())
    if name:
        recs[name] = "".join(buf)
    fq = out / "sim_reads.fq"
    with open(fq, "w") as fh:
        for c, p, r, x, w in windows:
            seq = recs[f"{c}:{w}-{w + L - 1}"].upper()
            off = p - w
            if seq[off] != r:
                continue
            for allele, base in (("R", r), ("A", x)):
                s2 = seq[:off] + base + seq[off + 1:]
                fh.write(f"@{c}|{p}|{r}|{x}|{w}|{allele}\n{s2}\n+\n{'I' * L}\n")
    subprocess.run([STAR, "--runThreadN", str(a.threads), "--genomeDir", a.index, "--readFilesIn", str(fq),
                    "--outSAMtype", "SAM", "--outFileNamePrefix", str(out / "sim."), "--outSAMattributes", "NH",
                    "--outSAMunmapped", "Within", "--outFilterMultimapNmax", "20"], check=True)
    hits = {}
    for line in open(out / "sim.Aligned.out.sam"):
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if int(f[1]) & 256:
            continue
        c, p, r, x, w, allele = f[0].split("|")
        nh = next((int(t[5:]) for t in f[11:] if t.startswith("NH:i:")), 0)
        ok = nh == 1 and f[2] == c and int(f[3]) == int(w)
        hits.setdefault((c, int(p), r, x, allele), []).append(ok)
    rows = []
    for (c, p, r, x) in {(k[0], k[1], k[2], k[3]) for k in hits}:
        ra, aa = hits.get((c, p, r, x, "R"), []), hits.get((c, p, r, x, "A"), [])
        rr, ar = (sum(ra) / len(ra) if ra else 0.0), (sum(aa) / len(aa) if aa else 0.0)
        keep = len(ra) >= 2 and len(aa) >= 2 and abs(rr - ar) <= 0.01 and min(rr, ar) >= 0.5
        rows.append((c, p, r, x, len(ra), rr, ar, keep))
    m = pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "n_windows", "return_ref", "return_alt", "keep"])
    m.to_csv(out / "tags_mappability.tsv", sep="\t", index=False)
    final = tags.merge(m[m["keep"]][["chrom", "pos", "ref", "alt"]], on=["chrom", "pos", "ref", "alt"])
    final.to_csv(out / "tags_final.tsv", sep="\t", index=False)
    print(f"tag SNVs {len(snv)}; with windows {len(m)}; kept {int(m['keep'].sum())}; "
          f"final tag rows {len(final)}; eGenes {final['gene_id'].nunique()}")


if __name__ == "__main__":
    main()
