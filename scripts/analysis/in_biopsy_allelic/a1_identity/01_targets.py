#!/usr/bin/env python3
"""A1 step 1: common exonic SNV targets for the RNA identity panel.

Selects biallelic dbSNP b157 SNVs with 1000 Genomes 30x global minor allele
frequency >= --min-maf inside exons of the --n-genes most expressed genes
(median CPM 20-500 in one reference cohort), drops the MHC, and thins to one site per
--thin-bp. Writes bcftools targets (chrom, pos, REF,ALT) and a regions file.
Reads no phenotype or outcome.
"""
import argparse
import gzip
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

MHC = ("chr6", 28_500_000, 33_500_000)


def top_genes(counts_path, n_genes, cpm_lo=20.0, cpm_hi=500.0):
    """Genes with median CPM in [cpm_lo, cpm_hi]: enough depth to genotype (about 30-750x
    at 30M reads) without the enormous read stacks of the most abundant liver genes,
    which made pileups I/O-bound (job 0931). Highest-CPM genes in the band first."""
    counts = pd.read_csv(counts_path, sep="\t", comment="#", index_col=0)
    mat = counts.iloc[:, 5:].to_numpy(dtype=float)  # featureCounts: 5 annotation columns after Geneid
    cpm = mat / mat.sum(axis=0, keepdims=True) * 1e6
    median_cpm = pd.Series(np.median(cpm, axis=1), index=counts.index)
    chrom = counts["Chr"].astype(str).str.split(";").str[0]
    keep = chrom.str.fullmatch(r"chr([0-9]+|X)") & median_cpm.between(cpm_lo, cpm_hi)
    return set(median_cpm[keep].sort_values(ascending=False).index[:n_genes].str.split(".").str[0])


def exon_intervals(gtf_path, genes):
    rows = []
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "exon":
                continue
            gid = f[8].split('gene_id "', 1)[1].split('"', 1)[0].split(".")[0]
            if gid in genes:
                rows.append((f[0], int(f[3]) - 1, int(f[4])))
    iv = pd.DataFrame(rows, columns=["chrom", "start", "end"]).sort_values(["chrom", "start"])
    merged = []
    for chrom, grp in iv.groupby("chrom", sort=False):
        cur_s, cur_e = None, None
        for s, e in zip(grp["start"], grp["end"]):
            if cur_s is None or s > cur_e:
                if cur_s is not None:
                    merged.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s, e
            else:
                cur_e = max(cur_e, e)
        merged.append((chrom, cur_s, cur_e))
    return pd.DataFrame(merged, columns=["chrom", "start", "end"])


def refseq_map(dbsnp):
    header = subprocess.run(["bcftools", "view", "-h", dbsnp], capture_output=True, text=True, check=True).stdout
    out = {}
    for line in header.splitlines():
        if line.startswith("##contig=<ID=NC_0000"):
            acc = line.split("ID=", 1)[1].split(">", 1)[0].split(",", 1)[0]
            num = int(acc.split(".")[0].replace("NC_", ""))
            name = {23: "chrX", 24: "chrY"}.get(num, f"chr{num}")
            out[name] = acc
    return out


def parse_freq(info_freq):
    for block in info_freq.split("|"):
        if block.startswith("1000Genomes_30X:"):
            return [0.0 if v == "." else float(v) for v in block.split(":", 1)[1].split(",")]
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--dbsnp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-genes", type=int, default=3000)
    ap.add_argument("--min-maf", type=float, default=0.2)
    ap.add_argument("--thin-bp", type=int, default=10_000)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    genes = top_genes(a.counts, a.n_genes)
    exons = exon_intervals(a.gtf, genes)
    acc = refseq_map(a.dbsnp)
    exons = exons[exons["chrom"].isin(acc)]
    reg = exons.assign(chrom=exons["chrom"].map(acc), start=exons["start"] + 1)
    reg_path = out / "dbsnp_query_regions.tsv"
    reg.to_csv(reg_path, sep="\t", header=False, index=False)

    q = subprocess.run(
        ["bcftools", "query", "-R", str(reg_path), "-i", 'INFO/VC="SNV"',
         "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/FREQ\n", a.dbsnp],
        capture_output=True, text=True, check=True).stdout
    back = {v: k for k, v in acc.items()}
    rows = []
    for line in q.splitlines():
        chrom, pos, ref, alts, freq = line.split("\t")
        if len(ref) != 1 or freq == ".":
            continue
        f = parse_freq(freq)
        alts = alts.split(",")
        if f is None or len(f) != len(alts) + 1:
            continue
        nonzero = [i for i, v in enumerate(f) if v > 0]
        if len(nonzero) != 2 or nonzero[0] != 0:
            continue
        alt = alts[nonzero[1] - 1]
        if len(alt) != 1:
            continue
        maf = min(f[0], f[nonzero[1]])
        if maf < a.min_maf:
            continue
        rows.append((back[chrom], int(pos), ref, alt, maf))
    t = pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "maf"]).drop_duplicates(["chrom", "pos"])
    in_mhc = (t["chrom"] == MHC[0]) & t["pos"].between(MHC[1], MHC[2])
    t = t[~in_mhc]
    t["bin"] = t["pos"] // a.thin_bp
    t = t.sort_values("maf", ascending=False).drop_duplicates(["chrom", "bin"])
    order = {c: i for i, c in enumerate([f"chr{i}" for i in range(1, 23)] + ["chrX"])}
    t = t.assign(o=t["chrom"].map(order)).sort_values(["o", "pos"])

    tgt = out / "targets.tsv"
    t.assign(alleles=t["ref"] + "," + t["alt"])[["chrom", "pos", "alleles"]].to_csv(tgt, sep="\t", header=False, index=False)
    subprocess.run(["bgzip", "-f", str(tgt)], check=True)
    subprocess.run(["tabix", "-f", "-s1", "-b2", "-e2", str(tgt) + ".gz"], check=True)
    t[["chrom", "pos"]].assign(end=t["pos"])[["chrom", "pos", "end"]].to_csv(out / "regions.tsv", sep="\t", header=False, index=False)
    t.drop(columns=["bin", "o"]).to_csv(out / "targets_annotated.tsv", sep="\t", index=False)
    print(f"genes={len(genes)} exon_intervals={len(exons)} targets={len(t)} chrX={int((t['chrom'] == 'chrX').sum())}")


if __name__ == "__main__":
    main()
