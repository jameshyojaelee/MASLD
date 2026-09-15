#!/usr/bin/env python3
"""F2/F3 step 1: the 1000 Genomes EUR common-SNV panel the corrected null draws from.

The P5 deposit's stamped prespecification says the matched null draws "both variants common in 1000G
EUR"; its code drew uniform reference positions with a fixed complementary alternate allele. This step
builds what the prespecification describes.

Per chromosome: PLINK 1.9 --freq on data/1kg_eur/chr<N>_eur.{bed,bim,fam} (379 EUR samples, hg19),
restricted to biallelic SNVs (both alleles one base of ACGT) with EUR MAF >= MIN_MAF; the hg19 positions
are lifted to hg38 with data/broadaway_eqtl/hg19ToHg38.over.chain (rtracklayer, first mapping, same rule
as 03_liftover_universeB.R); only single mappings that stay on the same chromosome are kept.

The FASTA check (the reference base at the lifted position must be one of the two 1000G alleles, and it
becomes the placed REF) happens at draw time in f2_02, not here, so the panel stays a pure population
object and the placement rule lives with the scorer that uses it.

Output: tables/null_panel/chr<N>.tsv.gz with pos_hg19, pos_hg38, a1, a2, maf.
"""

from __future__ import annotations

import gzip
import os
import pathlib
import subprocess
import sys

OUT = pathlib.Path(os.environ["F2_OUT_ROOT"])
PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
KG = PROJECT / "data/1kg_eur"
CHAIN = PROJECT / "data/broadaway_eqtl/hg19ToHg38.over.chain"
RSCRIPT = "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
LIFT_R = PROJECT / "scripts/analysis/alphagenome_program/f2_01_liftover.R"
SCRATCH = pathlib.Path(os.environ.get("F2_SCRATCH", "/tmp/f2_null_panel"))
MIN_MAF = 0.05
BASES = set("ACGT")


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def plink_freq(chrom: int) -> pathlib.Path:
    stem = SCRATCH / f"chr{chrom}_eur"
    frq = stem.with_suffix(".frq")
    if frq.exists():
        return frq
    cmd = ["plink", "--bfile", str(KG / f"chr{chrom}_eur"), "--freq", "--out", str(stem)]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    return frq


def build_chrom(chrom: int) -> int:
    out_path = OUT / "tables" / "null_panel" / f"chr{chrom}.tsv.gz"
    if out_path.exists():
        log(f"chr{chrom}: already built")
        return -1
    frq = plink_freq(chrom)
    # bim: chrom, id, cm, pos, a1, a2  (a1 = the allele PLINK counts, a2 the other)
    pos_by_id: dict[str, int] = {}
    with open(KG / f"chr{chrom}_eur.bim") as handle:
        for line in handle:
            p = line.split()
            pos_by_id[p[1]] = int(p[3])
    keep = []
    with open(frq) as handle:
        header = handle.readline().split()
        i_snp, i_a1, i_a2, i_maf = header.index("SNP"), header.index("A1"), header.index("A2"), header.index("MAF")
        for line in handle:
            p = line.split()
            a1, a2 = p[i_a1].upper(), p[i_a2].upper()
            if len(a1) != 1 or len(a2) != 1 or a1 not in BASES or a2 not in BASES or a1 == a2:
                continue
            try:
                maf = float(p[i_maf])
            except ValueError:
                continue
            if maf < MIN_MAF:
                continue
            pos = pos_by_id.get(p[i_snp])
            if pos is None:
                continue
            keep.append((pos, a1, a2, maf))
    keep.sort()
    log(f"chr{chrom}: {len(keep)} biallelic SNVs with EUR MAF >= {MIN_MAF}")
    lift_in = SCRATCH / f"chr{chrom}_lift_in.tsv"
    lift_out = SCRATCH / f"chr{chrom}_lift_out.tsv"
    with open(lift_in, "w") as handle:
        handle.write("chrom_hg19\tpos_hg19\n")
        for pos, _, _, _ in keep:
            handle.write(f"chr{chrom}\t{pos}\n")
    if lift_out.exists():
        lift_out.unlink()
    subprocess.run([RSCRIPT, str(LIFT_R), str(lift_in), str(lift_out), str(CHAIN)], check=True)
    lifted: list[tuple[int, str, int]] = []
    with open(lift_out) as handle:
        head = handle.readline().rstrip("\n").split("\t")
        j = {c: k for k, c in enumerate(head)}
        for line in handle:
            p = line.rstrip("\n").split("\t")
            lifted.append((int(p[j["n_liftover_mappings"]]), p[j["hg38_chrom"]], int(p[j["hg38_pos"]]) if p[j["hg38_pos"]] not in ("", "NA") else -1))
    if len(lifted) != len(keep):
        raise RuntimeError(f"chr{chrom}: liftover returned {len(lifted)} rows for {len(keep)} positions")
    n = 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(".tmp.gz")
    with gzip.open(tmp, "wt") as handle:
        handle.write("pos_hg19\tpos_hg38\ta1\ta2\tmaf\n")
        for (pos, a1, a2, maf), (nmap, hchrom, hpos) in zip(keep, lifted):
            if nmap != 1 or hchrom != f"chr{chrom}" or hpos < 1:
                continue
            handle.write(f"{pos}\t{hpos}\t{a1}\t{a2}\t{maf:.6f}\n")
            n += 1
    tmp.rename(out_path)
    log(f"chr{chrom}: {n} panel variants written (single-mapping, same chromosome)")
    return n


def main() -> None:
    SCRATCH.mkdir(parents=True, exist_ok=True)
    chroms = [int(c) for c in sys.argv[1:]] or list(range(1, 23))
    total = 0
    for c in chroms:
        n = build_chrom(c)
        if n > 0:
            total += n
    log(f"null panel: {total} variants over {len(chroms)} chromosomes")


if __name__ == "__main__":
    main()
