#!/usr/bin/env python3
"""Step 31 (P3a): build and retrieve the external same-variant benchmark sets.

usage: 31_query_benchmark_sets.py caqtl|sqtl|mpra

caqtl : every peak-lead variant of the Currin 2025 bulk-liver caQTL benchmark (`direction_features_caqtl.tsv`,
        label-frame ref/alt already FASTA-resolved by 79_score_panels.py; re-validated here; SNVs only, indels are not served).
sqtl  : GTEx v8 liver sGene lead variants (`Liver.v8.sgenes.txt.gz`) plus the significant pairs' variants for the
        same clusters capped per cluster (see SQTL_PER_CLUSTER); splice/RNA scorers are all retrieved (every
        scorer is archived, the selection happens downstream, fixed before any result is read).
mpra  : the Hu/Zhu 2026 MPRA library variants (built by 32_mpra_library.py into tables/mpra_library_variants.tsv).

Every set is written once as raw/query_<set>_variants.txt and archived under raw/atlas_<set>/ (chunked, resumable).
"""

from __future__ import annotations

import gzip
import sys
from collections import defaultdict

import pysam

import atlas_query as aq
import lib_atlas as la

ROOT = la.out_root()
RAW = ROOT / "raw"
TABLES = ROOT / "tables"
FASTA = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
CAQTL = la.PROJECT / "GWAS/finemapping/results/seqfunc/direction_features_caqtl.tsv"
SQTL_DIR = la.PROJECT / "data/external/gtex_v8_liver_sqtl/GTEx_Analysis_v8_sQTL"
SQTL_PER_CLUSTER = 10


def validated_uid(fasta, chrom: str, pos: int, ref: str, alt: str, stats: dict) -> str | None:
    """hg38 uid after checking the reference allele against the FASTA (first base for indels). Never repairs."""
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    if chrom not in fasta.references:
        stats["nonstandard_chrom"] += 1
        return None
    base = fasta.fetch(chrom, pos - 1, pos - 1 + len(ref)).upper()
    if base != ref.upper():
        stats["fasta_mismatch"] += 1
        return None
    if len(ref) != 1 or len(alt) != 1:
        # probed 2026-09-09: the Atlas point-query endpoint returns UNIMPLEMENTED for insertions and deletions (SDK 0.9.0)
        stats["indel_not_served"] += 1
        return None
    stats["snv"] += 1
    return la.variant_uid(chrom, pos, ref.upper(), alt.upper())


def caqtl_uids(fasta) -> list[str]:
    stats = defaultdict(int)
    uids = []
    rows = la.read_tsv(CAQTL)
    for r in rows:
        u = validated_uid(fasta, str(r["chr"]), int(float(r["pos_hg38"])), r["ref"], r["alt"], stats)
        if u:
            uids.append(u)
    if stats["fasta_mismatch"] > 0:
        raise la.ContractError(f"caQTL label-frame alleles disagree with the FASTA at {stats['fasta_mismatch']} rows")
    la.log(f"caqtl: {len(rows)} rows -> {len(set(uids))} uids ({dict(stats)})")
    return sorted(set(uids))


def sqtl_uids(fasta) -> list[str]:
    stats = defaultdict(int)
    uids, leads = [], set()
    with gzip.open(SQTL_DIR / "Liver.v8.sgenes.txt.gz", "rt") as h:
        header = h.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(header)}
        for line in h:
            p = line.rstrip("\n").split("\t")
            v = p[ix["variant_id"]]                       # chr1_902173_GGTGT_G_b38
            chrom, pos, ref, alt, _ = v.split("_")
            u = validated_uid(fasta, chrom, int(pos), ref, alt, stats)
            if u:
                leads.add(u)
    per_cluster = defaultdict(list)
    with gzip.open(SQTL_DIR / "Liver.v8.sqtl_signifpairs.txt.gz", "rt") as h:
        header = h.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(header)}
        for line in h:
            p = line.rstrip("\n").split("\t")
            cluster = p[ix["phenotype_id"]].split(":")[3]
            per_cluster[cluster].append((float(p[ix["pval_nominal"]]), p[ix["variant_id"]]))
    for cluster, items in per_cluster.items():
        for _, v in sorted(items)[:SQTL_PER_CLUSTER]:
            chrom, pos, ref, alt, _ = v.split("_")
            u = validated_uid(fasta, chrom, int(pos), ref, alt, stats)
            if u:
                uids.append(u)
    allu = sorted(leads | set(uids))
    if stats["fasta_mismatch"] > 0:
        raise la.ContractError(f"GTEx b38 alleles disagree with the FASTA at {stats['fasta_mismatch']} rows")
    la.log(f"sqtl: {len(leads)} sGene leads + top-{SQTL_PER_CLUSTER} per cluster over {len(per_cluster)} clusters -> {len(allu)} uids ({dict(stats)})")
    return allu


def mpra_uids(fasta) -> list[str]:
    stats = defaultdict(int)
    uids = []
    for r in la.read_tsv(TABLES / "mpra_library_variants.tsv"):
        u = validated_uid(fasta, r["chrom"], int(r["pos_hg38"]), r["ref"], r["alt"], stats)
        if u:
            uids.append(u)
    la.log(f"mpra: {len(set(uids))} uids ({dict(stats)})")
    return sorted(set(uids))


def main(which: str) -> None:
    fasta = pysam.FastaFile(FASTA)
    build = {"caqtl": caqtl_uids, "sqtl": sqtl_uids, "mpra": mpra_uids}[which]
    uids = build(fasta)
    listing = RAW / f"query_{which}_variants.txt"
    if listing.exists():
        prior = listing.read_text().split()
        if prior != uids:
            raise la.ContractError(f"{listing} exists with a different variant list; refusing to continue")
    else:
        listing.parent.mkdir(parents=True, exist_ok=True)
        listing.write_text("\n".join(uids) + "\n")
    client, sdk = aq.create_client()
    scorers = sorted(client.scorer_metadata())
    chunks = aq.query_archived(client, sdk, uids, scorers, RAW / f"atlas_{which}", extra={"stage": which})
    la.log(f"{which}: {len(chunks)} chunks archived")


if __name__ == "__main__":
    main(sys.argv[1])
