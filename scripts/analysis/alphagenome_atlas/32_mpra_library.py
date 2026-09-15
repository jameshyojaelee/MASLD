#!/usr/bin/env python3
"""Step 32 (P3a, MPRA): build the Hu/Zhu 2026 MPRA library table and retrieve interval saturation for every oligo.

The barcode map names oligos by their 126-bp hg38 interval only; the tested variant sits at offset 63 from the
interval start (verified on all 1,128 rows of mpra_substrate_truth.tsv). Differential-activity-variant (DAV)
calls per context come from Tables S2-S5 (converted to tables/mpra_dav_calls_S2_S5.tsv by readxl). Alleles are
known only where the interval overlaps the seqfunc substrate (ref_hg38/alt_hg38); for the rest the Atlas
returns every substitution at the centre position, and the analysis step uses the allele-agnostic maximum
(fixed rule, chosen before any score was read; MPRA outcomes never select tracks, windows or aggregation).

Outputs: tables/mpra_library_intervals.tsv, raw/atlas_mpra_intervals/<interval>/ (query_interval archives)
"""

from __future__ import annotations

import csv
import re
import time
from collections import defaultdict

import pysam
from alphagenome.data import genome

import atlas_archive as aa
import atlas_query as aq
import lib_atlas as la

ROOT = la.out_root()
RAW = ROOT / "raw"
TABLES = ROOT / "tables"
FASTA = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
MPRA = la.PROJECT / "GWAS/finemapping/data/seqfunc_external/hu2025_mpra"
TRUTH = la.PROJECT / "GWAS/finemapping/results/seqfunc/mpra_benchmark/mpra_substrate_truth.tsv"
CENTRE_OFFSET = 63
SCORERS = ["ATAC", "DNASE", "CHIP_HISTONE", "CHIP_TF", "CAGE", "RNA_SEQ", "AVI_SCORE"]
CONTEXTS = ["HepG2_ctrl", "HepG2_PAOA", "LX2_ctrl", "LX2_TGFB"]


def parse_interval(s: str) -> tuple[str, int, int]:
    chrom, rng_ = s.split(":")
    a, b = rng_.split("-")
    return chrom, int(a), int(b)


def build_table(fasta) -> list[dict]:
    intervals = [l.strip() for l in (MPRA / "tested_oligo_intervals.tsv").open() if re.match(r"^chr[0-9XY]+:\d+-\d+$", l.strip())]
    if len(intervals) != 3498:
        raise la.ContractError(f"expected 3,498 tested oligo intervals, parsed {len(intervals)}")
    dav = defaultdict(dict)
    for r in la.read_tsv(TABLES / "mpra_dav_calls_S2_S5.tsv"):
        dav[r["coordination"]][r["context"]] = (r["log2(Fold change)"], r["FDR"], r["rsid"])
    alleles = {}
    for r in la.read_tsv(TRUTH):
        if r.get("coord") and r.get("ref_hg38") and r.get("alt_hg38"):
            alleles[r["coord"]] = (r["ref_hg38"], r["alt_hg38"], r.get("rsid", ""))
    rows = []
    n_off = 0
    for iv in intervals:
        chrom, start, end = parse_interval(iv)
        centre = start + CENTRE_OFFSET
        base = fasta.fetch(chrom, centre - 1, centre).upper()
        ref, alt, rsid = alleles.get(iv, ("", "", ""))
        if ref and ref != base:
            n_off += 1
        row = {"interval": iv, "chrom": chrom, "start": start, "end": end, "centre_pos_hg38": centre, "fasta_base": base,
               "ref_hg38": ref, "alt_hg38": alt, "rsid": rsid or next((v[2] for v in dav.get(iv, {}).values()), ""),
               "allele_known": bool(ref and alt)}
        for ctx in CONTEXTS:
            l2, fdr, _ = dav.get(iv, {}).get(ctx, ("", "", ""))
            row[f"dav_{ctx}"] = bool(l2)
            row[f"log2fc_{ctx}"] = l2
            row[f"fdr_{ctx}"] = fdr
        rows.append(row)
    if n_off:
        raise la.ContractError(f"{n_off} substrate ref alleles disagree with the FASTA at the oligo centre")
    la.log(f"mpra: {len(rows)} intervals; alleles known {sum(r['allele_known'] for r in rows)}; DAV per context "
           + ", ".join(f"{c}={sum(r[f'dav_{c}'] for r in rows)}" for c in CONTEXTS))
    return rows


def main() -> None:
    fasta = pysam.FastaFile(FASTA)
    rows = build_table(fasta)
    out = TABLES / "mpra_library_intervals.tsv"
    if out.exists():
        la.log(f"{out.name} already written by an earlier run; continuing to retrieval")
    else:
        la.write_tsv_once(out, rows, list(rows[0].keys()))
    client, sdk = aq.create_client()
    out_root = RAW / "atlas_mpra_intervals"
    done = 0
    for r in rows:
        d = out_root / r["interval"].replace(":", "_")
        if (d / "request.json").exists():
            done += 1
            continue
        interval = genome.Interval(chromosome=r["chrom"], start=r["start"] - 1, end=r["end"])
        t0 = time.time()
        res = aq.call_with_quota_retry(lambda: client.query_interval(interval, requested_scorers=SCORERS, progress_bar=False, max_workers=4))
        aa.archive_scores(res, d, {"interval": r["interval"], "requested_scorers": SCORERS, "sdk_version": sdk, "elapsed_seconds": time.time() - t0,
                                   "rule": "126-bp MPRA oligo; variant at offset 63"})
        done += 1
        if done % 100 == 0:
            la.log(f"mpra intervals archived: {done}/{len(rows)}")
    la.log(f"mpra: {done} intervals archived")


if __name__ == "__main__":
    main()
