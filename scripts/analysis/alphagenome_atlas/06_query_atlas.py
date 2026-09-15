#!/usr/bin/env python3
"""Step 06: full Atlas retrieval for the query set, direct universes first, resumable.

usage: 06_query_atlas.py direct|enzyme
"""

from __future__ import annotations

import csv
import sys

import atlas_query as aq
import lib_atlas as la

ROOT = la.out_root()
TABLES = ROOT / "tables"
RAW = ROOT / "raw"
FASTA = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"


def main(stage: str) -> None:
    if stage not in ("direct", "enzyme", "gate"):
        raise SystemExit("stage must be direct, enzyme or gate")
    if stage == "gate":
        # the archived 211-variant AlphaGenome gate set (65b): hg38 alleles from the 65-series score file, FASTA-validated
        import pysam
        truth = la.read_tsv(la.PROJECT / "GWAS/finemapping/results/seqfunc/broadaway_benchmark_truth.tsv")
        ag = la.read_tsv(la.PROJECT / "GWAS/finemapping/results/seqfunc/alphagenome_eqtl_scores.tsv")
        fasta = pysam.FastaFile(FASTA)
        rows = la.gate_set(truth, ag, fasta.fetch)
        uids = sorted({r["variant_uid"] for r in rows})
        la.log(f"gate: {len(rows)} lead-gene rows, {len(uids)} unique SNV uids (FASTA-validated)")
        client, sdk = aq.create_client()
        scorers = sorted(client.scorer_metadata())
        chunks = aq.query_archived(client, sdk, uids, scorers, RAW / "atlas_gate", extra={"stage": stage, "n_lead_gene_rows": len(rows)})
        la.log(f"gate: {len(chunks)} chunks archived")
        return
    direct, enzyme = set(), set()
    with la.open_text(TABLES / "signal_variant_weights.tsv.gz") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            if r["in_query_set"] != "True":
                continue
            (enzyme if r["universe"] == "C_enzyme" else direct).add(r["variant_uid"])
    uids = sorted(direct) if stage == "direct" else sorted(enzyme - direct)
    la.log(f"{stage}: {len(uids)} unique variants to query (direct {len(direct)}, enzyme-only {len(enzyme - direct)})")
    client, sdk = aq.create_client()
    meta = client.scorer_metadata()
    known = {r["scorer"] for r in la.read_tsv(TABLES / "scorer_metadata.tsv")}
    if set(meta) != known:
        raise la.ContractError(f"scorer list changed since the pilot: {sorted(set(meta) ^ known)}")
    scorers = sorted(meta)
    with (RAW / f"query_{stage}_variants.txt").open("w") as handle:
        handle.write("\n".join(uids) + "\n")
    chunks = aq.query_archived(client, sdk, uids, scorers, RAW / f"atlas_{stage}", extra={"stage": stage})
    la.log(f"{stage}: {len(chunks)} chunks archived")


if __name__ == "__main__":
    main(sys.argv[1])
