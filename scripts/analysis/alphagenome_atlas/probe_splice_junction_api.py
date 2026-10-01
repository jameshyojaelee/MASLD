#!/usr/bin/env python3
"""Throwaway probe: what does the model API return for SPLICE_JUNCTIONS on one SNV?

Decides whether a signed, junction-specific splice effect can be built (step 36). Not a result.
usage: probe_splice_junction_api.py chr4 88231392
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import lib_atlas as la

FASTA = ("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/"
         "refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa")
LIVER = ["UBERON:0002107", "UBERON:0001114", "UBERON:0001115", "CL:0000182"]


def main() -> None:
    import pysam
    from alphagenome.data import genome
    from alphagenome.models import dna_client

    import atlas_query as aq

    chrom, pos1 = sys.argv[1], int(sys.argv[2])
    key = la.load_api_key()[0]
    print("api key length", len(key))
    model = dna_client.create(key)
    window = dna_client.SEQUENCE_LENGTH_1MB
    start0 = max(0, pos1 - 1 - window // 2)
    fasta = pysam.FastaFile(FASTA)
    ref_seq = fasta.fetch(chrom, start0, start0 + window).upper()
    base = ref_seq[pos1 - 1 - start0]
    alt_base = {"A": "G", "G": "A", "C": "T", "T": "C"}[base]
    alt_seq = ref_seq[:pos1 - 1 - start0] + alt_base + ref_seq[pos1 - start0:]
    print(f"{chrom}:{pos1} ref {base} -> probe alt {alt_base}; window {start0}-{start0+window}")
    interval = genome.Interval(chromosome=chrom, start=start0, end=start0 + window)
    outs = [dna_client.OutputType.SPLICE_JUNCTIONS, dna_client.OutputType.SPLICE_SITE_USAGE,
            dna_client.OutputType.SPLICE_SITES]
    got = {}
    for arm, seq in (("ref", ref_seq), ("alt", alt_seq)):
        got[arm] = aq.call_with_quota_retry(
            lambda s=seq: model.predict_sequence(sequence=s, requested_outputs=outs,
                                                 ontology_terms=LIVER, interval=interval),
            label=f"probe:{arm}")
    for arm in ("ref", "alt"):
        o = got[arm]
        print(f"\n--- {arm} ---")
        for attr in ("splice_junctions", "splice_site_usage", "splice_sites"):
            d = getattr(o, attr, None)
            if d is None:
                print(f"{attr}: None")
                continue
            v = np.asarray(d.values)
            print(f"{attr}: values {v.shape} dtype {v.dtype} min {np.nanmin(v):.4g} max {np.nanmax(v):.4g} "
                  f"share_negative {float(np.nanmean(v < 0)):.4f}")
            md = getattr(d, "metadata", None)
            if md is not None:
                print(f"  metadata cols {list(md.columns)[:10]} n {len(md)}")
            iv = getattr(d, "junctions", None)
            if iv is not None:
                print(f"  junctions type {type(iv)} n {len(iv)}")
                for j in list(iv)[:4]:
                    print("   ", j)
            if hasattr(d, "interval"):
                print("  interval", d.interval)
            if hasattr(d, "resolution"):
                print("  resolution", d.resolution)
    a, b = got["alt"].splice_junctions, got["ref"].splice_junctions
    if a is not None and b is not None:
        ja, jb = list(a.junctions), list(b.junctions)
        print(f"\njunction sets identical: {ja == jb}; n_alt {len(ja)} n_ref {len(jb)}")
        va, vb = np.asarray(a.values, float), np.asarray(b.values, float)
        if va.shape == vb.shape:
            d = np.nanmean(va - vb, axis=1)
            k = int(np.nanargmax(np.abs(d)))
            print(f"largest |mean delta| junction {ja[k]} delta {d[k]:.6g}; "
                  f"share nonzero delta {float(np.mean(np.abs(d) > 0)):.4f}")


if __name__ == "__main__":
    main()
