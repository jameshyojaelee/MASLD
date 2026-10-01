#!/usr/bin/env python3
"""Stage 3: reproduce the hosted-API HSD17B13 anchor with the LOCAL weights.

The hosted model API scored HSD17B13 rs72613567 (hg38 chr4:87310240 T>TA) at splice-site-usage
log2 -1.090159766332727, from p6f-indel-rescue-20260914T135052Z (and again, unchanged, in the
superseding 20260914T231617Z run).

Recipe, taken from scripts/analysis/alphagenome_atlas/63_indel_rescue.py and reproduced here without
importing or modifying anything in that directory:
  window   1,048,576 bp, start0 = pos1 - 1 - WINDOW//2
  alt      length-compensated at the far end (~500 kb from the readout)
  outputs  RNA_SEQ, ATAC, DNASE, CHIP_HISTONE, SPLICE_SITE_USAGE
  ontology UBERON:0002107, UBERON:0001114, UBERON:0001115, CL:0000182
  readout  mean |track| over +/-2,000 bp of the variant, H3K27ac columns only for CHIP_HISTONE
  effect   log2((alt + 1e-9) / (ref + 1e-9))

The local and hosted services are different implementations. The two values are reported side by side
and are never averaged, pooled, or substituted for one another.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import h1_runtime_common as C  # noqa: E402

CHROM = "chr4"
POS1 = 87_310_240
REF = "T"
ALT = "TA"
HOSTED = {
    "source": (
        "GWAS/finemapping/results/alphagenome_atlas/p6f-indel-rescue-20260914T135052Z/"
        "tables/indel_rescue_effects.tsv"
    ),
    "source_variant_id": "hg19:chr4:88231392:T:TA",
    "splice_log2": -1.090159766332727,
    "splice_ref": 1.774854851121615e-04,
    "splice_alt": 8.336606250876295e-05,
    "rna_log2": 0.35681125162957333,
    "atac_log2": 0.0807375345443907,
    "dnase_log2": 0.023492713888251197,
    "h3k27ac_log2": 0.040576708104114595,
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--length", type=int, default=1_048_576)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    import pysam
    from alphagenome.data import genome

    model, device, load_s = C.load_model()
    outs = C.requested_outputs()
    fasta = pysam.FastaFile(C.FASTA_PATH)

    window = args.length
    start0 = max(0, POS1 - 1 - window // 2)
    end = start0 + window
    ref_seq = fasta.fetch(CHROM, start0, end).upper()
    if len(ref_seq) != window:
        raise C.ContractError(f"fasta returned {len(ref_seq)} bp, wanted {window}")
    alt_seq = C.indel_alt_sequence(
        lambda c, a, b: fasta.fetch(c, a, b), CHROM, start0, end, POS1, REF, ALT
    )
    interval = genome.Interval(chromosome=CHROM, start=start0, end=end)
    mask = C.center_mask(window, POS1, start0, C.FLANK)

    rec = {
        "variant": f"{CHROM}:{POS1} {REF}>{ALT} (hg38, rs72613567)",
        "checkpoint": str(C.CHECKPOINT),
        "fasta": C.FASTA_PATH,
        "checkpoint_load_seconds": round(load_s, 2),
        "window_bp": window,
        "window_start0": start0,
        "window_end": end,
        "flank_bp": C.FLANK,
        "n_masked_positions": int(mask.sum()),
        "ontology_terms": C.LIVER_TERMS,
        "requested_outputs": C.OUTPUT_NAMES,
        "reference_base_check": ref_seq[POS1 - 1 - start0],
        "ref_alt_sequences_differ": bool(ref_seq != alt_seq),
        "hosted": HOSTED,
    }

    got = {}
    for arm, seq in (("ref", ref_seq), ("alt", alt_seq)):
        o = model.predict_sequence(
            sequence=seq, requested_outputs=outs, ontology_terms=C.LIVER_TERMS, interval=interval
        )
        got[arm] = C.summarise(o, mask)
        if arm == "ref":
            # Which tracks the pinned liver ontology actually selects, recorded once.
            sel = {}
            for name, attr in C.CHANNELS:
                td = getattr(o, attr, None)
                if td is None or td.values is None:
                    sel[name] = {"n_tracks": 0}
                    continue
                md = getattr(td, "metadata", None)
                cols = {}
                if md is not None:
                    for col in ("name", "ontology_curie", "biosample_name", "assay", "strand"):
                        if col in md:
                            cols[col] = [str(x) for x in md[col].astype(str).tolist()]
                sel[name] = {
                    "n_tracks": int(td.values.shape[1]),
                    "resolution_bp": int(td.resolution),
                    "metadata_columns": list(md.columns) if md is not None else None,
                    "tracks": cols,
                }
            rec["liver_track_selection"] = sel
        del o

    local = {}
    for ch in ("rna", "atac", "dnase", "h3k27ac", "splice"):
        a, b = got["alt"][ch], got["ref"][ch]
        local[f"{ch}_ref"] = b
        local[f"{ch}_alt"] = a
        local[f"{ch}_log2"] = C.log2_ratio(a, b)
    rec["local"] = local

    comp = {}
    for ch in ("rna", "atac", "dnase", "h3k27ac", "splice"):
        h = HOSTED.get(f"{ch}_log2")
        l = local[f"{ch}_log2"]
        comp[ch] = {
            "local_log2": l,
            "hosted_log2": h,
            "difference_local_minus_hosted": (l - h) if (h is not None and l == l) else None,
            "same_sign": (
                bool(np.sign(l) == np.sign(h)) if (h is not None and l == l and h != 0) else None
            ),
            "abs_ratio_local_over_hosted": (
                float(abs(l) / abs(h)) if (h not in (None, 0) and l == l) else None
            ),
        }
    rec["local_vs_hosted"] = comp
    rec["not_pooled"] = (
        "local and hosted values are different services and are reported side by side only; they are "
        "never averaged and the hosted value is not training-eligible"
    )
    rec["memory_peak"] = C.memory_stats(device)
    C.write_json(args.out, rec)


if __name__ == "__main__":
    main()
