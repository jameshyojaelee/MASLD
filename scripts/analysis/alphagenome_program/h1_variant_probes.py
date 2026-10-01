#!/usr/bin/env python3
"""Stage 2: prove the local runtime behaves as a variant-effect model before anyone uses it.

Four checks, each with its own pass criterion fixed here before the run:

(a) REF vs ALT on the SAME window give different outputs, and the difference is localized near the
    variant. Substitution only: an indel shifts every downstream base, so the whole downstream half of
    the window legitimately differs and localization cannot be read off it.
    PASS: max |alt-ref| > 0 in at least one channel AND the share of total |alt-ref| mass inside
    +/-2 kb of the variant exceeds the uniform-window expectation by at least 10x.

(b) Output coordinates align with input coordinates. Two independent checks.
    b1 metadata: every returned track's rows * resolution equals the input interval width, and the
       interval carried on the output equals the one passed in.
    b2 empirical: a 1,024-bp scramble is placed at three different offsets in the SAME window. The
       argmax of the per-position |alt-ref| response must MOVE with the offset and land inside the
       scrambled block. A test at one offset cannot separate "aligned" from "always reports the
       centre"; three offsets can.
    PASS: b1 exact for all channels, and b2 argmax within +/-1 output bin of the block at all three
    offsets for the strongest-responding 1-bp channel.

(c) Two runs of the same input agree to floating-point tolerance.
    PASS: max |run1 - run2| == 0 across every returned track. Reported as measured either way.

(d) An all-N or malformed input fails loudly rather than silently returning something usable.
    PASS requires an exception. Recorded honestly if it does not raise.
"""

from __future__ import annotations

import argparse
import sys
import traceback

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import h1_runtime_common as C  # noqa: E402

CHROM = "chr4"
POS1 = 87_310_240  # HSD17B13 rs72613567 in hg38; the probes reuse this locus
SCRAMBLE_OFFSETS = (300_000, 700_000, 900_000)
SCRAMBLE_WIDTH = 1_024
SCRAMBLE_SEED = 20260915


def per_row_response(out_ref, out_alt) -> dict:
    """Per-position summed |alt - ref| for every channel, with its resolution."""
    resp = {}
    for name, attr in C.CHANNELS:
        a = getattr(out_alt, attr, None)
        b = getattr(out_ref, attr, None)
        if a is None or b is None or a.values is None or b.values is None:
            continue
        if a.values.shape != b.values.shape:
            resp[name] = {"error": f"shape mismatch {a.values.shape} vs {b.values.shape}"}
            continue
        d = np.abs(np.asarray(a.values, dtype=np.float64) - np.asarray(b.values, dtype=np.float64))
        resp[name] = {
            "row_sum": d.sum(axis=1),
            "resolution": int(a.resolution),
            "total": float(d.sum()),
            "max": float(d.max()),
        }
    return resp


def scramble(seq: str, center_offset: int, width: int, seed: int) -> str:
    rng = np.random.default_rng(seed)
    lo = center_offset - width // 2
    hi = lo + width
    block = "".join(rng.choice(list("ACGT"), size=width))
    return seq[:lo] + block + seq[hi:]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--length", type=int, default=1_048_576)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    import pysam
    from alphagenome.data import genome

    rec = {"probe_length_bp": args.length, "locus": f"{CHROM}:{POS1}"}
    model, device, load_s = C.load_model()
    rec["checkpoint_load_seconds"] = round(load_s, 2)
    outs = C.requested_outputs()
    fasta = pysam.FastaFile(C.FASTA_PATH)

    start0 = max(0, POS1 - 1 - args.length // 2)
    end = start0 + args.length
    ref_seq = fasta.fetch(CHROM, start0, end).upper()
    interval = genome.Interval(chromosome=CHROM, start=start0, end=end)
    var_offset = POS1 - 1 - start0
    rec["window"] = {"start0": start0, "end": end, "variant_offset_in_window": var_offset}
    rec["reference_base_at_variant"] = ref_seq[var_offset]

    def call(seq, iv=interval):
        return model.predict_sequence(
            sequence=seq, requested_outputs=outs, ontology_terms=C.LIVER_TERMS, interval=iv
        )

    out_ref = call(ref_seq)

    # ---------------------------------------------------------------- (c) determinism, run first so a
    # later failure still leaves it recorded.
    out_ref2 = call(ref_seq)
    det = {}
    for name, attr in C.CHANNELS:
        a = getattr(out_ref, attr, None)
        b = getattr(out_ref2, attr, None)
        if a is None or b is None or a.values is None:
            det[name] = None
            continue
        d = np.abs(np.asarray(a.values, np.float64) - np.asarray(b.values, np.float64))
        det[name] = {
            "max_abs_diff": float(d.max()),
            "n_differing": int((d != 0).sum()),
            "n_values": int(d.size),
        }
    rec["check_c_determinism"] = {
        "per_channel": det,
        "max_abs_diff_overall": max(
            (v["max_abs_diff"] for v in det.values() if v), default=float("nan")
        ),
        "passed": all(v["max_abs_diff"] == 0.0 for v in det.values() if v),
    }
    del out_ref2

    # ---------------------------------------------------------------- (b1) metadata alignment
    b1 = {}
    for name, attr in C.CHANNELS:
        td = getattr(out_ref, attr, None)
        if td is None or td.values is None:
            b1[name] = None
            continue
        b1[name] = {
            "rows": int(td.values.shape[0]),
            "resolution_bp": int(td.resolution),
            "rows_times_resolution": int(td.values.shape[0]) * int(td.resolution),
            "input_width": args.length,
            "covers_input_exactly": int(td.values.shape[0]) * int(td.resolution) == args.length,
            "interval_on_output": str(td.interval),
            "interval_matches_input": (
                td.interval is not None
                and td.interval.chromosome == interval.chromosome
                and int(td.interval.start) == int(interval.start)
                and int(td.interval.end) == int(interval.end)
            ),
        }
    rec["check_b1_metadata_alignment"] = {
        "per_channel": b1,
        "passed": all(
            v["covers_input_exactly"] and v["interval_matches_input"] for v in b1.values() if v
        ),
    }

    # ---------------------------------------------------------------- (a) REF vs ALT, substitution
    alt_base = {"T": "G", "A": "C", "C": "A", "G": "T"}[ref_seq[var_offset]]
    snv_seq = ref_seq[:var_offset] + alt_base + ref_seq[var_offset + 1 :]
    out_snv = call(snv_seq)
    resp = per_row_response(out_ref, out_snv)
    a_rec = {
        "substitution": f"{ref_seq[var_offset]}>{alt_base} at {CHROM}:{POS1}",
        "note": "a substitution, not the rs72613567 insertion: an indel shifts all downstream bases",
        "per_channel": {},
    }
    for name, r in resp.items():
        if "row_sum" not in r:
            a_rec["per_channel"][name] = r
            continue
        rs = r["row_sum"]
        res = r["resolution"]
        vrow = var_offset // res
        for flank in (2_000, 10_000):
            half = max(1, flank // res)
            lo, hi = max(0, vrow - half), min(rs.size, vrow + half)
            share = float(rs[lo:hi].sum() / rs.sum()) if rs.sum() > 0 else float("nan")
            uniform = float((hi - lo) / rs.size)
            a_rec["per_channel"].setdefault(name, {})[f"share_within_{flank}bp"] = share
            a_rec["per_channel"][name][f"uniform_expectation_{flank}bp"] = uniform
            a_rec["per_channel"][name][f"enrichment_{flank}bp"] = (
                share / uniform if uniform > 0 else float("nan")
            )
        a_rec["per_channel"][name]["max_abs_diff"] = r["max"]
        a_rec["per_channel"][name]["total_abs_diff"] = r["total"]
        a_rec["per_channel"][name]["resolution_bp"] = res
    any_diff = any(r.get("max", 0) > 0 for r in resp.values() if "max" in r)
    enrich = [
        v["enrichment_2000bp"]
        for v in a_rec["per_channel"].values()
        if isinstance(v, dict) and v.get("enrichment_2000bp") == v.get("enrichment_2000bp")
    ]
    a_rec["outputs_differ"] = bool(any_diff)
    a_rec["min_enrichment_2000bp_across_channels"] = float(min(enrich)) if enrich else float("nan")
    a_rec["passed"] = bool(any_diff and enrich and min(enrich) >= 10.0)
    rec["check_a_ref_alt_localized"] = a_rec
    del out_snv

    # ---------------------------------------------------------------- (b2) empirical alignment
    b2 = {"scramble_width_bp": SCRAMBLE_WIDTH, "seed": SCRAMBLE_SEED, "offsets": {}}
    for off in SCRAMBLE_OFFSETS:
        if off + SCRAMBLE_WIDTH // 2 > args.length:
            continue
        s = scramble(ref_seq, off, SCRAMBLE_WIDTH, SCRAMBLE_SEED)
        o = call(s)
        r = per_row_response(out_ref, o)
        per = {}
        for name, rr in r.items():
            if "row_sum" not in rr:
                continue
            rs, res = rr["row_sum"], rr["resolution"]
            am = int(np.argmax(rs))
            per[name] = {
                "resolution_bp": res,
                "argmax_row": am,
                "argmax_offset_bp": am * res,
                "argmax_genomic_pos": start0 + am * res,
                "distance_to_block_center_bp": int(am * res - off),
                "inside_block": bool(abs(am * res - off) <= (SCRAMBLE_WIDTH // 2 + res)),
                "total_abs_diff": rr["total"],
                "max_abs_diff": rr["max"],
            }
        strongest = max(per, key=lambda k: per[k]["total_abs_diff"]) if per else None
        b2["offsets"][str(off)] = {
            "block_center_offset_bp": off,
            "block_center_genomic_pos": start0 + off,
            "per_channel": per,
            "strongest_channel": strongest,
            "strongest_inside_block": per[strongest]["inside_block"] if strongest else None,
        }
        del o
    b2["passed"] = bool(
        b2["offsets"]
        and all(v["strongest_inside_block"] for v in b2["offsets"].values())
    )
    rec["check_b2_empirical_alignment"] = b2

    # ---------------------------------------------------------------- (d) malformed inputs
    d_rec = {}
    cases = {
        # (sequence, interval passed alongside it)
        "all_N_full_length": ("N" * args.length, None),
        "invalid_character_Z_10000_positions": ("Z" * 10_000 + ref_seq[10_000:], None),
        "unsupported_length_1000bp": (ref_seq[:1000], None),
        "empty_sequence": ("", None),
        "length_off_by_one": (ref_seq[: args.length - 1], None),
        # The same off-by-one sequence but WITH the full-width interval, to record whether the
        # interval argument is itself a guard a caller can rely on.
        "length_off_by_one_with_interval": (ref_seq[: args.length - 1], interval),
        "all_N_with_interval": ("N" * args.length, interval),
    }
    for label, (seq, iv) in cases.items():
        entry = {"input_length": len(seq), "interval_passed": iv is not None}
        try:
            o = call(seq, iv=iv)
            entry["raised"] = False
            entry["outcome"] = "returned an Output without raising"
            summ = {}
            for name, attr in C.CHANNELS:
                td = getattr(o, attr, None)
                if td is None or td.values is None:
                    summ[name] = None
                    continue
                v = np.asarray(td.values, np.float64)
                ref_td = getattr(out_ref, attr, None)
                same_shape = ref_td is not None and ref_td.values.shape == v.shape
                summ[name] = {
                    "shape": list(v.shape),
                    "n_nan": int(np.isnan(v).sum()),
                    "all_zero": bool(np.all(v == 0)),
                    "mean_abs": float(np.nanmean(np.abs(v))),
                    "mean_abs_reference_run": (
                        float(np.nanmean(np.abs(np.asarray(ref_td.values, np.float64))))
                        if ref_td is not None and ref_td.values is not None
                        else None
                    ),
                    "identical_to_reference_run": (
                        bool(np.array_equal(v, np.asarray(ref_td.values, np.float64)))
                        if same_shape
                        else False
                    ),
                }
            entry["output_summary"] = summ
            del o
        except Exception as exc:  # noqa: BLE001
            entry["raised"] = True
            entry["error_type"] = type(exc).__name__
            entry["error"] = str(exc)[:1500]
            entry["traceback_tail"] = traceback.format_exc()[-800:]
        entry["passed"] = bool(entry["raised"])
        d_rec[label] = entry
    d_rec["passed"] = all(v["passed"] for v in d_rec.values() if isinstance(v, dict))
    rec["check_d_malformed_input"] = d_rec

    rec["memory_peak"] = C.memory_stats(device)
    rec["summary"] = {
        k: rec[k]["passed"]
        for k in (
            "check_a_ref_alt_localized",
            "check_b1_metadata_alignment",
            "check_b2_empirical_alignment",
            "check_c_determinism",
            "check_d_malformed_input",
        )
    }
    C.write_json(args.out, rec)


if __name__ == "__main__":
    main()
