#!/usr/bin/env python3
"""Stage 1: time and memory for ONE input length, in its own process.

One process per length is deliberate. `peak_bytes_in_use` is a high-water mark that never falls, so a
single process that walks 2 kb -> 1 Mb reports the 1 Mb peak for every length after it. Restarting gives
each length its own peak.

Usage: h1_bench_length.py --length 1048576 --out <dir>/bench_1048576.json [--repeats 3]
"""

from __future__ import annotations

import argparse
import sys
import time
import traceback

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import h1_runtime_common as C  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--length", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--repeats", type=int, default=3)
    args = ap.parse_args()

    import jax
    import pysam
    from alphagenome.data import genome

    rec = {
        "length_bp": args.length,
        "jax_version": jax.__version__,
        "jax_devices": [str(d) for d in jax.local_devices()],
        "requested_outputs": C.OUTPUT_NAMES,
        "ontology_terms": C.LIVER_TERMS,
    }

    model, device, load_s = C.load_model()
    rec["checkpoint_load_seconds"] = round(load_s, 2)
    rec["memory_after_load"] = C.memory_stats(device)

    outs = C.requested_outputs()
    fasta = pysam.FastaFile(C.FASTA_PATH)

    # Real genomic sequence anchored at the HSD17B13 locus so the benchmark runs on the same kind of
    # input the anchor check uses. chr4:87310240 is the variant position.
    pos1 = 87_310_240
    start0 = max(0, pos1 - 1 - args.length // 2)
    end = start0 + args.length
    seq = fasta.fetch("chr4", start0, end).upper()
    if len(seq) != args.length:
        rec["state"] = f"fasta returned {len(seq)} bp, wanted {args.length}"
        C.write_json(args.out, rec)
        return
    interval = genome.Interval(chromosome="chr4", start=start0, end=end)

    def call():
        return model.predict_sequence(
            sequence=seq,
            requested_outputs=outs,
            ontology_terms=C.LIVER_TERMS,
            interval=interval,
        )

    try:
        t0 = time.time()
        out = call()
        rec["first_call_seconds_includes_compile"] = round(time.time() - t0, 3)
        rec["memory_after_first_call"] = C.memory_stats(device)

        shapes = {}
        for name, attr in C.CHANNELS:
            td = getattr(out, attr, None)
            if td is None or td.values is None:
                shapes[name] = None
                continue
            shapes[name] = {
                "rows": int(td.values.shape[0]),
                "tracks": int(td.values.shape[1]),
                "resolution_bp": int(td.resolution),
                "rows_times_resolution": int(td.values.shape[0]) * int(td.resolution),
                "dtype": str(td.values.dtype),
            }
        rec["output_shapes"] = shapes
        rec["output_covers_full_input"] = {
            k: (v["rows_times_resolution"] == args.length) for k, v in shapes.items() if v
        }

        times = []
        for _ in range(args.repeats):
            t0 = time.time()
            o = call()
            del o
            times.append(time.time() - t0)
        rec["warm_call_seconds"] = [round(t, 3) for t in times]
        rec["warm_call_seconds_median"] = round(float(np.median(times)), 3)
        rec["memory_peak"] = C.memory_stats(device)
        rec["state"] = "ok"
    except Exception as exc:  # noqa: BLE001 - the point is to record the failure mode
        rec["state"] = "failed"
        rec["error_type"] = type(exc).__name__
        rec["error"] = str(exc)[:4000]
        rec["traceback_tail"] = traceback.format_exc()[-2000:]
        rec["memory_peak"] = C.memory_stats(device)

    C.write_json(args.out, rec)


if __name__ == "__main__":
    main()
