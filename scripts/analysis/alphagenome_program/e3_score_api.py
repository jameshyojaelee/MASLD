#!/usr/bin/env python3
"""E3 step 2: score the selected MPRA oligos through the hosted AlphaGenome model API.

Two padding arms with identical oligo placement:
  npad   - the 126-bp oligo at a fixed window offset, every other base N, no genomic interval.
  native - the same 16,384-bp hg38 window with native flanks; only the oligo's variant base differs
           between the REF and ALT sequences.

The call pattern (create the client from `lib_atlas.load_api_key`, wrap every request in
`atlas_query.call_with_quota_retry`, summarise a TrackData over a fixed base mask that is folded onto
the track's own row count) is taken from `scripts/analysis/alphagenome_atlas/63_indel_rescue.py`.
Those modules are imported read-only; nothing under `alphagenome_atlas/` is written or edited.

Rate: a fixed minimum interval between calls holds the request rate at or below 30 per minute,
because the saturation job shares this API key.

Outputs (raw/): api_scores.jsonl (one JSON object per element, appended on completion; a restart
skips elements already present), api_probe.json when --probe is given.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys
import time

import numpy as np
import pandas as pd

PROJ = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATLAS_CODE = PROJ / "scripts/analysis/alphagenome_atlas"
FASTA = ("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/"
         "refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa")

OUTPUTS = ["ATAC", "DNASE", "RNA_SEQ"]
CHANNELS = {"atac": "atac", "dnase": "dnase", "rna": "rna_seq"}
LIVER_TERMS = ["UBERON:0002107", "UBERON:0001114", "UBERON:0001115", "CL:0000182"]
WINDOW = 16_384
OLIGO_LEN = 126
READOUT_PRIMARY = 128     # bases centred on the oligo centre
READOUT_SECONDARY = 1024
MIN_SECONDS_PER_CALL = 2.0   # 30 calls / minute ceiling
EPS = 1e-6


def fold_mask(mask: np.ndarray, n_rows: int) -> np.ndarray:
    """Fold a base-resolution mask onto a track's own row count (binned tracks would silently mis-select)."""
    m = np.asarray(mask, dtype=bool)
    if m.size == n_rows:
        return m
    if n_rows <= 0 or m.size % n_rows:
        raise RuntimeError(f"mask of {m.size} does not fold onto {n_rows} rows")
    return m.reshape(n_rows, m.size // n_rows).any(axis=1)


def summarise(output, masks: dict) -> dict:
    """Mean track value inside each readout window, per channel, over all returned liver tracks."""
    out = {}
    for name, attr in CHANNELS.items():
        td = getattr(output, attr, None)
        if td is None or td.values is None or np.asarray(td.values).size == 0:
            for w in masks:
                out[f"{name}_{w}"] = math.nan
            out[f"{name}_n_tracks"] = 0
            continue
        v = np.asarray(td.values, dtype=float)
        out[f"{name}_n_tracks"] = int(v.shape[1])
        for w, mask in masks.items():
            m = fold_mask(mask, v.shape[0])
            out[f"{name}_{w}"] = float(np.nanmean(v[m, :])) if m.any() else math.nan
    return out


def track_names(output) -> dict:
    names = {}
    for name, attr in CHANNELS.items():
        td = getattr(output, attr, None)
        md = getattr(td, "metadata", None) if td is not None else None
        if md is None:
            names[name] = []
            continue
        cols = [c for c in ("name", "ontology_curie", "biosample_name") if c in md]
        names[name] = md[cols].astype(str).agg(" | ".join, axis=1).tolist() if cols else []
    return names


def build_sequences(fasta, row) -> dict:
    """REF/ALT sequences for both arms, with the oligo at one fixed window offset."""
    contig = str(row["contig"])
    start0, end0 = int(row["oligo_start0"]), int(row["oligo_end0"])
    if end0 - start0 != OLIGO_LEN:
        raise RuntimeError(f"{row['element_id']}: oligo is {end0 - start0} bp, expected {OLIGO_LEN}")
    centre = start0 + OLIGO_LEN // 2
    win_start = centre - WINDOW // 2
    win_end = win_start + WINDOW
    native_ref = fasta.fetch(contig, win_start, win_end).upper()
    if len(native_ref) != WINDOW:
        raise RuntimeError(f"{row['element_id']}: window off the end of {contig}")

    vi = int(row["variant_pos0"]) - win_start          # variant offset inside the window
    oi = start0 - win_start                            # oligo offset inside the window
    if native_ref[vi] != str(row["genomic_ref"]).upper():
        raise RuntimeError(f"{row['element_id']}: genome has {native_ref[vi]} at the variant, "
                           f"element says {row['genomic_ref']}")
    alt_base = str(row["genomic_alt"]).upper()
    native_alt = native_ref[:vi] + alt_base + native_ref[vi + 1:]

    oligo_ref = native_ref[oi:oi + OLIGO_LEN]
    oligo_alt = native_alt[oi:oi + OLIGO_LEN]
    pad = "N" * oi
    tail = "N" * (WINDOW - oi - OLIGO_LEN)
    npad_ref = pad + oligo_ref + tail
    npad_alt = pad + oligo_alt + tail
    for s in (native_ref, native_alt, npad_ref, npad_alt):
        if len(s) != WINDOW:
            raise RuntimeError(f"{row['element_id']}: built sequence is {len(s)} bp")

    centre_off = centre - win_start
    masks = {}
    for tag, width in (("primary", READOUT_PRIMARY), ("secondary", READOUT_SECONDARY)):
        m = np.zeros(WINDOW, bool)
        m[max(0, centre_off - width // 2):min(WINDOW, centre_off + width // 2)] = True
        masks[tag] = m
    return {"contig": contig, "win_start": win_start, "win_end": win_end,
            "variant_offset": vi, "oligo_offset": oi, "centre_offset": centre_off,
            "native": {"ref": native_ref, "alt": native_alt},
            "npad": {"ref": npad_ref, "alt": npad_alt}, "masks": masks}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("out_root")
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--max-calls", type=int, default=1200)
    args = ap.parse_args()

    sys.path.insert(0, str(ATLAS_CODE))
    import pysam
    from alphagenome.data import genome
    from alphagenome.models import dna_client

    import atlas_query as aq          # quota-retry wrapper, reused unmodified
    import lib_atlas as la            # API key loader

    out_root = pathlib.Path(args.out_root)
    raw = out_root / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    sel = pd.read_csv(out_root / "tables/selected_elements.tsv", sep="\t")

    if WINDOW != dna_client.SUPPORTED_SEQUENCE_LENGTHS["SEQUENCE_LENGTH_16KB"]:
        raise SystemExit("16 kb is not the SDK's 16 kb")

    key, src = la.load_api_key()
    print(f"[e3] API key from {src}", flush=True)
    model = dna_client.create(key)
    outs = [getattr(dna_client.OutputType, o) for o in OUTPUTS]
    fasta = pysam.FastaFile(FASTA)

    state = {"calls": 0, "last": 0.0}

    def call(seq: str, interval, label: str):
        wait = MIN_SECONDS_PER_CALL - (time.time() - state["last"])
        if wait > 0:
            time.sleep(wait)
        if state["calls"] >= args.max_calls:
            raise SystemExit(f"[e3] call budget {args.max_calls} exhausted at {label}")
        o = aq.call_with_quota_retry(
            lambda: model.predict_sequence(sequence=seq, requested_outputs=outs,
                                           ontology_terms=LIVER_TERMS, interval=interval),
            label=label)
        state["calls"] += 1
        state["last"] = time.time()
        return o

    if args.probe:
        row = sel.iloc[0]
        b = build_sequences(fasta, row)
        iv = genome.Interval(chromosome=b["contig"], start=b["win_start"], end=b["win_end"])
        probe = {"element_id": str(row["element_id"]), "offsets": {k: int(b[k]) for k in
                 ("variant_offset", "oligo_offset", "centre_offset")}}
        for arm, interval in (("native", iv), ("npad", None)):
            try:
                o = call(b[arm]["ref"], interval, f"probe:{arm}")
            except Exception as exc:                       # noqa: BLE001 - recorded, then re-raised
                probe[arm] = {"state": "refused", "error": f"{type(exc).__name__}: {exc}"}
                continue
            probe[arm] = {"state": "ok", "summary": summarise(o, b["masks"]),
                          "tracks": track_names(o),
                          "rows": {n: int(np.asarray(getattr(o, a).values).shape[0])
                                   for n, a in CHANNELS.items() if getattr(o, a, None) is not None}}
        json.dump(probe, (raw / "api_probe.json").open("w"), indent=1, default=float)
        print(json.dumps(probe, indent=1, default=float)[:4000], flush=True)
        return

    ck = raw / "api_scores.jsonl"
    done = set()
    if ck.exists():
        with ck.open() as h:
            for line in h:
                line = line.strip()
                if line:
                    done.add(json.loads(line)["element_id"])
    print(f"[e3] {len(done)} elements already scored; {sel.shape[0] - len(done)} to go", flush=True)

    t0 = time.time()
    with ck.open("a") as h:
        for i, row in sel.iterrows():
            eid = str(row["element_id"])
            if eid in done:
                continue
            b = build_sequences(fasta, row)
            iv = genome.Interval(chromosome=b["contig"], start=b["win_start"], end=b["win_end"])
            rec = {"element_id": eid, "contig": b["contig"], "win_start": b["win_start"],
                   "variant_offset": b["variant_offset"], "oligo_offset": b["oligo_offset"],
                   "in_random_stratum": bool(row["in_random_stratum"]),
                   "in_topeffect_stratum": bool(row["in_topeffect_stratum"]),
                   "long_range_block": str(row["long_range_block"]),
                   "mb_block": str(row["mb_block"])}
            for arm, interval in (("native", iv), ("npad", None)):
                for allele in ("ref", "alt"):
                    o = call(b[arm][allele], interval, f"{eid}:{arm}:{allele}")
                    s = summarise(o, b["masks"])
                    for k, v in s.items():
                        rec[f"{arm}_{allele}_{k}"] = v
            h.write(json.dumps(rec, default=float) + "\n")
            h.flush()
            if state["calls"] % 40 == 0:
                rate = state["calls"] / max(1e-9, (time.time() - t0) / 60.0)
                print(f"[e3] {state['calls']} calls, {rate:.1f}/min, element {i + 1}/{sel.shape[0]}",
                      flush=True)
    print(f"[e3] done: {state['calls']} calls in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
