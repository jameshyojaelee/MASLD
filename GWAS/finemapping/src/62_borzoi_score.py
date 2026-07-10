#!/usr/bin/env python3
"""
62_borzoi_score.py  --  Borzoi (johahi/borzoi-pytorch) variant ref-vs-alt scorer.

Stage-2 seqfunc infrastructure smoke test / scoring primitive for the MASLD
finemapping pipeline. Loads a local open-weights Borzoi regulatory-sequence
model, extracts the model's receptive-field window (524,288 bp) centered on a
variant from the hg38 FASTA, one-hot encodes ref and alt sequences, runs both
through the model on GPU, and returns per-track ref-vs-alt deltas -- focused on
liver RNA-seq and hepatocyte/HepG2 accessibility tracks.

ADDITIVE ONLY. New file (script number 62). Does not touch any existing script.

Borzoi geometry (calico/borzoi, Linder et al. 2025):
  * input length      = 524,288 bp (2^19)
  * output bins       = 6,144  (each 32 bp; central 196,608 bp of the input)
  * human head tracks = 7,611  (indices match examples/targets_human.txt)
  * the variant sits at the exact center of the input, i.e. output bin n_bins//2.

Usage (defaults score PNPLA3 rs738409 + SORT1 rs12740374):
  python 62_borzoi_score.py \
      --fasta   /path/genome.fa \
      --targets /path/targets_human.txt \
      --liver-tracks /path/liver_tracks_selected.csv \
      --out-dir /path/results/seqfunc

Variant spec format for --variants (repeatable): NAME,CHR,POS,REF,ALT  (1-based, hg38)
"""

import argparse
import os
import sys
import time
import json

import numpy as np

# ---- Borzoi geometry (fixed by the trained model) --------------------------
SEQ_LEN = 524_288          # model receptive field
N_HUMAN_TRACKS = 7_611     # human head
BASES = "ACGT"             # Borzoi/Basenji one-hot channel order
BASE_IDX = {b: i for i, b in enumerate(BASES)}

# Default smoke-test variants (hg38, 1-based).
#   PNPLA3 rs738409  chr22:43928847 C>G  (I148M; the flagship MASLD risk allele)
#   SORT1  rs12740374 chr1:109274968 G>T (canonical STRONG liver regulatory eQTL;
#                                         minor T allele creates a C/EBP site ->
#                                         drives SORT1/CELSR2/PSRC1 in liver,
#                                         Musunuru et al. Nature 2010)
DEFAULT_VARIANTS = [
    ("PNPLA3_rs738409", "chr22", 43_928_847, "C", "G"),
    ("SORT1_rs12740374", "chr1", 109_274_968, "G", "T"),
]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def parse_variants(spec_list):
    out = []
    for s in spec_list:
        parts = s.split(",")
        if len(parts) != 5:
            raise ValueError(f"--variants entry must be NAME,CHR,POS,REF,ALT : got {s!r}")
        name, chrom, pos, ref, alt = parts
        out.append((name, chrom, int(pos), ref.upper(), alt.upper()))
    return out


def one_hot(seq):
    """seq (str, len SEQ_LEN) -> float32 array (4, SEQ_LEN), ACGT channel order.
    N / any non-ACGT base -> all-zero column (Borzoi convention)."""
    arr = np.zeros((4, len(seq)), dtype=np.float32)
    idx = np.frombuffer(seq.upper().encode("ascii"), dtype=np.uint8)
    for b, i in BASE_IDX.items():
        arr[i] = (idx == ord(b)).astype(np.float32)
    return arr


def fetch_window(fasta, chrom, pos0, seq_len):
    """Fetch seq_len bp centered so that the variant (0-based pos0) lands at the
    exact center index seq_len//2. Pads with 'N' at chromosome boundaries."""
    half = seq_len // 2
    start = pos0 - half           # may be negative near chrom start
    end = start + seq_len         # may exceed chrom length
    clen = fasta.get_reference_length(chrom)
    fstart = max(0, start)
    fend = min(clen, end)
    core = fasta.fetch(chrom, fstart, fend).upper()
    left_pad = fstart - start     # >0 if we clipped at the start
    right_pad = end - fend        # >0 if we clipped at the end
    seq = ("N" * left_pad) + core + ("N" * right_pad)
    assert len(seq) == seq_len, f"window len {len(seq)} != {seq_len}"
    return seq


def load_targets(targets_file):
    """Return list of (index, identifier, description) from targets_human.txt."""
    rows = []
    with open(targets_file) as fh:
        header = fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            rows.append((int(f[0]), f[1], f[8]))
    return rows


def load_liver_tracks(path):
    """liver_tracks_selected.csv -> list of dicts."""
    import csv
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            out.append({
                "track_index": int(r["track_index"]),
                "identifier": r["identifier"],
                "assay": r["assay"],
                "category": r["category"],
                "description": r["description"],
            })
    return out


def find_axes(out_shape, n_targets):
    """Given a (batch, A, B) prediction, decide which axis is tracks (==n_targets)
    and which is bins. Returns (tracks_axis, bins_axis) in 0-based dims of the
    non-batch tensor (i.e. after squeezing batch)."""
    a, b = out_shape[1], out_shape[2]
    if a == n_targets:
        return ("tracks_first", a, b)   # (batch, tracks, bins)
    if b == n_targets:
        return ("bins_first", a, b)     # (batch, bins, tracks)
    # fall back: assume larger axis is bins (6144 > 7611 is false, so tracks>bins)
    # Borzoi: bins=6144, tracks=7611 -> tracks is the larger. Pick larger as tracks.
    if a > b:
        return ("tracks_first", a, b)
    return ("bins_first", a, b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--targets", required=True,
                    help="Borzoi examples/targets_human.txt (7611 tracks)")
    ap.add_argument("--liver-tracks", required=True,
                    help="liver_tracks_selected.csv from track manifest")
    ap.add_argument("--model", default="johahi/borzoi-replicate-0",
                    help="HuggingFace model id for borzoi-pytorch")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--variants", action="append", default=None,
                    help="NAME,CHR,POS,REF,ALT (repeatable). Omit for defaults.")
    ap.add_argument("--win-bins", type=int, default=16,
                    help="half-width (in 32bp bins) of the local sum window "
                         "around the center bin for the windowed delta")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    variants = parse_variants(args.variants) if args.variants else DEFAULT_VARIANTS

    import torch
    import pysam

    dev = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"torch {torch.__version__} | device={dev} | "
        f"cuda_avail={torch.cuda.is_available()}")
    if dev == "cuda":
        log(f"GPU: {torch.cuda.get_device_name(0)} | "
            f"{torch.cuda.get_device_properties(0).total_memory/1e9:.1f} GB")

    # ---- load model --------------------------------------------------------
    log(f"Loading Borzoi weights: {args.model}")
    from borzoi_pytorch import Borzoi
    model = Borzoi.from_pretrained(args.model)
    model = model.to(dev).eval()
    log("Model loaded.")

    fasta = pysam.FastaFile(args.fasta)
    targets = load_targets(args.targets)
    n_targets = len(targets)
    log(f"targets_human.txt: {n_targets} tracks")
    desc_by_idx = {t[0]: t[2] for t in targets}
    ident_by_idx = {t[0]: t[1] for t in targets}
    liver = load_liver_tracks(args.liver_tracks)
    liver_idx = [t["track_index"] for t in liver]
    log(f"liver-relevant tracks: {len(liver_idx)} "
        f"({sum(t['assay']=='accessibility' for t in liver)} accessibility / "
        f"{sum(t['assay']=='rna_liver' for t in liver)} rna_liver)")

    def predict(seq_str):
        """Return prediction as numpy (n_bins, n_tracks)."""
        oh = one_hot(seq_str)                      # (4, SEQ_LEN)
        x = torch.from_numpy(oh).unsqueeze(0).to(dev)   # (1,4,L)
        with torch.no_grad():
            try:
                out = model(x)
            except torch.cuda.OutOfMemoryError:
                log("  fp32 OOM -> retrying under bf16 autocast")
                torch.cuda.empty_cache()
                with torch.autocast("cuda", dtype=torch.bfloat16):
                    out = model(x)
        if isinstance(out, (tuple, list)):
            out = out[0]                            # human head
        out = out.float().cpu()
        # squeeze batch
        out = out[0]
        mode, a, b = find_axes((1,) + tuple(out.shape), n_targets)
        if mode == "bins_first":       # (bins, tracks)
            pred = out.numpy()
        else:                           # (tracks, bins) -> transpose
            pred = out.numpy().T
        return pred  # (n_bins, n_tracks)

    all_rows = []
    summary = {}
    for name, chrom, pos1, ref, alt in variants:
        log(f"=== {name}  {chrom}:{pos1} {ref}>{alt} ===")
        pos0 = pos1 - 1
        seq = fetch_window(fasta, chrom, pos0, SEQ_LEN)
        center = SEQ_LEN // 2
        fasta_ref = seq[center]
        note = ""
        # defensive: FASTA base is ground truth; guarantee a real substitution
        eff_ref = fasta_ref
        eff_alt = alt
        if fasta_ref != ref:
            note = f"WARN provided_ref={ref} != fasta_ref={fasta_ref}; using fasta base"
            log("  " + note)
        if eff_alt == eff_ref:
            eff_alt = next(b for b in BASES if b != eff_ref)
            note += f"; alt==ref -> using {eff_alt}"
            log(f"  alt equals ref base; substituting alt -> {eff_alt}")

        ref_seq = seq
        alt_seq = seq[:center] + eff_alt + seq[center + 1:]

        t0 = time.time()
        ref_pred = predict(ref_seq)
        alt_pred = predict(alt_seq)
        log(f"  predicted ref+alt in {time.time()-t0:.1f}s | "
            f"pred shape (bins,tracks)={ref_pred.shape}")

        n_bins = ref_pred.shape[0]
        cbin = n_bins // 2
        lo = max(0, cbin - args.win_bins)
        hi = min(n_bins, cbin + args.win_bins + 1)

        delta_center = alt_pred[cbin] - ref_pred[cbin]                 # (n_tracks,)
        delta_win = (alt_pred[lo:hi] - ref_pred[lo:hi]).sum(axis=0)    # (n_tracks,)

        max_abs = float(np.max(np.abs(delta_win)))
        top = np.argsort(-np.abs(delta_win))[:10]
        summary[name] = {
            "chr": chrom, "pos": pos1, "ref": ref, "alt": alt,
            "eff_ref": eff_ref, "eff_alt": eff_alt, "note": note,
            "n_bins": int(n_bins), "n_tracks": int(ref_pred.shape[1]),
            "center_bin": int(cbin), "win_bins": args.win_bins,
            "max_abs_windowed_delta": max_abs,
            "nonzero": bool(max_abs > 0),
            "top_tracks": [
                {"idx": int(i), "id": ident_by_idx.get(int(i), "?"),
                 "desc": desc_by_idx.get(int(i), "?"),
                 "delta_win": float(delta_win[i])}
                for i in top
            ],
        }

        # per liver track rows
        for t in liver:
            i = t["track_index"]
            all_rows.append({
                "variant": name, "chr": chrom, "pos": pos1,
                "ref": ref, "alt": alt, "eff_ref": eff_ref, "eff_alt": eff_alt,
                "track_index": i, "identifier": t["identifier"],
                "assay": t["assay"], "category": t["category"],
                "description": t["description"],
                "ref_center": float(ref_pred[cbin, i]),
                "alt_center": float(alt_pred[cbin, i]),
                "delta_center": float(delta_center[i]),
                "delta_win": float(delta_win[i]),
            })

        # console: liver accessibility + a couple RNA rows
        acc = [t for t in liver if t["assay"] == "accessibility"]
        log(f"  liver ACCESSIBILITY deltas (windowed, {2*args.win_bins+1} bins):")
        for t in acc:
            i = t["track_index"]
            log(f"    [{i:>4}] {t['description']:<42} "
                f"ref={ref_pred[cbin,i]:9.4f} alt={alt_pred[cbin,i]:9.4f} "
                f"dwin={delta_win[i]:+.5f}")
        # mean over liver RNA tracks
        rna_i = [t["track_index"] for t in liver if t["assay"] == "rna_liver"]
        if rna_i:
            log(f"  liver RNA mean windowed delta over {len(rna_i)} tracks: "
                f"{float(np.mean(delta_win[rna_i])):+.5f} "
                f"(max abs {float(np.max(np.abs(delta_win[rna_i]))):.5f})")
        log(f"  MAX |windowed delta| over all 7611 tracks = {max_abs:.5f}  "
            f"=> ref{'!=' if max_abs>0 else '=='}alt")

    # ---- write outputs -----------------------------------------------------
    import csv
    out_csv = os.path.join(args.out_dir, "borzoi_smoke_liver_deltas.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()))
        w.writeheader()
        w.writerows(all_rows)
    out_json = os.path.join(args.out_dir, "borzoi_smoke_summary.json")
    with open(out_json, "w") as fh:
        json.dump(summary, fh, indent=2)
    log(f"WROTE {out_csv} ({len(all_rows)} rows)")
    log(f"WROTE {out_json}")

    ok = all(v["nonzero"] for v in summary.values())
    log(f"SMOKE TEST {'PASS' if ok else 'FAIL'}: "
        f"all variants produced non-zero ref-vs-alt delta = {ok}")
    sys.exit(0 if ok else 2)


if __name__ == "__main__":
    main()
