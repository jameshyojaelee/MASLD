#!/usr/bin/env python3
"""C1 endpoint 2, step 2: Borzoi liver accessibility allele scores for the targets wave 1 missed.

Reuses the exact recipe of the archived src/79 Borzoi caQTL panel -- same port
(`johahi/borzoi-replicate-0`), same 524,288-bp window centred on the variant, same 2 ATAC and
5 DNase liver tracks from tracks/liver_tracks_selected.csv, same 16-bin-either-side summation,
same resolve_orient -> orient_mult handling that puts every delta in the label REF->ALT frame --
but writes into this package's directory instead of the read-only seqfunc results tree.

Output is appended in checkpoints so a truncated run still delivers the priority-0 and
priority-1 targets (the union allelic-imbalance sites and the Currin leads carrying an archived
AlphaGenome score, which are the ones the four-model matched set needs).
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import time

import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SRC = os.path.join(ROOT, "GWAS/finemapping/src")
SF = os.path.join(ROOT, "GWAS/finemapping/results/seqfunc")
FASTA = "/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
TRACKS = os.path.join(SF, "tracks/targets_human.txt")
LIVER_TRACKS = os.path.join(SF, "tracks/liver_tracks_selected.csv")

COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def load_module(num_name):
    path = os.path.join(SRC, num_name + ".py")
    spec = importlib.util.spec_from_file_location(num_name.replace("-", "_"), path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def resolve_orient(fasta_base, ref, alt):
    """(hg38_ref, hg38_alt, orient_mult, note); orient_mult maps a fasta ref->alt signed delta
    into the LABEL ref->alt frame. Identical to src/79_score_panels.resolve_orient."""
    ref, alt, b = ref.upper(), alt.upper(), fasta_base.upper()
    if b == ref:
        return ref, alt, 1, ""
    if b == alt:
        return alt, ref, -1, "hg19_ref_alt_swapped_on_hg38"
    if b == COMP.get(ref):
        return COMP[ref], COMP[alt], 1, "strand_flip"
    if b == COMP.get(alt):
        return COMP[alt], COMP[ref], -1, "strand_flip_swapped"
    return b, None, 0, f"UNRESOLVED_fasta={b}_alleles={ref}/{alt}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--model", default="johahi/borzoi-replicate-0")
    ap.add_argument("--win-bins", type=int, default=16)
    ap.add_argument("--checkpoint-every", type=int, default=500)
    ap.add_argument("--max-variants", type=int, default=0)
    args = ap.parse_args()

    m62 = load_module("62_borzoi_score")
    import torch
    import pysam
    # Hard GPU guard. The borzoi env ships torch cu130, which refuses to initialise against the
    # 12.4 driver on the L40S nodes and falls back to CPU; two 524,288-bp forward passes per
    # variant on CPU would burn the whole wall clock and produce nothing. Abort instead.
    if not torch.cuda.is_available():
        raise SystemExit(
            "FATAL: CUDA unavailable in this environment on this node. torch "
            f"{torch.__version__} needs a newer driver than this node carries; submit with "
            "--constraint=b6k. Refusing to run on CPU.")
    dev = "cuda"
    gpu = torch.cuda.get_device_name(0)
    log(f"torch {torch.__version__} device={dev} gpu={gpu} dtype=float32")
    from borzoi_pytorch import Borzoi
    model = Borzoi.from_pretrained(args.model).to(dev).eval()
    fasta = pysam.FastaFile(FASTA)
    targets = m62.load_targets(TRACKS)
    n_targets = len(targets)
    liver = m62.load_liver_tracks(LIVER_TRACKS)
    atac_idx = [t["track_index"] for t in liver if t["category"] == "ATAC"]
    dnase_idx = [t["track_index"] for t in liver if t["category"] == "DNASE"]
    log(f"liver accessibility tracks: {len(atac_idx)} ATAC {atac_idx}, "
        f"{len(dnase_idx)} DNase {dnase_idx}")

    def predict(seq):
        # float32 only. The archived src/79 panel ran float32; a silent bfloat16 fallback would
        # put part of the Borzoi column on a different numeric path from the archived part of the
        # same column, so an out-of-memory is raised instead of being absorbed.
        oh = m62.one_hot(seq)
        x = torch.from_numpy(oh).unsqueeze(0).to(dev)
        with torch.no_grad():
            try:
                out = model(x)
            except torch.cuda.OutOfMemoryError as e:
                torch.cuda.empty_cache()
                raise SystemExit(
                    "FATAL: out of memory in float32 on this GPU. The archived panel ran "
                    "float32; rerun on a larger GPU rather than silently switching precision. "
                    f"({e})")
        if isinstance(out, (tuple, list)):
            out = out[0]
        out = out.float().cpu()[0]
        mode, a, b = m62.find_axes((1,) + tuple(out.shape), n_targets)
        return out.numpy() if mode == "bins_first" else out.numpy().T

    df = pd.read_csv(args.targets, sep="\t", dtype={"chr": str})
    if args.max_variants > 0:
        df = df.head(args.max_variants)
    log(f"targets {len(df)}; priority counts "
        f"{df.priority.value_counts().sort_index().to_dict()}")

    SEQ_LEN = m62.SEQ_LEN
    rows, n_ok, n_unres = [], 0, 0
    t0 = time.time()
    header_written = False
    if os.path.exists(args.out):
        os.remove(args.out)

    def flush():
        nonlocal rows, header_written
        if not rows:
            return
        pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False,
                                  mode="a", header=not header_written)
        header_written = True
        rows = []

    for i, r in enumerate(df.itertuples()):
        chrom = "chr" + str(r.chr)
        pos1 = int(r.pos_hg38)
        pos0 = pos1 - 1
        ref, alt = str(r.ref).upper(), str(r.alt).upper()
        base = dict(canon=r.canon, chr=str(r.chr), pos_hg38=pos1, label_ref=ref, label_alt=alt,
                    target=r.target, priority=int(r.priority), orient_mult=1, note="",
                    fasta_base="", borzoi_atac_delta=np.nan, borzoi_dnase_delta=np.nan)
        try:
            seq = m62.fetch_window(fasta, chrom, pos0, SEQ_LEN)
        except Exception as e:  # noqa: BLE001
            base["note"] = f"fetch_fail:{e}"
            rows.append(base)
            continue
        center = SEQ_LEN // 2
        fbase = seq[center].upper()
        hg38_ref, hg38_alt, om, note = resolve_orient(fbase, ref, alt)
        base["orient_mult"] = om
        base["note"] = note
        base["fasta_base"] = fbase
        if hg38_alt is None:
            rows.append(base)
            n_unres += 1
            continue
        alt_seq = seq[:center] + hg38_alt + seq[center + 1:]
        rp = predict(seq)
        ap_ = predict(alt_seq)
        n_bins = rp.shape[0]
        cbin = n_bins // 2
        lo = max(0, cbin - args.win_bins)
        hi = min(n_bins, cbin + args.win_bins + 1)
        dwin = (ap_[lo:hi] - rp[lo:hi]).sum(axis=0)     # fasta ref -> fasta alt
        base["borzoi_atac_delta"] = float(np.nanmean([dwin[k] for k in atac_idx]))
        base["borzoi_dnase_delta"] = float(np.nanmean([dwin[k] for k in dnase_idx]))
        n_ok += 1
        rows.append(base)
        if (i + 1) % args.checkpoint_every == 0:
            flush()
            el = time.time() - t0
            log(f"  [{i+1}/{len(df)}] ok={n_ok} unres={n_unres} {el:.0f}s "
                f"({el/(i+1):.2f}s/variant)")
    flush()
    log(f"wrote {args.out}: ok={n_ok} unresolved={n_unres} of {len(df)}")


if __name__ == "__main__":
    main()
