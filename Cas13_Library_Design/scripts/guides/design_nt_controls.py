#!/usr/bin/env python
"""
design_nt_controls.py — design N non-targeting (NT) Cas13 control guides that
DO NOT target the mouse transcriptome.

Strategy (mirrors standard Cas13/CRISPR NT-control design):
  1. Generate random GC-balanced 23-mers with sequence guards (no homopolymer
     run >=4, GC in [0.35, 0.65], no BsmBI/BbsI sites, mutually unique).
  2. Screen EVERY candidate against the full GENCODE vM38 transcriptome (all
     transcripts: mRNA + lncRNA + everything) with Bowtie v1 `-v 3` (<=3 mismatch,
     both strands). RfxCas13d targets mature mRNA, so the transcriptome is the
     complete off-target space. A candidate is KEPT only if it has ZERO alignment
     -- i.e. its nearest transcript differs by >=4 mismatches, far beyond Cas13's
     cleavage tolerance. Random 23-mers are ~17-18 mismatches from any transcript,
     so survivors carry a large safety margin.
  3. Emit the first N survivors in the control-guide schema (control_type
     == 'non_targeting'), so they concatenate with cas13_control_guides_vM38.csv.

The Bowtie index of the transcriptome is built once by run_nt_design.sh and passed
via --idx. Bowtie must be on PATH (module load Bowtie). Deterministic (seed=42).
"""
from __future__ import annotations

import argparse
import random
import subprocess
import sys
from pathlib import Path

import pandas as pd

GUIDE_LEN = 23
GC_LO, GC_HI = 0.35, 0.65
# Disallowed: homopolymer run >=4 of any base; type-IIS cloning sites (both strands)
# commonly used for Cas13 library assembly (BsmBI/Esp3I CGTCTC, BbsI GAAGAC).
RESTRICTION = ("CGTCTC", "GAGACG", "GAAGAC", "GTCTTC")
RUNS = ("AAAA", "CCCC", "GGGG", "TTTT")

# Control-guide schema (must match cas13_control_guides_vM38.csv column order).
SCHEMA = [
    "guide_id", "gene_id_mouse", "gene_symbol_mouse", "biotype", "guide_seq",
    "target_seq", "region", "tiger_score", "cas13_score", "combined_score",
    "n_isoforms_targeted", "single_isoform", "tx_id_set", "tx_id_pos", "position",
    "n_target", "any_indel", "mismatch", "rank_within_gene", "n_available_pool",
    "human_symbol", "expected_direction", "control_type",
    "essentiality_chronos_liver", "n_liver_lines", "is_essential_liver",
]


def _ok(seq: str) -> bool:
    gc = (seq.count("G") + seq.count("C")) / len(seq)
    if not (GC_LO <= gc <= GC_HI):
        return False
    if any(r in seq for r in RUNS):
        return False
    if any(r in seq for r in RESTRICTION):
        return False
    return True


def generate_candidates(n_cand: int, rng: random.Random) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    tries = 0
    while len(out) < n_cand and tries < n_cand * 200:
        tries += 1
        s = "".join(rng.choice("ACGT") for _ in range(GUIDE_LEN))
        if s in seen or not _ok(s):
            continue
        seen.add(s)
        out.append(s)
    return out


def bowtie_aligned(cands: list[str], idx: str, workdir: Path, threads: int) -> set[str]:
    """Return the set of candidate ids that align to the transcriptome with <=3
    mismatches (both strands). Bowtie v1 prints a line only for aligned reads."""
    fa = workdir / "nt_candidates.fa"
    with fa.open("w") as fh:
        for i, s in enumerate(cands):
            fh.write(f">c{i}\n{s}\n")
    hits = workdir / "nt_hits.txt"
    # -v 3: <=3 mismatch end-to-end; -k 1: stop at first hit (presence is enough);
    # --suppress 2..8: print only the read name. Both strands searched by default.
    cmd = ["bowtie", "-v", "3", "-k", "1", "-p", str(threads), "-f",
           "--suppress", "2,3,4,5,6,7,8", idx, str(fa), str(hits)]
    print("[nt] bowtie:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    aligned = set()
    with hits.open() as fh:
        for line in fh:
            name = line.strip()
            if name:
                aligned.add(name)
    return aligned


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--idx", required=True, help="Bowtie transcriptome index prefix")
    ap.add_argument("--n", type=int, default=500, help="number of NT guides to emit")
    ap.add_argument("--out", required=True, help="output CSV (NT control guides)")
    ap.add_argument("--workdir", required=True, help="scratch dir for fasta/hits")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--oversample", type=int, default=8,
                    help="candidate multiple of --n to generate before screening")
    args = ap.parse_args()
    workdir = Path(args.workdir); workdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(42)

    n_cand = args.n * args.oversample
    cands = generate_candidates(n_cand, rng)
    print(f"[nt] generated {len(cands)} guarded candidates (GC {GC_LO}-{GC_HI}, "
          f"no run>=4, no BsmBI/BbsI)", flush=True)

    aligned = bowtie_aligned(cands, args.idx, workdir, args.threads)
    survivors = [s for i, s in enumerate(cands) if f"c{i}" not in aligned]
    print(f"[nt] transcriptome screen: {len(aligned)} candidates hit a transcript "
          f"(<=3 mm) -> rejected; {len(survivors)} non-targeting survivors", flush=True)
    if len(survivors) < args.n:
        sys.exit(f"ERROR: only {len(survivors)} survivors < requested {args.n}; "
                 f"raise --oversample.")
    keep = survivors[: args.n]

    rows = []
    for i, s in enumerate(keep, start=1):
        rows.append({
            "guide_id": f"NT_{i:04d}", "gene_id_mouse": "",
            "gene_symbol_mouse": "non-targeting", "biotype": "non_targeting",
            "guide_seq": s, "target_seq": "", "region": "non_targeting",
            "tiger_score": "", "cas13_score": "", "combined_score": "",
            "n_isoforms_targeted": 0, "single_isoform": "", "tx_id_set": "",
            "tx_id_pos": "", "position": "", "n_target": 0, "any_indel": "",
            "mismatch": "", "rank_within_gene": i, "n_available_pool": "",
            "human_symbol": "", "expected_direction": "neutral",
            "control_type": "non_targeting", "essentiality_chronos_liver": "",
            "n_liver_lines": "", "is_essential_liver": "",
        })
    df = pd.DataFrame(rows)[SCHEMA]
    df.to_csv(args.out, index=False)
    # invariants
    assert df["guide_seq"].str.len().eq(GUIDE_LEN).all(), "non-23bp NT guide"
    assert df["guide_seq"].is_unique, "duplicate NT guide_seq"
    print(f"[nt] wrote {len(df)} non-targeting guides -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
