#!/usr/bin/env python3
"""V2 target screen, step 2: mask the target references.

Two masks are applied, both fixed in PRESPEC_V2.md before any library is run.

1. Low-complexity mask. A 32-bp window slides with step 1 over each reference.
   Every base of the window is masked if any of:
     (a) the most frequent nucleotide occupies >= 0.80 of the window
     (b) the two most frequent nucleotides together occupy >= 0.92
     (c) the Shannon entropy of the window's 3-mer multiset is < 3.0 bits
   A window containing any non-ACGT base is masked outright.

2. Human-homology mask. Tiles of 50 bp at step 10 are cut from the unmasked
   reference and aligned end-to-end to the human genome with HISAT2 using the
   aligner's default end-to-end score floor. Every target base covered by a tile
   that aligned anywhere in the human genome is masked.

Masked bases are written as N, so they cannot carry an alignment and are
excluded from the breadth-of-coverage denominator.

Subcommands:
  tiles  <refs.fasta> <tiles.fasta>              write the tiles for step 2
  apply  <refs.fasta> <tiles.sam> <out.fasta> <report.tsv>
"""

from __future__ import annotations

import collections
import math
import sys

WINDOW = 32
STEP = 1
MAX1_FRAC = 0.80
MAX2_FRAC = 0.92
MIN_3MER_ENTROPY_BITS = 3.0
TILE_LEN = 50
TILE_STEP = 10


def read_fasta(path: str) -> list[tuple[str, str]]:
    recs: list[tuple[str, str]] = []
    name = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name = line[1:].split()[0]
                buf = []
            elif line:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def kmer_entropy_bits(window: str, k: int = 3) -> float:
    counts = collections.Counter(window[i:i + k] for i in range(len(window) - k + 1))
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return -sum((c / total) * math.log2(c / total) for c in counts.values())


def low_complexity_mask(seq: str) -> bytearray:
    n = len(seq)
    mask = bytearray(n)
    up = seq.upper()
    for start in range(0, max(1, n - WINDOW + 1), STEP):
        win = up[start:start + WINDOW]
        if len(win) < WINDOW:
            break
        counts = collections.Counter(win)
        bad_bases = sum(v for b, v in counts.items() if b not in "ACGT")
        if bad_bases:
            for i in range(start, start + WINDOW):
                mask[i] = 1
            continue
        freqs = sorted((v / WINDOW for v in counts.values()), reverse=True)
        top1 = freqs[0]
        top2 = top1 + (freqs[1] if len(freqs) > 1 else 0.0)
        if (top1 >= MAX1_FRAC or top2 >= MAX2_FRAC
                or kmer_entropy_bits(win) < MIN_3MER_ENTROPY_BITS):
            for i in range(start, start + WINDOW):
                mask[i] = 1
    return mask


def cmd_tiles(ref_path: str, out_path: str) -> int:
    with open(out_path, "w") as out:
        for name, seq in read_fasta(ref_path):
            for start in range(0, max(1, len(seq) - TILE_LEN + 1), TILE_STEP):
                tile = seq[start:start + TILE_LEN]
                if len(tile) < TILE_LEN:
                    break
                out.write(">%s|%d\n%s\n" % (name, start, tile))
    return 0


def cmd_apply(ref_path: str, sam_path: str, out_fasta: str,
              report_path: str) -> int:
    refs = read_fasta(ref_path)
    homology: dict[str, set[int]] = {name: set() for name, _ in refs}
    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            flag = int(f[1])
            if flag & 0x4:        # unmapped: no human homology from this tile
                continue
            qname = f[0]
            if "|" not in qname:
                continue
            contig, start_s = qname.rsplit("|", 1)
            if contig not in homology:
                continue
            start = int(start_s)
            homology[contig].update(range(start, start + TILE_LEN))

    rows = []
    with open(out_fasta, "w") as out:
        for name, seq in refs:
            lc = low_complexity_mask(seq)
            hm = homology[name]
            masked = bytearray(lc)
            for i in hm:
                if 0 <= i < len(masked):
                    masked[i] = 1
            n_lc = sum(lc)
            n_hm = len([i for i in hm if i < len(seq)])
            n_both = sum(1 for i in range(len(seq)) if lc[i] and i in hm)
            n_masked = sum(masked)
            out.write(">%s\n" % name)
            mseq = "".join("N" if masked[i] else seq[i] for i in range(len(seq)))
            for i in range(0, len(mseq), 70):
                out.write(mseq[i:i + 70] + "\n")
            rows.append({
                "contig": name, "length": len(seq),
                "masked_low_complexity": n_lc,
                "masked_human_homology": n_hm,
                "masked_both": n_both,
                "masked_total": n_masked,
                "usable_length": len(seq) - n_masked,
                "usable_fraction": round((len(seq) - n_masked) / len(seq), 4),
            })

    cols = ["contig", "length", "masked_low_complexity",
            "masked_human_homology", "masked_both", "masked_total",
            "usable_length", "usable_fraction"]
    with open(report_path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")
    return 0


def main() -> int:
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__ or "")
        return 2
    if sys.argv[1] == "tiles":
        return cmd_tiles(sys.argv[2], sys.argv[3])
    if sys.argv[1] == "apply":
        return cmd_apply(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    sys.stderr.write("unknown subcommand %r\n" % sys.argv[1])
    return 2


if __name__ == "__main__":
    sys.exit(main())
