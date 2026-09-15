#!/usr/bin/env python
"""10a: collapse GSE296875 ATAC peaks onto gene promoters.

A peak is assigned to a gene when it overlaps that gene's TSS +/- promoter_bp.
A peak overlapping several promoters is assigned to all of them; there is no
nearest-gene tie-break, which would silently invent a choice the data does not
make. This matches the rule the H3K27ac aggregator applies, so the two
chromatin columns of the panel mean the same thing.
"""
import argparse
import gzip
import json
import pathlib
import sys
from collections import defaultdict

import numpy as np

ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CONTRACT = json.loads((ROOT / "scripts/analysis/cross_assay_modules/00_contract.json").read_text())
GENE_TYPES = ("protein_coding", "lncRNA")


def load_tss(gtf: pathlib.Path):
    """TSS per gene from the GENCODE gene rows, strand-aware."""
    out = []
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attr = f[8]
            if not any(f'gene_type "{t}"' in attr for t in GENE_TYPES):
                continue
            name = attr.split('gene_name "', 1)
            if len(name) < 2:
                continue
            symbol = name[1].split('"', 1)[0]
            start, end, strand = int(f[3]) - 1, int(f[4]), f[6]
            tss = start if strand == "+" else end - 1
            out.append((f[0], tss, symbol))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-root", required=True, type=pathlib.Path)
    ap.add_argument("--promoter-bp", type=int, default=1000)
    args = ap.parse_args()

    out = args.out_root / "chromatin" / "atac_promoter"
    if (out / "gene_promoter_counts.npy").exists():
        print("[10a] already built; refusing to overwrite", flush=True)
        return 0
    out.mkdir(parents=True, exist_ok=True)

    donor_dir = ROOT / CONTRACT["inputs"]["gse296875_donor_dir"]
    peaks = [p.strip() for p in (donor_dir / "atac_peaks.txt").read_text().splitlines() if p.strip()]
    atac = np.load(donor_dir / "donor_atac_all.npy").astype(np.float64)
    print(f"[10a] ATAC {atac.shape[0]} donors x {atac.shape[1]} peaks", flush=True)
    assert atac.shape[1] == len(peaks), "peak axis and matrix disagree"

    # Peak labels are chr-start-end, 0-based half-open as deposited.
    by_chrom = defaultdict(list)
    for j, p in enumerate(peaks):
        chrom, start, end = p.rsplit("-", 2)
        by_chrom[chrom].append((int(start), int(end), j))
    for c in by_chrom:
        by_chrom[c].sort()
    starts = {c: np.asarray([v[0] for v in by_chrom[c]]) for c in by_chrom}

    tss = load_tss(pathlib.Path(CONTRACT["inputs"]["gencode_gtf"]))
    print(f"[10a] {len(tss)} gene TSS loaded", flush=True)

    gene_to_peaks = defaultdict(set)
    w = args.promoter_bp
    for chrom, pos, symbol in tss:
        if chrom not in by_chrom:
            continue
        lo, hi = pos - w, pos + w + 1
        arr = by_chrom[chrom]
        s = starts[chrom]
        # peaks starting before hi; walk back while they can still reach lo
        i = int(np.searchsorted(s, hi, side="right"))
        k = i - 1
        while k >= 0:
            ps, pe, idx = arr[k]
            if pe > lo:
                gene_to_peaks[symbol].add(idx)
            if ps < lo - 100000:
                break
            k -= 1

    symbols = sorted(gene_to_peaks)
    mat = np.zeros((atac.shape[0], len(symbols)), dtype=np.float64)
    for gi, sym in enumerate(symbols):
        cols = sorted(gene_to_peaks[sym])
        mat[:, gi] = atac[:, cols].sum(axis=1)
    lib = atac.sum(axis=1)

    np.save(out / "gene_promoter_counts.npy", mat)
    np.save(out / "library_size_all_peaks.npy", lib)
    (out / "gene_axis.txt").write_text("\n".join(symbols) + "\n")
    (out / "receipt.json").write_text(json.dumps({
        "promoter_bp": w,
        "n_donors": int(mat.shape[0]),
        "n_genes_with_promoter_peak": len(symbols),
        "n_peaks_total": len(peaks),
        "n_peaks_in_any_promoter": int(len({i for s in gene_to_peaks.values() for i in s})),
        "multi_gene_rule": "a peak is assigned to every promoter it overlaps; no nearest-gene tie-break",
        "library_size_rule": "sum over all peaks, never the promoter subset",
        "peak_coordinate_convention": "zero_based_half_open_as_deposited",
    }, indent=2))
    print(f"[10a] {len(symbols)} genes with a promoter peak", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
