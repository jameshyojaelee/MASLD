#!/usr/bin/env python3
"""Aggregate H3K27ac (or any interval-keyed chromatin signal) to GENCODE v49 gene-anchored windows.

Portable by construction: the gene-window axis is a deterministic function of (GENCODE release,
window sizes) and is hashed, so two laboratories that declare their coordinate convention and
genome build obtain the same ordered columns. The script REFUSES to run without a manifest that
declares both. All coordinates are converted once, at ingest, to 0-based half-open.

Windows per gene (protein_coding + lncRNA, primary contigs):
  promoter  = TSS +/- promoter_bp
  distal    = TSS +/- distal_bp, minus the promoter window
  body      = gene start..end (optional secondary axis)
A region overlapping several windows is assigned to EVERY one of them (no nearest-gene tie-break).
Library size is the sample's ORIGINAL region sum over all regions, never the window sum.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path

import numpy as np

REGION = re.compile(r"^(chr(?:[0-9]+|X|Y|M)):([0-9]+)-([0-9]+)$")
CONVENTIONS = ("zero_based_half_open", "one_based_inclusive")
BUILDS = ("GRCh38",)
GENE_TYPES = ("protein_coding", "lncRNA")


def to_half_open(start: int, end: int, convention: str) -> tuple[int, int]:
    """Convert one interval to 0-based half-open, once."""
    if convention == "one_based_inclusive":
        return start - 1, end
    if convention == "zero_based_half_open":
        return start, end
    raise ValueError(convention)


def parse_regions(labels: list[str], convention: str) -> tuple[list[str], np.ndarray, np.ndarray]:
    chroms, starts, ends = [], [], []
    for label in labels:
        m = REGION.match(label)
        if m is None:
            raise ValueError(f"region label not chr:start-end: {label}")
        s, e = to_half_open(int(m.group(2)), int(m.group(3)), convention)
        if e <= s:
            raise ValueError(f"empty or inverted interval after conversion: {label}")
        chroms.append(m.group(1)); starts.append(s); ends.append(e)
    return chroms, np.asarray(starts, np.int64), np.asarray(ends, np.int64)


def load_genes(gtf_gz: Path, gene_types=GENE_TYPES):
    """GENCODE gene records -> list of (gene_id, gene_name, gene_type, chrom, start0, end0, strand).
    GTF is 1-based inclusive; converted here, once."""
    attr_id = re.compile(r'gene_id "([^"]+)"')
    attr_name = re.compile(r'gene_name "([^"]+)"')
    attr_type = re.compile(r'gene_type "([^"]+)"')
    genes = []
    with gzip.open(gtf_gz, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            gt = attr_type.search(f[8])
            if gt is None or gt.group(1) not in gene_types:
                continue
            gid = attr_id.search(f[8]).group(1)
            gname_m = attr_name.search(f[8])
            gname = gname_m.group(1) if gname_m else gid
            s0, e0 = to_half_open(int(f[3]), int(f[4]), "one_based_inclusive")
            genes.append((gid, gname, gt.group(1), f[0], s0, e0, f[6]))
    genes.sort(key=lambda g: (g[3], g[4], g[0]))
    return genes


def build_windows(genes, promoter_bp: int, distal_bp: int):
    """Return ordered window records: (column_name, gene_id, kind, chrom, start0, end0)."""
    windows = []
    for gid, gname, gtype, chrom, s0, e0, strand in genes:
        tss = s0 if strand == "+" else e0 - 1          # 0-based position of the TSS base
        p_s, p_e = max(0, tss - promoter_bp), tss + promoter_bp + 1
        d_s, d_e = max(0, tss - distal_bp), tss + distal_bp + 1
        windows.append((f"{gid}|promoter", gid, "promoter", chrom, p_s, p_e, None))
        # distal = [d_s, d_e) minus [p_s, p_e): stored as the outer window with a hole
        windows.append((f"{gid}|distal", gid, "distal", chrom, d_s, d_e, (p_s, p_e)))
        windows.append((f"{gid}|body", gid, "body", chrom, s0, e0, None))
    return windows


def membership(chroms, starts, ends, windows):
    """Sparse membership: for each window, the region indices overlapping it.
    Returns dict column_name -> np.ndarray of region indices."""
    by_chrom: dict[str, list[int]] = {}
    for i, c in enumerate(chroms):
        by_chrom.setdefault(c, []).append(i)
    idx_by_chrom = {c: np.asarray(v) for c, v in by_chrom.items()}
    sorted_by_chrom = {}
    for c, idx in idx_by_chrom.items():
        order = np.argsort(starts[idx], kind="stable")
        sorted_by_chrom[c] = (idx[order], starts[idx][order], ends[idx][order])
    out = {}
    for name, gid, kind, chrom, ws, we, hole in windows:
        if chrom not in sorted_by_chrom:
            out[name] = np.empty(0, np.int64); continue
        idx, st, en = sorted_by_chrom[chrom]
        # overlap [ws,we) : region.start < we and region.end > ws
        hi = np.searchsorted(st, we, side="left")          # regions starting before we
        cand = idx[:hi]; cs, ce = st[:hi], en[:hi]
        keep = ce > ws
        if hole is not None:
            hs, he = hole
            # exclude regions lying ENTIRELY inside the promoter hole; regions straddling
            # the hole boundary count for both windows (assigned to every window overlapped)
            inside_hole = (cs >= hs) & (ce <= he)
            keep = keep & ~inside_hole
        out[name] = cand[keep]
    return out


def aggregate(counts: np.ndarray, member: dict[str, np.ndarray], columns: list[str],
              fractional: bool = False, region_multiplicity: np.ndarray | None = None):
    """G = counts @ W. counts is samples x regions. Optional fractional weighting splits a
    region's counts equally across the windows it belongs to."""
    n = counts.shape[0]
    G = np.zeros((n, len(columns)), np.float64)
    for j, name in enumerate(columns):
        r = member[name]
        if r.size == 0:
            continue
        if fractional and region_multiplicity is not None:
            G[:, j] = (counts[:, r] / region_multiplicity[r][None, :]).sum(1)
        else:
            G[:, j] = counts[:, r].sum(1)
    return G


def axis_sha256(columns: list[str]) -> str:
    return hashlib.sha256(("\n".join(columns) + "\n").encode("utf-8")).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, required=True,
                    help="JSON declaring coordinate_convention and genome_build; REQUIRED")
    ap.add_argument("--counts", type=Path, required=True, help="samples x regions .npy")
    ap.add_argument("--regions", type=Path, required=True, help="one chr:start-end per line, matrix order")
    ap.add_argument("--gtf", type=Path, required=True, help="GENCODE v49 primary GTF (gz)")
    ap.add_argument("--promoter-bp", type=int, default=1000)
    ap.add_argument("--distal-bp", type=int, default=50000)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    man = json.loads(args.manifest.read_text(encoding="utf-8"))
    conv, build = man.get("coordinate_convention"), man.get("genome_build")
    if conv not in CONVENTIONS or build not in BUILDS:
        raise SystemExit("REFUSED: manifest must declare coordinate_convention in "
                         f"{CONVENTIONS} and genome_build in {BUILDS}; got {conv!r}, {build!r}")
    args.output.mkdir(parents=True, exist_ok=False)

    labels = [l.strip() for l in args.regions.read_text().splitlines() if l.strip()]
    counts = np.load(args.counts).astype(np.float64)
    if counts.shape[1] != len(labels):
        raise SystemExit(f"counts has {counts.shape[1]} columns, regions has {len(labels)}")
    chroms, starts, ends = parse_regions(labels, conv)
    genes = load_genes(args.gtf)
    windows = build_windows(genes, args.promoter_bp, args.distal_bp)
    member = membership(chroms, starts, ends, windows)
    two_window_cols = [w[0] for w in windows if w[2] in ("promoter", "distal")]
    body_cols = [w[0] for w in windows if w[2] == "body"]
    mult = np.zeros(len(labels), np.int64)
    for name in two_window_cols:
        mult[member[name]] += 1
    lib = counts.sum(1)                                         # ORIGINAL region library size

    G2 = aggregate(counts, member, two_window_cols)
    G2f = aggregate(counts, member, two_window_cols, fractional=True, region_multiplicity=np.maximum(mult, 1))
    Gb = aggregate(counts, member, body_cols)
    np.save(args.output / "gene_two_window_counts.npy", G2)
    np.save(args.output / "gene_two_window_fractional_counts.npy", G2f)
    np.save(args.output / "gene_body_counts.npy", Gb)
    np.save(args.output / "library_size_original_region_sum.npy", lib)
    observed2 = np.array([member[c].size > 0 for c in two_window_cols])
    (args.output / "gene_two_window_axis.tsv").write_text(
        "column\tgene_id\tkind\twindow_observed\n" + "\n".join(
            f"{c}\t{c.split('|')[0]}\t{c.split('|')[1]}\t{bool(o)}" for c, o in zip(two_window_cols, observed2)) + "\n")
    (args.output / "gene_body_axis.tsv").write_text("column\n" + "\n".join(body_cols) + "\n")
    receipt = {
        "schema_version": "masld-bench-gene-window-aggregation-v1",
        "coordinate_convention_declared": conv, "genome_build_declared": build,
        "gtf": str(args.gtf), "gtf_sha256": hashlib.sha256(args.gtf.read_bytes()).hexdigest(),
        "gene_types": list(GENE_TYPES), "n_genes": len(genes),
        "promoter_bp": args.promoter_bp, "distal_bp": args.distal_bp,
        "declared_two_window_columns": len(two_window_cols),
        "nonempty_two_window_columns": int(observed2.sum()),
        "nonempty_promoter_columns": int(sum(member[c].size > 0 for c in two_window_cols if c.endswith("|promoter"))),
        "nonempty_distal_columns": int(sum(member[c].size > 0 for c in two_window_cols if c.endswith("|distal"))),
        "nonempty_body_columns": int(sum(member[c].size > 0 for c in body_cols)),
        "regions_total": len(labels),
        "regions_in_no_two_window": int((mult == 0).sum()),
        "regions_in_multiple_two_windows": int((mult > 1).sum()),
        "library_size_rule": "original region sum over all regions, never the window sum",
        "multi_window_rule": "a region is assigned to every window it overlaps; no nearest-gene tie-break",
        "two_window_axis_sha256": axis_sha256(two_window_cols),
        "body_axis_sha256": axis_sha256(body_cols),
        "n_samples": int(counts.shape[0]),
    }
    (args.output / "aggregation_receipt.json").write_text(json.dumps(receipt, indent=2))
    (args.output / "axis.sha256").write_text(receipt["two_window_axis_sha256"] + "\n")
    print(json.dumps({k: receipt[k] for k in ("n_genes", "declared_two_window_columns",
          "nonempty_two_window_columns", "nonempty_promoter_columns", "nonempty_distal_columns",
          "regions_in_no_two_window", "regions_in_multiple_two_windows", "two_window_axis_sha256")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
