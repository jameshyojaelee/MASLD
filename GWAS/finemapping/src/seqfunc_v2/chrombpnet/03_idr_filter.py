#!/usr/bin/env python3
"""IDR pseudoreplicates and retain pooled peaks by summit containment."""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import os
import subprocess
from pathlib import Path

ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
OUT = ROOT / "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2"
ENV = ROOT / ".mamba/seqfunc_chrombpnet"


def read_narrowpeak(path: Path) -> list[list[str]]:
    rows = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                raise SystemExit(f"malformed narrowPeak row: {path}")
            rows.append(fields)
    if not rows:
        raise SystemExit(f"empty narrowPeak: {path}")
    return rows


def merge_intervals(intervals: list[tuple[int, int]]) -> tuple[list[int], list[int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [x[0] for x in merged], [x[1] for x in merged]


def overlaps(index: tuple[list[int], list[int]], point: int) -> bool:
    starts, ends = index
    i = bisect.bisect_right(starts, point) - 1
    return i >= 0 and point < ends[i]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--idr", type=float, default=0.05)
    args = ap.parse_args()
    out = args.out.resolve()
    peak_root = out / "peaks"
    rep1 = peak_root / "pseudorep1/hepatocyte.pseudorep1.autosomal.blacklist_filtered.narrowPeak"
    rep2 = peak_root / "pseudorep2/hepatocyte.pseudorep2.autosomal.blacklist_filtered.narrowPeak"
    pooled = peak_root / "pooled/hepatocyte.pooled.autosomal.blacklist_filtered.narrowPeak"
    for path in (rep1, rep2, pooled):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"missing peak input: {path}")
    idr_dir = out / "idr"
    idr_dir.mkdir(exist_ok=True)
    sorted_inputs = []
    for i, path in enumerate((rep1, rep2), 1):
        rows = read_narrowpeak(path)
        # signalValue is narrowPeak column 7 (index 6); deterministic coordinate
        # tie-breaks make the IDR input byte-reproducible.
        rows.sort(key=lambda r: (-float(r[6]), r[0], int(r[1]), int(r[2]), r[3]))
        dest = idr_dir / f"pseudorep{i}.signal_ranked.narrowPeak"
        with dest.open("w") as handle:
            for row in rows:
                handle.write("\t".join(row[:10]) + "\n")
        sorted_inputs.append(dest)
    raw_idr = idr_dir / "pseudoreplicates.idr.raw.tsv"
    compat = Path(__file__).with_name("idr_numpy_compat.py")
    cmd = [
        str(ENV / "bin/python"), str(compat),
        "--samples", str(sorted_inputs[0]), str(sorted_inputs[1]),
        "--input-file-type", "narrowPeak",
        "--rank", "signal.value",
        "--soft-idr-threshold", str(args.idr),
        "--output-file", str(raw_idr),
        "--log-output-file", str(idr_dir / "idr.log"),
    ]
    subprocess.run(cmd, check=True)
    if not raw_idr.is_file() or raw_idr.stat().st_size == 0:
        raise SystemExit("IDR produced no output")

    cutoff = -math.log10(args.idr)
    reproducible = []
    with raw_idr.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise SystemExit("IDR output lacks the global-IDR column")
            global_idr = 10 ** (-float(fields[11]))
            if global_idr <= args.idr:
                reproducible.append((fields[0], int(fields[1]), int(fields[2]), global_idr, fields))
    if not reproducible:
        raise SystemExit("no pseudoreplicate intervals pass global IDR <= 0.05")
    reproducible.sort(key=lambda x: (x[0], x[1], x[2]))
    repro_bed = idr_dir / "pseudoreplicates.idr_0.05.bed"
    with repro_bed.open("w") as handle:
        for chrom, start, end, gidr, fields in reproducible:
            handle.write(f"{chrom}\t{start}\t{end}\t{fields[3]}\t{gidr:.8g}\n")

    by_chrom_raw: dict[str, list[tuple[int, int]]] = {}
    for chrom, start, end, _, _ in reproducible:
        by_chrom_raw.setdefault(chrom, []).append((start, end))
    by_chrom = {chrom: merge_intervals(intervals) for chrom, intervals in by_chrom_raw.items()}
    pooled_rows = read_narrowpeak(pooled)
    retained = []
    for row in pooled_rows:
        summit = int(row[1]) + int(row[9])
        if row[0] in by_chrom and overlaps(by_chrom[row[0]], summit):
            retained.append(row)
    if not retained:
        raise SystemExit("no pooled MACS3 peaks have summits in reproducible IDR intervals")
    retained.sort(key=lambda r: (r[0], int(r[1]), int(r[2])))
    final = out / "hepatocyte.idr_pooled_summit.narrowPeak"
    with final.open("w") as handle:
        for row in retained:
            handle.write("\t".join(row[:10]) + "\n")

    qc = {
        "pseudorep1_filtered_peaks": sum(1 for _ in rep1.open()),
        "pseudorep2_filtered_peaks": sum(1 for _ in rep2.open()),
        "idr_tested_pairs": sum(1 for _ in raw_idr.open()),
        "idr_reproducible_intervals": len(reproducible),
        "pooled_filtered_peaks": len(pooled_rows),
        "pooled_peaks_summit_in_idr": len(retained),
        "idr_threshold": args.idr,
        "idr_rank": "signal.value",
        "final_peak_rule": "pooled MACS3 summit contained in pseudoreplicate IDR interval",
    }
    with (out / "peak_reproducibility_qc.json").open("w") as handle:
        json.dump(qc, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(qc, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
