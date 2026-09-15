#!/usr/bin/env python
"""09a: put H3K27ac on a gene axis, by reusing the project's own aggregator.

The chromatin matrix is keyed by interval, not by gene, so it cannot be scored
with a gene module until the intervals are assigned to genes. That assignment
already exists in this repository as a portable, hashed, convention-declaring
script; this step prepares its three inputs and calls it rather than writing a
second, subtly different definition of "the promoter of a gene".

GSE267145 intervals are 1-based inclusive, recorded in
docs/decisions/2026-09-06-gse267145-coordinate-convention.md. The aggregator
refuses to run unless that convention is declared, which is the point.
"""
import argparse
import gzip
import json
import pathlib
import subprocess
import sys

import numpy as np

ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CONTRACT = json.loads((ROOT / "scripts/analysis/cross_assay_modules/00_contract.json").read_text())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-root", required=True, type=pathlib.Path)
    ap.add_argument("--promoter-bp", type=int, default=1000)
    args = ap.parse_args()

    stage = args.out_root / "chromatin" / "h3k27ac_input"
    stage.mkdir(parents=True, exist_ok=True)
    agg_out = args.out_root / "chromatin" / "h3k27ac_gene_windows"
    if agg_out.exists():
        print(f"[09a] {agg_out} exists; refusing to overwrite", flush=True)
        return 0

    src = ROOT / CONTRACT["inputs"]["h3k27ac_matrix"]
    with gzip.open(src, "rt") as fh:
        header = fh.readline().split()
        samples = [h.strip('"') for h in header]
        regions, rows = [], []
        for line in fh:
            parts = line.split()
            regions.append(parts[0].strip('"'))
            rows.append(np.asarray(parts[1:], dtype=np.float64))
    counts = np.vstack(rows).T  # samples x regions, the layout the aggregator wants
    print(f"[09a] {counts.shape[0]} samples x {counts.shape[1]} regions", flush=True)

    np.save(stage / "counts.npy", counts)
    (stage / "regions.txt").write_text("\n".join(regions) + "\n")
    (stage / "samples.txt").write_text("\n".join(samples) + "\n")
    (stage / "manifest.json").write_text(json.dumps({
        "coordinate_convention": "one_based_inclusive",
        "genome_build": "GRCh38",
        "source": str(src),
        "note": "convention recorded in docs/decisions/2026-09-06-gse267145-coordinate-convention.md",
    }, indent=2))

    cmd = [
        sys.executable,
        str(ROOT / "Analysis/MASLD_Model_Benchmark/scripts/aggregate_h3k27ac_to_gene_windows.py"),
        "--manifest", str(stage / "manifest.json"),
        "--counts", str(stage / "counts.npy"),
        "--regions", str(stage / "regions.txt"),
        "--gtf", CONTRACT["inputs"]["gencode_gtf"],
        "--promoter-bp", str(args.promoter_bp),
        "--output", str(agg_out),
    ]
    print("[09a] " + " ".join(cmd), flush=True)
    return subprocess.run(cmd, check=True).returncode


if __name__ == "__main__":
    sys.exit(main())
