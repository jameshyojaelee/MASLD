#!/usr/bin/env python
"""
build_guide_index.py — one-time converter: upstream genome-wide guide CSV
(~3.8 GB, ~12.6M rows) -> a slim, sorted parquet index for fast per-gene queries.

The upstream file (sfriedman track-cas13) is read-only. We never modify it; we
build a derived parquet under Cas13_Library_Design/data/guides/cache/ keyed by an
unversioned `gene_id_base` (ENSMUSG without the .N suffix) so it joins directly
against the library's `gene_id_mouse`. Sorting by gene_id_base lets parquet
row-group statistics prune the file on point queries (dynamic puller) while the
full-library build does one streaming `is_in` scan.

Run on a COMPUTE node (see build_guide_index.sbatch) — do NOT run on the login node.

Usage:
    python build_guide_index.py --release vM38
    python build_guide_index.py --release vM37 --streaming   # low-memory, unsorted
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import guides_config as cfg


def build(release: str, streaming: bool) -> Path:
    src = cfg.RELEASES[release]["genome"]
    if not src.exists():
        sys.exit(f"ERROR: source guide CSV not found: {src}")
    cfg.CACHE_DIR.mkdir(parents=True, exist_ok=True)
    out = cfg.index_path(release)

    print(f"[index] release={release}")
    print(f"[index] source : {src}")
    print(f"[index] output : {out}")
    print(f"[index] mode   : {'streaming (unsorted)' if streaming else 'eager (sorted)'}")
    t0 = time.time()

    # All columns read as Utf8 (infer_schema_length=0) — the upstream CSV mixes
    # numeric scores, sentinel -1.0 expression values, and ':'-delimited strings;
    # forcing string avoids inference errors. Cast the numeric columns we need.
    lazy = (
        pl.scan_csv(
            src,
            separator=",",
            infer_schema_length=0,
            quote_char='"',
        )
        .select(cfg.INDEX_COLS)
        .with_columns(
            pl.col("gene_id").str.split(".").list.first().alias("gene_id_base"),
            pl.col("tiger_score").cast(pl.Float64, strict=False),
            pl.col("cas13_score").cast(pl.Float64, strict=False),
            pl.col("combined_score").cast(pl.Float64, strict=False),
            pl.col("n_target").cast(pl.Int32, strict=False),
        )
    )

    if streaming:
        # Out-of-core sink, unsorted. Robust for low-memory nodes; point queries
        # still work via a full predicate-pushdown scan (a few seconds).
        lazy.sink_parquet(out, compression="zstd", row_group_size=131072)
    else:
        # Eager: sort by gene_id_base so row-group min/max stats prune lookups.
        df = lazy.collect()
        n = df.height
        df = df.sort("gene_id_base")
        df.write_parquet(out, compression="zstd", row_group_size=131072,
                         statistics=True)
        print(f"[index] rows={n:,}  genes={df['gene_id_base'].n_unique():,}")

    dt = time.time() - t0
    size_gb = out.stat().st_size / 1e9
    print(f"[index] DONE in {dt:.1f}s  ->  {out}  ({size_gb:.2f} GB)")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--release", default=cfg.DEFAULT_RELEASE,
                    choices=list(cfg.RELEASES.keys()))
    ap.add_argument("--streaming", action="store_true",
                    help="Low-memory out-of-core sink (unsorted parquet).")
    args = ap.parse_args()
    build(args.release, args.streaming)


if __name__ == "__main__":
    main()
