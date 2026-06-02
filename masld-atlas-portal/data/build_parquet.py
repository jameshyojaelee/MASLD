"""Build ``atlas.parquet`` from the source multi-evidence atlas CSV.

Usage (from the ``masld-atlas-portal`` directory):

    python -m data.build_parquet

The resulting file is written next to this package as
``masld-atlas-portal/atlas.parquet`` and takes precedence over the copy
served by the Next.js site (``masld-atlas-v2/public/data/atlas.parquet``).
Re-run whenever ``multi_evidence_atlas.csv`` changes.

Size target: ~50 MB (string columns dictionary-encoded, Snappy compression).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.csv as pcsv
import pyarrow.parquet as pq


DEFAULT_SOURCE_CSV = (
    Path(__file__).resolve().parent.parent.parent
    / "RNA-seq"
    / "results"
    / "multi_evidence"
    / "multi_evidence_atlas.csv"
)
DEFAULT_OUTPUT = Path(__file__).resolve().parent.parent / "atlas.parquet"


def build(source: Path, output: Path) -> Path:
    if not source.exists():
        raise FileNotFoundError(f"Source CSV not found: {source}")

    print(f"[build_parquet] reading {source}", file=sys.stderr)
    # Use pyarrow's native CSV reader (faster than pandas and robust to the
    # wide 323-column atlas that trips pandas 2.3).
    table = pcsv.read_csv(
        source,
        read_options=pcsv.ReadOptions(block_size=1 << 24),
        parse_options=pcsv.ParseOptions(newlines_in_values=False),
    )
    n_string_cols = sum(1 for f in table.schema if pa.types.is_string(f.type))
    print(
        f"[build_parquet] {table.num_rows:,} rows x {table.num_columns:,} columns "
        f"({n_string_cols} string cols)",
        file=sys.stderr,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(
        table,
        output,
        compression="snappy",
        use_dictionary=True,
        write_statistics=True,
    )
    size_mb = output.stat().st_size / 1024 / 1024
    print(f"[build_parquet] wrote {output} ({size_mb:.1f} MB)", file=sys.stderr)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description="Build atlas.parquet")
    parser.add_argument(
        "--source", type=Path, default=DEFAULT_SOURCE_CSV,
        help="Path to multi_evidence_atlas.csv",
    )
    parser.add_argument(
        "--output", type=Path, default=DEFAULT_OUTPUT,
        help="Output parquet path",
    )
    args = parser.parse_args()
    build(args.source, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
