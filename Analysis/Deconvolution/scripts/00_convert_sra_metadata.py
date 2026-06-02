#!/usr/bin/env python3
"""Convert SRA RunTable CSV into a TSV with sample_id for deconvolution."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input SRA RunTable CSV")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("--sample-col", default="Run", help="Column to use as sample_id")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    if args.sample_col not in df.columns:
        raise SystemExit(f"Missing column '{args.sample_col}' in {args.input}")

    df.insert(0, "sample_id", df[args.sample_col].astype(str))

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
