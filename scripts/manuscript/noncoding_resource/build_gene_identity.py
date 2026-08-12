#!/usr/bin/env python3
"""Command-line entry point for the GENCODE v49 identity candidate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from gene_identity import (
    ANNOTATION_RELEASE,
    build_release,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    census = build_release(
        args.gtf,
        args.output_dir,
        annotation_release=ANNOTATION_RELEASE,
    )
    print(json.dumps(census, sort_keys=True))


if __name__ == "__main__":
    main()
