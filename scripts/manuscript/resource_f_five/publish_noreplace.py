#!/usr/bin/env python3
"""Publish one completed directory with an atomic no-replace commit marker."""

from __future__ import annotations

import argparse
from pathlib import Path

from snapshot_io import publish_directory_noreplace


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--marker", required=True)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    publish_directory_noreplace(
        args.source,
        args.destination,
        marker=args.marker,
    )


if __name__ == "__main__":
    main()
