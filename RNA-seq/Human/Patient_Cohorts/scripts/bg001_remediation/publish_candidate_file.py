#!/usr/bin/env python3
"""Atomically publish one already-written same-parent candidate temporary file."""

from __future__ import annotations

import argparse
import os
import stat
from pathlib import Path

from safe_io import require_safe_parent


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--temporary", required=True, type=Path)
    parser.add_argument("--destination", required=True, type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve(strict=True)
    temporary = args.temporary
    destination = args.destination
    require_safe_parent(destination, root)
    require_safe_parent(temporary, root)
    if temporary.parent.resolve(strict=True) != destination.parent.resolve(strict=True):
        raise SystemExit("Candidate temporary and destination must have the same parent")
    mode = os.lstat(temporary).st_mode
    if stat.S_ISLNK(mode) or not stat.S_ISREG(mode):
        raise SystemExit("Candidate temporary is symlinked or non-regular")
    if not temporary.name.startswith(f".{destination.name}.") or not temporary.name.endswith(".tmp"):
        raise SystemExit("Candidate temporary name does not bind its destination")
    try:
        os.lstat(destination)
    except FileNotFoundError:
        pass
    else:
        raise SystemExit("Refusing existing candidate destination")
    descriptor = os.open(temporary, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.link(temporary, destination, follow_symlinks=False)
    os.unlink(temporary)
    directory_fd = os.open(destination.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)
    print(destination)


if __name__ == "__main__":
    main()
