#!/usr/bin/env python3
"""Validate visually corrected r3 of the full-composition panel."""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


def import_script(name: str, module_name: str):
    path = SCRIPT_DIR / name
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import script: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RENDERER = import_script("41_render_full_composition_figure_r3.py", "full_composition_figure_r3")
VALIDATOR = import_script("40_validate_full_composition_figure.py", "full_composition_figure_validator_base")
VALIDATOR.import_renderer = lambda _script_dir: RENDERER


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, default=None)
    parser.add_argument("--output-root", type=Path, default=None)
    args = parser.parse_args()
    try:
        VALIDATOR.validate(args.project_root, args.input_root, args.output_root)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
