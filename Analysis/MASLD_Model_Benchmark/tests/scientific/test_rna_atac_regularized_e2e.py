#!/usr/bin/env python
"""Exercise regularized RNA-to-ATAC baselines on one real donor fold."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, PACKAGE_ROOT.as_posix())

from tests.scientific import test_rna_atac_classical_e2e as harness


ADAPTER_PATH = (
    PACKAGE_ROOT / "src" / "masld_bench" / "adapters" / "rna_atac_regularized.py"
)
SPEC = importlib.util.spec_from_file_location(
    "masld_bench_standalone_rna_atac_regularized", ADAPTER_PATH
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load standalone scientific adapter: {ADAPTER_PATH}")
adapter = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = adapter
SPEC.loader.exec_module(adapter)
harness.adapter = adapter


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--environment-lock", required=True, type=Path)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    print(
        json.dumps(
            harness.run_test(
                environment_lock=arguments.environment_lock,
                data_path=arguments.data,
                output_root=arguments.output,
            ),
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
