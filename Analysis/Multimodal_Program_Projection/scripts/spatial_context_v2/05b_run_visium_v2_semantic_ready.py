#!/usr/bin/env python3
"""Run frozen v2 Visium projection after the validation-only v1 hotfix.

The original candidate producer remains byte-identical.  This versioned entry
point changes only how the v1 READY seal resolves its validation report; it
then delegates the complete scientific computation to the frozen producer and
byte-pinned v1 engine.  Setting the delegated module's ``__file__`` makes the
new entry point, rather than the old producer, appear in output provenance.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

from visium_rerun_lib import (
    RELEASE_ID,
    V1_ENGINE_SHA256,
    build_paths,
    read_tsv,
    sha256_file,
)


ORIGINAL_PRODUCER_SHA256 = (
    "5172f92d0ae0edc4d2f5cb1c2ef35faf813ca87f5b8e2f05529f583666870052"
)


def load_original_producer(script_root: Path):
    path = script_root / "05_run_visium_rerun.py"
    observed = sha256_file(path)
    if observed != ORIGINAL_PRODUCER_SHA256:
        raise RuntimeError(
            f"original candidate producer drift: expected {ORIGINAL_PRODUCER_SHA256}, "
            f"observed {observed}"
        )
    spec = importlib.util.spec_from_file_location("_frozen_visium_v2_producer", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not import frozen v2 producer")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def semantic_ready_check(paths) -> str:
    root = paths.native_root / "v1_regression"
    ready_path = root / "READY"
    rows = read_tsv(ready_path)
    if len(rows) != 1:
        raise RuntimeError("v1 regression READY must contain exactly one row")
    row = rows[0]
    expected = {
        "release_id": RELEASE_ID,
        "registry_version": "v1",
        "status": "pass_v1_regression",
        "v1_engine_sha256": V1_ENGINE_SHA256,
    }
    for key, value in expected.items():
        if row.get(key) != value:
            raise RuntimeError(f"v1 semantic READY {key} mismatch")

    report_name = row.get("validation_report_path", "validation_report.tsv")
    report_path = (root / report_name).resolve()
    if report_path.parent != root.resolve() or not report_path.is_file():
        raise RuntimeError("v1 semantic READY validation-report path is unsafe or missing")
    if row.get("validation_report_sha256") != sha256_file(report_path):
        raise RuntimeError("v1 semantic READY validation-report hash mismatch")

    strict_report = root / "validation_report.tsv"
    if row.get("strict_failed_report_sha256") != sha256_file(strict_report):
        raise RuntimeError("v1 semantic READY does not pin the original failed report")
    hotfix_manifest = root / "validator_hotfix_manifest.tsv"
    if row.get("validator_hotfix_manifest_sha256") != sha256_file(hotfix_manifest):
        raise RuntimeError("v1 semantic READY hotfix-manifest hash mismatch")
    hotfix_rows = read_tsv(hotfix_manifest)
    if (
        len(hotfix_rows) != 1
        or hotfix_rows[0].get("scope") != "validation_only"
        or hotfix_rows[0].get("results_recomputed") != "FALSE"
        or hotfix_rows[0].get("semantic_report_sha256") != sha256_file(report_path)
    ):
        raise RuntimeError("v1 semantic validation hotfix contract is invalid")
    return sha256_file(ready_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=None)
    args = parser.parse_args()
    try:
        paths = build_paths(args.project_root)
        semantic_ready_check(paths)
        script = Path(__file__).resolve()
        producer = load_original_producer(script.parent)
        producer.require_v1_regression_ready = semantic_ready_check
        # The delegated run uses module.__file__ for the producer manifest.
        producer.__file__ = str(script)
        producer.run(paths, "v2")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
